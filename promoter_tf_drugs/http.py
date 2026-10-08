"""Tiny cached HTTP client. Stdlib only, so the package needs only numpy.

Every external resource used here is a public API with rate limits. Caching to
disk means a re-run of the same gene costs nothing and you stop hammering
Ensembl/ChEMBL while you iterate on thresholds.
"""

from __future__ import annotations

import gzip
import hashlib
import random
import http.client
import json
import logging
import os
import socket
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

log = logging.getLogger(__name__)

DEFAULT_CACHE = Path(os.environ.get("PTD_CACHE", Path.home() / ".cache" / "promoter_tf_drugs"))
USER_AGENT = "promoter-tf-drugs/0.1 (research pipeline)"

# Retries, and why they are this patient. Ensembl, ChEMBL and the
# maayanlab services all return HTTP 500 when they are merely busy, so a
# 500 is a backoff signal rather than a verdict. Six attempts with capped
# jittered backoff spans roughly a minute, which rides out the outages
# seen in practice; the old three attempts over 7s did not, and giving up
# there discards every lookup the run had already completed.
DEFAULT_RETRIES = int(os.environ.get("PTD_RETRIES", 6))
RETRY_CAP = float(os.environ.get("PTD_RETRY_CAP", 20.0))


class HttpError(RuntimeError):
    pass


class Client:
    def __init__(self, cache_dir: Path | str = DEFAULT_CACHE, ttl_days: float = 30,
                 offline: bool = False, min_interval: float = 0.15,
                 max_requests: int | None = None):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.ttl = ttl_days * 86400
        self.offline = offline
        self.min_interval = min_interval
        self._last_call = 0.0
        # How many requests this run may send, or None for no ceiling.
        #
        # It counts requests that LEAVE. A cache hit is free and stays
        # free, or a resumed run would behave differently from a fresh
        # one for no reason the caller could see.
        #
        # Retries do not count: one attempt sequence against one URL is
        # one question, however many times the connection had to be
        # made before an answer came back. A redirect DOES count, and
        # the difference is deliberate - a 307 to a new host is a second
        # round trip, and a server that redirects repeatedly inside
        # max_redirects produces real traffic on every hop. A budget
        # that declined to count that would be measuring intent rather
        # than cost, which is not what a ceiling is for.
        #
        # Running out raises HttpError - the same exception a dead host
        # raises - so every `_optional` and `except HttpError` path
        # degrades exactly as it already does for an outage, and the
        # budget introduces no new failure mode. `refused` is what makes
        # that honest: a ceiling that silently truncates evidence is
        # worse than no ceiling, because the report would show failed
        # sources with nothing to say the caller caused it.
        if max_requests is not None and max_requests < 1:
            raise ValueError("max_requests must be at least 1")
        self.max_requests = max_requests
        self.requests_made = 0
        self.requests_refused = 0
        # Hosts that have already exhausted their retries this run. Patience
        # is the right response to a BUSY service and the wrong one to a
        # DEAD service: six jittered attempts span ~34s, which rides out an
        # Ensembl hiccup but is pure waiting against a host that has been
        # down for hours. SigCom was down for an entire session, costing
        # ~68s per run across its two endpoints, and ChEMBL is called
        # hundreds of times. So a host gets the full patient treatment once,
        # and fails fast afterwards.
        self._down: dict = {}

    # ---------------------------------------------------------------- cache
    def _path(self, url: str, body: str | None, *, method: str = "GET",
              parse_json: bool = True, headers: dict | None = None) -> Path:
        identity = json.dumps([url, method.upper(), body, parse_json,
                               sorted((k.lower(), v) for k, v in
                                      (headers or {}).items())], sort_keys=True)
        key = hashlib.sha256(identity.encode()).hexdigest()[:32]
        host = urllib.parse.urlparse(url).netloc.replace(":", "_") or "misc"
        d = self.cache_dir / host
        d.mkdir(parents=True, exist_ok=True)
        return d / f"{key}.json.gz"

    def _read_cache(self, path: Path):
        if not path.exists():
            return None
        if self.ttl and (time.time() - path.stat().st_mtime) > self.ttl and not self.offline:
            return None
        try:
            with gzip.open(path, "rt", encoding="utf-8") as fh:
                cached = json.load(fh)
            return cached if isinstance(cached, dict) and "payload" in cached else None
        except Exception:
            return None

    def _write_cache(self, path: Path, payload) -> None:
        # A distinct file per writer prevents concurrent runs from replacing
        # or deleting one another's unfinished response.
        tmp = None
        try:
            with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".tmp",
                                             delete=False) as fh:
                tmp = Path(fh.name)
            with gzip.open(tmp, "wt", encoding="utf-8") as fh:
                json.dump(payload, fh)
            tmp.replace(path)
        except OSError as exc:
            # Cache persistence is optional. In particular, Windows may deny
            # a replacement while another process has this entry open.
            log.warning("could not cache response at %s: %s", path, exc)
        finally:
            if tmp is not None:
                try:
                    tmp.unlink(missing_ok=True)
                except OSError:
                    log.debug("could not remove cache temporary file %s", tmp)

    # ------------------------------------------------------------- requests
    def _throttle(self) -> None:
        gap = time.monotonic() - self._last_call
        if gap < self.min_interval:
            time.sleep(self.min_interval - gap)
        self._last_call = time.monotonic()

    def request(self, url: str, *, method: str = "GET", data: dict | None = None,
                headers: dict | None = None, parse_json: bool = True,
                retries: int | None = None, timeout: float = 60,
                max_redirects: int = 5):
        retries = DEFAULT_RETRIES if retries is None else retries
        if retries < 1:
            raise ValueError("retries must be at least 1")
        if max_redirects < 0:
            raise ValueError("max_redirects must be non-negative")
        method = method.upper()
        body = json.dumps(data, sort_keys=True) if data is not None else None
        cache_path = self._path(url, body, method=method,
                                parse_json=parse_json, headers=headers)

        cached = self._read_cache(cache_path)
        if cached is not None:
            return cached["payload"]

        host = urllib.parse.urlsplit(url).netloc
        if host in self._down:
            # Already proven unreachable this run. The caller still gets an
            # HttpError, so every _optional / except HttpError path behaves
            # exactly as before — it just arrives immediately.
            raise HttpError(
                f"{host} already failed this run ({self._down[host]}); "
                f"not retrying {url}")

        if self.offline:
            raise HttpError(
                f"offline mode and no cached response for {url}. "
                "Run once with network access, or use --demo."
            )

        if self.max_requests is not None \
                and self.requests_made >= self.max_requests:
            self.requests_refused += 1
            raise HttpError(
                f"request budget of {self.max_requests} spent; not "
                f"fetching {url}. Raise --max-requests, or rerun - "
                f"everything already fetched is cached, so a second run "
                f"resumes rather than repeating.")
        self.requests_made += 1

        hdrs = {"User-Agent": USER_AGENT, "Accept": "application/json"}
        if body is not None:
            hdrs["Content-Type"] = "application/json"
        hdrs.update(headers or {})

        last_exc: Exception | None = None
        for attempt in range(retries):
            self._throttle()
            req = urllib.request.Request(
                url, method=method, headers=hdrs,
                data=body.encode() if body is not None else None,
            )
            try:
                with urllib.request.urlopen(req, timeout=timeout) as resp:
                    raw = resp.read()
                    if resp.headers.get("Content-Encoding") == "gzip":
                        raw = gzip.decompress(raw)
                text = raw.decode("utf-8", "replace")
                payload = json.loads(text) if parse_json else text
                self._write_cache(cache_path, {"url": url, "payload": payload})
                return payload
            except urllib.error.HTTPError as exc:
                last_exc = exc
                # urllib will not replay a POST body across 307/308, so it
                # surfaces the redirect as an error. Follow it manually.
                if exc.code in (307, 308) and method == "POST":
                    target = exc.headers.get("Location")
                    if target and target != url:
                        if max_redirects == 0:
                            raise HttpError(f"too many redirects for {url}") from exc
                        log.info("following %s redirect to %s", exc.code, target)
                        return self.request(
                            urllib.parse.urljoin(url, target), method=method,
                            data=data, headers=headers, parse_json=parse_json,
                            retries=retries, timeout=timeout,
                            max_redirects=max_redirects - 1)
                if exc.code == 404:
                    raise HttpError(f"404 for {url}") from exc
                if exc.code in (429, 500, 502, 503, 504):
                    # Public bioinformatics endpoints 500 under load as
                    # often as they 429, so a 500 is a backoff signal rather
                    # than a verdict. _retry_delay honours Retry-After in
                    # either the seconds or the HTTP-date form; the cap and
                    # jitter stop parallel callers resynchronising onto the
                    # same retry instant, and there is no point sleeping
                    # after the final attempt.
                    if attempt + 1 == retries:
                        break
                    wait = _retry_delay(exc.headers.get("Retry-After"), attempt)
                    if not exc.headers.get("Retry-After"):
                        wait = min(RETRY_CAP, wait) + random.uniform(0, 1.0)
                    log.warning("HTTP %s from %s, retrying in %.1fs "
                                "(attempt %d/%d)", exc.code, url, wait,
                                attempt + 1, retries)
                    time.sleep(wait)
                    continue
                raise HttpError(f"HTTP {exc.code} for {url}") from exc
            except (urllib.error.URLError, http.client.HTTPException,
                    socket.error, OSError, TimeoutError,
                    json.JSONDecodeError) as exc:
                # RemoteDisconnected subclasses ConnectionResetError -> OSError,
                # not URLError, so a narrower clause lets it escape the retry
                # loop and abort the run. Catch the whole family.
                last_exc = exc
                if attempt + 1 < retries:
                    time.sleep(2 ** attempt)

        reason = f"{type(last_exc).__name__}: {str(last_exc)[:80]}"
        if host not in self._down:
            self._down[host] = reason
            log.warning(
                "%s is unreachable after %d attempts (%s). Treating it as "
                "down for the rest of this run: further requests to it fail "
                "immediately instead of waiting ~%.0fs each. Everything not "
                "depending on it is unaffected.",
                host, retries, reason, sum(min(RETRY_CAP, 2.0 ** a)
                                           for a in range(retries - 1)))
        raise HttpError(f"request failed after {retries} attempts: {url} ({last_exc})")

    def spend(self) -> dict:
        """What this run cost, and whether the ceiling bound it.

        `requests` counts round trips that left the machine: cache hits
        are not among them, a retried request counts once, and each
        redirect hop counts separately.

        `refused` above zero means evidence is missing because of the
        budget rather than because a source had nothing - a distinction
        the report has to draw, since the two look identical in every
        table underneath it.
        """
        return {"requests": self.requests_made,
                "refused": self.requests_refused,
                "budget": self.max_requests}

    def get(self, url: str, params: dict | None = None, **kw):
        if params:
            parts = urllib.parse.urlsplit(url)
            query = urllib.parse.urlencode(params, doseq=True)
            url = urllib.parse.urlunsplit(parts._replace(
                query=f"{parts.query}&{query}" if parts.query else query))
        return self.request(url, **kw)

    def post(self, url: str, data: dict, **kw):
        return self.request(url, method="POST", data=data, **kw)


def _retry_delay(value: str | None, attempt: int) -> float:
    """Accept Retry-After seconds or an HTTP date; tolerate malformed values."""
    import math
    from email.utils import parsedate_to_datetime

    if value:
        try:
            delay = float(value)
        except ValueError:
            try:
                delay = parsedate_to_datetime(value).timestamp() - time.time()
            except (ValueError, TypeError, OverflowError):
                delay = float("nan")
        if math.isfinite(delay):
            return max(0.0, delay)
    return float(2 ** attempt)
