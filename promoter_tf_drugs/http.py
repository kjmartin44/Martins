"""Tiny cached HTTP client. Stdlib only, so the package needs only numpy.

Every external resource used here is a public API with rate limits. Caching to
disk means a re-run of the same gene costs nothing and you stop hammering
Ensembl/ChEMBL while you iterate on thresholds.
"""

from __future__ import annotations

import gzip
import hashlib
import http.client
import json
import logging
import os
import socket
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

log = logging.getLogger(__name__)

DEFAULT_CACHE = Path(os.environ.get("PTD_CACHE", Path.home() / ".cache" / "promoter_tf_drugs"))
USER_AGENT = "promoter-tf-drugs/0.1 (research pipeline)"


class HttpError(RuntimeError):
    pass


class Client:
    def __init__(self, cache_dir: Path | str = DEFAULT_CACHE, ttl_days: float = 30,
                 offline: bool = False, min_interval: float = 0.15):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.ttl = ttl_days * 86400
        self.offline = offline
        self.min_interval = min_interval
        self._last_call = 0.0

    # ---------------------------------------------------------------- cache
    def _path(self, url: str, body: str | None) -> Path:
        key = hashlib.sha256((url + "\0" + (body or "")).encode()).hexdigest()[:32]
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
                return json.load(fh)
        except Exception:
            return None

    def _write_cache(self, path: Path, payload) -> None:
        tmp = path.with_suffix(".tmp")
        with gzip.open(tmp, "wt", encoding="utf-8") as fh:
            json.dump(payload, fh)
        tmp.replace(path)

    # ------------------------------------------------------------- requests
    def _throttle(self) -> None:
        gap = time.monotonic() - self._last_call
        if gap < self.min_interval:
            time.sleep(self.min_interval - gap)
        self._last_call = time.monotonic()

    def request(self, url: str, *, method: str = "GET", data: dict | None = None,
                headers: dict | None = None, parse_json: bool = True,
                retries: int = 3, timeout: float = 60):
        body = json.dumps(data, sort_keys=True) if data is not None else None
        cache_path = self._path(url, body)

        cached = self._read_cache(cache_path)
        if cached is not None:
            return cached["payload"]

        if self.offline:
            raise HttpError(
                f"offline mode and no cached response for {url}. "
                "Run once with network access, or use --demo."
            )

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
                        log.info("following %s redirect to %s", exc.code, target)
                        return self.request(
                            urllib.parse.urljoin(url, target), method=method,
                            data=data, headers=headers, parse_json=parse_json,
                            retries=retries, timeout=timeout)
                if exc.code == 404:
                    raise HttpError(f"404 for {url}") from exc
                if exc.code in (429, 500, 502, 503, 504):
                    wait = float(exc.headers.get("Retry-After", 2 ** attempt))
                    log.warning("HTTP %s from %s, retrying in %.1fs", exc.code, url, wait)
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
                time.sleep(2 ** attempt)

        raise HttpError(f"request failed after {retries} attempts: {url} ({last_exc})")

    def get(self, url: str, params: dict | None = None, **kw):
        if params:
            url = f"{url}?{urllib.parse.urlencode(params)}"
        return self.request(url, **kw)

    def post(self, url: str, data: dict, **kw):
        return self.request(url, method="POST", data=data, **kw)
