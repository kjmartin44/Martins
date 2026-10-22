"""Telling the user what the run is doing while it does it.

A real gene takes somewhere between two minutes and a quarter of an
hour, and until now it printed nothing at all until it finished. A
working run and a wedged one look identical for that long, which is the
question that actually gets asked ten minutes in, and the answer people
reach for - run it again with `-v` - replaces silence with thousands of
lines of per-factor lookups and is no easier to read.

So this is a separate, curated stream: one line per stage, with the
elapsed time, and a short detail saying what the stage found.

Three decisions worth knowing.

**It writes to stderr.** The report goes to stdout, and a run whose
output is piped into a file must not have progress lines spliced through
the middle of the markdown.

**It is off by default.** `pipeline.run` is a library function and a
library that prints is a nuisance. The commands turn it on, so a person
at a terminal gets it and a program embedding the package does not.

**It never raises.** A closed pipe or an encoding a console cannot
render must not take a fifteen-minute run down at the last step, so
everything here is wrapped and failure means silence rather than a
traceback.
"""

from __future__ import annotations

import sys
import time

__all__ = ["enable", "disable", "enabled", "start", "step", "finish",
           "stage_count"]

# The fixed stages of `pipeline.run`, in order, so a line can say how far
# along it is. Stages that a run's arguments skip still count: "6/10,
# skipped" is more informative than a denominator that changes per run.
STAGES = 10

_state: dict = {"on": False, "stream": None, "t0": None, "prefix": ""}


def enable(stream=None, *, prefix: str = "") -> None:
    """Start reporting. `prefix` labels the gene in a multi-gene run."""
    _state["on"] = True
    _state["stream"] = stream if stream is not None else sys.stderr
    _state["prefix"] = prefix


def disable() -> None:
    _state["on"] = False


def enabled() -> bool:
    return bool(_state["on"])


def stage_count() -> int:
    return STAGES


def _write(text: str) -> None:
    if not _state["on"]:
        return
    stream = _state["stream"] or sys.stderr
    try:
        stream.write(text + "\n")
        stream.flush()
    except Exception:                                       # noqa: BLE001
        # A closed pipe, or a console that cannot encode the text. Not a
        # reason to lose the run.
        pass


def _elapsed() -> str:
    t0 = _state["t0"]
    if t0 is None:
        return "  0:00"
    s = int(time.time() - t0)
    return f"{s // 60:3d}:{s % 60:02d}"


def start(what: str) -> None:
    """Begin timing, and say what is being worked on."""
    _state["t0"] = time.time()
    _write(f"[{_elapsed()}] {_state['prefix']}{what}")


def step(n: int, label: str, detail: str = "") -> None:
    """One stage finished. `detail` is what it found, not what it did."""
    line = f"[{_elapsed()}] {_state['prefix']}{n:2d}/{STAGES}  {label}"
    if detail:
        line += f" - {detail}"
    _write(line)


def finish(detail: str = "") -> None:
    _write(f"[{_elapsed()}] {_state['prefix']}done"
           + (f" - {detail}" if detail else ""))


def note(text: str) -> None:
    """Something worth saying between stages, indented under them."""
    _write(f"[{_elapsed()}] {_state['prefix']}        {text}")
