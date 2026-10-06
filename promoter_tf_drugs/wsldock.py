"""Run a Linux docking binary from Windows, through WSL.

Why this exists
---------------
The docking programs worth using are Linux software. smina has a Windows
build and AutoDock Vina has one; gnina, which is the reason to bother,
does not. So on Windows the choice was a weaker program or nothing, and
this module removes it: WSL runs the Linux binary and the package keeps
passing Windows paths around as it always has.

What makes gnina the better program
-----------------------------------
Vina and smina score a pose with an empirical function - a weighted sum
of distance terms fitted to binding data. gnina docks the same way and
then RESCORES the poses with a convolutional network trained on
crystallographic complexes, which is a different kind of judgement about
whether a pose looks like a real one. On the standard benchmarks that
moves top-1 pose accuracy up substantially over Vina.

It gives two numbers per pose and they answer different questions:

  CNNscore     how pose-like this geometry is, 0 to 1. The network's
               confidence that it is looking at a real binding mode.
  CNNaffinity  predicted binding affinity in pK units, so larger is
               tighter. Comparable across ligands in a way the Vina
               score is not.

The Vina-style `minimizedAffinity` is still produced and still reported,
because the three disagree often enough to be worth seeing.

What it does not change
-----------------------
**A pocket still has to validate.** A better scoring function ranks a
validated pocket better; it does not make an unvalidated one meaningful,
and `screen_by_docking` still refuses below its AUC and enrichment
floors.

**Direction is still unaddressed.** A CNN that is confident a molecule
binds says nothing about whether binding activates or inhibits, which is
what decides the direction your gene moves.

**A GPU changes the speed, not the answer.** With CUDA the CNN rescoring
costs little; without it gnina falls back to CPU and the same numbers
arrive slower.
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
from pathlib import Path

log = logging.getLogger(__name__)

# Preference order. gnina first because of the CNN rescoring, then smina
# for its scoring flexibility, then vina.
LINUX_BINARIES = ("gnina", "smina", "vina")

# Where a binary installed by hand usually lands, searched when it is not
# on the non-interactive PATH - `wsl bash -c` does not read a login
# profile, so ~/bin is routinely invisible to it.
EXTRA_DIRS = ("$HOME/bin", "$HOME/.local/bin", "/usr/local/bin", "/opt/bin")

_CACHE: dict = {}


def available() -> bool:
    """Is WSL usable at all from here?"""
    if not shutil.which("wsl.exe"):
        return False
    try:
        p = subprocess.run(["wsl.exe", "--", "true"],
                           capture_output=True, timeout=60)
        return p.returncode == 0
    except Exception:                                       # noqa: BLE001
        return False


def to_wsl_path(path) -> str:
    """`C:\\Users\\x\\y` -> `/mnt/c/Users/x/y`.

    Done here rather than with `wslpath` because that is one subprocess
    per path and a screen converts thousands.
    """
    p = Path(path).resolve()
    s = str(p).replace("\\\\", "/").replace("\\", "/")
    if len(s) > 1 and s[1] == ":":
        return "/mnt/" + s[0].lower() + s[2:]
    return s


def find_binary(prefer: str | None = None) -> tuple:
    """`(name, path)` of the best Linux docking binary, or `(None, None)`.

    Looks on the PATH a non-interactive shell sees, then in the places a
    hand-installed binary lands. The second step matters: `wsl bash -c`
    skips the login profile, so a binary in ~/bin is on no PATH at all
    and the obvious check reports nothing installed.
    """
    key = ("find", prefer)
    if key in _CACHE:
        return _CACHE[key]
    if not available():
        _CACHE[key] = (None, None)
        return _CACHE[key]

    order = ([prefer] if prefer else []) + [
        b for b in LINUX_BINARIES if b != prefer]
    dirs = " ".join(EXTRA_DIRS)
    for name in order:
        script = (
            f'command -v {name} 2>/dev/null && exit 0; '
            f'for d in {dirs}; do '
            f'  if [ -x "$d/{name}" ]; then echo "$d/{name}"; exit 0; fi; '
            f'done; exit 1')
        try:
            p = subprocess.run(["wsl.exe", "--", "bash", "-c", script],
                               capture_output=True, text=True, timeout=120)
        except Exception:                                   # noqa: BLE001
            continue
        out = (p.stdout or "").replace("\\x00", "").strip().splitlines()
        if p.returncode == 0 and out:
            found = out[0].strip()
            if found:
                _CACHE[key] = (name, found)
                log.info("docking binary in WSL: %s at %s", name, found)
                return _CACHE[key]
    _CACHE[key] = (None, None)
    return _CACHE[key]


def has_gpu() -> bool:
    """Does WSL see a CUDA device? Decides speed, never the answer."""
    if "gpu" in _CACHE:
        return _CACHE["gpu"]
    ok = False
    if available():
        try:
            p = subprocess.run(
                ["wsl.exe", "--", "bash", "-c",
                 "nvidia-smi --query-gpu=name --format=csv,noheader"],
                capture_output=True, text=True, timeout=120)
            ok = p.returncode == 0 and bool((p.stdout or "").strip())
        except Exception:                                   # noqa: BLE001
            ok = False
    _CACHE["gpu"] = ok
    return ok


def gpu_name() -> str | None:
    if not has_gpu():
        return None
    p = subprocess.run(
        ["wsl.exe", "--", "bash", "-c",
         "nvidia-smi --query-gpu=name,memory.total --format=csv,noheader"],
        capture_output=True, text=True, timeout=120)
    return (p.stdout or "").replace("\x00", "").strip() or None


def run(args: list, *, timeout: float = 1800) -> subprocess.CompletedProcess:
    """Run a command inside WSL, with Windows paths already translated."""
    return subprocess.run(["wsl.exe", "--"] + args,
                          capture_output=True, text=True, timeout=timeout)


def describe() -> str:
    """One paragraph on what is set up, for the report and the CLI."""
    if not available():
        return ("WSL is not available from here, so the Linux docking "
                "binaries cannot be reached. A Windows smina or vina on "
                "PATH is still used if one is installed.")
    name, path = find_binary()
    if not name:
        return ("WSL is available and no docking binary is installed in "
                "it. gnina is the one worth having - it rescores poses "
                "with a convolutional network rather than an empirical "
                "function. Install it with `python -m "
                "promoter_tf_drugs.wsldock --install`.")
    gpu = gpu_name()
    bits = [f"{name} at {path}, through WSL."]
    if name == "gnina":
        bits.append(
            "Poses are rescored by a convolutional network, which gives a "
            "CNNscore (how pose-like the geometry is) and a CNNaffinity "
            "(predicted pK, larger is tighter) alongside the Vina-style "
            "score. The three disagree often enough to be worth seeing.")
    else:
        bits.append(
            "This is an empirical scoring function. gnina would add CNN "
            "rescoring, which is the one upgrade that changes pose "
            "ranking rather than speed.")
    bits.append(f"GPU: {gpu}." if gpu else
                "No CUDA device visible to WSL, so gnina would rescore on "
                "CPU - slower, same numbers.")
    return " ".join(bits)
