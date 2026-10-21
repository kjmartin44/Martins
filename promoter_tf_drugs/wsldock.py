"""Run a Linux docking binary from Windows, through WSL.

Why this exists
---------------
The docking programs worth using are Linux software. smina has a Windows
build and AutoDock Vina has one; gnina, which is the reason to bother,
does not. So on Windows the choice was a weaker program or nothing, and
this removes it: WSL runs the Linux binary while the package keeps
passing Windows paths around as it always has.

What makes gnina the better program
-----------------------------------
Vina and smina score a pose with an empirical function - a weighted sum
of distance terms fitted to binding data. gnina docks the same way and
then RESCORES the poses with a convolutional network trained on
crystallographic complexes, which is a different kind of judgement about
whether a pose looks like a real one.

It gives two numbers per pose, answering different questions:

  CNNscore     how pose-like the geometry is, 0 to 1 - the network's
               confidence that it is looking at a real binding mode.
  CNNaffinity  predicted affinity in pK units, larger being tighter, and
               comparable across ligands in a way a Vina score is not.

The Vina-style `minimizedAffinity` is still produced and still reported,
because the three disagree often enough to be worth seeing.

Two things that bite, both found the hard way
---------------------------------------------
**Scripts go over stdin, never as `bash -c`.** The Windows command line
is rebuilt on the way through the interop layer and a script passed as an
argument does not survive it: `wsl.exe -- bash -c "echo $PATH"` comes
back with a bash syntax error. There is no command line to mangle when
the script arrives on stdin.

**Output carries interleaved NUL bytes.** Every comparison against it
fails in a way that reads like the command having produced nothing.

What it does not change
-----------------------
**A pocket still has to validate.** A better scoring function ranks a
validated pocket better; it does not make an unvalidated one meaningful,
and `screen_by_docking` still refuses below its AUC and enrichment
floors.

**Direction is still unaddressed.** A network confident that a molecule
binds says nothing about whether binding activates or inhibits, which is
what decides the direction your gene moves.

**A GPU changes the speed, not the answer.** With CUDA the rescoring
costs little; without it gnina falls back to CPU and the same numbers
arrive slower.
"""

from __future__ import annotations

import logging
import os
import shlex
import shutil
import subprocess
from pathlib import Path, PureWindowsPath

log = logging.getLogger(__name__)

NUL = chr(0)

# Preference order. gnina first for the CNN rescoring, then smina for its
# scoring flexibility, then vina.
LINUX_BINARIES = ("gnina", "smina", "vina")

# Where a hand-installed binary lands. Searched because `bash` started
# this way reads no login profile, so ~/bin is on no PATH at all and the
# obvious `command -v` check reports nothing installed.
EXTRA_DIRS = ("$HOME/bin", "$HOME/.local/bin", "/usr/local/bin", "/opt/bin")

# Where the CUDA libraries gnina needs may have been unpacked.
LIB_DIRS = ("$HOME/cuda-libs", "$HOME/.local/cuda-libs",
            "/usr/local/cuda/lib64")

_CACHE: dict = {}


def run_script(script: str, *, timeout: float = 1800
               ) -> subprocess.CompletedProcess:
    """Run a shell script inside WSL, passed over STDIN.

    Not `bash -c`: a script passed as an argument is mangled by the
    interop layer's command-line rebuilding, and
    `wsl.exe -- bash -c "echo $PATH"` comes back with a bash syntax
    error. There is no command line to mangle over stdin.

    The input is sent as BYTES. With text=True Python translates
    every newline to CRLF on Windows and bash then fails with
    "syntax error near unexpected token $'do\r'" on the first
    multi-line construct - which reads like a bug in the script
    rather than in how it was delivered.

    NUL bytes are stripped from the output. The interop layer
    interleaves them, and they make every comparison against the
    output fail in a way that reads like the command having produced
    nothing.
    """
    clean = script.replace("\r\n", "\n")
    proc = subprocess.run(["wsl.exe", "--", "bash", "-s"],
                          input=clean.encode("utf-8"),
                          capture_output=True, timeout=timeout)

    def _text(raw) -> str:
        return (raw or b"").decode("utf-8", "replace").replace(NUL, "")

    return subprocess.CompletedProcess(
        proc.args, proc.returncode, _text(proc.stdout),
        _text(proc.stderr))


def available() -> bool:
    """Is WSL usable at all from here?"""
    if "available" in _CACHE:
        return _CACHE["available"]
    ok = False
    if shutil.which("wsl.exe"):
        try:
            ok = run_script("true", timeout=60).returncode == 0
        except Exception:                                   # noqa: BLE001
            ok = False
    _CACHE["available"] = ok
    return ok


def _from_windows_path(path) -> str | None:
    """`C:\\Users\\x\\y` -> `/mnt/c/Users/x/y`, parsed as Windows syntax.

    `PureWindowsPath` reads a Windows path on any host, which `Path`
    does not: on POSIX a backslash is an ordinary filename character,
    so `Path("C:\\Users\\x")` is a single relative name. Returns None
    for anything that is not drive-absolute - a POSIX path, a relative
    path, `C:file` (relative to a drive's current directory) or a UNC
    path - so the caller can handle those as it always did.
    """
    win = PureWindowsPath(str(path))
    if win.is_absolute() and len(win.drive) == 2 and win.drive[1] == ":":
        return "/".join(["/mnt", win.drive[0].lower(), *win.parts[1:]])
    return None


def to_wsl_path(path) -> str:
    """`C:\\Users\\x\\y` -> `/mnt/c/Users/x/y`.

    Done here rather than by calling `wslpath`, which is one
    subprocess per path and a screen converts thousands.

    A path that is already POSIX-absolute comes back untouched. It
    has to: `Path("/home/u/x").resolve()` on Windows resolves the
    leading slash against the current drive and gives C:\\home\\u\\x,
    which this would then hand to WSL as /mnt/c/home/u/x - a
    location that does not exist, reported as a missing file rather
    than as a path that was converted when it should not have been.
    """
    raw = str(path).replace("\\", "/")
    if raw.startswith("/"):
        return raw
    if os.name != "nt":
        # `Path` cannot read a Windows path here, and `resolve()` would
        # prefix the current directory to it. On Windows itself the
        # code below stays as it was: `resolve()` collapses `..` and
        # follows junctions, which parsing the text does not.
        converted = _from_windows_path(path)
        if converted is not None:
            return converted
    s = str(Path(path).resolve()).replace("\\", "/")
    if len(s) > 1 and s[1] == ":":
        return "/mnt/" + s[0].lower() + s[2:]
    return s


def _first_existing(candidates: str, test: str) -> str | None:
    """First of `candidates` satisfying `test`, as an absolute path."""
    script = (f"for d in {candidates}; do\n"
              f"  if {test}; then echo \"$d\"; exit 0; fi\n"
              f"done\nexit 1\n")
    try:
        r = run_script(script, timeout=120)
    except Exception as exc:                                # noqa: BLE001
        log.warning("could not run a search inside WSL (%s: %s), so "
                    "nothing was found there", type(exc).__name__,
                    str(exc)[:100])
        return None
    if r.returncode == 0 and r.stdout.strip():
        return r.stdout.strip().splitlines()[0].strip()
    return None


def find_binary(prefer: str | None = None) -> tuple:
    """`(name, path)` of the best Linux docking binary, or `(None, None)`."""
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
        script = (f"command -v {name} 2>/dev/null && exit 0\n"
                  f"for d in {dirs}; do\n"
                  f'  if [ -x "$d/{name}" ]; then echo "$d/{name}"; '
                  f"exit 0; fi\n"
                  f"done\nexit 1\n")
        try:
            r = run_script(script, timeout=120)
        except Exception as exc:                            # noqa: BLE001
            log.warning("could not look for %s inside WSL (%s: %s)", name,
                        type(exc).__name__, str(exc)[:100])
            continue
        if r.returncode == 0 and r.stdout.strip():
            found = r.stdout.strip().splitlines()[0].strip()
            if found:
                _CACHE[key] = (name, found)
                log.info("docking binary in WSL: %s at %s", name, found)
                return _CACHE[key]
    _CACHE[key] = (None, None)
    return _CACHE[key]


def library_path() -> str | None:
    """Directory holding the CUDA libraries gnina needs, if there is one.

    gnina's "static" release is static in everything except CUDA: it
    still links against seven shared libraries that ship with the
    toolkit. Where no toolkit is installed,
    `scripts/install_gnina_libs.py` unpacks them from NVIDIA's PyPI
    wheels, and this finds that directory again.
    """
    if "libdir" in _CACHE:
        return _CACHE["libdir"]
    found = None
    if available():
        found = _first_existing(" ".join(LIB_DIRS),
                                '[ -f "$d/libcudnn.so.9" ]')
    _CACHE["libdir"] = found
    return found


def has_gpu() -> bool:
    """Does WSL see a CUDA device? Decides speed, never the answer."""
    if "gpu" in _CACHE:
        return _CACHE["gpu"]
    ok = False
    if available():
        try:
            r = run_script(
                "nvidia-smi --query-gpu=name --format=csv,noheader",
                timeout=120)
            ok = r.returncode == 0 and bool(r.stdout.strip())
        except Exception:                                   # noqa: BLE001
            ok = False
    _CACHE["gpu"] = ok
    return ok


def gpu_name() -> str | None:
    if not has_gpu():
        return None
    try:
        r = run_script(
            "nvidia-smi --query-gpu=name,memory.total --format=csv,noheader",
            timeout=120)
        return r.stdout.strip().splitlines()[0].strip() or None
    except Exception as exc:                                # noqa: BLE001
        log.debug("could not read the GPU name: %s", type(exc).__name__)
        return None


def dock_script(binary: str, args: list) -> str:
    """A shell script running `binary` with the CUDA libraries in scope.

    Built here so no caller has to know about LD_LIBRARY_PATH, and so a
    machine with a real CUDA toolkit simply gets no prefix.

    Every word is quoted for the shell. Joining them raw split any path
    holding a space - `/mnt/c/Users/Ada Lovelace/rec.pdb` arrived as two
    arguments - and the binary then exited non-zero, which the caller
    records as a compound that did not score. A whole screen came back
    empty rather than broken, for a username.
    """
    libdir = library_path()
    prefix = (f'export LD_LIBRARY_PATH={shlex.quote(f"{libdir}:")}'
              f'"$LD_LIBRARY_PATH"\n' if libdir else "")
    return prefix + " ".join(
        shlex.quote(str(a)) for a in [binary] + list(args)) + "\n"


def reset_cache() -> None:
    """Forget what was discovered. For tests, and after an install."""
    _CACHE.clear()


def describe() -> str:
    """One paragraph on what is set up, for the report and the CLI."""
    if not available():
        return ("WSL is not available from here, so the Linux docking "
                "binaries cannot be reached. A Windows smina or vina on "
                "PATH is still used if one is installed.")
    name, path = find_binary()
    if not name:
        return ("WSL is available and no docking binary is installed in "
                "it. gnina is the one worth having: it rescores poses "
                "with a convolutional network rather than an empirical "
                "function. See `python -m promoter_tf_drugs.wsldock`.")
    bits = [f"{name} at {path}, through WSL."]
    if name == "gnina":
        bits.append(
            "Poses are rescored by a convolutional network, giving a "
            "CNNscore - how pose-like the geometry is - and a CNNaffinity "
            "in pK units alongside the Vina-style score. The three "
            "disagree often enough to be worth seeing.")
        if not library_path():
            bits.append(
                "No CUDA library directory was found, and gnina's release "
                "needs one: run scripts/install_gnina_libs.py inside WSL "
                "if it fails to start.")
    else:
        bits.append(
            "This is an empirical scoring function. gnina would add CNN "
            "rescoring, the one upgrade that changes pose ranking rather "
            "than speed.")
    gpu = gpu_name()
    bits.append(f"GPU: {gpu}." if gpu else
                "No CUDA device visible to WSL, so rescoring runs on CPU - "
                "slower, same numbers.")
    return " ".join(bits)


def main(argv=None) -> int:
    """`python -m promoter_tf_drugs.wsldock` - report what is set up."""
    import argparse

    p = argparse.ArgumentParser(
        prog="python -m promoter_tf_drugs.wsldock",
        description="What docking program is reachable through WSL.")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(name)s: %(message)s")

    print(describe())
    name, _ = find_binary()
    if name:
        print()
        print("Docking will use it automatically. The pocket still has to "
              "validate: a better scoring function ranks a validated "
              "pocket better and does not make an unvalidated one mean "
              "anything.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
