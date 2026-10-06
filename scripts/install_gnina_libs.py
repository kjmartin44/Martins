"""Fetch the CUDA runtime libraries gnina needs, without root.

gnina's "static" release is static in everything except CUDA: it still
wants seven shared libraries that only ship with the CUDA toolkit. On a
WSL install with no passwordless sudo and no NVIDIA apt repository,
installing that toolkit is not available.

NVIDIA publishes the same libraries as PyPI wheels, and a wheel is a zip.
This downloads the ones gnina asks for and unpacks their .so files into a
directory of its own, which gnina is then pointed at with
LD_LIBRARY_PATH. Nothing is installed system-wide and removing the
directory undoes it.

Run inside WSL:
    python3 scripts/install_gnina_libs.py
"""

import json
import os
import sys
import urllib.request
import zipfile
from pathlib import Path

# Each wheel, and which of gnina's missing libraries it carries.
WHEELS = {
    "nvidia-cuda-runtime-cu12": ["libcudart.so.12"],
    "nvidia-cublas-cu12": ["libcublas.so.12", "libcublasLt.so.12"],
    "nvidia-cudnn-cu12": ["libcudnn.so.9"],
    "nvidia-cufft-cu12": ["libcufft.so.11"],
    "nvidia-cusolver-cu12": ["libcusolver.so.11"],
    "nvidia-cusparse-cu12": ["libcusparse.so.12"],
    # Transitive: cusolver and cusparse link against it, so it does not
    # appear in gnina's own `ldd` until the two above are resolved.
    "nvidia-nvjitlink-cu12": ["libnvJitLink.so.12"],
    # dlopen'd at runtime by the PyTorch gnina embeds, so it appears in
    # no `ldd` output at all - only in a crash once docking starts.
    "nvidia-cuda-nvrtc-cu12": ["libnvrtc.so.12"],
}

DEST = Path(os.path.expanduser("~/cuda-libs"))
UA = {"User-Agent": "promoter-tf-drugs"}


def newest_linux_wheel(pkg: str) -> tuple:
    """(url, size) of the newest manylinux x86_64 wheel for `pkg`."""
    req = urllib.request.Request(
        f"https://pypi.org/pypi/{pkg}/json", headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        data = json.loads(r.read())

    def key(v):
        return [int(x) if x.isdigit() else 0 for x in v.split(".")]

    for version in sorted(data["releases"], key=key, reverse=True):
        for f in data["releases"][version]:
            name = f["filename"]
            if (name.endswith(".whl") and "manylinux" in name
                    and "x86_64" in name and "aarch" not in name):
                return f["url"], f["size"], version
    raise SystemExit(f"no linux wheel found for {pkg}")


def main() -> int:
    DEST.mkdir(parents=True, exist_ok=True)
    total = 0
    for pkg, libs in WHEELS.items():
        if all((DEST / lib).exists() for lib in libs):
            print(f"  {pkg:28s} already present")
            continue
        url, size, version = newest_linux_wheel(pkg)
        print(f"  {pkg:28s} {version:12s} {size / 1e6:6.0f} MB", flush=True)
        tmp = DEST / "wheel.tmp"
        with urllib.request.urlopen(
                urllib.request.Request(url, headers=UA), timeout=600) as r, \
                open(tmp, "wb") as fh:
            while True:
                chunk = r.read(1 << 22)
                if not chunk:
                    break
                fh.write(chunk)
        with zipfile.ZipFile(tmp) as z:
            for info in z.infolist():
                base = Path(info.filename).name
                # Keep every .so, not only the exact names asked for: a
                # versioned library usually needs its own soname links,
                # and the extras cost disk rather than correctness.
                if ".so" not in base:
                    continue
                out = DEST / base
                with z.open(info) as src, open(out, "wb") as dst:
                    while True:
                        chunk = src.read(1 << 22)
                        if not chunk:
                            break
                        dst.write(chunk)
                total += 1
        tmp.unlink()

    print(f"\n{total} library file(s) in {DEST}")
    missing = [lib for libs in WHEELS.values() for lib in libs
               if not (DEST / lib).exists()]
    if missing:
        # A wheel can ship `libcudnn.so.9.1.0` and no `.so.9` symlink,
        # which is the name the binary actually asks for.
        for lib in list(missing):
            stem = lib.split(".so")[0]
            cands = sorted(p for p in DEST.glob(stem + ".so*")
                           if p.name != lib)
            if cands:
                target = DEST / lib
                try:
                    target.symlink_to(cands[-1].name)
                    print(f"  linked {lib} -> {cands[-1].name}")
                    missing.remove(lib)
                except OSError as exc:
                    print(f"  could not link {lib}: {exc}")
    if missing:
        print("STILL MISSING:", ", ".join(missing))
        return 1
    print("every library gnina asked for is present")
    print(f"\nRun gnina with:  LD_LIBRARY_PATH={DEST} ~/bin/gnina")
    return 0


if __name__ == "__main__":
    sys.exit(main())
