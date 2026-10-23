"""Declarative argument checking for the entry points.

Why this exists
---------------
`pipeline.run` takes fifty-eight arguments. Fifty-seven carry a type
annotation and about fifteen were actually enforced, by a block of
hand-written `if ... raise ValueError` lines. An annotation Python never
reads is a comment, so `upstream="1000"` reached the window arithmetic
and failed much later, somewhere that looked unrelated.

This makes the contract a table. Every argument gets a rule, a test
asserts the table and the signature still agree, and a new argument
cannot be added without saying how it is checked.

Why not pydantic
----------------
It was considered, and it is the obvious tool for this shape of problem.
Two reasons against it here.

The dependency. This package keeps NumPy as its only required runtime
dependency, and pydantic's core is a compiled extension. That is a real
cost for a package whose point is that it runs wherever Python does, on
machines where installing anything is the hard part.

More to the point, it would not have caught the bugs this code actually
produces. Those have all been well-typed values carrying the wrong
meaning: a POSIX path where a Windows one was needed, a modified amino
acid selected as a ligand, a variance parsed as an affinity, a failed
lookup recorded as a measured zero. Every one passes a type check. What
catches those is a check that knows what the value MEANS - which is
`why` below, and the static guards in the test suite - not one that
knows it is a float.

So: the declarative half is worth having and is here. The dependency is
not, because the half it would add is the half that was already fine.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

__all__ = ["Spec", "check", "check_all", "coverage"]


@dataclass(frozen=True)
class Spec:
    """How one argument is checked, and what a bad value would do.

    `why` is not decoration. A message saying "pvalue must be between 0
    and 1" tells you the rule; one saying it "is a per-base motif
    threshold, and 1.0 would call a site at every position" tells you
    what you got wrong, which is usually the part that was not obvious.
    """

    kind: str                                  # see _KINDS
    optional: bool = False                     # None is allowed
    minimum: float | None = None
    maximum: float | None = None
    exclusive_min: bool = False
    exclusive_max: bool = False
    choices: tuple = ()
    must_exist: bool = False                   # for kind="path"
    length: int | None = None                  # exact, for list/tuple
    why: str = ""


_KINDS = ("int", "float", "str", "bool", "path", "list", "tuple", "any")


def _type_ok(value, kind: str) -> bool:
    if kind == "any":
        return True
    if kind == "bool":
        return isinstance(value, bool)
    if kind == "int":
        # bool is a subclass of int, and letting True through as 1 is how
        # a flag silently becomes a window size.
        return isinstance(value, int) and not isinstance(value, bool)
    if kind == "float":
        return (isinstance(value, (int, float))
                and not isinstance(value, bool))
    if kind == "str":
        return isinstance(value, str)
    if kind == "path":
        return isinstance(value, (str, Path))
    if kind == "list":
        return isinstance(value, (list, tuple))
    if kind == "tuple":
        return isinstance(value, (list, tuple))
    raise ValueError(f"unknown kind {kind!r}")


def check(name: str, value, spec: Spec) -> None:
    """Raise if `value` does not satisfy `spec`. Returns None otherwise."""
    tail = f" {spec.why}" if spec.why else ""

    if value is None:
        if spec.optional:
            return
        raise ValueError(f"{name} must not be None.{tail}")

    if not _type_ok(value, spec.kind):
        raise ValueError(
            f"{name} must be {spec.kind}, got "
            f"{type(value).__name__}.{tail}")

    if spec.kind == "str" and not spec.choices and not value.strip():
        raise ValueError(f"{name} must not be empty.{tail}")

    if spec.choices and value not in spec.choices:
        shown = ", ".join(repr(c) for c in spec.choices)
        raise ValueError(f"{name} must be one of {shown}, got "
                         f"{value!r}.{tail}")

    if spec.kind in ("int", "float"):
        if not math.isfinite(value):
            raise ValueError(f"{name} must be finite, got {value!r}.{tail}")
        if spec.minimum is not None:
            bad = (value <= spec.minimum if spec.exclusive_min
                   else value < spec.minimum)
            if bad:
                rel = "greater than" if spec.exclusive_min else "at least"
                raise ValueError(f"{name} must be {rel} {spec.minimum}, "
                                 f"got {value!r}.{tail}")
        if spec.maximum is not None:
            bad = (value >= spec.maximum if spec.exclusive_max
                   else value > spec.maximum)
            if bad:
                rel = "less than" if spec.exclusive_max else "at most"
                raise ValueError(f"{name} must be {rel} {spec.maximum}, "
                                 f"got {value!r}.{tail}")

    if spec.length is not None and spec.kind in ("list", "tuple"):
        if len(value) != spec.length:
            raise ValueError(
                f"{name} must have exactly {spec.length} item(s), got "
                f"{len(value)}.{tail}")

    if spec.kind == "path" and spec.must_exist:
        # FileNotFoundError rather than ValueError: a missing file is a
        # different problem from a malformed argument, and the caller
        # usually wants to tell them apart.
        if not Path(value).is_file():
            raise FileNotFoundError(f"{name}: file not found: {value}.{tail}")


def check_all(specs: dict, values: dict, *, where: str = "") -> None:
    """Check every value that has a spec.

    Values without a spec are ignored here; `coverage` is what stops a
    new argument from quietly arriving unchecked, because that belongs
    in a test rather than at runtime.
    """
    for name, spec in specs.items():
        if name in values:
            try:
                check(name, values[name], spec)
            except (ValueError, FileNotFoundError) as exc:
                if where:
                    raise type(exc)(f"{where}: {exc}") from None
                raise


def coverage(specs: dict, func) -> tuple:
    """`(missing, extra)` between a spec table and a function signature.

    The point of the table is that it cannot drift from the thing it
    describes. A new argument added to `run` with no rule shows up in
    `missing`, and a rule for an argument that no longer exists shows up
    in `extra`.
    """
    import inspect

    params = set(inspect.signature(func).parameters)
    params.discard("self")
    params.discard("cls")
    return tuple(sorted(params - set(specs))), tuple(
        sorted(set(specs) - params))
