"""Which filter of the upstream search removes a drug that is in reach.

The replay is only useful if it applies the search's own rules: paths of three
steps or fewer, signs that multiply to the effect the goal needs, and a drug
whose action at the node matches what that path asks for.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

from generoutes.network import Edge, build_graph


def _script():
    path = Path(__file__).resolve().parents[1] / "scripts" / "upstream_misses.py"
    spec = importlib.util.spec_from_file_location("_upstream_misses", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


um = _script()


def graph(*edges):
    return build_graph([Edge(a, b, s, 1, "") for a, b, s in edges])


def test_signs_multiply_along_a_path_and_depth_is_limited():
    incoming, _, _ = graph(("K", "M", -1), ("M", "F", -1), ("X", "K", 1),
                           ("Y", "X", 1), ("Z", "Y", 1))
    reach = um.signed_reach(incoming, "F")
    assert reach["M"] == {-1}
    assert reach["K"] == {1}, "K inhibits M and M inhibits F, so K raises F"
    assert reach["X"] == {1} and "Y" not in reach, "four steps is too far"


def test_a_node_reached_by_both_signs_keeps_both():
    incoming, _, _ = graph(("A", "B", 1), ("A", "C", -1), ("B", "F", 1),
                           ("C", "F", 1))
    assert um.signed_reach(incoming, "F")["A"] == {1, -1}


def _classify(drug, targets, want=1, dirs=None, searched=True, ranked=True):
    incoming, _, _ = graph(("N", "F", 1))
    sr = um.signed_reach(incoming, "F")
    return um.classify(
        drug, set(targets), {"F": (want, sr)} if searched else {},
        {"N": 1} if ranked else {},
        {um.breadth_key(drug): {"N": dirs}} if dirs else {})


def test_a_drug_passing_every_filter_is_left_to_the_limits():
    assert _classify("D", ["N"], want=1, dirs={"up"}) == um.SURVIVES
    assert _classify("D", ["N"], want=1, dirs={"unclear"}) == um.SURVIVES


def test_the_wrong_sign_the_wrong_way_and_the_missing_index_are_told_apart():
    assert _classify("D", ["N"], want=-1, dirs={"up"}) == um.WRONG_SIGN
    assert _classify("D", ["N"], want=1, dirs={"down"}) == um.WRONG_WAY
    assert _classify("D", ["N"], want=1, dirs=None) == um.NOT_INDEXED


def test_a_factor_that_was_not_searched_is_named_as_such():
    assert _classify("D", ["N"], searched=False) == um.OUTSIDE
    assert _classify("D", ["N"], searched=False, ranked=False) == um.UNREACHED
