"""Multi-hop upstream search: promoter -> TF -> kinase -> receptor -> drug.

Why this exists
---------------
The rest of the package is one hop: a drug that binds your TF. That ceiling is
severe, because most of the ~1,600 human TFs have no direct pharmacology at
all. CREB1 has no approved ligand — but PKA does, adenylyl cyclase does, and
the GPCRs above them do. Walking upstream turns "no drugs for this factor" into
a ranked list of real options.

It also introduces four ways to be confidently wrong. All four are handled
here, and it is worth knowing which, because network-based repurposing papers
routinely skip them.

**1. Signs must propagate.** Inhibiting an inhibitor is activation. A path's
net direction is the *product* of its edge signs. A 4-hop path with two
inhibitory edges is net activating. Unsigned networks - plain PPI, most
STRING-style resources - cannot support this analysis at all, which is why
OmniPath/SIGNOR signed causal edges are required rather than optional.

**2. Specificity collapses with every hop.** This is the uncomfortable one,
because it works directly against the specificity axis the rest of the package
optimises. A receptor four steps upstream reaches your gene through one branch
and several thousand other genes through the others. Distance buys
druggability and spends selectivity. Hop count therefore carries an explicit
penalty and out-degree fan-out is reported per path, so the trade is visible
instead of hidden inside a rank.

**3. Random walks find hubs, not biology.** Node degree in curated networks
tracks how much a protein has been *studied*. AKT1, TP53, EGFR and MAPK1 have
enormous edge counts, so an unnormalised walk lands on them from any seed. The
walk here is degree-normalised, and hub nodes are flagged in the output. If
your top hit is AKT1 from every gene you try, that is the artefact, not a
finding.

**4. Confidence decays multiplicatively.** Each curated edge is uncertain.
Four hops at 0.75 per edge is 0.32 before the drug is even considered. That
decay is applied, not waved at.

Data
----
OmniPath (https://omnipathdb.org) aggregates signed, directed, literature-
curated causal interactions - kinase-substrate, TF-target, receptor-effector -
with provenance. No key required. SIGNOR is included within it and is the
highest-quality signed causal source.

Paths are reported alongside every ranking. A rank with no mechanism attached
is not actionable; you should be able to read the chain and judge whether it
is plausible in your cell type.
"""

from __future__ import annotations

import logging
from collections import defaultdict, deque
from dataclasses import dataclass, asdict, field

import numpy as np

from .http import Client, HttpError

log = logging.getLogger(__name__)

OMNIPATH = "https://omnipathdb.org/interactions"

# Degree above which a node is treated as a study-bias hub rather than a
# specific route. These are the nodes an unnormalised walk always finds.
HUB_DEGREE = 150
PER_HOP_CONFIDENCE = 0.75
MAX_DEPTH = 4


@dataclass
class Edge:
    source: str
    target: str
    sign: int                 # +1 stimulation, -1 inhibition, 0 unsigned
    n_references: int = 0
    sources: str = ""


@dataclass
class Path:
    nodes: list                # drug target -> ... -> TF
    signs: list                # per-edge signs
    net_sign: int
    hops: int
    confidence: float
    fan_out: int              # downstream breadth of the entry node
    specificity: float
    hub_nodes: list = field(default_factory=list)

    def describe(self) -> str:
        parts = []
        for i, node in enumerate(self.nodes):
            parts.append(node)
            if i < len(self.signs):
                parts.append("—|" if self.signs[i] < 0 else "-->")
        return " ".join(parts)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["description"] = self.describe()
        return d


@dataclass
class UpstreamHit:
    node: str                 # the druggable protein upstream
    tf_symbol: str
    best_path: dict
    n_paths: int
    score: float
    rwr_score: float
    net_sign: int
    # Which way the TF itself must move. Derived from the run's goal and
    # the factor's role, and recorded here because `attach_drugs` cannot
    # recover it: `upstream_paths` keeps only paths whose net sign equals
    # this, so net_sign alone carries no information about the goal.
    want_tf_sign: int = 1
    drugs: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


# ------------------------------------------------------------------ loading
def fetch_network(client: Client, organism: int = 9606,
                  datasets: str = "omnipath,kinaseextra,tfregulons",
                  min_references: int = 1) -> list[Edge]:
    """Signed, directed, literature-curated causal edges from OmniPath."""
    try:
        res = client.get(OMNIPATH, {
            "organisms": organism,
            "datasets": datasets,
            "fields": "sources,references,curation_effort",
            "directed": 1,
            "signed": 1,
            "genesymbols": 1,
            "format": "json",
        })
    except HttpError as exc:
        log.warning("OmniPath fetch failed: %s", exc)
        return []

    rows = res if isinstance(res, list) else res.get("data", []) or []
    edges: list[Edge] = []
    for r in rows:
        if not isinstance(r, dict):
            continue
        src = r.get("source_genesymbol") or r.get("source")
        tgt = r.get("target_genesymbol") or r.get("target")
        if not src or not tgt or src == tgt:
            continue
        stim = bool(int(r.get("is_stimulation") or 0))
        inhib = bool(int(r.get("is_inhibition") or 0))
        if stim == inhib:          # unsigned or contradictory: unusable here
            continue
        refs = str(r.get("references") or "")
        n_refs = len([x for x in refs.split(";") if x.strip()])
        if n_refs < min_references:
            continue
        edges.append(Edge(str(src).upper(), str(tgt).upper(),
                          1 if stim else -1, n_refs,
                          str(r.get("sources") or "")))
    log.info("OmniPath: %d signed directed edges", len(edges))
    return edges


def build_graph(edges: list[Edge]):
    """Return (incoming, outgoing, degree) adjacency maps."""
    incoming: dict[str, list[Edge]] = defaultdict(list)
    outgoing: dict[str, list[Edge]] = defaultdict(list)
    for e in edges:
        incoming[e.target].append(e)
        outgoing[e.source].append(e)
    degree = {n: len(incoming[n]) + len(outgoing[n])
              for n in set(incoming) | set(outgoing)}
    return dict(incoming), dict(outgoing), degree


# ------------------------------------------------------- upstream traversal
def upstream_paths(incoming: dict, outgoing: dict, degree: dict, tf: str, *,
                   max_depth: int = MAX_DEPTH, want_sign: int = 1,
                   max_paths: int = 2000) -> list[Path]:
    """Enumerate signed paths INTO `tf`, walking against edge direction.

    Only paths whose product of signs matches `want_sign` are kept: those are
    the ones that would move the TF in the direction you asked for.
    """
    tf = tf.upper()
    results: list[Path] = []
    queue = deque([(tf, [tf], [])])
    seen_states = set()

    while queue and len(results) < max_paths:
        node, chain, signs = queue.popleft()
        if len(chain) - 1 >= max_depth:
            continue
        for e in incoming.get(node, []):
            up = e.source
            if up in chain:                      # no cycles
                continue
            new_chain = [up] + chain
            new_signs = [e.sign] + signs
            state = (up, len(new_chain))
            if state in seen_states and len(results) > max_paths // 2:
                continue
            seen_states.add(state)

            net = int(np.prod(new_signs))
            hops = len(new_chain) - 1
            if net == want_sign:
                fan = len(outgoing.get(up, []))
                conf = PER_HOP_CONFIDENCE ** hops
                # Specificity falls with both distance and branching.
                spec = (1.0 / (1.0 + 0.6 * (hops - 1))) * (1.0 / (1.0 + fan / 40.0))
                hubs = [n for n in new_chain if degree.get(n, 0) >= HUB_DEGREE]
                results.append(Path(new_chain, new_signs, net, hops,
                                    round(conf, 3), fan, round(spec, 3), hubs))
            queue.append((up, new_chain, new_signs))

    results.sort(key=lambda p: (-p.confidence * p.specificity, p.hops))
    return results


# ------------------------------------------------- random walk with restart
def random_walk_with_restart(incoming: dict, outgoing: dict, degree: dict,
                             seeds: list[str], *, restart: float = 0.3,
                             iterations: int = 60,
                             degree_normalise: bool = True) -> dict[str, float]:
    """RWR over the reversed graph (upstream propagation from the TF).

    Degree normalisation is on by default and should stay on. Without it the
    walk concentrates on whichever nodes have been studied most, and every
    query returns the same handful of hubs.
    """
    nodes = sorted(set(incoming) | set(outgoing))
    if not nodes:
        return {}
    index = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)

    # Transition matrix over reversed edges: probability flows TF -> upstream.
    T = np.zeros((n, n), dtype=np.float32)
    for target, edges in incoming.items():
        j = index[target]
        for e in edges:
            i = index[e.source]
            w = 1.0
            if degree_normalise:
                # Penalise arriving at well-connected (well-studied) nodes.
                w /= (1.0 + np.log1p(degree.get(e.source, 1)))
            T[i, j] += w

    col_sums = T.sum(axis=0, keepdims=True)
    col_sums[col_sums == 0] = 1.0
    T /= col_sums

    p0 = np.zeros(n, dtype=np.float32)
    seed_idx = [index[s.upper()] for s in seeds if s.upper() in index]
    if not seed_idx:
        return {}
    p0[seed_idx] = 1.0 / len(seed_idx)

    p = p0.copy()
    for _ in range(iterations):
        p_new = (1 - restart) * (T @ p) + restart * p0
        if np.abs(p_new - p).sum() < 1e-9:
            p = p_new
            break
        p = p_new

    return {nodes[i]: float(p[i]) for i in np.nonzero(p)[0]}


# ------------------------------------------------------------------ driver
def find_upstream_targets(client: Client, tf: str, *, goal: str = "up",
                          tf_role: str = "activator", max_depth: int = MAX_DEPTH,
                          edges: list[Edge] | None = None,
                          top_n: int = 25) -> list[UpstreamHit]:
    """Rank druggable proteins upstream of a TF, with the path that explains it."""
    # To raise the gene via an activating TF you must raise the TF; via a
    # repressor you must lower it. Then invert again for a downregulation goal.
    want_tf_direction = 1 if (tf_role in ("activator", "dual")) == (goal == "up") else -1

    edges = edges if edges is not None else fetch_network(client)
    if not edges:
        return []
    incoming, outgoing, degree = build_graph(edges)
    if tf.upper() not in incoming:
        log.info("%s has no incoming signed edges in the network", tf)
        return []

    paths = upstream_paths(incoming, outgoing, degree, tf,
                           max_depth=max_depth, want_sign=want_tf_direction)
    rwr = random_walk_with_restart(incoming, outgoing, degree, [tf])

    by_entry: dict[str, list[Path]] = defaultdict(list)
    for p in paths:
        by_entry[p.nodes[0]].append(p)

    hits: list[UpstreamHit] = []
    for node, node_paths in by_entry.items():
        best = max(node_paths, key=lambda p: p.confidence * p.specificity)
        r = rwr.get(node, 0.0)
        score = 100.0 * best.confidence * best.specificity * (1.0 + min(r * 50, 1.0))

        warnings = []
        if degree.get(node, 0) >= HUB_DEGREE:
            warnings.append(
                f"{node} is a network hub (degree {degree[node]}); high rank here "
                "may reflect how well studied it is rather than a specific route")
        if best.hops >= 3:
            warnings.append(
                f"{best.hops} hops upstream: druggable, but it reaches your gene "
                f"through one of {best.fan_out} downstream branches")
        if best.hub_nodes:
            warnings.append("path passes through hub(s): "
                            + ", ".join(best.hub_nodes))

        hits.append(UpstreamHit(
            node=node, tf_symbol=tf.upper(), best_path=best.to_dict(),
            n_paths=len(node_paths), score=round(score, 1),
            rwr_score=round(r, 6), net_sign=best.net_sign,
            want_tf_sign=want_tf_direction, warnings=warnings))

    hits.sort(key=lambda h: -h.score)
    return hits[:top_n]


def attach_drugs(client: Client, hits: list[UpstreamHit], *, goal_on_node: str = "up",
                 approved_only: bool = True) -> list[UpstreamHit]:
    """Find approved drugs for each upstream node, respecting the path sign."""
    from .drugs import chembl_drugs

    for h in hits:
        # Two signs decide this, and only one of them was being read.
        #
        #   n * net_sign = want_tf_sign,  so  n = want_tf_sign * net_sign
        #
        # where n is the direction the DRUG must move the upstream node.
        # The previous version used `want = "up" if net_sign == 1`, which
        # silently assumes want_tf_sign is always +1. It is not:
        # `find_upstream_targets` flips it for a repressor, and again for
        # a downregulation goal. Worse, `upstream_paths` keeps only paths
        # whose net sign EQUALS want_tf_sign, so the two were always equal
        # and the old expression reduced to "follow want_tf_sign" - giving
        # the opposite of the right answer whenever it was -1.
        #
        # Concretely, for KINASE -| REPRESSOR with the gene to be raised:
        # you lower the repressor, which means RAISING the kinase. The old
        # code asked for a kinase inhibitor, which raises the repressor and
        # lowers the gene.
        want_tf = getattr(h, "want_tf_sign", None)
        if want_tf is None:
            want_tf = 1 if goal_on_node == "up" else -1
        want = "up" if want_tf * h.net_sign == 1 else "down"
        try:
            found = chembl_drugs(client, h.node, None, approved_only)
        except Exception as exc:                             # noqa: BLE001
            log.debug("drug lookup failed for %s: %s", h.node, exc)
            continue
        h.drugs = [d.to_dict() for d in found if d.direction in (want, "unclear")][:5]
    return hits
