"""Adapters for the GESTALT zebrafish lineage-tracing benchmark (McKenna 2016)."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from b01_clock.model.observations import CharacterMatrix, EDITED, MISSING, UNEDITED
from b01_clock.model.tree import LineageTree

GESTALT_GEO = "GSE81713"
GESTALT_DRYAD = None
GESTALT_ROLE = "real_data_consistency"

# GESTALT cells were harvested at 5 days post-fertilization (5 dpf)
GESTALT_TOTAL_TIME_DP = 5.0


def load_mix_matrix(path: str | Path) -> "tuple[CharacterMatrix, dict]":
    """Parse a PHYLIP MIX input matrix (name + binary string per taxon).

    Returns (matrix, meta) where meta records the header vs. observed site
    counts. Some GESTALT MIX files declare 1154 but carry 1156 consistent
    columns; we trust the self-consistent observed length.
    """
    path = Path(path)
    import gzip as _gzip
    opener = _gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        header = f.readline()
        declared = int(header.strip().split()[1])
        leaf_states = {}
        for line in f:
            tok = line.split()
            if len(tok) < 2:
                continue
            name, vec = tok[0], "".join(tok[1:])
            row = np.zeros(len(vec), dtype=int)
            for j, ch in enumerate(vec):
                if ch == "1":
                    row[j] = EDITED
                elif ch == "0":
                    row[j] = UNEDITED
                else:
                    row[j] = MISSING
            leaf_states[name] = row
    lens = {len(v) for v in leaf_states.values()}
    n_sites = int(max(lens))
    matrix = CharacterMatrix(leaf_states=leaf_states, n_sites=n_sites)
    meta = {"declared_n_sites": declared, "observed_n_sites": n_sites,
            "site_count_consistent": len(lens) == 1}
    return matrix, meta


def parse_newick(newick: str, total_time: float) -> LineageTree:
    """Minimal Newick parser for labeled trees (no branch lengths required).

    PHYLIP MIX output is topology-only. Unnamed internal nodes get N<k> names;
    edges default to total_time/2 for a sane age init (the optimizer rescales
    ages anyway). Raises on malformed input instead of silently falling back.
    """
    import re
    # strip PHYLIP bracket annotations ([0.0100] axis-length tail etc.) and whitespace
    s = re.sub(r"\[[^\]]*\]", "", newick)
    s = re.sub(r"\s+", "", s)
    pos = [0]
    counter = {"n": 0, "i": 0}

    def name_of(label):
        if label:
            return label
        n = f"int{counter['n']}"
        counter["n"] += 1
        return n

    # Returns node record: {"name", "children": [(name, brlen), ...]}
    def parse_clade():
        if s[pos[0]] == "(":
            pos[0] += 1
            child_recs = []
            while True:
                child_recs.append(parse_clade())
                if s[pos[0]] == ",":
                    pos[0] += 1
                    continue
                break
            if s[pos[0]] != ")":
                raise ValueError(f"expected ')' at {s[pos[0]:pos[0]+12]!r}")
            pos[0] += 1
            name = parse_label()
            rec = {"name": name, "children": child_recs}
            return rec
        name = parse_label()
        return {"name": name, "children": ()}

    def parse_label():
        if s[pos[0]] in "(),;":
            return ""  # unnamed internal node
        m = re.match(r"[^():,;]+", s[pos[0]:])
        if not m:
            raise ValueError(f"parse error at {s[pos[0]:pos[0]+12]!r}")
        label = m.group(0)
        pos[0] += len(label)
        if s[pos[0]:pos[0] + 1] == ":":
            pos[0] += 1
            m2 = re.match(r"[\d.eE+-]+", s[pos[0]:])
            pos[0] += len(m2.group(0)) if m2 else 0
        return label if label != "root" else "root"

    root_clades = []
    while pos[0] < len(s) and s[pos[0]] != ";":
        root_clades.append(parse_clade())
    if not root_clades:
        raise ValueError("empty newick")

    tree = LineageTree("root")

    def attach(rec, parent):
        if rec["name"] and rec["name"] != "root":
            name = rec["name"]
            if name in tree.nodes:  # duplicate label in this dialect
                name = f"{name}_dup{counter['n']}"
                counter["n"] += 1
        else:
            name = f"int{counter['n']}"
            counter["n"] += 1
        if parent in tree.nodes:
            tree.add_child(parent, name, total_time / 2.0)
        for ch in rec["children"]:
            attach(ch, name)
        return name

    for rec in root_clades:
        attach(rec, "root")
    return tree


def load_gestalt_newick(path: str | Path, total_time: float = GESTALT_TOTAL_TIME_DP) -> LineageTree:
    """Load a PHYLIP MIX newick (topology-only) as a LineageTree."""
    path = Path(path)
    import gzip as _gzip
    opener = _gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        newick = f.read()
    return parse_newick(newick, total_time)


def provenance_record(paths: dict, total_time: float = GESTALT_TOTAL_TIME_DP) -> dict:
    """Return a reproducible provenance record for an auditable run."""
    import hashlib, json
    shas = {}
    for key, p in paths.items():
        if p is None or not Path(p).exists():
            shas[key] = None
            continue
        h = hashlib.sha256()
        with Path(p).open("rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        shas[key] = h.hexdigest()
    return {
        "accession": GESTALT_GEO,
        "role": GESTALT_ROLE,
        "total_time_dp": total_time,
        "files": {k: str(v) for k, v in paths.items()},
        "sha256": shas,
    }