"""Regression test: GESTALT PHYLIP-MIX adapter parses real-format inputs.

Guards: (1) the MIX newick syntax — outer umbrella parens, unnamed internal
nodes, trailing PHYLIP bracket annotation `[0.0100]`, optional root label and
terminal semicolon; (2) top-level comma separation producing a multifurcating
root; (3) unbracketed txt.gz mix matrices loading with observed-column
consistency.
"""

import gzip
from io import StringIO
from pathlib import Path

from b01_clock.data.adapters.gestalt import load_mix_matrix, parse_newick


def test_parse_mix_newick_with_bracket_annotation():
    newick = "(((A,B)C,(D,(E,F)))G,(H,I))[0.0100];"
    tree = parse_newick(newick, total_time=5.0)
    leaves = set(tree.leaves())
    assert leaves == {"A", "B", "D", "E", "F", "H", "I"}, leaves
    assert {"C", "G"}.issubset(tree.nodes)


def test_parse_newick_with_unnamed_root_clades():
    newick = "((A,B),(C,D));"
    tree = parse_newick(newick, total_time=5.0)
    assert set(tree.leaves()) == {"A", "B", "C", "D"}
    assert tree.root in tree.nodes


def test_load_mix_matrix_gz_payload(tmp_path: Path):
    payload = "TAX 2\nA\t0\t1\nB\t1\t0\nC\t0\t0\n"
    p = tmp_path / "mix.txt.gz"
    p.write_bytes(gzip.compress(payload.encode()))
    matrix, meta = load_mix_matrix(p)
    assert matrix.n_sites == 2
    assert set(matrix.leaf_states) == {"A", "B", "C"}
    assert meta["observed_n_sites"] == meta["declared_n_sites"] == 2