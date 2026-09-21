"""Lineage tree data structures for timing certificates."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np


@dataclass
class Node:
    name: str
    parent: Optional[str] = None
    children: List[str] = field(default_factory=list)
    time_from_parent: float = 0.0  # branch duration into this node
    age: float = 0.0  # time from root
    is_leaf: bool = False


class LineageTree:
    """Rooted tree with named nodes and branch durations."""

    def __init__(self, root: str = "root"):
        self.root = root
        self.nodes: Dict[str, Node] = {root: Node(name=root, age=0.0)}

    def add_child(self, parent: str, child: str, time_from_parent: float = 0.0) -> None:
        if parent not in self.nodes:
            raise KeyError(f"Unknown parent {parent}")
        if child in self.nodes:
            raise ValueError(f"Node {child} already exists")
        self.nodes[parent].children.append(child)
        self.nodes[parent].is_leaf = False
        age = self.nodes[parent].age + time_from_parent
        self.nodes[child] = Node(
            name=child,
            parent=parent,
            time_from_parent=time_from_parent,
            age=age,
            is_leaf=True,
        )

    def edges(self) -> List[Tuple[str, str]]:
        # Deterministic by child name: copy() rebuilds self.nodes in BFS order, so a
        # dict-insertion ordering would mis-pair branch-times between an original tree
        # and its copy (observed ages > T on the ctree->copy boundary).
        return sorted(
            ((n.parent, n.name) for n in self.nodes.values() if n.parent is not None),
            key=lambda e: e[1],
        )

    def edge_index(self) -> Dict[Tuple[str, str], int]:
        return {e: i for i, e in enumerate(self.edges())}

    def leaves(self) -> List[str]:
        return [n for n, node in self.nodes.items() if node.is_leaf or not node.children]

    def internal_nodes(self) -> List[str]:
        return [n for n, node in self.nodes.items() if node.children]

    def path_to_root(self, node: str) -> List[str]:
        path = [node]
        while self.nodes[node].parent is not None:
            node = self.nodes[node].parent  # type: ignore[assignment]
            path.append(node)
        return path

    def path_edges(self, node: str) -> List[Tuple[str, str]]:
        path = self.path_to_root(node)
        return list(zip(path[1:], path[:-1]))  # parent → child toward leaf

    def branch_times_vector(self) -> np.ndarray:
        return np.array([self.nodes[c].time_from_parent for _, c in self.edges()], dtype=float)

    def set_branch_times(self, times: Sequence[float]) -> None:
        edges = self.edges()
        if len(times) != len(edges):
            raise ValueError("times length must match number of edges")
        for (p, c), t in zip(edges, times):
            self.nodes[c].time_from_parent = float(t)
        self._recompute_ages()

    def _recompute_ages(self) -> None:
        def dfs(u: str) -> None:
            for v in self.nodes[u].children:
                self.nodes[v].age = self.nodes[u].age + self.nodes[v].time_from_parent
                dfs(v)

        self.nodes[self.root].age = 0.0
        dfs(self.root)

    def total_depth(self) -> float:
        leaves = self.leaves()
        if not leaves:
            return 0.0
        return float(max(self.nodes[ℓ].age for ℓ in leaves))

    def is_ultrametric(self, tol: float = 1e-8) -> bool:
        ages = [self.nodes[ℓ].age for ℓ in self.leaves()]
        return max(ages) - min(ages) <= tol

    def node_age(self, node: str) -> float:
        return float(self.nodes[node].age)

    def basal_split(self) -> Tuple[str, str, str]:
        """First resolved (>=2-child) split: (split_node, child_a, child_b).

        Walks past single-child chains near the root — real inferred trees often
        have an unresolved basal collapse (root -> one child -> ...).
        """
        u = self.root
        while len(self.nodes[u].children) < 2:
            if not self.nodes[u].children:
                raise ValueError("single-leaf tree has no split")
            u = self.nodes[u].children[0]
        a, b = sorted(self.nodes[u].children[:2])  # deterministic target naming
        return u, a, b

    def induced_subtree(
        self, keep_leaves: Sequence[str], root: str = "root"
    ) -> "LineageTree":
        """Induced subtree over keep_leaves, degree-2 (non-root) nodes suppressed.

        Branch lengths are summed through suppressed nodes so leaf ages are
        preserved. Used to keep the profile solver's free-parameter count tractable
        while using the official Dryad edit states.
        """
        keep = set(keep_leaves)
        nodes = set(keep)
        for ℓ in keep:
            path = self.path_to_root(ℓ)
            nodes.update(path)

        def distance(a: str, b: str) -> float:
            return abs(self.nodes[a].age - self.nodes[b].age)

        edges = [(p, c) for c in nodes for p in [self.nodes[c].parent] if p in nodes]
        # suppress degree-2 non-root nodes
        degree = {n: sum(1 for (p, c) in edges if c == n) + sum(1 for (p, c) in edges if p == n) for n in nodes}
        kept = set(nodes)
        for n in list(nodes):
            if n == root or n in keep:
                continue
            if degree[n] == 2:
                kept.discard(n)
        new_edges = []
        for (p, c) in edges:
            if c not in kept:
                continue
            up = p
            while up not in kept:
                up = self.nodes[up].parent  # type: ignore[assignment]
            new_edges.append((up, c, distance(up, c)))
        return LineageTree.from_edges(new_edges, root=root)

    def copy(self) -> "LineageTree":
        t = LineageTree(self.root)
        # rebuild in preorder
        order = [self.root]
        seen = {self.root}
        while order:
            u = order.pop(0)
            for v in self.nodes[u].children:
                if v not in seen:
                    t.add_child(u, v, self.nodes[v].time_from_parent)
                    seen.add(v)
                    order.append(v)
        return t

    def newick(self) -> str:
        def fmt(u: str) -> str:
            node = self.nodes[u]
            if not node.children:
                return f"{u}:{node.time_from_parent:.6g}"
            inner = ",".join(fmt(c) for c in node.children)
            if u == self.root:
                return f"({inner}){u};"
            return f"({inner}){u}:{node.time_from_parent:.6g}"

        return fmt(self.root)

    @staticmethod
    def balanced_binary(n_leaves: int, total_time: float, names: Optional[Sequence[str]] = None) -> "LineageTree":
        """Build a balanced ultrametric binary tree with equal branch times per level."""
        if n_leaves < 2 or (n_leaves & (n_leaves - 1)) != 0:
            raise ValueError("n_leaves must be a power of 2, >= 2")
        depth = int(np.log2(n_leaves))
        t_level = total_time / depth
        tree = LineageTree("root")
        leaf_names = list(names) if names is not None else [f"L{i}" for i in range(n_leaves)]
        # BFS create internal structure
        frontier = ["root"]
        next_id = 0
        for level in range(depth - 1):
            new_frontier = []
            for u in frontier:
                for _ in range(2):
                    name = f"N{next_id}"
                    next_id += 1
                    tree.add_child(u, name, t_level)
                    new_frontier.append(name)
            frontier = new_frontier
        # attach leaves
        li = 0
        for u in frontier:
            for _ in range(2):
                tree.add_child(u, leaf_names[li], t_level)
                li += 1
        return tree

    @staticmethod
    def from_edges(
        edges: Iterable[Tuple[str, str, float]],
        root: str = "root",
    ) -> "LineageTree":
        tree = LineageTree(root)
        remaining = list(edges)
        safety = 0
        while remaining and safety < 10000:
            safety += 1
            progress = False
            still = []
            for p, c, t in remaining:
                if p in tree.nodes and c not in tree.nodes:
                    tree.add_child(p, c, t)
                    progress = True
                else:
                    still.append((p, c, t))
            remaining = still
            if not progress:
                break
        if remaining:
            raise ValueError(f"Could not attach edges: {remaining}")
        return tree


def small_scout_tree(total_time: float = 1.0) -> LineageTree:
    """Canonical 4-leaf tree for the decisive scout.

          root
         /    \\
       A        B
      / \\      / \\
     L0 L1   L2 L3
    """
    t = total_time / 2.0
    return LineageTree.from_edges(
        [
            ("root", "A", t),
            ("root", "B", t),
            ("A", "L0", t),
            ("A", "L1", t),
            ("B", "L2", t),
            ("B", "L3", t),
        ],
        root="root",
    )
