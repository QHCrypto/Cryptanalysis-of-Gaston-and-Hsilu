"""Exact first-round nonce quadratic forms and reproducible cube subspaces."""
import json
import random
from pathlib import Path
from permutations import gaston_linear, gaston_effect, chi

ROOT = Path(__file__).resolve().parents[1]


def state_from_nonce(n):
    return [0, 0, 0, n & ((1 << 64) - 1), n >> 64]


def packed(s):
    return sum(v << (64 * i) for i, v in enumerate(s))


def rank(values):
    basis = {}
    for v in values:
        while v:
            j = v.bit_length() - 1
            if j in basis:
                v ^= basis[j]
            else:
                basis[j] = v
                break
    return len(basis)


def nullspace_from_columns(cols):
    pivots = {}
    null = []
    for j, col in enumerate(cols):
        v, x = col, 1 << j
        while v:
            p = v.bit_length() - 1
            if p not in pivots:
                pivots[p] = (v, x)
                break
            v ^= pivots[p][0]
            x ^= pivots[p][1]
        if not v:
            null.append(x)
    return null


def polar(a, b):
    # Second derivative of chi, independent of the affine base and the key.
    return packed([((a[(i + 1) % 5] & b[(i + 2) % 5]) ^
                    (b[(i + 1) % 5] & a[(i + 2) % 5])) for i in range(5)])


def coordinate_graph(name):
    inputs = [state_from_nonce(1 << j) for j in range(128)]
    if name == "gaston":
        inputs = [gaston_linear(s) for s in inputs]
    adjacent = [set() for _ in range(128)]
    for i in range(128):
        for j in range(i):
            if polar(inputs[i], inputs[j]):
                adjacent[i].add(j)
                adjacent[j].add(i)
    return adjacent


def independent_greedy(graph, trials=2000):
    rng = random.Random(20260928)
    best = []
    for trial in range(trials):
        remaining = set(range(len(graph)))
        selected = []
        while remaining:
            options = sorted(remaining, key=lambda v: (len(graph[v] & remaining), v))
            v = options[rng.randrange(min(4, len(options)))] if trial else options[0]
            selected.append(v)
            remaining.difference_update(graph[v] | {v})
        if len(selected) > len(best):
            best = sorted(selected)
    assert all(j not in graph[i] for i in best for j in best)
    return best


def maximum_independent_exact(graph):
    """Exhaustive branch-and-bound on complement cliques; no heuristic pruning."""
    n = len(graph)
    complement = [sum(1 << j for j in range(n) if j != i and j not in graph[i])
                  for i in range(n)]
    best = []
    nodes = 0

    def visit(chosen, candidates):
        nonlocal best, nodes
        nodes += 1
        if len(chosen) + candidates.bit_count() <= len(best):
            return
        if len(chosen) > len(best):
            best = chosen
        while candidates:
            if len(chosen) + candidates.bit_count() <= len(best):
                return
            bit = candidates & -candidates
            v = bit.bit_length() - 1
            candidates ^= bit
            visit(chosen + [v], candidates & complement[v])

    visit([], (1 << n) - 1)
    return best, nodes


def run():
    report = {"method": "exact polar forms; exhaustive complement clique for Gaston; disjoint-edge proof for Hsilu"}
    from permutations import hsilu_linear
    report["linear_layer_ranks"] = {}
    for name, mapping in (("gaston",gaston_linear),("hsilu",hsilu_linear)):
        columns = []
        for j in range(320):
            s = [0]*5
            s[j//64] = 1 << (j%64)
            columns.append(packed(mapping(s)))
        report["linear_layer_ranks"][name] = rank(columns)
        assert rank(columns) == 320
    for name in ("gaston", "hsilu"):
        graph = coordinate_graph(name)
        if name == "gaston":
            chosen, nodes = maximum_independent_exact(graph)
        else:
            assert all(graph[i] == {i ^ 64} for i in range(128))
            chosen, nodes = list(range(64)), 0
        report[name] = {
            "coordinate_conflict_edges": sum(map(len, graph)) // 2,
            "degree_range": [min(map(len, graph)), max(map(len, graph))],
            "certified_affine_nonce_coordinate_cube": chosen,
            "dimension": len(chosen),
            "optimal_in_coordinate_cube_class": True,
            "exhaustive_search_nodes": nodes,
        }
    columns = [gaston_effect(state_from_nonce(1 << j)) for j in range(128)]
    kernel = nullspace_from_columns(columns)
    assert all(gaston_effect(state_from_nonce(v)) == 0 for v in kernel)
    report["gaston_nonce_effect_kernel"] = {
        "rank": rank(columns), "dimension": len(kernel),
        "basis_hex": [f"{x:032x}" for x in kernel],
    }
    kernel_images = [gaston_linear(state_from_nonce(v)) for v in kernel]
    for i, a in enumerate(kernel_images):
        for j, b in enumerate(kernel_images[:i]):
            if polar(a, b):
                report["gaston_nonce_effect_kernel"]["non_affine_witness"] = {
                    "basis_indices": [i,j], "chi_polar_hex": f"{polar(a,b):080x}"}
                break
        if "non_affine_witness" in report["gaston_nonce_effect_kernel"]:
            break
    # Intersection with equal nonce forms at the input of chi:
    # ROL(N0, 9) = ROL(N1, 47), hence N1 = ROL(N0,26).
    from permutations import rol
    equal_cols = [gaston_effect(state_from_nonce((1 << j) | (rol(1 << j, 26) << 64)))
                  for j in range(64)]
    equal_kernel = nullspace_from_columns(equal_cols)
    equal_basis = [v | (rol(v, 26) << 64) for v in equal_kernel]
    images = [gaston_linear(state_from_nonce(v)) for v in equal_basis]
    assert all(polar(a, b) == 0 for a in images for b in images)
    report["gaston_effect_zero_and_equal_chi_nonce_rows"] = {
        "dimension": len(equal_basis),
        "basis_hex": [f"{x:032x}" for x in equal_basis],
        "relation": "N1=ROL(N0,26)",
    }
    (ROOT / "results" / "first_round_structure.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({k: v if k == "method" else {a:b for a,b in v.items() if a != "basis_hex"}
                      for k,v in report.items()}, indent=2))


if __name__ == "__main__":
    run()
