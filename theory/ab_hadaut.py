r"""Automorphism group of a Hadamard matrix (pairs of signed permutations P H Q^T = H) via the standard colored
graph (rows r_i^+-, columns c_j^+-; r_i^s ~ c_j^t iff s t H_ij = +1) and BLISS (python-igraph) [M].
usage (library): aut(H) -> (group order, generators as vertex permutations)"""
import numpy as np
import igraph as ig


def graph(H):
    n = H.shape[0]
    edges = []
    for i in range(n):
        for j in range(n):
            for s in (0, 1):
                for t in (0, 1):
                    if (1 - 2 * s) * (1 - 2 * t) * H[i, j] == 1:
                        edges.append((2 * i + s, 2 * n + 2 * j + t))
    g = ig.Graph(n=4 * n, edges=edges)
    colors = [0] * (2 * n) + [1] * (2 * n)
    return g, colors


def aut(H):
    g, colors = graph(np.asarray(H))
    gens = g.automorphism_group(color=colors)
    order = g.count_automorphisms(color=colors)
    return int(order), gens
