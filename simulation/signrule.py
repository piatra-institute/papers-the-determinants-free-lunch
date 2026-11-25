"""Sign-rule faults in the Leibniz expansion.

For a set E of row pairs, D_E(A) = sum_sigma sgn(sigma) (-1)^{inv_E(sigma)} prod_i A[i, sigma(i)],
where inv_E(sigma) counts inversions {i < j, sigma(i) > sigma(j)} with {i, j} in E. E = {} gives det,
E = all pairs gives perm. D_E is a determinant after entry sign changes iff
(-1)^{inv_E(sigma)} = c prod_i s[i, sigma(i)] for all sigma, i.e. iff the parity vector
f_E(sigma) = inv_E(sigma) mod 2 lies in the GF(2) span of the indicators [sigma(i) = j] and the
constant. Because f_E is linear in E, the admissible E form a GF(2) subspace.
"""
from itertools import combinations, permutations

import numpy as np


def gf2_rank_basis_nullspace(M):
    """Return a basis (rows) of the left null space {y : y M = 0} over GF(2)."""
    A = (M.T % 2).astype(np.uint8)  # we need the null space of M^T
    rows, cols = A.shape
    A = A.copy()
    piv_cols = []
    r = 0
    for c in range(cols):
        k = next((i for i in range(r, rows) if A[i, c]), None)
        if k is None:
            continue
        A[[r, k]] = A[[k, r]]
        for i in range(rows):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        piv_cols.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in piv_cols]
    basis = []
    for f in free:
        y = np.zeros(cols, dtype=np.uint8)
        y[f] = 1
        for i, pc in enumerate(piv_cols):
            if A[i, f]:
                y[pc] = 1
        basis.append(y)
    return np.array(basis, dtype=np.uint8).reshape(len(basis), cols)


def gf2_kernel(B):
    """Kernel {z : B z = 0} over GF(2), returned as a list of 0/1 vectors (exhaustive span)."""
    N = gf2_rank_basis_nullspace(B.T)  # left null space of B^T = kernel of B
    k = N.shape[0]
    out = []
    for mask in range(1 << k):
        z = np.zeros(B.shape[1], dtype=np.uint8)
        for t in range(k):
            if mask >> t & 1:
                z ^= N[t]
        out.append(z)
    return out


def cut_space(n):
    pairs = list(combinations(range(n), 2))
    cuts = set()
    for mask in range(1 << n):
        S = {i for i in range(n) if mask >> i & 1}
        cuts.add(tuple(int((i in S) != (j in S)) for i, j in pairs))
    return cuts


def census(n):
    pairs = list(combinations(range(n), 2))
    perms = list(permutations(range(n)))
    M = np.zeros((len(perms), n * n + 1), dtype=np.uint8)
    F = np.zeros((len(perms), len(pairs)), dtype=np.uint8)
    for a, s in enumerate(perms):
        for i in range(n):
            M[a, i * n + s[i]] = 1
        M[a, -1] = 1
        for e, (i, j) in enumerate(pairs):
            F[a, e] = s[i] > s[j]
    N = gf2_rank_basis_nullspace(M)          # y with y M = 0
    B = (N.astype(np.int64) @ F.astype(np.int64)) % 2  # f_E in span iff B 1_E = 0
    ker = {tuple(int(v) for v in z) for z in gf2_kernel(B.astype(np.uint8))}
    cuts = cut_space(n)
    full = tuple([1] * len(pairs))
    perm_type = {tuple(1 - v for v in z) for z in ker}
    return {"n": n, "graphs": 2 ** len(pairs), "det_type": len(ker), "perm_type": len(perm_type),
            "neither": 2 ** len(pairs) - len(ker | perm_type),
            "det_type_equals_cuts": ker == cuts, "complete_graph_det_type": full in ker}


def D_E(A, E, n):
    tot = 0
    for s in permutations(range(n)):
        sign = 1
        for i, j in combinations(range(n), 2):
            if s[i] > s[j] and (i, j) not in E:
                sign = -sign
        t = sign
        for i in range(n):
            t *= A[i][s[i]]
        tot += t
    return tot


def clique_projection_check(n, T, extra_edges, rng):
    """Rows outside T fixed to unit rows; rows in T supported on T. Then D_E equals the permanent of
    the T-submatrix with explicit entry signs, whenever T is a clique of E."""
    E = set(combinations(sorted(T), 2)) | set(extra_edges)
    A = [[0] * n for _ in range(n)]
    for k in range(n):
        if k not in T:
            A[k][k] = 1
    for i in T:
        for j in T:
            A[i][j] = int(rng.integers(-5, 6))
    lhs = D_E(A, E, n)
    Tl = sorted(T)
    S = {}
    for i in Tl:
        for j in Tl:
            cnt = sum(1 for k in range(n) if k not in T and (min(i, k), max(i, k)) not in E
                      and (i < k < j or j < k < i))
            S[(i, j)] = -1 if cnt % 2 else 1
    rhs = 0
    for s in permutations(Tl):
        t = 1
        for i, j in zip(Tl, s):
            t *= S[(i, j)] * A[i][j]
        rhs += t
    return lhs == rhs
