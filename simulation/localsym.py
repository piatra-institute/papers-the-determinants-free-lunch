"""Pair-local symmetries of det and perm over F_p.

A pair of rows (x, y) in F_p^n x F_p^n is stored as one vector of length 2n.
By Laplace expansion along the two rows, a map of the pair preserves det (perm) of every
n x n matrix containing it exactly when it preserves all 2 x 2 minors
x_s y_t - x_t y_s (all 2 x 2 permanents x_s y_t + x_t y_s), s < t.
"""
from itertools import combinations, product

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def pair_forms(V, n, p, kind):
    """V: (..., 2n) int array. Returns (..., n(n-1)/2) array of 2x2 minors or permanents mod p."""
    X, Y = V[..., :n], V[..., n:]
    sgn = -1 if kind == "det" else 1
    cols = [(X[..., s] * Y[..., t] + sgn * X[..., t] * Y[..., s]) % p for s, t in combinations(range(n), 2)]
    return np.stack(cols, axis=-1)


def polar(U, W, n, p, kind):
    return (pair_forms(U + W, n, p, kind) - pair_forms(U, n, p, kind) - pair_forms(W, n, p, kind)) % p


def all_vectors(m, p):
    return np.array(list(product(range(p), repeat=m)), dtype=np.int64)


def enumerate_linear_preservers(n, p, kind):
    """All linear maps Phi of F_p^{2n} (images of the 2n basis vectors) that preserve the pair
    forms. A quadratic map is preserved iff its values on the basis and its polar form on basis
    pairs are preserved, so this backtracking search is exhaustive and exact."""
    m = 2 * n
    allv = all_vectors(m, p)
    cand = allv[(pair_forms(allv, n, p, kind) == 0).all(axis=1)]  # every basis vector has zero forms
    basis = np.eye(m, dtype=np.int64)
    order = []
    for s in range(n):
        order += [s, n + s]  # interleave x_s, y_s for early pruning
    target = {(a, b): polar(basis[a], basis[b], n, p, kind) for a in range(m) for b in range(m)}
    sols = []

    def rec(k, chosen):
        if k == m:
            img = np.zeros((m, m), dtype=np.int64)
            for a, v in chosen.items():
                img[:, a] = v
            sols.append(img)
            return
        a = order[k]
        mask = np.ones(len(cand), dtype=bool)
        for b, v in chosen.items():
            pol = polar(cand, v[None, :], n, p, kind)
            mask &= (pol == target[(a, b)][None, :]).all(axis=1)
        for w in cand[mask]:
            chosen[a] = w
            rec(k + 1, chosen)
            del chosen[a]

    rec(0, {})
    return sols


def is_scalar_pair_map(M, n, p):
    """True iff Phi(x, y) = (a x + b y, c x + d y) for scalars a, b, c, d."""
    a, b, c, d = M[0, 0], M[0, n], M[n, 0], M[n, n]
    G = np.zeros_like(M)
    I = np.eye(n, dtype=np.int64)
    G[:n, :n], G[:n, n:], G[n:, :n], G[n:, n:] = a * I, b * I, c * I, d * I
    return bool(((M - G) % p == 0).all()), (int(a), int(b), int(c), int(d))


def classify_preservers(n, p, kind):
    sols = enumerate_linear_preservers(n, p, kind)
    scalar = 0
    det_one = 0
    monomial = 0
    for M in sols:
        ok, (a, b, c, d) = is_scalar_pair_map(M, n, p)
        scalar += ok
        if ok:
            det_one += (a * d - b * c) % p == 1
            monomial += ((b == 0 and c == 0 and (a * d) % p == 1) or (a == 0 and d == 0 and (b * c) % p == 1))
    return {"n": n, "p": p, "kind": kind, "count": len(sols), "scalar_pair_maps": scalar,
            "det_one": int(det_one), "monomial_perm_one": int(monomial),
            "sl2_order": p * (p * p - 1), "monomial_order": 2 * (p - 1)}


# ---------------------------------------------------------------- state-dependent moves: fibres

def keys_of(F, p):
    w = p ** np.arange(F.shape[-1], dtype=np.int64)
    return (F * w).sum(axis=-1)


def fibre_census(n, p):
    """For every pair (x, y), the number of pairs with the same 2x2 permanents (fibre of the
    permanental pair map) compared with the size of its monomial orbit, and the same for det."""
    V = all_vectors(2 * n, p)
    X, Y = V[:, :n], V[:, n:]
    wedge = keys_of(pair_forms(V, n, p, "det"), p)
    indep = wedge != 0
    out = {"n": n, "p": p, "pairs": int(len(V)), "independent_pairs": int(indep.sum())}
    _, inv_d, cnt_d = np.unique(wedge, return_inverse=True, return_counts=True)
    fd = cnt_d[inv_d][indep]
    out["det_fibre_min"], out["det_fibre_max"] = int(fd.min()), int(fd.max())
    pk = keys_of(pair_forms(V, n, p, "perm"), p)
    _, inv_p, cnt_p = np.unique(pk, return_inverse=True, return_counts=True)
    fp = cnt_p[inv_p]
    orbit = 2 * (p - 1)
    dense = (X != 0).all(axis=1) & (Y != 0).all(axis=1)
    out["perm_orbit_size"] = orbit
    out["perm_fibre_equals_orbit_independent"] = int((fp[indep] == orbit).sum())
    out["perm_fibre_equals_orbit_dense"] = int((fp[dense & indep] == orbit).sum())
    out["dense_independent_pairs"] = int((dense & indep).sum())
    out["perm_fibre_min_independent"] = int(fp[indep].min())
    out["perm_fibre_max_independent"] = int(fp[indep].max())
    out["perm_fibre_mean_independent"] = float(fp[indep].mean())
    # can a permanent-preserving pair move change the number of zero entries?
    zeros = (V == 0).sum(axis=1)
    zmin = np.full(len(cnt_p), 10 ** 9)
    zmax = np.full(len(cnt_p), -1)
    np.minimum.at(zmin, inv_p, zeros)
    np.maximum.at(zmax, inv_p, zeros)
    changes = (zmin[inv_p] != zmax[inv_p])
    out["perm_zero_count_changeable_independent"] = int(changes[indep].sum())
    lookup = {"perm_key": pk, "perm_fibre": fp, "indep": indep}
    return out, lookup


def frozen_fraction(n, p, lookup, samples, seed):
    """Fraction of uniformly random n x n matrices whose every row pair has a permanental fibre
    equal to its monomial orbit (no permanent-preserving pair move beyond scaling and swapping)."""
    rng = np.random.default_rng(seed)
    A = rng.integers(0, p, size=(samples, n, n))
    w = p ** np.arange(2 * n - 1, -1, -1, dtype=np.int64)  # big-endian, matching all_vectors
    frozen = np.ones(samples, dtype=bool)
    allindep = np.ones(samples, dtype=bool)
    for i, j in combinations(range(n), 2):
        V = np.concatenate([A[:, i, :], A[:, j, :]], axis=1)
        idx = (V * w).sum(axis=1)
        frozen &= lookup["perm_fibre"][idx] == 2 * (p - 1)
        allindep &= lookup["indep"][idx]
    return {"n": n, "p": p, "samples": samples, "seed": seed, "frozen": int(frozen.sum()),
            "all_pairs_independent": int(allindep.sum()), "frozen_and_independent": int((frozen & allindep).sum())}


# ---------------------------------------------------------------- reachability with state-dependent moves, n = 3

def leads_np(A):
    """A: (N, n, n). Lead column per row, n for a zero row."""
    nz = A != 0
    n = A.shape[-1]
    first = np.where(nz.any(axis=-1), nz.argmax(axis=-1), n)
    return first


def readable_mask(A):
    L = leads_np(A)
    n = A.shape[-1]
    ok = (L < n).all(axis=1)
    Ls = np.sort(L, axis=1)
    ok &= (np.diff(Ls, axis=1) != 0).all(axis=1)
    return ok


def perm3(A, p):
    return (A[:, 0, 0] * (A[:, 1, 1] * A[:, 2, 2] + A[:, 1, 2] * A[:, 2, 1])
            + A[:, 0, 1] * (A[:, 1, 0] * A[:, 2, 2] + A[:, 1, 2] * A[:, 2, 0])
            + A[:, 0, 2] * (A[:, 1, 0] * A[:, 2, 1] + A[:, 1, 1] * A[:, 2, 0])) % p


def det3(A, p):
    return (A[:, 0, 0] * (A[:, 1, 1] * A[:, 2, 2] - A[:, 1, 2] * A[:, 2, 1])
            - A[:, 0, 1] * (A[:, 1, 0] * A[:, 2, 2] - A[:, 1, 2] * A[:, 2, 0])
            + A[:, 0, 2] * (A[:, 1, 0] * A[:, 2, 1] - A[:, 1, 1] * A[:, 2, 0])) % p


def state_dependent_reachability(p, edges=((0, 1), (0, 2), (1, 2))):
    """All 3 x 3 matrices over F_p. A move on edge {i, j} replaces rows (i, j) by ANY pair that
    preserves det (resp. perm) of every matrix containing it. Report how many matrices with
    nonzero det (perm) can reach a readable (distinct-lead) form."""
    n = 3
    V = all_vectors(n * n, p)
    A = V.reshape(-1, n, n)
    N = len(A)
    rd = readable_mask(A)
    res = {"p": p, "edges": [list(e) for e in edges], "matrices": int(N)}
    for kind, val in (("det", det3(A, p)), ("perm", perm3(A, p))):
        rows_src, rows_dst = [], []
        offset = N
        for (i, j) in edges:
            k = 3 - i - j
            pair = np.concatenate([A[:, i, :], A[:, j, :]], axis=1)
            fk = keys_of(pair_forms(pair, n, p, kind), p)
            ok = keys_of(A[:, k, :], p)
            key = fk * (p ** n) + ok
            _, g = np.unique(key, return_inverse=True)
            rows_src.append(np.arange(N))
            rows_dst.append(offset + g)
            offset += g.max() + 1
        src = np.concatenate(rows_src)
        dst = np.concatenate(rows_dst)
        G = coo_matrix((np.ones(len(src), dtype=np.int8), (src, dst)), shape=(offset, offset))
        _, lab = connected_components(G, directed=False)
        lab = lab[:N]
        comp_has_readable = np.zeros(lab.max() + 1, dtype=bool)
        comp_has_readable[lab[rd]] = True
        nonzero = val != 0
        res[kind] = {"nonzero_value": int(nonzero.sum()),
                     "already_readable": int((rd & nonzero).sum()),
                     "can_reach_readable": int((comp_has_readable[lab] & nonzero).sum()),
                     "components_nonzero": int(len(np.unique(lab[nonzero])))}
    return res


# ---------------------------------------------------------------- greedy agents with oracle access to all pair moves

class FibreOracle:
    """All pairs over F_p^n x F_p^n grouped by their pair-form vector (det minors or perm minors)."""

    def __init__(self, n, p, kind):
        self.n, self.p, self.kind = n, p, kind
        V = all_vectors(2 * n, p).astype(np.int8)
        key = keys_of(pair_forms(V.astype(np.int64), n, p, kind), p)
        order = np.argsort(key, kind="stable")
        self.sorted_key = key[order]
        self.order = order
        self.key = key
        self.V = V
        L = leads_np(V.reshape(-1, 2, n))
        self.leadsum = L.sum(axis=1)
        self.w = p ** np.arange(2 * n - 1, -1, -1, dtype=np.int64)

    def best_moves(self, x, y):
        idx = int((np.concatenate([x, y]) * self.w).sum())
        k = self.key[idx]
        lo = np.searchsorted(self.sorted_key, k, side="left")
        hi = np.searchsorted(self.sorted_key, k, side="right")
        members = self.order[lo:hi]
        ls = self.leadsum[members]
        return members[ls == ls.max()], hi - lo


def greedy_trials(n, p, kind, trials, budget, seed, oracle=None):
    """Agents on the complete graph. On each contact the pair moves to a uniformly random pair,
    among ALL pairs preserving det (perm) of every completion, that maximises the pair's lead sum."""
    orc = oracle or FibreOracle(n, p, kind)
    rng = np.random.default_rng(seed)
    pairs = list(combinations(range(n), 2))
    done, contacts, started = 0, [], 0
    value = (lambda A: int(round(np.linalg.det(A))) % p) if False else None
    from fieldlin import det_gauss, perm_leibniz
    while started < trials:
        A = rng.integers(0, p, size=(n, n))
        v0 = det_gauss(A.tolist(), p) if kind == "det" else perm_leibniz(A.tolist(), p)
        if v0 == 0:
            continue
        started += 1
        t_done = None
        if readable_mask(A[None])[0]:
            t_done = 0
        for t in range(1, budget + 1):
            if t_done is not None:
                break
            i, j = pairs[rng.integers(len(pairs))]
            best, _ = orc.best_moves(A[i], A[j])
            m = orc.V[best[rng.integers(len(best))]].astype(np.int64)
            A[i], A[j] = m[:n], m[n:]
            if readable_mask(A[None])[0]:
                t_done = t
        v1 = det_gauss(A.tolist(), p) if kind == "det" else perm_leibniz(A.tolist(), p)
        assert v1 == v0, "value not preserved"
        if t_done is not None:
            done += 1
            contacts.append(t_done)
    return {"n": n, "p": p, "kind": kind, "trials": trials, "budget": budget, "seed": seed,
            "completed": done, "median_contacts": float(np.median(contacts)) if contacts else None}
