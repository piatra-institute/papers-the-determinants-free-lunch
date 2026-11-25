"""Unit tests: run with `uv run python -m pytest -q`."""
import random
from itertools import combinations

import numpy as np

import collective
import localsym
import signrule
from fieldlin import det_gauss, det_leibniz, perm_leibniz, pivot_set, readout


def test_moves_have_determinant_one():
    p = 101
    rng = random.Random(1)
    for _ in range(200):
        rows = collective.random_invertible(5, p, rng)
        d0 = det_gauss(rows, p)
        i, j = sorted(rng.sample(range(5), 2))
        collective.contact(rows, i, j, p, coin=True)
        assert det_gauss(rows, p) == d0


def test_unsigned_swap_flips_sign():
    p = 7
    rows = [[1, 2, 0], [0, 1, 3], [2, 0, 1]]
    d = det_gauss(rows, p)
    rows[0], rows[1] = rows[1], rows[0]
    assert det_gauss(rows, p) == (-d) % p


def test_determinants_agree_and_readout():
    p = 5
    rng = random.Random(2)
    for _ in range(50):
        rows = [[rng.randrange(p) for _ in range(4)] for _ in range(4)]
        assert det_gauss(rows, p) == det_leibniz(rows, p)
    upper = [[0, 2, 1, 4], [3, 1, 0, 0], [0, 0, 0, 2], [0, 0, 4, 1]]
    assert readout(upper, p) == det_gauss(upper, p)


def test_pivot_profile_is_set_of_leads():
    p = 3
    rows = [(1, 0, 2, 0), (1, 1, 0, 0)]
    S = pivot_set(rows, p)
    leads = set()
    for a in range(p):
        for b in range(p):
            v = [(a * x + b * y) % p for x, y in zip(*rows)]
            nz = [k for k, t in enumerate(v) if t]
            if nz:
                leads.add(nz[0])
    assert S == leads


def test_block_rows_have_disjoint_profiles():
    rng = random.Random(4)
    rows, profiles = collective.block_rows(6, [[0, 2, 4], [1, 3, 5]], 101, rng)
    from fieldlin import deficit_q
    edges = [(0, 2), (2, 4), (1, 3), (3, 5)]
    assert deficit_q(rows, edges, 101) == 0


def test_linear_preservers_small():
    r = localsym.classify_preservers(3, 3, "perm")
    assert r["count"] == 4 == r["monomial_perm_one"]
    r = localsym.classify_preservers(3, 3, "det")
    assert r["count"] == 24 == r["det_one"]


def test_pair_forms_laplace():
    rng = np.random.default_rng(0)
    n, p = 4, 7
    for _ in range(20):
        A = rng.integers(0, p, size=(n, n))
        for kind in ("det", "perm"):
            forms = localsym.pair_forms(np.concatenate([A[0], A[1]])[None, :], n, p, kind)[0]
            tot = 0
            for k, (s, t) in enumerate(combinations(range(n), 2)):
                rest = [c for c in range(n) if c not in (s, t)]
                sub = A[2:, rest].tolist()
                if kind == "det":
                    sign = (-1) ** (s + t + 1)  # rows 0,1 and columns s,t (0-based)
                    tot += sign * int(forms[k]) * det_leibniz(sub, p)
                else:
                    tot += int(forms[k]) * perm_leibniz(sub, p)
            ref = det_leibniz(A.tolist(), p) if kind == "det" else perm_leibniz(A.tolist(), p)
            assert tot % p == ref


def test_cut_rules_are_signed_determinants():
    n = 4
    rng = np.random.default_rng(1)
    S = {0, 2}
    E = {(i, j) for i, j in combinations(range(n), 2) if (i in S) != (j in S)}
    for _ in range(5):
        A = rng.integers(-4, 5, size=(n, n)).tolist()
        B = [[A[i][j] * ((-1) ** (i + j) if i in S else 1) for j in range(n)] for i in range(n)]
        lhs = signrule.D_E(A, E, n)
        d = det_leibniz(B, 10 ** 9 + 7)
        assert lhs % (10 ** 9 + 7) in (d, (-d) % (10 ** 9 + 7))


def test_census_n4():
    r = signrule.census(4)
    assert r["det_type_equals_cuts"] and r["det_type"] == 8
