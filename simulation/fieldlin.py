"""Exact linear algebra over a prime field F_p (Python integers only)."""
from itertools import combinations, permutations


def lead(row):
    """Index of the first nonzero entry, or None for a zero row."""
    for k, v in enumerate(row):
        if v:
            return k
    return None


def inv(a, p):
    return pow(a % p, p - 2, p)


def echelon(rows, p):
    """Row-reduce a list of rows over F_p; return (reduced rows, pivot columns)."""
    m = [list(r) for r in rows]
    piv = []
    r = 0
    ncol = len(m[0]) if m else 0
    for c in range(ncol):
        k = next((i for i in range(r, len(m)) if m[i][c] % p), None)
        if k is None:
            continue
        m[r], m[k] = m[k], m[r]
        f = inv(m[r][c], p)
        m[r] = [(v * f) % p for v in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] % p:
                g = m[i][c]
                m[i] = [(a - g * b) % p for a, b in zip(m[i], m[r])]
        piv.append(c)
        r += 1
        if r == len(m):
            break
    return m[:r], piv


def rank(rows, p):
    return len(echelon(rows, p)[1])


def pivot_set(rows, p):
    """Pivot columns of the row space spanned by rows: the possible leads of its nonzero vectors."""
    return frozenset(echelon(rows, p)[1])


def det_gauss(rows, p):
    """Determinant over F_p by Gaussian elimination with sign tracking."""
    m = [list(r) for r in rows]
    n = len(m)
    d = 1
    for c in range(n):
        k = next((i for i in range(c, n) if m[i][c] % p), None)
        if k is None:
            return 0
        if k != c:
            m[c], m[k] = m[k], m[c]
            d = -d
        d = d * m[c][c] % p
        f = inv(m[c][c], p)
        for i in range(c + 1, n):
            if m[i][c] % p:
                g = m[i][c] * f % p
                m[i] = [(a - g * b) % p for a, b in zip(m[i], m[c])]
    return d % p


def perm_sign(seq):
    s = 1
    seq = list(seq)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def det_leibniz(rows, p):
    """Independent determinant by permutation expansion (small n only)."""
    n = len(rows)
    tot = 0
    for s in permutations(range(n)):
        t = perm_sign(s)
        for i in range(n):
            t *= rows[i][s[i]]
        tot += t
    return tot % p


def perm_leibniz(rows, p=None):
    n = len(rows)
    tot = 0
    for s in permutations(range(n)):
        t = 1
        for i in range(n):
            t *= rows[i][s[i]]
        tot += t
    return tot % p if p else tot


def plucker(rows, p):
    """All maximal minors of a k x n block (rows in the given order), as a tuple."""
    k = len(rows)
    n = len(rows[0])
    return tuple(det_gauss([[r[c] for c in cols] for r in rows], p) for cols in combinations(range(n), k))


def components(n, edges):
    adj = {i: set() for i in range(n)}
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    seen, comps = set(), []
    for s in range(n):
        if s in seen:
            continue
        stack, comp = [s], []
        seen.add(s)
        while stack:
            u = stack.pop()
            comp.append(u)
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        comps.append(sorted(comp))
    return comps


def deficit_q(rows, edges, p):
    """q = n - |union of component pivot profiles| (the obstruction index)."""
    n = len(rows)
    cover = set()
    for comp in components(n, edges):
        cover |= pivot_set([rows[i] for i in comp], p)
    return n - len(cover)


def readout(rows, p):
    """Determinant from (lead, pivot value) summaries, valid when all leads are distinct."""
    n = len(rows)
    leads = [lead(r) for r in rows]
    if None in leads or len(set(leads)) != n:
        return None
    order = sorted(range(n), key=lambda i: leads[i])
    d = perm_sign(order)
    for i in order:
        d = d * rows[i][leads[i]]
    return d % p
