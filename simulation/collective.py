"""Row-holding agents that preserve the determinant with pair-local SL_2 moves.

Agent i holds row i. On a contact along edge {i, j} (i < j) the pair sees only rows i and j:
  * equal leads c: r_j <- r_j - (r_j[c] / r_i[c]) r_i          (a transvection, determinant 1)
  * different leads, mobility on, coin true: (r_i, r_j) <- (r_j, -r_i)   (signed swap, determinant 1)
  * otherwise: no change.
"""
import random
from collections import deque
from itertools import combinations, product

from fieldlin import (components, deficit_q, det_gauss, det_leibniz, inv, lead,
                      pivot_set, plucker, rank, readout)


def contact(rows, i, j, p, mobility=True, coin=True):
    li, lj = lead(rows[i]), lead(rows[j])
    if li is not None and li == lj:
        c = rows[j][li] * inv(rows[i][li], p) % p
        rows[j] = [(b - c * a) % p for a, b in zip(rows[i], rows[j])]
        return "elim"
    if mobility and coin:
        rows[i], rows[j] = rows[j], [(-a) % p for a in rows[i]]
        return "swap"
    return "noop"


def distinct_leads(rows):
    ls = [lead(r) for r in rows]
    return len(set(ls)) - (1 if None in ls else 0)


# ---------------------------------------------------------------- exhaustive small cases

def all_invertible(n, p):
    vecs = list(product(range(p), repeat=n))
    out = []
    for rows in product(vecs, repeat=n):
        if rank(rows, p) == n:
            out.append(tuple(rows))
    return out


def all_graphs(n):
    pairs = list(combinations(range(n), 2))
    for mask in range(1 << len(pairs)):
        yield tuple(e for k, e in enumerate(pairs) if mask >> k & 1)


def reverse_reach(targets, preds):
    seen = set(targets)
    dq = deque(targets)
    while dq:
        u = dq.popleft()
        for v in preds.get(u, ()):
            if v not in seen:
                seen.add(v)
                dq.append(v)
    return seen


def exhaustive_protocol(n, p):
    """For every invertible n x n matrix over F_p and every labelled graph on n agents,
    compare the predicted obstruction index q with exact reachability under the protocol
    and with orbits of the full pair-local SL_2 action."""
    states = all_invertible(n, p)
    index = {s: k for k, s in enumerate(states)}
    N = len(states)
    pairs = list(combinations(range(n), 2))
    move = {}  # edge -> list of (kind, successor index)
    gen_by_edge = {}
    for (i, j) in pairs:
        mv, g1, g2 = [], [], []
        for s in states:
            rows = [list(r) for r in s]
            kind = contact(rows, i, j, p, mobility=True, coin=True)
            mv.append((kind, index[tuple(tuple(r) for r in rows)]))
            a = [list(r) for r in s]
            a[i] = [(x + y) % p for x, y in zip(s[i], s[j])]
            b = [list(r) for r in s]
            b[j] = [(x + y) % p for x, y in zip(s[j], s[i])]
            g1.append(index[tuple(tuple(r) for r in a)])
            g2.append(index[tuple(tuple(r) for r in b)])
        move[(i, j)] = mv
        gen_by_edge[(i, j)] = (g1, g2)
    nd = [distinct_leads(s) for s in states]
    readable = [k for k in range(N) if nd[k] == n]
    dets = [det_gauss(s, p) for s in states]
    piv_cache, pl_cache = {}, {}

    def piv(block):
        if block not in piv_cache:
            piv_cache[block] = pivot_set(block, p)
        return piv_cache[block]

    def pl(block):
        if block not in pl_cache:
            pl_cache[block] = plucker(block, p)
        return pl_cache[block]

    tot = {"pairs": 0, "q0": 0, "reach_mismatch": 0, "as_mismatch": 0, "deficit_not_attained": 0,
           "profile_not_invariant": 0, "orbit_mismatch_graphs": 0, "det_violations": 0,
           "transitions_checked": 0, "readout_errors": 0, "no_mobility_success": 0}
    per_graph = []
    for G in all_graphs(n):
        comps = components(n, G)
        prof = []
        q = []
        for s in states:
            ps = tuple(piv(tuple(s[i] for i in c)) for c in comps)
            prof.append(ps)
            q.append(n - len(frozenset().union(*ps)))
        preds, preds_nm = {}, {}
        for e in G:
            for k, (kind, t) in enumerate(move[e]):
                if t == k:
                    continue
                tot["transitions_checked"] += 1
                preds.setdefault(t, []).append(k)
                if kind == "elim":
                    preds_nm.setdefault(t, []).append(k)
                tot["det_violations"] += dets[t] != dets[k]
                tot["profile_not_invariant"] += prof[t] != prof[k]
        can = reverse_reach(readable, preds)
        reach_bad = reverse_reach([k for k in range(N) if k not in can], preds)
        reach_attain = reverse_reach([k for k in range(N) if n - nd[k] == q[k]], preds)
        can_nm = reverse_reach(readable, preds_nm)
        for k in range(N):
            pred = q[k] == 0
            tot["pairs"] += 1
            tot["q0"] += pred
            tot["reach_mismatch"] += (k in can) != pred
            tot["as_mismatch"] += (k not in reach_bad) != pred
            tot["deficit_not_attained"] += k not in reach_attain
        tot["no_mobility_success"] += len(can_nm)
        parent = list(range(N))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for e in G:
            for g in gen_by_edge[e]:
                for k, t in enumerate(g):
                    a, b = find(k), find(t)
                    if a != b:
                        parent[a] = b
        inv_key = {}
        for k, s in enumerate(states):
            key = tuple(pl(tuple(s[i] for i in c)) for c in comps)
            inv_key.setdefault(key, set()).add(find(k))
        roots = {find(k) for k in range(N)}
        orbit_ok = all(len(v) == 1 for v in inv_key.values()) and len(inv_key) == len(roots)
        tot["orbit_mismatch_graphs"] += not orbit_ok
        per_graph.append({"edges": [list(e) for e in G], "q0": sum(1 for x in q if x == 0),
                          "reachable": len(can), "almost_sure": N - len(reach_bad),
                          "no_mobility": len(can_nm), "orbits": len(roots), "invariant_classes": len(inv_key)})
    for k in readable:
        if readout(states[k], p) != dets[k] or det_leibniz(states[k], p) != dets[k]:
            tot["readout_errors"] += 1
    return {"n": n, "p": p, "states": N, "graphs": len(per_graph), **tot, "per_graph": per_graph}


# ---------------------------------------------------------------- randomised trials over F_101

def random_invertible(n, p, rng):
    while True:
        rows = [[rng.randrange(p) for _ in range(n)] for _ in range(n)]
        if rank(rows, p) == n:
            return rows


def connected(n, edges):
    return len(components(n, edges)) == 1


def ring_plus_chords(vertices, rng, chords):
    v = list(vertices)
    E = set()
    if len(v) > 1:
        for a, b in zip(v, v[1:] + v[:1]):
            if a != b:
                E.add((min(a, b), max(a, b)))
    for _ in range(chords):
        a, b = rng.sample(v, 2) if len(v) > 1 else (v[0], v[0])
        if a != b:
            E.add((min(a, b), max(a, b)))
    return sorted(E)


def block_rows(n, comp_sets, p, rng):
    """Rows for which each component's row space has a prescribed, disjoint pivot profile."""
    cols = list(range(n))
    rng.shuffle(cols)
    rows = [None] * n
    pos = 0
    profiles = []
    for comp in comp_sets:
        S = sorted(cols[pos:pos + len(comp)])
        pos += len(comp)
        profiles.append(S)
        basis = []
        for c in S:
            r = [0] * n
            r[c] = rng.randrange(1, p)
            for d in range(c + 1, n):
                r[d] = rng.randrange(p)
            basis.append(r)
        k = len(comp)
        while True:  # random invertible mixing inside the component hides the profile
            M = [[rng.randrange(p) for _ in range(k)] for _ in range(k)]
            if rank(M, p) == k:
                break
        for a, i in enumerate(comp):
            rows[i] = [sum(M[a][b] * basis[b][d] for b in range(k)) % p for d in range(n)]
    return rows, profiles


def make_case(condition, n, p, rng):
    verts = list(range(n))
    if condition == "intact":
        edges = list(combinations(verts, 2))
        rows = random_invertible(n, p, rng)
    elif condition == "damaged_connected":
        while True:
            base = ring_plus_chords(verts, rng, n)
            kept = [e for e in base if rng.random() > 0.3]
            if connected(n, kept):
                edges = kept
                break
        rows = random_invertible(n, p, rng)
    else:
        k = rng.randrange(2, n - 1)
        perm = verts[:]
        rng.shuffle(perm)
        comp_sets = [sorted(perm[:k]), sorted(perm[k:])]
        edges = sorted(set(ring_plus_chords(comp_sets[0], rng, 1)) | set(ring_plus_chords(comp_sets[1], rng, 1)))
        if condition == "split_compatible":
            rows, _ = block_rows(n, comp_sets, p, rng)
        else:
            rows = random_invertible(n, p, rng)
    return rows, edges


def run_trial(rows, edges, p, rng, budget, mobility=True, audit_every=500):
    rows = [list(r) for r in rows]
    n = len(rows)
    d0 = det_gauss(rows, p)
    q = deficit_q(rows, edges, p)
    audits = 0
    violations = 0
    t_done = None
    elim = swaps = 0
    if distinct_leads(rows) == n:
        t_done = 0
    for t in range(1, budget + 1):
        if t_done is not None:
            break
        i, j = rng.choice(edges)
        kind = contact(rows, i, j, p, mobility=mobility, coin=rng.random() < 0.5)
        elim += kind == "elim"
        swaps += kind == "swap"
        if t % audit_every == 0:
            audits += 1
            violations += det_gauss(rows, p) != d0
        if kind != "noop" and distinct_leads(rows) == n:
            t_done = t
    audits += 1
    violations += det_gauss(rows, p) != d0
    ro = readout(rows, p)
    return {"n": n, "q": q, "completed": t_done is not None, "contacts": t_done, "eliminations": elim,
            "swaps": swaps, "final_deficit": n - distinct_leads(rows), "audits": audits,
            "det_violations": violations, "readout_ok": (ro == d0) if ro is not None else None}


def randomized_trials(seed=20261125, p=101, sizes=(4, 6, 8, 10, 12), per_cell=10, budget_factor=50):
    rng = random.Random(seed)
    conditions = ["intact", "damaged_connected", "split_compatible", "split_dense"]
    trials = []
    for cond in conditions:
        for n in sizes:
            for rep in range(per_cell):
                rows, edges = make_case(cond, n, p, rng)
                budget = budget_factor * n ** 3
                trial_seed = rng.randrange(2 ** 31)
                live = run_trial(rows, edges, p, random.Random(trial_seed), budget)
                rec = {"condition": cond, "rep": rep, "edges": len(edges), "budget": budget, "live": live}
                if cond in ("intact", "damaged_connected"):
                    rec["no_mobility"] = run_trial(rows, edges, p, random.Random(trial_seed), budget, mobility=False)
                trials.append(rec)
    return trials
