"""Recorded run: every number in the paper comes from output/results.json written here."""
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import collective
import localsym
import signrule

OUT = Path(__file__).parent / "output"
OUT.mkdir(exist_ok=True)


def summarise_randomized(trials):
    s = {"cases": len(trials), "predicted_feasible": 0, "feasible_completed": 0, "obstructed": 0,
         "obstructed_completed": 0, "final_deficit_equals_q": 0, "det_audits": 0, "det_violations": 0,
         "readouts": 0, "readout_correct": 0, "by_condition": {}, "no_mobility": {}}
    max_c = (0, None)
    for t in trials:
        l = t["live"]
        c = s["by_condition"].setdefault(t["condition"], {"cases": 0, "q0": 0, "completed": 0, "contacts": []})
        c["cases"] += 1
        c["q0"] += l["q"] == 0
        c["completed"] += l["completed"]
        if l["completed"]:
            c["contacts"].append([l["n"], l["contacts"]])
            if l["contacts"] > max_c[0]:
                max_c = (l["contacts"], l["n"], t["budget"])
        if l["q"] == 0:
            s["predicted_feasible"] += 1
            s["feasible_completed"] += l["completed"]
        else:
            s["obstructed"] += 1
            s["obstructed_completed"] += l["completed"]
        s["final_deficit_equals_q"] += l["final_deficit"] == l["q"]
        s["det_audits"] += l["audits"]
        s["det_violations"] += l["det_violations"]
        if l["readout_ok"] is not None:
            s["readouts"] += 1
            s["readout_correct"] += bool(l["readout_ok"])
        if "no_mobility" in t:
            nm = t["no_mobility"]
            d = s["no_mobility"].setdefault(t["condition"], {"cases": 0, "completed": 0, "audits": 0, "det_violations": 0})
            d["cases"] += 1
            d["completed"] += nm["completed"]
            d["audits"] += nm["audits"]
            d["det_violations"] += nm["det_violations"]
    s["slowest_contacts"], s["slowest_n"], s["slowest_budget"] = max_c
    s["obstructed_q_min"] = min((t["live"]["q"] for t in trials if t["live"]["q"] > 0), default=None)
    s["obstructed_q_max"] = max((t["live"]["q"] for t in trials if t["live"]["q"] > 0), default=None)
    return s


def figure(res, path):
    plt.rcParams.update({"font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14,
                         "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 12})
    fig, ax = plt.subplots(2, 2, figsize=(12, 9.6))
    R = res["randomized"]
    conds = ["intact", "damaged_connected", "split_compatible", "split_dense"]
    labels = ["intact\ncomplete", "damaged,\nconnected", "split,\ncompatible", "split,\ndense"]
    live = [R["by_condition"][c]["completed"] / R["by_condition"][c]["cases"] for c in conds]
    nm = [R["no_mobility"][c]["completed"] / R["no_mobility"][c]["cases"] if c in R["no_mobility"] else np.nan for c in conds]
    pred = [R["by_condition"][c]["q0"] / R["by_condition"][c]["cases"] for c in conds]
    x = np.arange(4)
    a = ax[0, 0]
    a.bar(x - 0.27, pred, 0.27, color="#bbbbbb", label="predicted feasible ($q=0$)")
    a.bar(x, live, 0.27, color="#1f5fa8", label="completed, with swaps")
    a.bar(x + 0.27, nm, 0.27, color="#d08c2b", label="completed, no swaps")
    a.set_xticks(x, labels)
    a.set_ylabel("fraction of cases")
    a.set_ylim(0, 1.55)
    a.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    a.legend(loc="upper center", frameon=False)
    a.set_title("(a) Determinant agents over $\\mathbb{F}_{101}$")
    a = ax[0, 1]
    mk = {"intact": "o", "damaged_connected": "s", "split_compatible": "^"}
    col = {"intact": "#1f5fa8", "damaged_connected": "#c0392b", "split_compatible": "#2e8b57"}
    for c in ["intact", "damaged_connected", "split_compatible"]:
        pts = np.array(R["by_condition"][c]["contacts"], dtype=float)
        jitter = (np.random.default_rng(0).random(len(pts)) - 0.5) * 0.6
        a.scatter(pts[:, 0] + jitter, np.maximum(pts[:, 1], 1), marker=mk[c], s=36, alpha=0.75, color=col[c],
                  label=c.replace("_", ", "))
    a.set_yscale("log")
    a.set_xlabel("matrix size $n$")
    a.set_ylabel("contacts to readable form")
    a.legend(frameon=False, loc="upper left")
    a.set_title("(b) Time to completion")
    a = ax[1, 0]
    F = res["fibres"]
    C = res["certificate"]
    for p, c in ((3, "#1f5fa8"), (5, "#c0392b")):
        ks = sorted([k for k in F if F[k]["p"] == p], key=lambda k: F[k]["n"])
        ns = [F[k]["n"] for k in ks]
        fr = [F[k]["perm_fibre_equals_orbit_independent"] / F[k]["independent_pairs"] for k in ks]
        a.plot(ns, fr, "o-", color=c, lw=2, ms=8, label=f"frozen, $p={p}$")
        ce = [C[f"exhaustive_n{F[k]['n']}_p{p}"]["certified"] / F[k]["independent_pairs"]
              if f"exhaustive_n{F[k]['n']}_p{p}" in C else 0.0 for k in ks]
        a.plot(ns, ce, "s--", color=c, lw=1.5, ms=7, label=f"certified, $p={p}$")
    a.axhline(1.0, color="#555555", ls=":", lw=1.5, label="determinant (all pairs)")
    a.set_xticks([3, 4, 5, 6])
    a.set_ylim(-0.05, 1.75)
    a.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    a.set_xlabel("matrix size $n$")
    a.set_ylabel("share of independent pairs")
    a.set_title("(c) Permanent pairs with no extra moves")
    a.legend(frameon=False, loc="upper left", ncol=2, fontsize=11)
    a = ax[1, 1]
    G = res["greedy_paired"]
    keys = [("3", 5), ("4", 5), ("5", 5), ("5", 3), ("6", 3)]
    xl, dv, pv = [], [], []
    for n, p in keys:
        g = G[f"n{n}_p{p}"]
        xl.append(f"$n={n}$\n$p={p}$")
        dv.append(g["det"]["completed"] / g["trials"])
        pv.append(g["perm"]["completed"] / g["trials"])
    x = np.arange(len(keys))
    a.bar(x - 0.2, dv, 0.4, color="#1f5fa8", label="determinant-preserving")
    a.bar(x + 0.2, pv, 0.4, color="#d08c2b", label="permanent-preserving")
    a.set_xticks(x, xl)
    a.set_ylim(0, 1.3)
    a.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    a.set_ylabel("fraction reaching readable form")
    a.legend(frameon=False, loc="upper right")
    a.set_title("(d) Greedy agents, same matrices")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main():
    t0 = time.time()
    import pytest
    rc = int(pytest.main(["-q", str(Path(__file__).parent / "test_study.py")]))
    if rc != 0:
        raise SystemExit("unit tests failed")
    res = {"tests": {"exit_code": rc, "file": "test_study.py"},
           "environment": {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__,
                           "matplotlib": matplotlib.__version__, "platform": platform.platform()}}
    # 1. exhaustive protocol and orbit checks
    ex = {}
    for n, p in ((3, 2), (3, 3), (4, 2)):
        r = collective.exhaustive_protocol(n, p)
        if n == 4:
            r["per_graph"] = r["per_graph"][:0]
        ex[f"n{n}_p{p}"] = r
    ex["total_pairs"] = sum(ex[k]["pairs"] for k in ("n3_p2", "n3_p3", "n4_p2"))
    ex["total_mismatches"] = sum(ex[k]["reach_mismatch"] + ex[k]["as_mismatch"] + ex[k]["deficit_not_attained"]
                                 + ex[k]["profile_not_invariant"] + ex[k]["orbit_mismatch_graphs"]
                                 + ex[k]["det_violations"] + ex[k]["readout_errors"] for k in ("n3_p2", "n3_p3", "n4_p2"))
    ex["total_transitions"] = sum(ex[k]["transitions_checked"] for k in ("n3_p2", "n3_p3", "n4_p2"))
    res["exhaustive_protocol"] = ex
    # 2. randomised trials over F_101
    trials = collective.randomized_trials(seed=20261125, p=101, sizes=(4, 6, 8, 12, 16), per_cell=20, budget_factor=50)
    res["randomized"] = summarise_randomized(trials)
    (OUT / "trials.json").write_text(json.dumps(trials))
    # 3. linear pair-local preservers, exhaustive
    lp = {}
    for n, p in ((2, 3), (3, 2), (3, 3), (3, 5), (4, 2), (4, 3)):
        for kind in ("det", "perm"):
            lp[f"{kind}_n{n}_p{p}"] = localsym.classify_preservers(n, p, kind)
    res["linear_preservers"] = lp
    # 4. fibres of the pair maps (state-dependent moves)
    fib, frozen = {}, {}
    for n, p in ((3, 3), (3, 5), (4, 3), (4, 5), (5, 3), (5, 5), (6, 3)):
        r, lk = localsym.fibre_census(n, p)
        fib[f"n{n}_p{p}"] = r
        if n == 5:
            frozen[f"n{n}_p{p}"] = localsym.frozen_fraction(n, p, lk, 20000, seed=7)
    res["fibres"] = fib
    res["frozen"] = frozen
    # 5. exhaustive reachability with state-dependent moves at n = 3
    reach = {}
    for p in (3, 5):
        for name, E in (("complete", ((0, 1), (0, 2), (1, 2))), ("path", ((0, 1), (1, 2)))):
            reach[f"p{p}_{name}"] = localsym.state_dependent_reachability(p, E)
    res["reach_n3"] = reach
    # 6. greedy agents with oracle access to every value-preserving pair move, paired design
    gr = {}
    for n, p in ((3, 5), (4, 5), (5, 5), (5, 3), (6, 3)):
        gr[f"n{n}_p{p}"] = localsym.greedy_paired(n, p, trials=200, budget=3000, seed=11 + n + p)
    res["greedy_paired"] = gr
    # 6b. rank-two rigidity certificate: exhaustive soundness check and seeded large-field shares
    cert = {}
    for n, p in ((5, 3), (5, 5), (6, 3)):
        cert[f"exhaustive_n{n}_p{p}"] = localsym.certificate_exhaustive(n, p)
    cert["sampled_n5_p101"] = localsym.certificate_sampled(5, 101, 200000, seed=3)
    cert["sampled_n6_p1009"] = localsym.certificate_sampled(6, 1009, 200000, seed=4)
    cert["sampled_n5_p1009"] = localsym.certificate_sampled(5, 1009, 100000, seed=5, matrices=20000)
    res["certificate"] = cert
    # 7. sign-rule census
    rng = np.random.default_rng(5)
    res["sign_census"] = {f"n{n}": signrule.census(n) for n in range(3, 8)}
    checks = [signrule.clique_projection_check(6, {0, 2, 3, 5}, [(1, 4), (0, 1)], rng) for _ in range(10)]
    checks += [signrule.clique_projection_check(7, {1, 2, 4, 6}, [(0, 5), (3, 4), (2, 5)], rng) for _ in range(10)]
    res["clique_projection"] = {"checks": len(checks), "passed": int(sum(checks))}
    res["runtime_seconds"] = round(time.time() - t0, 1)
    (OUT / "results.json").write_text(json.dumps(res, indent=1, default=int))
    figure(res, OUT / "figure_main.png")
    paper_fig = Path(__file__).parent.parent / "paper" / "figure_main.png"
    figure(res, paper_fig)
    print(json.dumps({k: res[k] for k in ("randomized",)}, indent=1, default=int)[:3000])


if __name__ == "__main__":
    main()
