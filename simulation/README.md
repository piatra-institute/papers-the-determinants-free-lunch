# Simulation for "The Determinant's Free Lunch"

Exact finite-field computations behind every number in the paper. All arithmetic is over a prime field with Python integers or NumPy integer arrays; every random process is seeded.

- `fieldlin.py`: echelon forms, pivot profiles, determinants (Gaussian elimination and permutation expansion), the readout from (lead, pivot) summaries.
- `collective.py`: row-holding agents with determinant-one pair moves; exhaustive reachability and orbit checks for every invertible matrix and every labelled communication graph at (n, p) = (3, 2), (3, 3), (4, 2); seeded randomised trials over F_101.
- `localsym.py`: exhaustive enumeration of the linear pair-local maps that preserve det or perm for every completion; fibre census of the pair maps; exhaustive n = 3 reachability with state-dependent moves; greedy agents with oracle access to every value-preserving pair move.
- `signrule.py`: GF(2) census of sign-rule faults in the Leibniz expansion and the clique projection check.
- `test_study.py`: unit tests, executed at the start of `run_all.py`.

Run from this directory:

```bash
uv lock
uv run python run_all.py
```

Outputs: `output/results.json` (all reported numbers), `output/trials.json` (per-trial records of the randomised experiment), `output/figure_main.png` (also written to `../paper/figure_main.png`). The full run takes about two minutes on a laptop.
