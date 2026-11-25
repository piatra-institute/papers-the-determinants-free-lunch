# The Determinant's Free Lunch: Local Symmetry and the Permanent's Missing Moves

Agents that each hold one row of a matrix over a finite field can preserve its determinant using only pairwise contacts. The paper proves which matrices such agents can reach on a communication graph, gives an exact criterion for reaching a form whose determinant is readable from per-row summaries, and proves that the permanent's pair-local linear symmetries are only scalings and exchanges. Exhaustive enumeration over small fields and seeded trials over F_101 check the theorems, and further computations measure how far the permanent's local rigidity survives when moves may depend on the state.

[Read the paper](https://piatra.institute/papers/the-determinants-free-lunch.pdf)

## Build

```bash
uv run build.py          # -> paper/PAPER.pdf
```

Requires `pandoc` and `xelatex` on PATH. From the workspace you can also run `papers build the-determinants-free-lunch`.

## Reproduce

```bash
cd simulation
uv lock
uv run python run_all.py
```

The run writes `simulation/output/results.json`, per-trial records and the figure. `python3 scripts/numeral_audit.py` lists every numeral in the manuscript and whether `claims.yaml` binds it.
