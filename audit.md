# Audit

## 2026-11-26 — v1 authored

Scope: new paper written from seed leads; all computations reimplemented in `simulation/` (no seed code executed).
Changes:
  - Theorems 1-4, Corollary 5, Proposition 6 with proofs in paper/PAPER.md.
  - Recorded run `main` (`papers run`, receipt verification/main.json): unit tests, exhaustive enumeration, seeded trials, census computations, figure.
  - claims.yaml: 45 computation claims, 6 source claims (all verified: false), definitions, assumptions and interpretations; `scripts/numeral_audit.py` lists every numeral and confirms unbound prose numerals are labels, years or design values.
Verification: see `papers check` output in the final report; review records pending an independent reviewer.
