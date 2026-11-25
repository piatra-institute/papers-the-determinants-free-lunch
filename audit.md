# Audit

## 2026-11-26 — revision round 1 (Codex referee)

Scope: findings codex:F1 to codex:F5.
Changes:
  - Theorem 6 (rank-two rigidity) with proof; certificate_exhaustive and certificate_sampled in simulation/localsym.py; dimension-count remark removed.
  - Paired greedy experiment (greedy_paired) replacing separate populations; Figure 1 panels (c) and (d) updated.
  - Abstract adds the characteristic restriction and the rigidity result; discussion no longer claims quadratic solving is necessary.
  - References added: anderson1956, drton2007, schwartz1980; marcus1962 now read in full and bound.
  - claims.yaml regenerated (55 computation, 8 source claims, all source claims verified: false); src-mulmuley rebound.
Verification:
  - recorded run `main` re-recorded (verification/main.json); numeral audit: unbound prose numerals are labels, years or design values only.
  - build: 14 pages, no missing glyphs; all pages rendered and inspected.
  - standards PASS (248 abstract words); refs --verify no mismatch; check => NEEDS_REVIEW (no FAIL).

## 2026-11-26 — v1 authored

Scope: new paper written from seed leads; all computations reimplemented in `simulation/` (no seed code executed).
Changes:
  - Theorems 1-4, Corollary 5, Proposition 6 with proofs in paper/PAPER.md.
  - Recorded run `main` (`papers run`, receipt verification/main.json): unit tests, exhaustive enumeration, seeded trials, census computations, figure.
  - claims.yaml: 45 computation claims, 6 source claims (all verified: false), definitions, assumptions and interpretations; `scripts/numeral_audit.py` lists every numeral and confirms unbound prose numerals are labels, years or design values.
Verification: see `papers check` output in the final report; review records pending an independent reviewer.
