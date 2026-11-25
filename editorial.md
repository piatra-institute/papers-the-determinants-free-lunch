# Editorial record

Status: author's v1; independent review pending.

## Contribution and reader

For readers in algebraic complexity, distributed computing and the "basal intelligence" discussion: an exact, checkable account of why locally blind agents can expose a determinant (pair-local SL_2 generates SL on each component; disjoint pivot profiles decide readability; a randomised protocol attains it almost surely) and a proof that the permanent's pair-local linear symmetries are only scalings and exchanges. A measured negative result limits the no-go: state-dependent permanent-preserving moves reach readable forms for all 3x3 matrices over F_5 and often at n = 4, 5, 6 over small fields.

## Reverse outline

1. Opening: the question, the sorting analogy and Levin's framing as framing, the three-part answer, prior work.
2. Model: moves, readable form, readout, pivot profiles, q, protocol.
3. Theorems 1-3 with proofs: orbits, readability, almost-sure completion.
4. Checks: exhaustive enumeration and randomised trials; Figure 1.
5. Theorem 4 (pair-local preservers, SL_2 vs monomial) with proof, exhaustive search, Corollary 5.
6. State-dependent moves: det unchanged; permanent fibres, frozen pairs, exhaustive n = 3 reachability, greedy agents.
7. Sign-rule faults: Proposition 6 (cuts), relation to Pólya-Szegő and Pfaffian orientations, clique projection.
8. Discussion: what the free lunch is, its costs, the limit of the no-go, VP vs VNP not addressed, open questions.

## Strongest rival

The results restate Gaussian elimination and the classical global preserver theorems (Frobenius, Dieudonné, Marcus-May, Botta). The paper's answer: the pair-local statements and the graph-and-data criterion are not consequences of the global theorems as stated, and the state-dependent measurements show where the local contrast fails. Novelty is claimed only for these localised statements and measurements, which are elementary.

## Decisions and remaining checks

- Seed results were reimplemented; the seed's F_2 counts (168 matrices, 8 graphs) were reproduced but are not cited in the paper.
- Sign-defect section kept short (owner instruction) and credits Pólya-Szegő, Kasteleyn, Little, Vazirani-Yannakakis, Robertson-Seymour-Thomas and McCuaig.
- Open: Botta and Marcus-May full texts not inspected; the dimension-count remark on frozen pairs is unproved and labelled so.
- Pending independent review: claims coverage, source verdicts, editorial and visual review records.
