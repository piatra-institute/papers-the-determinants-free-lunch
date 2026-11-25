# Brief

## Question

Agents that each hold one row of an invertible matrix, and act only on their own row and one neighbour's, can preserve the determinant without seeing it. What exactly supplies this capacity, when does it let them reach a form from which the determinant can be read, and does the permanent offer the same local moves?

## Claim

On any communication graph, determinant-one moves of adjacent row pairs generate the full special linear group on each connected component; a readable form is reachable exactly when the components' pivot profiles are disjoint, and a randomised local protocol reaches it almost surely under that condition. For the permanent, the linear pair-local preservers (n at least 3, characteristic not two) are only scalings and exchanges, so linear permanent-preserving agents cannot create zeros. Whether state-dependent permanent-preserving moves are equally rigid is a separate question, to be measured.

## Kind

formal-model with proofs and a simulation (has_simulation: true): theorems, exhaustive enumeration over small fields, seeded randomised trials over F_101.

## Strongest serious alternative

The determinant's local capacity is just Gaussian elimination, known for two centuries, and the permanent's lack of it is the classical Marcus-May-Botta preserver theorem; nothing new follows. Response: the global theorems do not say which symmetries two rows acting alone possess, nor what a network of such pairs can reach; the paper proves the pair-local versions, the exact graph-and-data feasibility criterion, and measures (with a negative result) how far permanent rigidity survives state-dependent moves.

## Cornerstone literature

Valiant 1979 (both papers); Bürgisser 2000; Frobenius and Dieudonné on determinant preservers (via Li and Pierce 2001); Marcus and May 1962; Botta 1967; Mulmuley and Sohoni 2001; Pólya and Szegő 1913 (via McCuaig 2004); Kasteleyn; Little 1975; Vazirani and Yannakakis 1989; Robertson, Seymour and Thomas 1999; McCuaig 2004; Dijkstra 1974; Boyd et al. 2006; Deb, Médard and Choute 2006; Angluin et al. 2006; Tiskin 2007; Zhang, Goldstein and Levin 2024; Levin 2025.
