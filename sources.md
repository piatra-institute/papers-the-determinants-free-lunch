# Sources

Every DOI was checked against Crossref (or OpenAlex) before entry; `papers refs --verify` records the identity cache. Depth of inspection per source:

- anderson1956. Full text (Berkeley digital-collection scan of the Third Berkeley Symposium Proceedings, vol. 5, pp. 111-150; identity via OpenAlex W1748168094 and Project Euclid euclid.bsmsp/1200511860; no DOI). Lemma 5.1 and Theorem 5.1 with proof on pp. 118-119: identification of a low-rank factor structure from off-diagonal entries under a row-deletion rank condition. Greek symbols in the quote restored from the OCR text layer. Bound in claims.yaml (src-anderson-rubin).
- angluin2006. Identity (Crossref). Cited for the general scope of population protocols; no specific claim bound.
- botta1967. Identity (Crossref, DOI 10.1090/S0002-9939-1967-0213376-6). The AMS PDF was blocked by a bot check, so the statement of its theorem (fields, hypotheses) was not inspected; cited only for the well-known permanent-preserver result and not bound in claims.yaml. The paper's own Theorem 4 is proved independently.
- boyd2006. Abstract (OpenAlex).
- burgisser2000. Identity (Crossref). Cited for VNP-completeness of the permanent in characteristic other than two, well-known content; not bound.
- deb2006. Abstract (OpenAlex).
- drton2007. Identity (Crossref, DOI 10.1007/s00440-006-0033-2; Crossref issued year 2006, print volume 2007). Cited for the algebraic setting of low-rank symmetric matrices with unobserved diagonal; full text not inspected; not bound.
- dieudonne1949. Identity (Crossref gives issued year 1948, Archiv der Mathematik 1(4); often cited as 1949). French title marked language fr. Not inspected; not bound.
- dijkstra1974. Abstract (ACM via OpenAlex). Bound in claims.yaml.
- kasteleyn1963. Identity (Crossref). Not inspected; credited as the origin of Pfaffian orientations following McCuaig's and Robertson-Seymour-Thomas's attributions.
- levin2025. Full text of PsyArXiv version 3 (OSF file 7k32r), Section 4, p. 16, quoted in claims.yaml. Preprint, not peer-reviewed.
- li2001. Identity (Crossref). Survey; cited for the Frobenius determinant-preserver theorem. Not inspected; not bound.
- little1975. Identity (Crossref). Content stated via Robertson-Seymour-Thomas statement 1.2.
- marcus1962. Full text (Cambridge publisher PDF), main Theorem on p. 179: linear maps of m x n matrices preserving permanental compounds of order r, field with at least r elements, characteristic not 2; the permanent half of Theorem 4 is its case of two rows and r = 2. Bound in claims.yaml (src-marcus-may).
- mccuaig2004. Full text (open-access EJC PDF), p. 11 and abstract. Bound in claims.yaml.
- mulmuley2001. Abstract (OpenAlex). Bound in claims.yaml.
- robertson1999. Full text of arXiv version math/9911268 (Section 1); the Annals version of record was not compared page by page. Bound in claims.yaml.
- schwartz1980. Identity (Crossref, DOI 10.1145/322217.322225). Cited for the standard bound that a nonzero polynomial of degree d vanishes at a uniform point of F_p^m with probability at most d/p; full text not inspected; not bound.
- tiskin2007. Identity (Crossref); title supports "communication-efficient parallel generic pairwise elimination". Not bound.
- valiant1979classes, valiant1979permanent. Identity (Crossref). Well-known content (VNP-completeness; #P-completeness); not bound.
- vazirani1989. Abstract (OpenAlex). The equivalence statement used in the paper is taken from Robertson-Seymour-Thomas statement 1.1.
- zhang2024. Version-of-record abstract (OpenAlex) and full text of arXiv 2401.05375v1; the Sage PDF could not be retrieved. Bound in claims.yaml at abstract depth.

Not cited because they could not be inspected: Frobenius (1897), Pólya (1913), Szegő (1913); their results are reported through McCuaig (2004) and Li and Pierce (2001). No reference title contains an arrow or other glyph missing from the build font.
