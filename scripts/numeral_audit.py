"""List every numeral token in paper/PAPER.md and report whether claims.yaml binds it.

A token is bound when a computation claim reports exactly that token and the claim's passage
covers the token's position in its paragraph. Tokens inside inline math are listed separately:
they are formula constants, field sizes or labels unless a claim binds them. Unbound prose tokens
are printed with context so the author can confirm each is a label, year, locator or design value.
Run from the paper repository root: python3 scripts/numeral_audit.py
"""
import re
from pathlib import Path

import yaml

root = Path(__file__).resolve().parent.parent
body = (root / "paper" / "PAPER.md").read_text().split("---", 2)[2].split("## References")[0]
claims = yaml.safe_load((root / "claims.yaml").read_text())["claims"]


def norm(t):
    return " ".join(t.split())


bound = [(norm(c["passage"]), str(c["evidence"]["reported"])) for c in claims if c["kind"] in ("computation", "empirical")]
tok = re.compile(r"(?<![\w.,])\d[\d,]*(?:\.\d+)?%?(?![\w])")
total = math_unbound = 0
prose_unbound = []
for para in body.split("\n"):
    para = norm(para)
    if not para:
        continue
    spans = []
    for pas, rep in bound:
        for m in re.finditer(re.escape(pas), para):
            spans.append((m.start(), m.end(), rep))
    math = [(m.start(), m.end()) for m in re.finditer(r"\$[^$]*\$", para)]
    cites = [(m.start(), m.end()) for m in re.finditer(r"\[@[^\]]*\]", para)]
    for m in tok.finditer(para):
        t = m.group(0).rstrip(",")
        if any(a <= m.start() < b for a, b in cites):
            continue
        total += 1
        if any(a <= m.start() < b and rep == t for a, b, rep in spans):
            continue
        if any(a <= m.start() < b for a, b in math):
            math_unbound += 1
            continue
        prose_unbound.append((t, para[max(0, m.start() - 50): m.end() + 30]))
print(f"numeral tokens (outside citations): {total}; unbound inside math: {math_unbound}; unbound in prose: {len(prose_unbound)}")
for t, ctx in prose_unbound:
    print(f"  {t!r:>12}  ...{ctx}...")
