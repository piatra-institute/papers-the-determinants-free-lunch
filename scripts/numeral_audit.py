"""Exhaustive coverage audit of numerals and number words in paper/PAPER.md.

Scope: from the abstract to the end of the last section before "Code and data availability".
Tokens: every digit token (including grouped numbers and percentages) and every number word
(zero to twelve, twice, single, double, half, both, first, second, third), case-insensitive.
An occurrence is covered when some claims.yaml passage contains it (the passage occurs in the
paragraph and spans the token). Digit tokens in prose (outside inline math and citations) are
additionally listed when no computation claim reports that exact token over its position; such
tokens must be labels, years or design values bound by definition claims.

Run from the paper repository root: python3 scripts/numeral_audit.py
Exit status 1 when any occurrence is uncovered.
"""
import re
import sys
from pathlib import Path

import yaml

WORDS = ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten",
         "eleven", "twelve", "twice", "single", "double", "half", "both", "first", "second", "third"]
WORD_RE = re.compile(r"\b(" + "|".join(WORDS) + r")\b", re.IGNORECASE)
DIGIT_RE = re.compile(r"(?<![\w.,])\d[\d,]*(?:\.\d+)?%?(?![\w])")


def norm(t):
    return " ".join(t.split())


def scope_paragraphs(text):
    body = text.split("---", 2)[2]
    start = body.index("## Abstract")
    end = body.index("## Code and data availability")
    out = []
    for line in body[start:end].split("\n"):
        line = norm(line)
        if line and not line.startswith("## "):
            out.append(line)
    return out


def audit(text, claims):
    passages = [(norm(c["passage"]), c) for c in claims]
    report = {"tokens": 0, "uncovered": [], "prose_digits_without_computation": []}
    for para in scope_paragraphs(text):
        spans = []
        for pas, c in passages:
            for m in re.finditer(re.escape(pas), para):
                spans.append((m.start(), m.end(), c))
        cites = [(m.start(), m.end()) for m in re.finditer(r"\[@[^\]]*\]", para)]
        math = [(m.start(), m.end()) for m in re.finditer(r"\$[^$]*\$", para)]
        occ = [(m.start(), m.group(0), "word") for m in WORD_RE.finditer(para)]
        occ += [(m.start(), m.group(0).rstrip(","), "digit") for m in DIGIT_RE.finditer(para)
                if not any(a <= m.start() < b for a, b in cites)]
        for pos, tok, kind in occ:
            report["tokens"] += 1
            cov = [c for a, b, c in spans if a <= pos < b]
            if not cov:
                report["uncovered"].append((tok, para[max(0, pos - 50): pos + 30], para[:60]))
            elif kind == "digit" and not any(a <= pos < b for a, b in math):
                if not any(c["kind"] in ("computation", "empirical") and str(c["evidence"]["reported"]) == tok for c in cov):
                    report["prose_digits_without_computation"].append((tok, para[max(0, pos - 40): pos + 20]))
    return report


def main():
    root = Path(__file__).resolve().parent.parent
    text = (root / "paper" / "PAPER.md").read_text()
    claims = yaml.safe_load((root / "claims.yaml").read_text())["claims"]
    r = audit(text, claims)
    print(f"tokens audited: {r['tokens']}; uncovered: {len(r['uncovered'])}; "
          f"prose digits covered only by definition/assumption/interpretation claims: {len(r['prose_digits_without_computation'])}")
    for tok, ctx, head in r["uncovered"]:
        print(f"  UNCOVERED {tok!r}: ...{ctx}...   [paragraph: {head}]")
    for tok, ctx in r["prose_digits_without_computation"]:
        print(f"  label/design digit {tok!r}: ...{ctx}...")
    return 1 if r["uncovered"] else 0


if __name__ == "__main__":
    sys.exit(main())
