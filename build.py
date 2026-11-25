#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["pyyaml"]
# ///
"""Canonical Piatra paper build.

Renders <paper_dir>/PAPER.md -> <paper_dir>/PAPER.pdf via pandoc + xelatex,
using a shared LaTeX preamble (inlined below). Title and subtitle are separated
at the first colon, followed by the institute name and publication date.

Two modes, one file:
  * vendored   — copied to a paper repo root; `uv run build.py` builds that repo
                 (REPO = the directory this file sits in).
  * canonical  — `uv run build/build.py <repo_path>` builds an arbitrary repo
                 (REPO = argv[1]); used when a paper has no vendored recipe.

Configuration comes from <REPO>/metadata.yaml when present:
  paper_dir : "paper" (default) | "docs"   — where PAPER.md lives
  header    : running-header text          — falls back to `title`, then folder name
Everything else is convention. pandoc + xelatex are the only external tools;
the inlined preamble means the build has zero repo-relative dependencies.
"""

import hashlib
import json
import os
import re
import sys
import subprocess
import tempfile
import pathlib
from datetime import datetime, timezone

import yaml

# Shared preamble for standalone builds: consistent typography and running
# headers, with table, math, and theorem support available to every paper.
# @@HEADER@@ is replaced with the running-header text at build time.
PREAMBLE = r"""
% Avoid isolated first/last lines in paragraphs and reference entries.
\widowpenalty=10000
\clubpenalty=10000
\displaywidowpenalty=10000

% Running headers
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small\textit{@@HEADER@@}}
\fancyhead[R]{\small\thepage}
\renewcommand{\headrulewidth}{0.4pt}

% Title page: no header, just centered page number at bottom
\fancypagestyle{plain}{\fancyhf{}\fancyfoot[C]{\thepage}\renewcommand{\headrulewidth}{0pt}}

% Title page styling
\usepackage{titling}
\pretitle{\vspace*{\fill}\begin{center}\LARGE\bfseries}
\posttitle{\end{center}\vskip 3em}
\preauthor{\begin{center}\large}
\postauthor{\end{center}\vskip 1.5em}
\predate{\begin{center}\large}
\postdate{\end{center}\vspace*{\fill}\thispagestyle{empty}\newpage}

% Math packages
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{mathtools}

% Tables
\usepackage{booktabs}

% Keep references readable without the body text's expanded leading.
\usepackage{etoolbox}
\AtBeginEnvironment{CSLReferences}{\setstretch{1.0}}

% Theorem environments
\usepackage{amsthm}
\newtheorem{definition}{Definition}
\newtheorem{proposition}{Proposition}
\newtheorem{lemma}{Lemma}
""".lstrip()


def tex_escape(s: str) -> str:
    replacements = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%",
                    "#": r"\#", "_": r"\_", "$": r"\$", "{": r"\{",
                    "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(replacements.get(c, c) for c in s)


def title_for_pdf(title: str) -> str:
    """Give the subtitle its own line without changing the stored listing title.

    Pandoc reads backslash-newline as a hard line break. Existing breaks within
    a subtitle remain intact. Titles without a colon retain their punctuation
    and any deliberate break (for example, a title ending in a question mark).
    """
    return re.sub(r":(?:[ \t]*\\[ \t]*\r?\n|\s+)",
                  lambda _: ":\\\n", title, count=1).strip()


def title_filter(title: str) -> str:
    # CLI -M values are literal strings, not Markdown. Parse the title through
    # Pandoc so hard breaks, emphasis and math retain their native inline types.
    value = json.dumps(title_for_pdf(title), ensure_ascii=False)
    return ("function Meta(meta)\n"
            f"  local blocks = pandoc.read({value}, 'markdown').blocks\n"
            "  if #blocks > 0 then meta.title = pandoc.MetaInlines(blocks[1].content) end\n"
            "  return meta\nend\n")


def load_meta(repo: pathlib.Path) -> dict:
    p = repo / "metadata.yaml"
    if p.exists():
        return yaml.safe_load(p.read_text()) or {}
    return {}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_inputs(repo, meta):
    """Inputs include all local render resources, not only manuscript mtime."""
    import os
    paper = repo / meta.get("paper_dir", "paper")
    files = {repo / "metadata.yaml", repo / "references.yaml", repo / "references.json"}
    extensions = {".md", ".tex", ".sty", ".cls", ".csl", ".bib", ".yaml", ".json",
                  ".png", ".jpg", ".jpeg", ".svg", ".eps", ".pdf"}
    for root, dirs, names in os.walk(paper):
        dirs[:] = sorted(d for d in dirs if not d.startswith(".") and d != "__pycache__")
        for name in names:
            path = pathlib.Path(root) / name
            if name not in {"PAPER.pdf", "build-manifest.json"} and path.suffix.lower() in extensions:
                files.add(path)
    for name in meta.get("build_resources", []):
        path = (repo / name).resolve()
        if not path.is_relative_to(repo.resolve()) or not path.is_file():
            raise ValueError(f"invalid or missing build resource: {name}")
        files.add(path)
    return {str(p.relative_to(repo)): digest(p) for p in sorted(files) if p.is_file()}


def versions():
    return {tool: subprocess.run([tool, "--version"], check=True, capture_output=True,
                                text=True).stdout.splitlines()[0]
            for tool in ("pandoc", "xelatex")}


def publication_date(value):
    """UTC midnight for the publication date, without changing its display text."""
    if value is None or value == "":
        return datetime(1970, 1, 1, tzinfo=timezone.utc)
    raw = str(value).strip()
    if re.fullmatch(r"\d{4}-\d{2}(?:-\d{2})?", raw):
        parts = [int(part) for part in raw.split("-")]
        if len(parts) == 2:
            parts.append(1)
        return datetime(*parts, tzinfo=timezone.utc)
    months = "January February March April May June July August September October November December".split()
    match = re.fullmatch(r"([A-Za-z]+) (\d{4})", raw)
    if match and match[1].casefold() in [month.casefold() for month in months]:
        month = [name.casefold() for name in months].index(match[1].casefold()) + 1
        return datetime(int(match[2]), month, 1, tzinfo=timezone.utc)
    raise ValueError(f"invalid publication date {value!r}; expected YYYY-MM-DD, YYYY-MM or Month YYYY")


def render(repo, md=None, pdf=None, *, record=True):
    repo = repo.resolve()
    meta = load_meta(repo)
    paper_dir = repo / meta.get("paper_dir", "paper")
    header = meta.get("header") or meta.get("title") or repo.name
    md = md or paper_dir / "PAPER.md"
    pdf = pdf or paper_dir / "PAPER.pdf"
    if not md.exists():
        print(f"error: {md} not found", file=sys.stderr)
        return 1

    inputs = build_inputs(repo, meta)
    tool_versions = versions()
    preamble = PREAMBLE.replace("@@HEADER@@", tex_escape(str(header)))
    # Temporary paths can affect the trailer ID even with SOURCE_DATE_EPOCH.
    # Pandoc emits these two 16-byte IDs through xdvipdfmx's first-page special.
    document_id = hashlib.sha256(json.dumps(inputs, sort_keys=True).encode("utf-8")).hexdigest()[:32]
    with tempfile.TemporaryDirectory(prefix="piatra-build-") as temp:
        temp_dir = pathlib.Path(temp)
        preamble_path = temp_dir / "preamble.tex"
        preamble_path.write_text(preamble)
        built = temp_dir / "paper.pdf"
        cmd = [
            "pandoc", str(md),
            "-o", str(built),
            "--pdf-engine=xelatex",
            # xdvipdfmx hashes the XDV/PDF filenames into font subset tags.
            # Stable relative names keep Pandoc's random temp path out of them.
            "--pdf-engine-opt=-output-directory=.",
            "--resource-path", f"{paper_dir}:{repo}",
            "-H", str(preamble_path),
            "-V", f"pdf-trailer-id=<{document_id}><{document_id}>",
            "-V", "mainfont=Palatino",
            "-V", "mathfont=Palatino",
            "-V", "monofont=Menlo",
            "-V", "fontsize=11pt",
            "-V", "geometry:margin=1in",
            "-V", "linestretch=1.15",
        ]
        front = re.match(r"\A---\s*\n(.*?)\n(?:---|\.\.\.)\s*(?:\n|$)", md.read_text(), re.S)
        metadata = (yaml.safe_load(front.group(1)) or {}) if front else {}
        title = metadata.get("title", meta.get("title", ""))
        if title:
            filter_path = temp_dir / "title.lua"
            filter_path.write_text(title_filter(str(title)))
            cmd.extend(["--lua-filter", str(filter_path)])
        cmd.extend(["--metadata", "author=PIATRA . INSTITUTE"])
        date = metadata.get("date") or meta.get("date")
        epoch = int(publication_date(date).timestamp())
        environment = dict(os.environ, SOURCE_DATE_EPOCH=str(epoch), FORCE_SOURCE_DATE="1")
        environment["TEXINPUTS"] = str(repo) + os.pathsep + environment.get("TEXINPUTS", "")
        if date:
            cmd.extend(["--metadata", "date=" + str(date)])
        bibliography = metadata.get("bibliography", meta.get("bibliography"))
        if not bibliography:
            bibliography = [name for name in ("references.yaml", "references.json") if (repo / name).exists()]
        if bibliography:
            cmd.append("--citeproc")
            for name in bibliography if isinstance(bibliography, list) else [bibliography]:
                source = md.parent / name
                if not source.exists():
                    source = repo / name
                cmd.extend(["--bibliography", str(source.resolve())])
        proc = subprocess.run(cmd, cwd=temp_dir, capture_output=True, text=True, env=environment)
        log = proc.stdout + proc.stderr
        log_path = pdf.with_suffix(".build.log")
        log_path.write_text(log)
        if proc.returncode or re.search(r"Missing character:|[Cc]itation .+ not found", log):
            print(log, file=sys.stderr)
            print(f"build failed; previous PDF preserved. Log: {log_path}", file=sys.stderr)
            return proc.returncode or 1
        if inputs != build_inputs(repo, meta):
            print("build inputs changed during rendering; previous PDF preserved", file=sys.stderr)
            return 1
        import shutil
        with tempfile.NamedTemporaryFile(dir=pdf.parent, prefix=".paper-", suffix=".pdf", delete=False) as staged:
            staged_path = pathlib.Path(staged.name)
        try:
            shutil.copyfile(built, staged_path)
            staged_path.replace(pdf)
        finally:
            staged_path.unlink(missing_ok=True)
    if record:
        manifest = {"schema": 1, "built": datetime.now(timezone.utc).isoformat(),
                    "source_date_epoch": epoch,
                    "inputs": inputs, "recipe_sha256": digest(pathlib.Path(__file__)),
                    "versions": tool_versions, "pdf_sha256": digest(pdf)}
        (pdf.parent / "build-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pdf)
    return 0


def main() -> int:
    repo = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else pathlib.Path(__file__).resolve().parent
    try:
        return render(repo)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"build failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
