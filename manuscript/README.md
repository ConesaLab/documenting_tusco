# Manuscript Layout

This directory hosts the LaTeX sources and supporting assets for the TUSCO paper.

## Structure

- `paper/`
  - `src/` – primary LaTeX sources for the manuscript and supplementary PDFs.
  - `assets/` – final figure PDFs consumed by the paper build.
  - `styles/` – Springer Nature class and bibliography styles required by the SN template.
  - `scripts/` – helpers used to refresh numerical values in the manuscript.
  - `notes/` – planning and conversion notes provided by the authors.
  - `build/` – LaTeX build artefacts (ignored in git).
- `supplementary-note/`
  - `src/` – LaTeX source for the standalone supplementary note.
  - `assets/` – Illustrator sources and derived figure panels referenced by the note.
  - `build/` – LaTeX build artefacts (ignored in git).

## Building

From `manuscript/paper/src`:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error tusco_paper.tex
```

Repeat with `tusco_paper_supplement_figures.tex` and `tusco_paper_supplement_tables.tex` as needed.

From repo root, refresh the staged paper figures before compiling:

```bash
python scripts/publish_assets.py --target paper
```

From `manuscript/supplementary-note/src`:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error -outdir=../build supplementary_note.tex
```

All commands assume `latexmk` is available on your PATH. Clean builds can be produced with `latexmk -C -outdir=../build <file.tex>`.
