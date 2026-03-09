# Figure Script I/O Conventions

This repository now provides a shared `figure_context()` helper in `scripts/figure_utils.R`
to standardise argument parsing, data resolution, logging, and output directories across
all figure scripts.

## Using the context
- Construct the context near the top of each script:
  ```r
  ctx <- figure_context(defaults = list(out_dir = "..", width = 7.09, height = 3.55))
  params <- ctx$params
  ```
- Write outputs via `ctx$plot_dir`, `ctx$table_dir`, or the convenience wrappers
  `ctx$save_plot()` and `ctx$write_table()`.
- Resolve inputs with `ctx$resolve_input(...)`; it searches the figure-local `data/`,
  shared `figs/data/`, and repository `data/raw`/`data/processed` trees.
- Call `ctx$start_log()` when a run log is needed; the helper returns a closure that
  should be invoked on exit to close the sinks.

## Migration notes
- **figure-03** scripts: drop the bespoke `resolve_path()` and `read_tsv_safe()`
  definitions in favour of the shared helpers. Replace manual `dir.create()` logic
  with the context outputs.
- **figure-04** scripts: replace the nested `first_existing()` searches for TUSCO and
  evaluation assets by `ctx$resolve_input()`; minimise duplicated logging setup once
  `ctx$start_log()` is in place.
- **figure-01** scripts: swap the inline project-root traversal for
  `figure_context()` and `ctx$resolve_input()` so GTF/TPM files resolve consistently.
- **supplementary figures** (`supp-fig-02/03/04`): remove per-script data-root
  discovery; use the context search list so LOCAL_ONLY runs keep working without
  absolute paths.
- **Standalone utilities** (`tissue_grid_plot.R`, analysis helpers): create the
  context with `defaults = list(out_dir = ".")` when outputs should remain beside
  the script, and re-use the shared readers.

As additional figures migrate, prefer augmenting `figure_context()` rather than
introducing new ad-hoc I/O helpers.
