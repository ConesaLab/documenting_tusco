# tusco-paper

Code for generating the plots in [the associated bioRxiv paper](https://www.biorxiv.org/content/10.1101/2025.08.23.671926v1), including the `tusco-novel` and `tusco-selector` modules, lives under `src/`.

Project structure (reorganized for clarity and reproducibility):

- `data/`
  - `raw/`: External inputs and large downloads (reference, lrgasp, expression, nih, spike-ins)
  - `processed/`: Derived, versioned outputs (e.g., `processed/tusco/{hsa,mmu}`)
- `figs/`
  - `figure-0N/` and `supp-fig-0N/`: Each with `code/`, `plots/`, `tables/`
- `src/`: Reusable Python code (`tusco_selector`, `tusco_novel_simulator`)
- `R/`: Shared R helpers (e.g., `R/paths.R`)
- `envs/`: Conda environments (e.g., `envs/tusco_selector.yml`)
- `config/`: Central configuration (e.g., `config/project.yml`)
- `workflows/`: Pipelines and SLURM jobs (optional)

Notes
- All scripts now read from `data/raw` and write reusable results to `data/processed`.
- Run Python modules from repo root with `export PYTHONPATH=src`.
- Download the Tusco dataset archive from https://tusco-paper-data.s3.eu-north-1.amazonaws.com/data.zip and extract it into the repository root so a `data/` directory is available (the helper script below will do this automatically if missing).
- Run `/Users/tianyuan/Desktop/github_dev/documenting_tusco/scripts/run_all_figs.sh` from the repository root to regenerate the figure outputs; it downloads the dataset on-demand and executes every figure script.
