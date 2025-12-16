#!/usr/bin/env python3
"""
Run TUSCO pipeline for all sensitivity analysis configurations.

This script:
1. Loads configurations from sensitivity_configs.py
2. Runs tusco_selector for each variant
3. Collects output gene sets for evaluation
4. Generates summary statistics

Usage:
    python run_sensitivity_pipeline.py [--dry-run] [--parallel N] [--config NAME]
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

# Add src to path for tusco_selector imports
REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
sys.path.insert(0, str(SRC_DIR))

from sensitivity_configs import get_unique_configs, TUSCOConfig


def run_single_config(config: TUSCOConfig, dry_run: bool = False) -> dict:
    """Run TUSCO pipeline for a single configuration.

    Args:
        config: Configuration to run
        dry_run: If True, just print the command without running

    Returns:
        Dictionary with run results
    """
    output_dir = REPO_ROOT / "reviewer_response" / "sensitivity_runs" / config.name
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build command
    cmd = [
        sys.executable, "-m", "tusco_selector",
        *config.to_cli_args(),
    ]

    # Override output directory in CLI args
    cmd_str = " ".join(cmd)

    result = {
        "name": config.name,
        "species": config.species,
        "sweep_type": config.sweep_type,
        "is_default": config.is_default,
        "output_dir": str(output_dir),
        "command": cmd_str,
        "success": False,
        "error": None,
    }

    if dry_run:
        print(f"[DRY RUN] {config.name}: {cmd_str}")
        result["success"] = True
        return result

    print(f"[RUNNING] {config.name}...")

    try:
        # Set up environment
        env = os.environ.copy()
        env["PYTHONPATH"] = f"{SRC_DIR}:{env.get('PYTHONPATH', '')}"

        # Run the pipeline
        proc = subprocess.run(
            cmd,
            cwd=str(REPO_ROOT),
            env=env,
            capture_output=True,
            text=True,
            timeout=3600,  # 1 hour timeout
        )

        if proc.returncode == 0:
            result["success"] = True
            print(f"[SUCCESS] {config.name}")
        else:
            result["error"] = proc.stderr or "Unknown error"
            print(f"[FAILED] {config.name}: {proc.stderr[:200] if proc.stderr else 'Unknown error'}")

        # Save stdout/stderr to log files
        log_dir = output_dir / "logs"
        log_dir.mkdir(exist_ok=True)
        (log_dir / "stdout.log").write_text(proc.stdout or "")
        (log_dir / "stderr.log").write_text(proc.stderr or "")

    except subprocess.TimeoutExpired:
        result["error"] = "Timeout (1 hour)"
        print(f"[TIMEOUT] {config.name}")
    except Exception as e:
        result["error"] = str(e)
        print(f"[ERROR] {config.name}: {e}")

    return result


def count_genes_in_output(output_dir: Path, species: str) -> dict:
    """Count genes in the output TUSCO files.

    Args:
        output_dir: Output directory from pipeline run
        species: Species code ("hsa" or "mmu")

    Returns:
        Dictionary with gene counts
    """
    species_name = "human" if species == "hsa" else "mouse"

    counts = {
        "total_genes": 0,
        "single_exon_genes": 0,
        "multi_exon_genes": 0,
    }

    # Main TUSCO file
    main_file = output_dir / f"tusco_{species_name}.tsv"
    if main_file.exists():
        with main_file.open() as f:
            # Skip header, count non-empty lines
            lines = [l for l in f if l.strip() and not l.startswith("#")]
            counts["total_genes"] = len(lines)

    # Single exon file
    single_file = output_dir / f"tusco_{species_name}_single_exon.tsv"
    if single_file.exists():
        with single_file.open() as f:
            lines = [l for l in f if l.strip() and not l.startswith("#")]
            counts["single_exon_genes"] = len(lines)

    # Multi-exon = total - single
    counts["multi_exon_genes"] = counts["total_genes"] - counts["single_exon_genes"]

    return counts


def run_all_configs(
    configs: list[TUSCOConfig],
    dry_run: bool = False,
    parallel: int = 1,
) -> list[dict]:
    """Run all configurations.

    Args:
        configs: List of configurations to run
        dry_run: If True, just print commands
        parallel: Number of parallel processes (1 = sequential)

    Returns:
        List of result dictionaries
    """
    results = []

    if parallel > 1 and not dry_run:
        # Parallel execution
        with ProcessPoolExecutor(max_workers=parallel) as executor:
            futures = {
                executor.submit(run_single_config, cfg, dry_run): cfg
                for cfg in configs
            }
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
    else:
        # Sequential execution
        for cfg in configs:
            result = run_single_config(cfg, dry_run)
            results.append(result)

    return results


def generate_summary(results: list[dict]) -> dict:
    """Generate summary statistics from run results."""
    summary = {
        "total_runs": len(results),
        "successful": sum(1 for r in results if r["success"]),
        "failed": sum(1 for r in results if not r["success"]),
        "by_sweep_type": {},
        "by_species": {},
        "gene_counts": {},
    }

    for result in results:
        sweep = result["sweep_type"]
        species = result["species"]

        # Count by sweep type
        if sweep not in summary["by_sweep_type"]:
            summary["by_sweep_type"][sweep] = {"total": 0, "success": 0}
        summary["by_sweep_type"][sweep]["total"] += 1
        if result["success"]:
            summary["by_sweep_type"][sweep]["success"] += 1

        # Count by species
        if species not in summary["by_species"]:
            summary["by_species"][species] = {"total": 0, "success": 0}
        summary["by_species"][species]["total"] += 1
        if result["success"]:
            summary["by_species"][species]["success"] += 1

        # Get gene counts for successful runs
        if result["success"]:
            output_dir = Path(result["output_dir"])
            if output_dir.exists():
                counts = count_genes_in_output(output_dir, species)
                summary["gene_counts"][result["name"]] = counts

    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Run TUSCO pipeline for sensitivity analysis configurations"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running",
    )
    parser.add_argument(
        "--parallel",
        type=int,
        default=1,
        help="Number of parallel processes (default: 1)",
    )
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Run only specific config (by name)",
    )
    parser.add_argument(
        "--sweep",
        type=str,
        choices=["splice", "tss", "expression"],
        default=None,
        help="Run only configs from specific sweep type",
    )
    parser.add_argument(
        "--species",
        type=str,
        choices=["hsa", "mmu"],
        default=None,
        help="Run only configs for specific species",
    )
    parser.add_argument(
        "--skip-defaults",
        action="store_true",
        help="Skip default configurations (use existing data)",
    )
    args = parser.parse_args()

    # Get configurations
    configs = get_unique_configs()

    # Filter by name
    if args.config:
        configs = [c for c in configs if c.name == args.config]
        if not configs:
            print(f"Error: Config '{args.config}' not found")
            sys.exit(1)

    # Filter by sweep type
    if args.sweep:
        configs = [c for c in configs if c.sweep_type == args.sweep]

    # Filter by species
    if args.species:
        configs = [c for c in configs if c.species == args.species]

    # Skip defaults
    if args.skip_defaults:
        configs = [c for c in configs if not c.is_default]

    if not configs:
        print("No configurations to run after filtering")
        sys.exit(0)

    print(f"Running {len(configs)} configurations...")
    if args.dry_run:
        print("(DRY RUN - commands will be printed but not executed)")
    print()

    # Run configurations
    results = run_all_configs(configs, dry_run=args.dry_run, parallel=args.parallel)

    # Generate and save summary
    summary = generate_summary(results)

    output_dir = REPO_ROOT / "reviewer_response" / "sensitivity_runs"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save results
    results_file = output_dir / "run_results.json"
    with results_file.open("w") as f:
        json.dump(results, f, indent=2)

    # Save summary
    summary_file = output_dir / "summary.json"
    with summary_file.open("w") as f:
        json.dump(summary, f, indent=2)

    # Print summary
    print()
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total runs: {summary['total_runs']}")
    print(f"Successful: {summary['successful']}")
    print(f"Failed: {summary['failed']}")
    print()

    if summary["gene_counts"]:
        print("Gene counts by configuration:")
        for name, counts in sorted(summary["gene_counts"].items()):
            print(f"  {name}: {counts['total_genes']} genes "
                  f"({counts['single_exon_genes']} single-exon, "
                  f"{counts['multi_exon_genes']} multi-exon)")

    print()
    print(f"Results saved to: {results_file}")
    print(f"Summary saved to: {summary_file}")


if __name__ == "__main__":
    main()
