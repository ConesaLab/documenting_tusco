#!/usr/bin/env python3
"""Publish figure assets consumed by manuscript and reviewer-response TeX sources."""

from __future__ import annotations

import argparse
import filecmp
import json
import shutil
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = REPO_ROOT / "config" / "submission_assets.json"


def load_manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text())


def is_same_file(source: Path, dest: Path) -> bool:
    if not dest.exists():
        return False
    return filecmp.cmp(source, dest, shallow=False)


def publish_target(name: str, target: dict, check: bool, prune: bool) -> int:
    output_root = REPO_ROOT / target["output_root"]
    output_root.mkdir(parents=True, exist_ok=True)

    errors = 0
    expected_names = set()

    for entry in target["files"]:
        source = REPO_ROOT / entry["source"]
        dest = output_root / entry["name"]
        expected_names.add(dest.name)

        if not source.exists():
            print(f"[missing] {name}: {source.relative_to(REPO_ROOT)}", file=sys.stderr)
            errors += 1
            continue

        if check:
            if not is_same_file(source, dest):
                print(f"[outdated] {name}: {dest.relative_to(REPO_ROOT)}", file=sys.stderr)
                errors += 1
            continue

        dest.parent.mkdir(parents=True, exist_ok=True)
        if not is_same_file(source, dest):
            shutil.copy2(source, dest)
            print(f"[copied] {source.relative_to(REPO_ROOT)} -> {dest.relative_to(REPO_ROOT)}")
        else:
            print(f"[ok] {dest.relative_to(REPO_ROOT)}")

    if prune:
        for path in output_root.iterdir():
            if path.is_file() and path.name not in expected_names:
                if check:
                    print(f"[orphan] {name}: {path.relative_to(REPO_ROOT)}", file=sys.stderr)
                    errors += 1
                else:
                    path.unlink()
                    print(f"[pruned] {path.relative_to(REPO_ROOT)}")

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description="Publish TeX-consumed figure assets from the manifest.")
    parser.add_argument(
        "--target",
        choices=["paper", "reviewer_round_1", "all"],
        default="all",
        help="Which asset target to publish or check.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Validate that staged assets match the manifest instead of copying.",
    )
    parser.add_argument(
        "--prune",
        action="store_true",
        help="Remove files in each output_root that are not declared in the manifest.",
    )
    args = parser.parse_args()

    manifest = load_manifest()
    targets = manifest.keys() if args.target == "all" else [args.target]

    errors = 0
    for name in targets:
        errors += publish_target(name, manifest[name], check=args.check, prune=args.prune)

    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
