#!/usr/bin/env python3
"""Main entry point for running tusco_selector as a module.

Usage:
    python -m tusco_selector hsa --output-dir output_dir [options]
    python -m tusco_selector mmu --output-dir output_dir [options]
"""

from tusco_selector.cli import main

if __name__ == "__main__":
    main()
