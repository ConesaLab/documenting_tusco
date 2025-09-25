"""
Utilities to download external resources declared in species resource modules.

Behaviour is unchanged; this refactor adds type hints and clarifying docstrings
and keeps network semantics identical to the original implementation.
"""

import os
from typing import Dict, Iterator, List, Tuple

import requests

# Optional dependencies
from tusco_selector.optional_deps import HAS_COLORS, HAS_TQDM, Fore, Style, tqdm


def fetch(url: str, output: str, cache_dir: str) -> str:
    """Download a single file to ``cache_dir/output``.

    If the destination file already exists, the download is skipped and the
    existing path is returned. Progress is shown when ``tqdm`` is available.
    """
    # Ensure the cache directory exists
    os.makedirs(cache_dir, exist_ok=True)
    dest = os.path.join(cache_dir, output)

    # Skip download if already present
    if os.path.exists(dest):
        if HAS_COLORS:
            print(f"{Fore.YELLOW}Skipping download, exists:{Style.RESET_ALL} {dest}")
        else:
            print(f"Skipping download, exists: {dest}")
        return dest

    # Stream download to file with progress bar
    response = requests.get(url, stream=True)
    response.raise_for_status()

    # Gracefully cope with lightweight mock objects (tests) that might not
    # provide a ``headers`` attribute.
    headers = getattr(response, "headers", {}) or {}
    total_size = int(headers.get("content-length", 0))
    block_size = 1024 * 1024  # 1MB chunks

    if HAS_COLORS:
        print(f"{Fore.CYAN}Downloading:{Style.RESET_ALL} {output}")
    else:
        print(f"Downloading: {output}")

    if HAS_TQDM:
        # Use tqdm even when file size is unknown
        with open(dest, "wb") as f:
            # If total_size is known, show that, otherwise show an indeterminate progress bar
            if total_size > 0:
                progress_bar = tqdm(
                    total=total_size, unit="B", unit_scale=True, desc=output, ncols=80
                )
                for chunk in response.iter_content(chunk_size=block_size):
                    if chunk:
                        f.write(chunk)
                        progress_bar.update(len(chunk))
                progress_bar.close()
            else:
                # For unknown file sizes, use tqdm with unit_scale but no total
                for chunk in tqdm(
                    response.iter_content(chunk_size=block_size),
                    desc=output,
                    unit="B",
                    unit_scale=True,
                    ncols=80,
                ):
                    if chunk:
                        f.write(chunk)
    else:
        # Fallback to basic output without progress bar
        with open(dest, "wb") as f:
            for chunk in response.iter_content(chunk_size=block_size):
                if chunk:
                    f.write(chunk)
                    print(".", end="", flush=True)
            print()  # New line after download completes

    if HAS_COLORS:
        print(f"{Fore.GREEN}Downloaded:{Style.RESET_ALL} {dest}")
    else:
        print(f"Downloaded: {dest}")
    return dest


def fetch_all(resources: dict, base_dir: str) -> List[str]:
    """Download all leaf resources described in a nested mapping.

    The nested ``resources`` should contain leaves with ``{"url", "output"}``
    keys (as used by species ``EXTERNAL_RESOURCES``). The output directory tree
    mirrors the nested keys up to the leaf. Returns the list of downloaded (or
    pre-existing) absolute file paths.
    """

    def _walk(node: Dict, key_path: List[str]) -> Iterator[Tuple[List[str], str, Dict]]:
        """Yield (path_parts, leaf_key, meta) for every downloadable leaf."""
        for leaf_key, subnode in node.items():
            if not isinstance(subnode, dict):
                continue
            if {"url", "output"}.issubset(subnode.keys()):
                yield key_path, leaf_key, subnode
            else:
                yield from _walk(subnode, key_path + [leaf_key])

    downloaded: List[str] = []

    # Count the total number of files to download for better reporting
    resources_to_download = list(_walk(resources, []))
    total_count = len(resources_to_download)

    if HAS_COLORS:
        print(f"{Fore.BLUE}Found {total_count} resources to process{Style.RESET_ALL}")
    else:
        print(f"Found {total_count} resources to process")

    for idx, (path_parts, leaf_key, meta) in enumerate(resources_to_download, 1):
        rel_dir = os.path.join(*path_parts) if path_parts else ""
        dest_dir = os.path.join(base_dir, rel_dir)

        url = meta["url"]
        output_name = meta["output"]

        pretty_name = "/".join(path_parts + [leaf_key])

        # Show progress info with colors if available
        if HAS_COLORS:
            print(f"\n{Fore.CYAN}[{idx}/{total_count}]{Style.RESET_ALL} {pretty_name}")
            if "description" in meta:
                print(
                    f"{Fore.WHITE}Description:{Style.RESET_ALL} {meta['description']}"
                )
        else:
            print(f"\n[{idx}/{total_count}] {pretty_name}")
            if "description" in meta:
                print(f"Description: {meta['description']}")

        downloaded.append(fetch(url, output_name, dest_dir))

    return downloaded
