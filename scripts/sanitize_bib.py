#!/usr/bin/env python3
"""Create a BibTeX-safe copy of a .bib file by replacing non-Latin-1 characters
with their LaTeX equivalents.

BibTeX is a legacy 8-bit tool that cannot handle characters outside Latin-1
(ISO 8859-1).  This script reads a UTF-8 .bib file and writes a copy where
problematic characters are replaced with standard LaTeX accent commands.

Usage:
    python sanitize_bib.py <input.bib> <output.bib>
"""

import sys

# Mapping of non-Latin-1 Unicode characters to LaTeX equivalents.
# Add new entries here as needed.
REPLACEMENTS = {
    "Ž": r"{\v{Z}}",
    "ž": r"{\v{z}}",
    "ś": r"{\'s}",
    "Ś": r"{\'S}",
    "ł": r"{\l}",
    "Ł": r"{\L}",
}


def sanitize(text: str) -> str:
    for char, latex in REPLACEMENTS.items():
        text = text.replace(char, latex)
    return text


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.bib> <output.bib>", file=sys.stderr)
        sys.exit(1)

    src, dst = sys.argv[1], sys.argv[2]
    with open(src, encoding="utf-8") as f:
        content = f.read()

    sanitized = sanitize(content)

    with open(dst, "w", encoding="utf-8") as f:
        f.write(sanitized)


if __name__ == "__main__":
    main()
