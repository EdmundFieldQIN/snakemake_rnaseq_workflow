#!/usr/bin/env python3
"""
Add a 'gene' column header to the first column of tab-separated files
when it is unnamed. Detects missing header by comparing field counts
between header and first data row, not by string matching.
"""

import sys
import os


def needs_gene_header(filepath: str) -> bool:
    """Return True if the first column lacks a header name."""
    with open(filepath, "r") as f:
        header = f.readline().rstrip("\n")
        data_line = f.readline().rstrip("\n")

    if not header or not data_line:
        return False

    header_fields = header.split("\t")
    data_fields = data_line.split("\t")

    # If header has one fewer field than data, the first column is unnamed
    return len(header_fields) == len(data_fields) - 1


def add_gene_header(filepath: str) -> bool:
    """Add 'gene' as first column header if missing. Returns True if modified."""
    if not needs_gene_header(filepath):
        return False

    with open(filepath, "r") as f:
        lines = f.readlines()

    if not lines:
        return False

    lines[0] = "gene\t" + lines[0]
    with open(filepath, "w") as f:
        f.writelines(lines)
    return True


def main():
    if len(sys.argv) < 2:
        sys.exit(0)

    for filepath in sys.argv[1:]:
        if os.path.isfile(filepath):
            modified = add_gene_header(filepath)
            if modified:
                print(f"Added gene header: {filepath}")


if __name__ == "__main__":
    main()
