#!/usr/bin/env python3
"""
dist_to_nj.py - Convert distance matrix to NJ tree (Newick)

Supports multiple input formats:
  - phylip: PHYLIP lower-triangular distance matrix
  - csv: CSV distance matrix (first row = headers)
  - tsv: TSV distance matrix

Usage:
    python3 dist_to_nj.py <input_dist> <output.newick> <format>

Requires: dendropy
"""

import sys
import os
import csv
import io

def read_phylip_distance(filepath):
    """Read PHYLIP lower-triangular distance matrix."""
    import dendropy

    with open(filepath) as f:
        lines = [l.strip() for l in f if l.strip()]

    if not lines:
        return None, None

    # First line is count (or first taxon)
    try:
        n = int(lines[0])
        data_lines = lines[1:]
    except ValueError:
        data_lines = lines
        n = len(data_lines)

    taxa = []
    matrix = []

    for line in data_lines:
        parts = line.split()
        if not parts:
            continue
        taxa.append(parts[0])
        vals = [float(x) for x in parts[1:]]
        matrix.append(vals)

    # Handle both lower-triangular and full matrix
    # Build full symmetric matrix
    full_matrix = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(len(matrix[i])):
            full_matrix[i][j] = matrix[i][j]
            if j < n:
                full_matrix[j][i] = matrix[i][j]

    return taxa, full_matrix


def read_csv_distance(filepath, delimiter=','):
    """Read CSV/TSV distance matrix."""
    with open(filepath) as f:
        reader = csv.reader(f, delimiter=delimiter)
        header = next(reader)
        # Remove empty strings
        taxa = [h.strip() for h in header if h.strip()]

        # If first column is taxa names, skip it
        matrix = []
        for row in reader:
            if not row:
                continue
            # Determine if first element is a name or number
            try:
                float(row[0])
                start_idx = 0
            except ValueError:
                start_idx = 1

            vals = []
            for v in row[start_idx:]:
                try:
                    vals.append(float(v.strip()))
                except (ValueError, IndexError):
                    vals.append(0.0)
            matrix.append(vals)

    n = len(taxa)
    # Ensure matrix is square
    full_matrix = [[0.0]*n for _ in range(n)]
    for i in range(min(n, len(matrix))):
        for j in range(min(n, len(matrix[i]))):
            full_matrix[i][j] = matrix[i][j]
            full_matrix[j][i] = matrix[i][j]

    return taxa, full_matrix


def build_nj_tree(taxa, distance_matrix):
    """Build NJ tree from distance matrix using dendropy or biopython."""
    # Try dendropy first
    try:
        import dendropy

        n = len(taxa)
        # Build a PHYLIP-format string in memory
        lines = [str(n)]
        for i in range(n):
            row_vals = " ".join(f"{distance_matrix[i][j]:.10f}" for j in range(n))
            lines.append(f"{taxa[i]}  {row_vals}")
        phylip_str = "\n".join(lines)

        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
            src=io.StringIO(phylip_str),
            delimiter=" ",
            is_allow_new_taxa=True
        )
        nj_tree = pdm.nj_tree()
        return nj_tree
    except Exception as e1:
        pass

    # Fallback: try biopython NJ
    try:
        from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
        from Bio import Phylo as BioPhylo

        dm_list = []
        for i in range(len(taxa)):
            row = [distance_matrix[i][j] for j in range(i + 1)]
            dm_list.append(row)

        dm = DistanceMatrix(taxa, dm_list)
        constructor = DistanceTreeConstructor()
        nj_tree = constructor.nj(dm)
        return nj_tree
    except Exception as e2:
        raise RuntimeError(f"Both dendropy ({e1}) and biopython ({e2}) NJ failed")


def main():
    if len(sys.argv) < 4:
        print("Usage: python3 dist_to_nj.py <input_dist> <output.newick> <format>")
        print("  format: phylip, csv, tsv")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    fmt = sys.argv[3].lower()

    if not os.path.exists(input_file):
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)

    # Read distance matrix
    if fmt == "phylip":
        taxa, matrix = read_phylip_distance(input_file)
    elif fmt == "csv":
        taxa, matrix = read_csv_distance(input_file, delimiter=',')
    elif fmt == "tsv":
        taxa, matrix = read_csv_distance(input_file, delimiter='\t')
    else:
        print(f"Error: Unknown format '{fmt}'. Use: phylip, csv, tsv")
        sys.exit(1)

    if taxa is None or len(taxa) < 3:
        print(f"Error: Need at least 3 taxa, got {len(taxa) if taxa else 0}")
        sys.exit(1)

    # Build NJ tree
    try:
        tree = build_nj_tree(taxa, matrix)
        # Write output — handle both dendropy and biopython tree objects
        if hasattr(tree, 'write'):
            # dendropy tree
            tree.write(path=output_file, schema="newick")
        else:
            # biopython tree
            from Bio import Phylo
            Phylo.write(tree, output_file, "newick")
        print(f"  NJ tree saved: {output_file} ({len(taxa)} taxa)")
    except Exception as e:
        print(f"  Error building NJ tree: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
