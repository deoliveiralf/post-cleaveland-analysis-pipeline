#!/usr/bin/env python3

import csv
import sys
import re

def strip_version(gene_id):
    # Removes last .NUMBER if present, otherwise returns as is
    if re.match(r".+\.\d+$", gene_id):
        return gene_id.rsplit('.', 1)[0]
    else:
        return gene_id

def process_files(mirna_file, annotation_file, output_file):
    mirna_data = []

    with open(mirna_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        print("miRNA file headers:", header)

        mirna_id_index = None
        gene_index = None

        for i, col in enumerate(header):
            colname = col.strip().lower()
            if colname == "mirna id":
                mirna_id_index = i
            elif colname == "gene":
                gene_index = i

        if mirna_id_index is None or gene_index is None:
            print("Error: Could not find 'miRNA ID' or 'Gene' columns in miRNA file")
            print("Available columns:", header)
            return

        for row in reader:
            if len(row) > max(mirna_id_index, gene_index):
                mirna_id = row[mirna_id_index].strip()
                gene = row[gene_index].strip()
                gene_base = strip_version(gene)
                print(f"DEBUG: original_gene='{gene}', gene_base='{gene_base}'")  # Debug print
                mirna_data.append((mirna_id, gene, gene_base))

    print(f"Loaded {len(mirna_data)} miRNA entries")
    print("Sample gene IDs from miRNA file:", [g[2] for g in mirna_data[:5]])

    # Read annotation data
    annotation_data = {}

    with open(annotation_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        print("Annotation file headers:", header)

        setaria_id_index = None
        function_index = None

        for i, col in enumerate(header):
            colname = col.strip().lower()
            if colname == "setaria viridis id":
                setaria_id_index = i
            elif colname == "function":
                function_index = i

        if setaria_id_index is None or function_index is None:
            print("Error: Could not find 'Setaria viridis ID' or 'Function' columns in annotation file")
            print("Available columns:", header)
            return

        for row in reader:
            if len(row) > max(setaria_id_index, function_index):
                setaria_id = row[setaria_id_index].strip()
                setaria_id_base = strip_version(setaria_id)
                function = row[function_index].strip()
                annotation_data[setaria_id_base] = function

    print(f"Loaded {len(annotation_data)} annotation entries")
    print("Sample Setaria viridis IDs from annotation file:", list(annotation_data.keys())[:5])

    # Generate output
    output_rows = []
    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["miRNA ID", "Gene", "Function"])

        for mirna_id, original_gene, gene_base in mirna_data:
            function = annotation_data.get(gene_base, "")
            print(f"Trying to match gene_base '{gene_base}' -> function: '{function}'")
            if not function and '.' in original_gene:
                base_gene = strip_version(original_gene)
                function = annotation_data.get(base_gene, "")
                print(f"Trying base gene: {base_gene}, found function: '{function}'")

            if function:
                writer.writerow([mirna_id, gene_base, function])
                output_rows.append([mirna_id, gene_base, function])

    print(f"Generated {len(output_rows)} output rows")
    return output_rows

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_script.py <mirna_file.tsv> <annotation_file.tsv> <output_file.tsv>")
        sys.exit(1)

    mirna_file = sys.argv[1]
    annotation_file = sys.argv[2]
    output_file = sys.argv[3]

    print(f"Processing files: {mirna_file}, {annotation_file}, {output_file}")
    results = process_files(mirna_file, annotation_file, output_file)

    if not results or len(results) == 0:
        print("Warning: No matching entries found!")
        print("Please check that:")
        print("1. The column names exactly match 'miRNA ID', 'Gene', 'Setaria viridis ID', and 'Function'")
        print("2. The gene IDs in the miRNA file match the Setaria viridis IDs in the annotation file")
        print("3. There are no extra spaces or special characters in the data")
