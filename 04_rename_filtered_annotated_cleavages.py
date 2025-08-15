#!/usr/bin/env python3
import csv
import argparse

def parse_fasta(fasta_file):
    """Simple FASTA parser without biopython"""
    name_map = {}
    current_id = None
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                full_id = line[1:].strip()
                parts = full_id.split('_')
                base_id = '_'.join(parts[-2:])  # Gets ChrXX_XXXXX part
                name_map[base_id] = full_id
    return name_map

def rename_miRNAs(table_file, fasta_file, output_file):
    miRNA_name_map = parse_fasta(fasta_file)
    
    with open(table_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
        writer.writeheader()
        
        for row in reader:
            original_name = row['miRNA']
            row['miRNA'] = miRNA_name_map.get(original_name, f'novel404_{original_name}')
            writer.writerow(row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rename miRNA IDs in a table using FASTA headers.")
    parser.add_argument("--table", required=True, help="Input CSV table with miRNA names to rename")
    parser.add_argument("--fasta", required=True, help="FASTA file containing miRNA IDs")
    parser.add_argument("--output", default="renamed_table.csv", help="Output CSV file")
    args = parser.parse_args()
    
    rename_miRNAs(args.table, args.fasta, args.output)
    print(f"Completed! Output saved to {args.output}")
