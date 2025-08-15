#!/usr/bin/env python3
"""
Script to rename miRNA IDs in a table based on matching patterns in a FASTA file.
"""
# USAGE: python rename_mirnas_table.py input_table.txt input_fasta.fa -o renamed_table.txt

import pandas as pd
import argparse
import sys
from pathlib import Path

def parse_fasta(fasta_file):
    """
    Parse FASTA file and extract miRNA base names.
    Returns a dictionary mapping the trailing ID (e.g., Chr01_688) to full name (e.g., miR168_Chr01_688).
    """
    id_mapping = {}
    seq_count = 0
    
    with open(fasta_file, 'r') as f:
        current_header = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                seq_count += 1
                full_name = line[1:]  # Remove '>'
                # Extract the trailing ID (last part after last underscore or similar)
                parts = full_name.split('_')
                if len(parts) >= 2:
                    trailing_id = '_'.join(parts[-2:])  # Gets ChrXX_XXXXX part
                    id_mapping[trailing_id] = full_name
                current_header = full_name
            elif current_header:
                # This is sequence data, we can skip or process if needed
                pass
    
    print(f"Processed {seq_count} sequences from FASTA")
    return id_mapping

def create_renaming_map(fasta_mapping, table_ids):
    """
    Create a renaming map based on FASTA IDs and table miRNA IDs.
    """
    rename_map = {}
    unmatched_ids = []
    
    for gene_id in table_ids:
        # Check if the gene_id exists in FASTA mapping
        if gene_id in fasta_mapping:
            rename_map[gene_id] = fasta_mapping[gene_id]
        else:
            # Try to find partial matches
            found_match = False
            
            # Look for matches where the FASTA key ends with our ID
            for fasta_id, fasta_name in fasta_mapping.items():
                if fasta_id.endswith(gene_id):
                    rename_map[gene_id] = fasta_name
                    found_match = True
                    break
            
            # If no match found, keep original name
            if not found_match:
                rename_map[gene_id] = gene_id
                unmatched_ids.append(gene_id)
    
    if unmatched_ids:
        print(f"\nWarning: No FASTA match found for {len(unmatched_ids)} miRNA IDs. First 10 examples:")
        for id in unmatched_ids[:10]:
            print(f"  {id}")
        if len(unmatched_ids) > 10:
            print(f"  ... and {len(unmatched_ids) - 10} more")
    
    return rename_map

def rename_table(fasta_file, input_file, output_file=None):
    """
    Main function to rename miRNAs in table based on FASTA file.
    """
    
    # Parse FASTA file
    print("\nParsing FASTA file...")
    fasta_mapping = parse_fasta(fasta_file)
    if not fasta_mapping:
        print("Error: No valid miRNA sequences found in FASTA file")
        print("Please check that your FASTA file contains entries like: '>miR168_Chr01_688'")
        sys.exit(1)
    print(f"Found {len(fasta_mapping)} miRNA patterns in FASTA file")
    
    # Read input table
    print("\nReading input table...")
    try:
        df = pd.read_csv(input_file, sep='\t')
        print(f"Found {len(df)} rows in input table")
        
        # Check for required columns
        if 'miRNA ID' not in df.columns:
            print("\nError: Input table must contain a 'miRNA ID' column")
            print("Available columns:", df.columns.tolist())
            sys.exit(1)
            
    except Exception as e:
        print(f"\nError reading input table: {str(e)}")
        sys.exit(1)
    
    # Get unique miRNA IDs from table
    mirna_ids = df['miRNA ID'].unique().tolist()
    print(f"Found {len(mirna_ids)} unique miRNA IDs in table")
    
    # Create renaming map
    print("\nCreating renaming map...")
    rename_map = create_renaming_map(fasta_mapping, mirna_ids)
    
    # Count how many will be renamed
    renamed_count = sum(1 for old, new in rename_map.items() if old != new)
    print(f"\nWill rename {renamed_count} out of {len(rename_map)} miRNA IDs")
    
    # Apply renaming
    df_renamed = df.copy()
    df_renamed['miRNA ID'] = df_renamed['miRNA ID'].map(rename_map)
    
    # Save results
    if output_file is None:
        output_file = input_file.replace('.txt', '_renamed.txt').replace('.tsv', '_renamed.tsv')
        if output_file == input_file:  # If no extension was replaced
            output_file = input_file + '_renamed'
    
    print(f"\nSaving renamed table to {output_file}")
    df_renamed.to_csv(output_file, sep='\t', index=False)
    
    # Print summary
    print("\nRenaming Summary:")
    print("-" * 50)
    changes_made = [(old, new) for old, new in rename_map.items() if old != new]
    
    if changes_made:
        print(f"Successfully renamed {len(changes_made)} miRNA IDs:")
        for old, new in changes_made[:10]:  # Show first 10 changes
            print(f"  {old} → {new}")
        if len(changes_made) > 10:
            print(f"  ... and {len(changes_made) - 10} more")
    else:
        print("No miRNA IDs were renamed (no matches found between table and FASTA)")
    
    return df_renamed, rename_map

def main():
    parser = argparse.ArgumentParser(
        description="Rename miRNA IDs in table based on patterns in FASTA file"
    )
    parser.add_argument("input_file", help="Input table file (tab-separated)")
    parser.add_argument("fasta_file", help="Input FASTA file with miRNA sequences")
    parser.add_argument("-o", "--output", help="Output file name (optional)")
    
    args = parser.parse_args()
    
    # Check if files exist
    if not Path(args.fasta_file).exists():
        print(f"\nError: FASTA file '{args.fasta_file}' not found")
        sys.exit(1)
    
    if not Path(args.input_file).exists():
        print(f"\nError: Input file '{args.input_file}' not found")
        sys.exit(1)
    
    try:
        renamed_df, rename_map = rename_table(args.fasta_file, args.input_file, args.output)
        print("\nScript completed successfully!")
        
    except Exception as e:
        print(f"\nError: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
