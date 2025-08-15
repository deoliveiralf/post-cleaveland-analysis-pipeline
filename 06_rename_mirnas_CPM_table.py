#!/usr/bin/env python3
"""
Script to rename miRNAs in CPM expression table based on FASTA file identifiers.
"""
#USAGE: python rename_mirnas.py your_fasta_file.fa your_cpm_table.txt -o renamed_cpm_table.txt

import pandas as pd
import argparse
import sys
from pathlib import Path

def parse_fasta(fasta_file):
    """
    Parse FASTA file and create a mapping dictionary.
    Returns a dictionary mapping sequence IDs to their headers.
    """
    id_mapping = {}
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Extract the sequence ID (everything after '>')
                seq_id = line[1:]
                # The mapping will be from seq_id to seq_id (same name)
                # But we can extract meaningful names if needed
                id_mapping[seq_id] = seq_id
    
    return id_mapping

def create_renaming_map(fasta_mapping, cpm_gene_ids):
    """
    Create a renaming map based on FASTA IDs and CPM table gene IDs.
    This function tries to match CPM gene IDs with FASTA sequence IDs.
    """
    rename_map = {}
    
    # Get all FASTA sequence IDs
    fasta_ids = set(fasta_mapping.keys())
    
    for gene_id in cpm_gene_ids:
        # Check if the gene_id exists in FASTA
        if gene_id in fasta_ids:
            rename_map[gene_id] = gene_id  # Keep the same name if found
        else:
            # Try to find partial matches or similar patterns
            # This handles cases where naming might be slightly different
            found_match = False
            
            # Extract the base name (before the first underscore if present)
            if '_' in gene_id:
                base_name = gene_id.split('_')[0]
                
                # Look for matches that start with the base name
                for fasta_id in fasta_ids:
                    if fasta_id.startswith(base_name):
                        rename_map[gene_id] = fasta_id
                        found_match = True
                        break
            
            # If no match found, keep original name
            if not found_match:
                rename_map[gene_id] = gene_id
                print(f"Warning: No FASTA match found for {gene_id}")
    
    return rename_map

def rename_cpm_table(cpm_file, fasta_file, output_file=None):
    """
    Main function to rename miRNAs in CPM table based on FASTA file.
    """
    
    # Parse FASTA file
    print("Parsing FASTA file...")
    fasta_mapping = parse_fasta(fasta_file)
    print(f"Found {len(fasta_mapping)} sequences in FASTA file")
    
    # Read CPM table
    print("Reading CPM expression table...")
    cpm_df = pd.read_csv(cpm_file, sep='\t', index_col=0)
    print(f"Found {len(cpm_df)} genes in CPM table")
    
    # Create renaming map
    print("Creating renaming map...")
    rename_map = create_renaming_map(fasta_mapping, cpm_df.index.tolist())
    
    # Count how many will be renamed
    renamed_count = sum(1 for old, new in rename_map.items() if old != new)
    print(f"Will rename {renamed_count} out of {len(rename_map)} gene IDs")
    
    # Apply renaming
    cpm_df_renamed = cpm_df.copy()
    cpm_df_renamed.index = cpm_df_renamed.index.map(rename_map)
    
    # Save results
    if output_file is None:
        output_file = cpm_file.replace('.txt', '_renamed.txt').replace('.tsv', '_renamed.tsv')
        if output_file == cpm_file:  # If no extension was replaced
            output_file = cpm_file + '_renamed'
    
    print(f"Saving renamed CPM table to {output_file}")
    cpm_df_renamed.to_csv(output_file, sep='\t')
    
    # Print summary
    print("\nRenaming Summary:")
    print("-" * 50)
    changes_made = [(old, new) for old, new in rename_map.items() if old != new]
    
    if changes_made:
        print(f"Successfully renamed {len(changes_made)} gene IDs:")
        for old, new in changes_made[:10]:  # Show first 10 changes
            print(f"  {old} → {new}")
        if len(changes_made) > 10:
            print(f"  ... and {len(changes_made) - 10} more")
    else:
        print("No gene IDs were renamed (all names already match)")
    
    return cpm_df_renamed, rename_map

def main():
    parser = argparse.ArgumentParser(
        description="Rename miRNAs in CPM expression table based on FASTA file identifiers"
    )
    parser.add_argument("fasta_file", help="Input FASTA file with miRNA sequences")
    parser.add_argument("cpm_file", help="Input CPM expression table (tab-separated)")
    parser.add_argument("-o", "--output", help="Output file name (optional)")
    
    args = parser.parse_args()
    
    # Check if files exist
    if not Path(args.fasta_file).exists():
        print(f"Error: FASTA file '{args.fasta_file}' not found")
        sys.exit(1)
    
    if not Path(args.cpm_file).exists():
        print(f"Error: CPM file '{args.cpm_file}' not found")
        sys.exit(1)
    
    try:
        renamed_df, rename_map = rename_cpm_table(args.cpm_file, args.fasta_file, args.output)
        print("\nScript completed successfully!")
        
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
