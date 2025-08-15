#!/usr/bin/env python3

# This script adds miRNA information to a table based on a fasta file
import sys

def parse_fasta(fasta_file):
    """Parse a fasta file and return a dictionary mapping chromosome IDs to full miRNA names and miRNA types"""
    mapping = {}
    current_name = None

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Process header line
                full_name = line[1:]  # Remove '>' character

                # Extract the chromosome ID (part after first '_')
                parts = full_name.split('_', 1)
                if len(parts) > 1:
                    miRNA_type = parts[0]  # e.g., miR166 or novel1
                    chr_id = parts[1]      # e.g., Chr09_45260

                    # Store the mapping
                    mapping[chr_id] = {
                        'full_name': full_name,
                        'miRNA_type': miRNA_type
                    }

    return mapping

def process_table(table_file, output_file, miRNA_mapping):
    """Process the table file and add the new columns"""
    with open(table_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Write header
        outfile.write("Original miRNA ID\tmiRNA_Chr ID\tmiRNA ID\tTarget ID\n")

        # Skip the original header line
        next(infile)

        # Process data lines
        for line in infile:
            parts = line.strip().split()
            if len(parts) >= 2:
                chr_id = parts[0]  # e.g., Chr09_45260
                target_id = parts[1]  # e.g., Sevir.1G015900.1

                # Look up the miRNA information
                if chr_id in miRNA_mapping:
                    full_name = miRNA_mapping[chr_id]['full_name']
                    miRNA_type = miRNA_mapping[chr_id]['miRNA_type']
                else:
                    # If no mapping found, use placeholders
                    full_name = "Unknown"
                    miRNA_type = "Unknown"

                # Write the new line
                outfile.write(f"{chr_id}\t{full_name}\t{miRNA_type}\t{target_id}\n")

# Main execution
def main():
    import sys

    # Check if correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python miRNA_mapping_script.py miRNA.fasta miRNA_table.txt > miRNA_output.txt")
        sys.exit(1)

    fasta_file = sys.argv[1]
    table_file = sys.argv[2]

    # Parse the fasta file
    miRNA_mapping = parse_fasta(fasta_file)

    # Process the table and output to stdout
    process_table_to_stdout(table_file, miRNA_mapping)

def process_table_to_stdout(table_file, miRNA_mapping):
    """Process the table file and output to stdout"""
    import sys

    # Write header
    sys.stdout.write("Original miRNA ID\tmiRNA_Chr ID\tmiRNA ID\tTarget ID\tGene\n")

    with open(table_file, 'r') as infile:
        # Skip the original header line
        next(infile)

        # Process data lines
        for line in infile:
            parts = line.strip().split()
            if len(parts) >= 2:
                chr_id = parts[0]  # e.g., Chr09_45260
                target_id = parts[1]  # e.g., Sevir.1G015900.1

                # Extract gene name (everything before the last dot)
                gene_parts = target_id.rsplit('.', 1)
                gene = gene_parts[0] if len(gene_parts) > 1 else target_id

                # Look up the miRNA information
                if chr_id in miRNA_mapping:
                    full_name = miRNA_mapping[chr_id]['full_name']
                    miRNA_type = miRNA_mapping[chr_id]['miRNA_type']
                else:
                    # If no mapping found, use placeholders
                    full_name = "Unknown"
                    miRNA_type = "Unknown"

                # Write the new line to stdout
                sys.stdout.write(f"{chr_id}\t{full_name}\t{miRNA_type}\t{target_id}\t{gene}\n")

if __name__ == "__main__":
    main()
