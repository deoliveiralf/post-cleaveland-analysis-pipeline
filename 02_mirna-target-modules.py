#!/usr/bin/env python3

#USAGE: python mirna-target-modules.py filtered_results.txt > mirna-target-modules-table.txt

import sys

def extract_and_print_ids(file_path):
    """
    Extract miRNA ID and Target ID from filtered_results.txt file
    based on the T-Plot file lines and print results.
    """
    # Print the header first
    print("miRNA ID\tTarget ID")
    found_count = 0

    try:
        with open(file_path, 'r') as file:
            for line in file:
                if "T-Plot file:" in line:
                    # Extract the filename part from the path
                    # Example: T-Plot file: cleaveland_results/Chr09_39038_Sevir.1G015900.1_1032_TPlot.pdf
                    filename = line.strip().split('/')[-1]

                    # Split by underscore to get components
                    parts = filename.split('_')

                    if len(parts) >= 4:
                        # Extract miRNA ID (first two parts) and Target ID (third part)
                        mirna_id = f"{parts[0]}_{parts[1]}"
                        target_id = parts[2]

                        # Print as tab-separated values
                        print(f"{mirna_id}\t{target_id}")
                        found_count += 1

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)

    # Print summary to stderr
    print(f"Processed file: {file_path}", file=sys.stderr)
    print(f"Found {found_count} entries", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python extract_ids.py filtered_results.txt")
        sys.exit(1)

    extract_and_print_ids(sys.argv[1])
