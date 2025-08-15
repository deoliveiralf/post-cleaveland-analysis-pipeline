#!/usr/bin/env python3

##USAGE: python filter-cleaveland-results.py <input_file> --pdf_dir <pdf_directory> [options]
#*Example: python filter-cleaveland-results.py full_results.txt --pdf_dir ./pdf_files --min-mfe-ratio 0.7 --max-allen-score 5 --max-pvalue 0.05 --category 0 1 2

# Required Arguments:

# input_file: Path to the original CleaveLand full_results.txt file
# --pdf_dir: Directory containing PDF files to be filtered

# Optional Arguments:

# --min-mfe-ratio: Minimum MFE ratio to keep (e.g., 0.7)
# --max-allen-score: Maximum Allen et al. score to keep (e.g., 5)
# --category: Degradome categories to keep (e.g., 0 1 2)
# --max-pvalue: Maximum p-value to keep (e.g., 0.05)
# --output_dir: Directory to save text files (default: current directory)
# --pdf_output_dir: Directory to copy matching PDFs to (default: matched_pdfs)

import os
import re
import shutil
import argparse
from pathlib import Path

def filter_cleaveland_output(input_file, output_file, min_mfe_ratio=None, max_allen_score=None,
                             category=None, max_pvalue=None):
    """
    Filter CleaveLand output based on specified criteria.

    This function integrates functionality from filter_cleaveland.py
    """
    def parse_cleaveland_output(file_path):
        """Parse the CleaveLand output file and return a list of records."""
        records = []
        current_record = {}

        with open(file_path, 'r') as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            line = lines[i].strip()

            # Start of a new record
            if line.startswith("SiteID:"):
                if current_record:
                    records.append(current_record)
                    current_record = {}

                # Parse site ID
                current_record['site_id'] = line.split("SiteID:")[1].strip()

            # Parse MFE values
            elif line.startswith("MFE of perfect match:"):
                current_record['mfe_perfect'] = float(line.split(":")[1].strip())
            elif line.startswith("MFE of this site:"):
                current_record['mfe_site'] = float(line.split(":")[1].strip())
            elif line.startswith("MFEratio:"):
                current_record['mfe_ratio'] = float(line.split(":")[1].strip())

            # Parse Allen score
            elif line.startswith("Allen et al. score:"):
                current_record['allen_score'] = float(line.split(":")[1].strip())

            # Parse paired regions
            elif line.startswith("Paired Regions"):
                current_record['paired_regions'] = []
                i += 1
                while i < len(lines) and not lines[i].strip().startswith("Unpaired Regions"):
                    region = lines[i].strip()
                    if region and not region.startswith("Paired Regions"):
                        current_record['paired_regions'].append(region)
                    i += 1
                continue

            # Parse unpaired regions
            elif line.startswith("Unpaired Regions"):
                current_record['unpaired_regions'] = []
                i += 1
                while i < len(lines) and not (
                    lines[i].strip().startswith("Degradome") or
                    lines[i].strip().startswith("Degardome")
                ):
                    region = lines[i].strip()
                    if region and not region.startswith("Unpaired Regions"):
                        current_record['unpaired_regions'].append(region)
                    i += 1
                continue

            # Parse degradome data
            elif line.startswith("Degradome data file:") or line.startswith("Degardome data file:"):
                current_record['degradome_file'] = line.split(":")[1].strip()
            elif line.startswith("Degradome Category:") or line.startswith("Degardome Category:"):
                try:
                    current_record['degradome_category'] = int(line.split(":")[1].strip())
                except ValueError:
                    # In case category isn't an integer
                    current_record['degradome_category'] = line.split(":")[1].strip()
            elif line.startswith("Degradome p-value:") or line.startswith("Degardome p-value:"):
                current_record['degradome_pvalue'] = float(line.split(":")[1].strip())
            elif line.startswith("T-Plot file:"):
                current_record['tplot_file'] = line.split(":")[1].strip()

            # Parse position data
            elif re.match(r"^\d+\s+\d+\s+\d+", line):
                if 'positions' not in current_record:
                    current_record['positions'] = []

                parts = line.split()
                if len(parts) >= 3:
                    pos_data = {
                        'position': int(parts[0]),
                        'reads': int(parts[1]),
                        'category': int(parts[2])
                    }
                    current_record['positions'].append(pos_data)

            i += 1

        # Add the last record
        if current_record:
            records.append(current_record)

        return records


    def filter_records(records):
        """Filter records based on the provided criteria."""
        filtered = []

        for record in records:
            # Check if all required fields exist for filtering
            pass_filter = True

            # Apply MFE ratio filter if specified
            if min_mfe_ratio is not None:
                if 'mfe_ratio' not in record or record['mfe_ratio'] < min_mfe_ratio:
                    pass_filter = False

            # Apply Allen score filter if specified
            if max_allen_score is not None and pass_filter:
                if 'allen_score' not in record or record['allen_score'] > max_allen_score:
                    pass_filter = False

            # Apply category filter if specified
            if category is not None and pass_filter:
                if 'degradome_category' not in record or record['degradome_category'] not in category:
                    pass_filter = False

            # Apply p-value filter if specified
            if max_pvalue is not None and pass_filter:
                if 'degradome_pvalue' not in record or record['degradome_pvalue'] > max_pvalue:
                    pass_filter = False

            if pass_filter:
                filtered.append(record)

        return filtered

    def write_output(records, output_file):
        """Write the filtered records to the output file."""
        with open(output_file, 'w') as f:
            for record in records:
                f.write(f"SiteID: {record.get('site_id', 'Unknown')}\n")
                f.write(f"MFE of perfect match: {record.get('mfe_perfect', 'N/A')}\n")
                f.write(f"MFE of this site: {record.get('mfe_site', 'N/A')}\n")
                f.write(f"MFEratio: {record.get('mfe_ratio', 'N/A')}\n")
                f.write(f"Allen et al. score: {record.get('allen_score', 'N/A')}\n")

                # Write paired regions
                f.write("Paired Regions\n")
                for region in record.get('paired_regions', []):
                    f.write(f"    {region}\n")

                # Write unpaired regions
                f.write("Unpaired Regions\n")
                for region in record.get('unpaired_regions', []):
                    f.write(f"    {region}\n")

                # Write degradome data
                f.write(f"Degradome data file: {record.get('degradome_file', 'N/A')}\n")
                f.write(f"Degradome Category: {record.get('degradome_category', 'N/A')}\n")
                f.write(f"Degradome p-value: {record.get('degradome_pvalue', 'N/A')}\n")
                if 'tplot_file' in record:
                    f.write(f"T-Plot file: {record.get('tplot_file', 'N/A')}\n")

                # Write position data
                if 'positions' in record and record['positions']:
                    f.write("\nPosition\tReads\tCategory\n")
                    for pos in record['positions']:
                        f.write(f"{pos['position']}\t{pos['reads']}\t{pos['category']}\n")

                f.write("\n" + "-"*50 + "\n\n")

    # Check if input file exists
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' does not exist")

    # Parse the input file
    records = parse_cleaveland_output(input_file)
    print(f"Parsed {len(records)} records from the input file")

    # Filter the records
    filtered_records = filter_records(records)
    print(f"Filtered to {len(filtered_records)} records")

    # Write the filtered records to the output file
    write_output(filtered_records, output_file)
    print(f"Filtered results written to '{output_file}'")

    return len(filtered_records)

def extract_ids_from_cleaveland(file_path):
    """
    Extract all SiteID values from CleaveLand output file.
    Remove any characters from ':' and afterward.

    This function integrates functionality from extract_ids.py
    """
    ids = []

    with open(file_path, 'r') as file:
        content = file.read()

        # Split the content by potential record separators to process each record
        records = re.split(r'\n\n+', content)

        for record in records:
            # Look for SiteID pattern
            match = re.search(r'SiteID:\s*(\S+)', record)
            if match:
                full_id = match.group(1)
                # Remove characters from ':' and afterward
                clean_id = full_id.split(':')[0]
                ids.append(clean_id)

    return ids

def filter_and_copy_pdfs(ids, pdf_dir, pdf_output_dir):
    """
    Filter PDF files based on IDs and copy to output directory.

    This function integrates functionality from filter_copy_pdfs.py

    Args:
        ids: List of extracted IDs to match
        pdf_dir: Directory containing PDF files
        pdf_output_dir: Directory to copy matching PDFs to

    Returns:
        List of copied files
    """
    # Create output directory if it doesn't exist
    os.makedirs(pdf_output_dir, exist_ok=True)

    copied_files = []

    # Get all PDF files in the directory
    pdf_files = [f for f in os.listdir(pdf_dir) if f.endswith('.pdf')]
    print(f"Found {len(pdf_files)} PDF files in {pdf_dir}")

    for pdf_file in pdf_files:
        # Extract the Sevir ID from the PDF filename
        # Pattern looks for: Sevir.XXXXXXXX.X in the filename
        match = re.search(r'Sevir\.[\w\d]+\.\d+', pdf_file)
        if match:
            pdf_id = match.group(0)

            # Check if this ID is in our list of extracted IDs
            if pdf_id in ids:
                src_path = os.path.join(pdf_dir, pdf_file)
                dst_path = os.path.join(pdf_output_dir, pdf_file)

                shutil.copy2(src_path, dst_path)
                copied_files.append(pdf_file)
                print(f"Copied: {pdf_file}")

    return copied_files

def main():
    parser = argparse.ArgumentParser(description='Integrated CleaveLand workflow: filter results, extract IDs, and copy PDFs')

    # Filter CleaveLand arguments
    parser.add_argument('input_file', help='Path to the CleaveLand full_results.txt file')
    parser.add_argument('--min-mfe-ratio', type=float, help='Minimum MFE ratio to keep')
    parser.add_argument('--max-allen-score', type=float, help='Maximum Allen et al. score to keep')
    parser.add_argument('--category', type=int, nargs='+', help='Degradome categories to keep (e.g., 0 1 2)')
    parser.add_argument('--max-pvalue', type=float, help='Maximum p-value to keep')

    # Output directory arguments
    parser.add_argument('--pdf_dir', required=True, help='Directory containing PDF files')
    parser.add_argument('--output_dir', default='.', help='Directory to save text files (default: current directory)')
    parser.add_argument('--pdf_output_dir', default='matched_pdfs', help='Directory to copy matching PDFs to (default: matched_pdfs)')

    args = parser.parse_args()

    try:
        # Create output directory for text files if it doesn't exist
        os.makedirs(args.output_dir, exist_ok=True)

        # Step 1: Filter CleaveLand results
        filtered_output_file = os.path.join(args.output_dir, "filtered_cleaveland_results.txt")
        records_count = filter_cleaveland_output(
            args.input_file,
            filtered_output_file,
            min_mfe_ratio=args.min_mfe_ratio,
            max_allen_score=args.max_allen_score,
            category=args.category,
            max_pvalue=args.max_pvalue
        )

        # Step 2: Extract IDs from filtered results
        print("\nExtracting IDs from filtered results...")
        extracted_ids = extract_ids_from_cleaveland(filtered_output_file)

        # Save extracted IDs
        ids_file = os.path.join(args.output_dir, "extracted_ids.txt")
        with open(ids_file, "w") as outfile:
            for id_value in extracted_ids:
                outfile.write(f"{id_value}\n")
        print(f"Extracted {len(extracted_ids)} IDs saved to {ids_file}")

        # Step 3: Copy matching PDFs
        print("\nCopying matching PDF files...")
        copied_files = filter_and_copy_pdfs(extracted_ids, args.pdf_dir, args.pdf_output_dir)

        # Save list of copied files
        copied_files_path = os.path.join(args.output_dir, "copied_files.txt")
        with open(copied_files_path, "w") as outfile:
            for file in copied_files:
                outfile.write(f"{file}\n")

        # Print summary
        print(f"\nWorkflow complete!")
        print(f"- Filtered {records_count} CleaveLand records")
        print(f"- Extracted {len(extracted_ids)} unique IDs")
        print(f"- Text files saved to: {args.output_dir}/")
        print(f"- Copied {len(copied_files)} matching PDF files to: {args.pdf_output_dir}/")
        print(f"- Summary saved to {copied_files_path}")

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
