#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
import re
import argparse
import logging
from datetime import datetime
from matplotlib import rcParams

# Set font for better visualization
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
rcParams['font.size'] = 11

def setup_logging(output_dir):
    """Set up comprehensive logging system"""
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, "cleavage_plots.log")
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler()
        ]
    )
    logging.info(f"Log file created at: {os.path.abspath(log_file)}")
    return log_file

def extract_mirna_base(name):
    """Extract consistent miRNA identifier from different naming formats"""
    if pd.isna(name):
        return name
    
    name = str(name).strip()
    parts = name.split('_')
    
    # Handle cases like: novel123_Chr01_12345, miR160_Chr01_12345, Chr01_12345
    if name.startswith(('novel', 'miR')):
        if len(parts) >= 3:
            return '_'.join(parts[1:3])  # Gets Chr01_12345
        else:
            return '_'.join(parts[1:])   # Fallback for shorter names
    
    # For names that start with Chr
    if name.startswith('Chr') and len(parts) >= 2:
        return '_'.join(parts[:2])  # Gets Chr01_12345
    
    # Fallback: return original name
    return name

def create_mirna_mapping(filtered_df, renamed_df):
    """Create robust miRNA name mapping with validation"""
    mapping = {}
    
    # Create mapping from original to renamed
    for _, row in renamed_df.iterrows():
        orig_name = row['miRNA']
        new_name = row['miRNA']  # Assuming renamed file has the new names
        
        # Extract base identifiers
        orig_base = extract_mirna_base(orig_name)
        
        # Store mapping
        mapping[orig_base] = new_name
    
    # If no mapping found, create identity mapping
    if not mapping:
        logging.warning("No mapping found in renamed file, using original names")
        for mirna in filtered_df['miRNA'].unique():
            base = extract_mirna_base(mirna)
            mapping[base] = mirna
    
    # Validate mapping
    sample_mappings = list(mapping.items())[:3]
    logging.info(f"Created miRNA name mapping with {len(mapping)} entries")
    logging.info(f"Sample mappings: {sample_mappings}")
    return mapping

def parse_cleaveland_output(cleaveland_file):
    """Enhanced parser for cleaveland output with better sequence extraction"""
    alignments = {}
    line_count = 0
    skipped_lines = 0
    
    try:
        with open(cleaveland_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        logging.error(f"Cleaveland file not found: {cleaveland_file}")
        return {}
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Look for alignment start (line starting with "5'")
        if line.startswith("5'"):
            try:
                # Parse the 5' line: 5' SEQUENCE 3' Transcript: NAME:START-END Slice_Site:POSITION
                parts = line.split()
                if len(parts) < 5:
                    raise ValueError("Insufficient parts in alignment line")
                
                # Extract target sequence (between 5' and 3')
                target_seq = parts[1]
                
                # Find transcript information and slice site
                transcript_name = None
                start_pos = None
                end_pos = None
                slice_site = None
                
                # Parse transcript info from the 5' line
                # Format: Transcript: NAME:START-END Slice_Site:POSITION
                transcript_match = re.search(r'Transcript:\s*([^:]+):(\d+)-(\d+)', line)
                if transcript_match:
                    transcript_name = transcript_match.group(1).strip()
                    start_pos = int(transcript_match.group(2))
                    end_pos = int(transcript_match.group(3))
                
                # Extract slice site
                slice_match = re.search(r'Slice\s*Site:\s*(\d+)', line)
                if slice_match:
                    slice_site = int(slice_match.group(1))
                
                # If we couldn't parse the expected format, try alternative parsing
                if not transcript_name or start_pos is None or end_pos is None or slice_site is None:
                    # Try to extract from the full line using different patterns
                    # Look for patterns like "Sevir.9G545800.1:663-686"
                    alt_match = re.search(r'([^:\s]+):(\d+)-(\d+)', line)
                    if alt_match:
                        transcript_name = alt_match.group(1)
                        start_pos = int(alt_match.group(2))
                        end_pos = int(alt_match.group(3))
                    
                    # If still no slice site, try to find any number after "Slice"
                    if slice_site is None:
                        slice_alt = re.search(r'Slice[^:]*:?\s*(\d+)', line)
                        if slice_alt:
                            slice_site = int(slice_alt.group(1))
                        else:
                            # Last resort: use numbers from the line
                            numbers = re.findall(r'\d+', line)
                            if len(numbers) >= 3:
                                slice_site = int(numbers[-1])  # Assume last number is slice site
                
                # Look for the query line (3' line) to get miRNA info
                mirna_seq = None
                mirna_name = None
                alignment_line = None
                
                # Check the next few lines for the query information
                for j in range(1, min(5, len(lines) - i)):
                    if i + j >= len(lines):
                        break
                    check_line = lines[i + j].strip()
                    
                    # Look for alignment line (usually contains | symbols)
                    if '|' in check_line and not check_line.startswith(('5\'', '3\'')):
                        alignment_line = check_line
                    
                    # Look for query line (starts with 3')
                    elif check_line.startswith('3\''):
                        query_parts = check_line.split()
                        if len(query_parts) >= 2:
                            mirna_seq = query_parts[1]
                        
                        # Extract miRNA name from Query: part
                        query_match = re.search(r'Query:\s*([^\s]+)', check_line)
                        if query_match:
                            mirna_name = query_match.group(1)
                        break
                
                # Validate that we have all necessary information
                if not all([transcript_name, start_pos is not None, end_pos is not None, 
                           slice_site is not None, mirna_name]):
                    logging.debug(f"Missing required info at line {i+1}: "
                                f"transcript={transcript_name}, start={start_pos}, end={end_pos}, "
                                f"slice={slice_site}, mirna={mirna_name}")
                    i += 1
                    continue
                
                # Create miRNA key
                mirna_key = extract_mirna_base(mirna_name)
                
                # Store alignment with all the extracted information
                alignment_key = (mirna_key, transcript_name, slice_site)
                alignments[alignment_key] = {
                    'mirna_seq': mirna_seq,
                    'target_seq': target_seq,
                    'start': start_pos,
                    'end': end_pos,
                    'slice_site': slice_site,
                    'mirna_name': mirna_name,
                    'transcript': transcript_name,
                    'alignment_line': alignment_line,
                    'alignment_length': end_pos - start_pos + 1,
                    'full_5prime_line': line,
                    'full_3prime_line': lines[i + j].strip() if j < len(lines) - i else ""
                }
                
                line_count += 1
                logging.debug(f"Parsed alignment: {mirna_key} -> {transcript_name} "
                            f"({start_pos}-{end_pos}, slice:{slice_site})")
                i += 5  # Skip the processed lines
                
            except Exception as e:
                skipped_lines += 1
                logging.debug(f"Skipped line {i+1}: {str(e)} - {line[:100]}...")
                i += 1
        else:
            i += 1
    
    logging.info(f"Parsed {line_count} valid alignments (skipped {skipped_lines} lines)")
    if line_count > 0:
        sample_keys = list(alignments.keys())[:3]
        logging.info(f"Sample alignment keys: {sample_keys}")
        # Show sample alignment details
        for key in sample_keys[:1]:
            alignment = alignments[key]
            logging.info(f"Sample alignment details: {key} -> "
                        f"start={alignment['start']}, end={alignment['end']}, "
                        f"slice={alignment['slice_site']}")
    
    return alignments

def find_alignment_match(mirna, target, cleavage_pos, alignment_info):
    """Find matching alignment with flexible key matching"""
    
    # Generate possible miRNA variants
    mirna_variants = [
        mirna,
        extract_mirna_base(mirna),
        mirna.replace('miR', '').replace('novel', ''),
        '_'.join(mirna.split('_')[-2:]) if '_' in mirna else mirna,
        '_'.join(mirna.split('_')[1:]) if '_' in mirna else mirna
    ]
    
    # Remove duplicates and empty strings
    mirna_variants = list(set([v for v in mirna_variants if v]))
    
    # Also generate target variants (in case of naming inconsistencies)
    target_variants = [target]
    if '.' in target:
        target_variants.append(target.split('.')[0])
    
    # Try exact matches first
    for mirna_variant in mirna_variants:
        for target_variant in target_variants:
            key = (mirna_variant, target_variant, cleavage_pos)
            if key in alignment_info:
                logging.debug(f"Found exact match with key: {key}")
                return alignment_info[key]
    
    # Try fuzzy matches (ignore cleavage position)
    for mirna_variant in mirna_variants:
        for target_variant in target_variants:
            for key, alignment in alignment_info.items():
                if key[0] == mirna_variant and key[1] == target_variant:
                    logging.debug(f"Found fuzzy match with key: {key}")
                    return alignment
    
    # Try even more flexible matching
    for mirna_variant in mirna_variants:
        for target_variant in target_variants:
            for key, alignment in alignment_info.items():
                if (mirna_variant in key[0] or key[0] in mirna_variant) and \
                   (target_variant in key[1] or key[1] in target_variant):
                    logging.debug(f"Found flexible match with key: {key}")
                    return alignment
    
    logging.warning(f"No alignment found for {mirna} - {target}")
    return None

def plot_cleavage_site(mirna, target, cleavage_data, alignment_info, output_dir):
    """Generate clean cleavage site plot matching the reference layout"""
    try:
        # Get region information for this miRNA-target pair
        regions = cleavage_data[(cleavage_data['miRNA'] == mirna) & 
                              (cleavage_data['Transcript'] == target)]
        
        if regions.empty:
            logging.warning(f"No region data for {mirna} - {target}")
            return False
        
        # Get the first region info
        region_row = regions.iloc[0]
        region_type = region_row.get('Region', 'unknown')
        cleavage_pos = region_row.get('Cleavage_Pos', 'unknown')
        
        # Find alignment for this miRNA-target pair
        alignment = find_alignment_match(mirna, target, cleavage_pos, alignment_info)
        
        # Create figure with clean layout
        fig, ax = plt.subplots(1, 1, figsize=(12, 4))
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        
        # Set up the plot area
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 4)
        ax.axis('off')
        
        # ==================== TITLE SECTION ====================
        # miRNA name
        ax.text(5, 3.5, f"miRNA: {mirna}", 
                ha='center', va='center', fontsize=16, fontweight='bold')
        
        # Transcript name
        ax.text(5, 3.1, f"Transcript: {target}", 
                ha='center', va='center', fontsize=14, fontweight='bold')
        
        # Region info
        region_label = region_type
        if region_type == 'three_prime_UTR':
            region_label = "3' UTR"
        elif region_type == 'five_prime_UTR':
            region_label = "5' UTR"
        elif region_type == '3UTR':
            region_label = "3' UTR"
        elif region_type == '5UTR':
            region_label = "5' UTR"
        
        ax.text(5, 2.7, f"Region: {region_label}", 
                ha='center', va='center', fontsize=12)
        
        # ==================== GENOMIC DIAGRAM ====================
        # Define region colors matching the screenshot
        region_colors = {
            'CDS': '#4F81BD',           # Blue
            'three_prime_UTR': '#F79646',  # Orange
            'five_prime_UTR': '#9BBB59',   # Green
            '3UTR': '#F79646',          # Orange
            '5UTR': '#9BBB59',          # Green
            'intergenic': '#C5504B',    # Red
            'exon': '#4F81BD',          # Blue
            'intron': '#8064A2',        # Purple
            'unknown': '#757575'        # Gray
        }
        
        # Get the color for this region
        region_color = region_colors.get(region_type, '#757575')
        
        # Draw the main genomic region bar
        bar_y = 1.5
        bar_height = 0.3
        bar_start = 1.5
        bar_end = 8.5
        bar_width = bar_end - bar_start
        
        # Create the region rectangle
        rect = patches.Rectangle(
            (bar_start, bar_y - bar_height/2), bar_width, bar_height,
            facecolor=region_color,
            edgecolor='black',
            linewidth=1.5
        )
        ax.add_patch(rect)
        
        # Add region label inside the bar
        ax.text(5, bar_y, region_label, 
                ha='center', va='center', fontsize=11, fontweight='bold', color='white')
        
        # ==================== POSITION MARKERS ====================
        # Get position information from alignment data
        if alignment:
            start_pos = alignment['start']
            end_pos = alignment['end']
            slice_site = alignment['slice_site']
            logging.debug(f"Using alignment positions for {mirna}-{target}: "
                         f"start={start_pos}, end={end_pos}, slice={slice_site}")
        else:
            # If no alignment found, we cannot create accurate plot
            logging.error(f"No alignment data found for {mirna} - {target}, cannot determine positions")
            plt.close()
            return False
        
        # Position markers
        positions = [
            (bar_start, start_pos, 'Start', '#4F81BD'),      # Blue
            (5, slice_site, 'Slice site', '#C5504B'),        # Red
            (bar_end, end_pos, 'End', '#4F81BD')             # Blue
        ]
        
        # Draw position markers
        for x_pos, pos_value, label, color in positions:
            # Draw vertical line
            ax.plot([x_pos, x_pos], [bar_y - bar_height/2 - 0.1, bar_y + bar_height/2 + 0.1], 
                   color=color, linewidth=2)
            
            # Add position number below
            ax.text(x_pos, bar_y - bar_height/2 - 0.25, str(pos_value), 
                   ha='center', va='top', fontsize=10, fontweight='bold', color=color)
            
            # Add label below position number
            ax.text(x_pos, bar_y - bar_height/2 - 0.45, label, 
                   ha='center', va='top', fontsize=10, fontweight='bold', 
                   color=color, style='italic')
        
        # Special highlighting for slice site
        # Draw a thicker red line for the cleavage site
        ax.plot([5, 5], [bar_y - bar_height/2 - 0.15, bar_y + bar_height/2 + 0.15], 
               color='#C5504B', linewidth=4)
        
        # Add horizontal line extensions
        line_y = bar_y - bar_height/2 - 0.1
        ax.plot([0.5, bar_start], [line_y, line_y], 'k-', linewidth=1)
        ax.plot([bar_end, 9.5], [line_y, line_y], 'k-', linewidth=1)
        
        # Adjust layout and save
        plt.tight_layout()
        
        # Save figure
        safe_mirna = re.sub(r'[^\w\-_\.]', '_', mirna)
        safe_target = re.sub(r'[^\w\-_\.]', '_', target)
        output_file = os.path.join(output_dir, f"{safe_mirna}_{safe_target}_cleavage.png")
        
        plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        logging.info(f"Created plot: {output_file} (start={start_pos}, end={end_pos}, slice={slice_site})")
        return True
        
    except Exception as e:
        logging.error(f"Error plotting {mirna}/{target}: {str(e)}", exc_info=True)
        plt.close()
        return False

def main():
    parser = argparse.ArgumentParser(
        description='Generate clean PNG plots of miRNA cleavage sites with region annotations',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--filtered', required=True,
                        help='Input filtered_annotated_cleavages.csv file')
    parser.add_argument('-r', '--renamed', 
                        help='Input renamed_filtered_annotated_cleavages.csv file (optional)')
    parser.add_argument('-c', '--cleaveland', required=True,
                        help='Input cleaveland_output.txt file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory for PNG files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose debug logging')
    
    args = parser.parse_args()
    
    # Setup logging
    log_file = setup_logging(args.output)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    logging.info(f"Starting processing with arguments: {vars(args)}")
    
    try:
        # Load and validate input files
        logging.info("Loading input files...")
        
        # Check if files exist
        if not os.path.exists(args.filtered):
            raise FileNotFoundError(f"Filtered file not found: {args.filtered}")
        if not os.path.exists(args.cleaveland):
            raise FileNotFoundError(f"Cleaveland file not found: {args.cleaveland}")
        
        filtered_df = pd.read_csv(args.filtered)
        logging.info(f"Loaded {len(filtered_df)} cleavage records from {args.filtered}")
        logging.info(f"Columns in filtered data: {list(filtered_df.columns)}")
        
        # Load renamed file if provided
        if args.renamed and os.path.exists(args.renamed):
            renamed_df = pd.read_csv(args.renamed)
            logging.info(f"Loaded renamed data from {args.renamed}")
        else:
            logging.info("No renamed file provided or file not found, using original names")
            renamed_df = filtered_df.copy()
        
        # Create miRNA mapping
        mirna_mapping = create_mirna_mapping(filtered_df, renamed_df)
        
        # Apply mapping to filtered data
        filtered_df['original_mirna'] = filtered_df['miRNA'].copy()
        filtered_df['mirna_base'] = filtered_df['miRNA'].apply(extract_mirna_base)
        filtered_df['miRNA'] = filtered_df['mirna_base'].map(mirna_mapping).fillna(filtered_df['miRNA'])
        
        logging.info(f"After mapping, {filtered_df['miRNA'].nunique()} unique miRNAs")
        
        # Parse alignments
        logging.info(f"Parsing cleaveland output from {args.cleaveland}...")
        alignment_info = parse_cleaveland_output(args.cleaveland)
        if not alignment_info:
            logging.error("No valid alignments found in cleaveland output!")
            logging.error("Cannot generate plots without position information from cleaveland file.")
            return
        
        # Process each miRNA-target pair
        grouped = filtered_df.groupby(['miRNA', 'Transcript'])
        total_pairs = len(grouped)
        logging.info(f"Found {total_pairs} miRNA-target pairs to process")
        
        successful_plots = 0
        failed_plots = 0
        
        for i, ((mirna, target), group) in enumerate(grouped, 1):
            logging.debug(f"Processing pair {i}/{total_pairs}: {mirna} - {target}")
            
            if plot_cleavage_site(mirna, target, filtered_df, alignment_info, args.output):
                successful_plots += 1
            else:
                failed_plots += 1
            
            # Progress update every 10 plots
            if i % 10 == 0:
                logging.info(f"Progress: {i}/{total_pairs} pairs processed")
        
        logging.info(f"Processing complete!")
        logging.info(f"Successfully created: {successful_plots} plots")
        logging.info(f"Failed to create: {failed_plots} plots")
        
        if failed_plots > 0:
            logging.warning(f"Some plots failed - check that miRNA/transcript names match between "
                          f"CSV and cleaveland files, and that cleaveland output format is correct")
        
        logging.info(f"Log file: {os.path.abspath(log_file)}")
        logging.info(f"Output directory: {os.path.abspath(args.output)}")
        
    except Exception as e:
        logging.critical(f"Fatal error: {str(e)}", exc_info=True)
        raise

if __name__ == "__main__":
    main()
