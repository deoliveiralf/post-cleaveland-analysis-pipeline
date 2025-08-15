#!/usr/bin/env python3
import sys
import os
import glob
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import Counter

def parse_gff3(gff3_file):
    """Parse GFF3 files with space/tab delimiters"""
    transcripts = {}
    current_transcript = None
    
    with open(gff3_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
                
            parts = line.strip().split()
            if len(parts) < 9:
                continue
                
            feature = parts[2]
            strand = parts[6]
            attrs = ' '.join(parts[8:])
            attr_dict = {}
            
            for attr in attrs.split(';'):
                attr = attr.strip()
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attr_dict[key.strip()] = val.strip()
            
            if feature == 'mRNA' and 'Name' in attr_dict:
                transcript_id = attr_dict['Name']
                transcripts[transcript_id] = {
                    'utr5': [], 'cds': [], 'utr3': [],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'strand': strand,
                    'total_cds': 0
                }
                current_transcript = transcript_id
                
            elif 'Parent' in attr_dict and current_transcript:
                parent = attr_dict['Parent']
                if parent.startswith(current_transcript):
                    start, end = int(parts[3]), int(parts[4])
                    if feature == 'five_prime_UTR':
                        transcripts[current_transcript]['utr5'].append((start, end))
                    elif feature == 'CDS':
                        phase = int(parts[7]) if parts[7] != '.' else 0
                        transcripts[current_transcript]['cds'].append((start, end, phase))
                        transcripts[current_transcript]['total_cds'] += 1
                    elif feature == 'three_prime_UTR':
                        transcripts[current_transcript]['utr3'].append((start, end))
    
    return transcripts

def parse_mirna_targets(targets_file):
    """Parse miRNA to target mappings - handle multiple miRNAs per target"""
    mirna_targets = {}
    with open(targets_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                target = parts[1]
                mirna = parts[0]
                if target not in mirna_targets:
                    mirna_targets[target] = []
                mirna_targets[target].append(mirna)
    return mirna_targets

def parse_cleaveland_batch(cleaveland_file, mirna_targets):
    """Parse cleaveland output file with multiple miRNA binding sites"""
    results = []
    
    # Build reverse mapping for quick lookup: {transcript: [mirnas]}
    transcript_to_mirnas = {}
    for mirna_list in mirna_targets.values():
        for mirna in mirna_list:
            transcript_to_mirnas[mirna] = mirna_list
    
    current_mirna = None
    current_binding_start = None
    current_binding_end = None
    
    with open(cleaveland_file, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Reset miRNA context on new sections or when we see separators
        if line.startswith("Query:") or line.startswith("---") or line.startswith("==="):
            current_mirna = None
            current_binding_start = None
            current_binding_end = None
        
        # Try multiple strategies to identify miRNA
        if line.startswith("Query:"):
            # Extract miRNA name from query line
            try:
                query_parts = line.split()
                if len(query_parts) > 1:
                    current_mirna = query_parts[1]
            except:
                current_mirna = None
        
        # Alternative: look for miRNA in other common formats
        elif "miRNA:" in line or "miR:" in line:
            try:
                # Extract miRNA ID from various formats
                if "miRNA:" in line:
                    current_mirna = line.split("miRNA:")[1].split()[0]
                elif "miR:" in line:
                    current_mirna = line.split("miR:")[1].split()[0]
            except:
                current_mirna = None
        
        # Check if line contains a miRNA ID that matches our targets
        elif any(mirna in line for target_mirnas in mirna_targets.values() for mirna in target_mirnas):
            # Find which miRNA from our targets is mentioned in this line
            for target_mirnas in mirna_targets.values():
                for mirna in target_mirnas:
                    if mirna in line:
                        current_mirna = mirna
                        break
                if current_mirna:
                    break
        
        # Get binding positions from paired regions
        elif line.startswith("Paired Regions"):
            if i + 1 < len(lines):
                next_line = lines[i + 1].strip()
                parts = next_line.split(',')
                if len(parts) >= 2:
                    try:
                        current_binding_start, current_binding_end = map(int, parts[0].split('-'))
                    except ValueError:
                        current_binding_start = current_binding_end = None
        
        # Get transcript and cleavage position
        elif line.startswith("SiteID:"):
            try:
                site_info = line.split()[1]
                transcript_id, pos = site_info.split(':')
                
                # Get the list of miRNAs that target this transcript
                mirnas = mirna_targets.get(transcript_id, [])
                
                # Determine which miRNA to use
                if current_mirna and current_mirna in mirnas:
                    # Use the current miRNA if it's valid for this transcript
                    miRNA = current_mirna
                elif len(mirnas) == 1:
                    # If only one miRNA targets this transcript, use it
                    miRNA = mirnas[0]
                elif len(mirnas) > 1:
                    # Multiple miRNAs target this transcript, but we can't determine which
                    # For now, let's try to use the first one and warn
                    miRNA = mirnas[0]
                    print(f"Warning: Multiple miRNAs target {transcript_id}, using {miRNA}")
                    print(f"  All possible miRNAs: {', '.join(mirnas)}")
                else:
                    # No miRNA found for this transcript
                    miRNA = "unknown_miRNA"
                    print(f"Warning: No miRNA found for transcript {transcript_id}")
                
                results.append({
                    'miRNA': miRNA,
                    'transcript_id': transcript_id,
                    'slice_pos': int(pos),
                    'binding_start': current_binding_start,
                    'binding_end': current_binding_end
                })
            except (IndexError, ValueError) as e:
                print(f"Warning: Could not parse SiteID line: {line}")
                continue
        
        i += 1
    
    return results

def map_position(transcript_data, mrna_pos):
    """Map mRNA position to genomic coordinates and CDS segment"""
    current_pos = 1
    
    # Check 5'UTR
    for start, end in transcript_data['utr5']:
        length = end - start + 1
        if mrna_pos <= current_pos + length - 1:
            genomic_pos = start + (mrna_pos - current_pos)
            return genomic_pos, "5'UTR", ".", "."
        current_pos += length
    
    # Check CDS
    for i, (start, end, _) in enumerate(transcript_data['cds'], 1):
        length = end - start + 1
        if mrna_pos <= current_pos + length - 1:
            genomic_pos = start + (mrna_pos - current_pos)
            return genomic_pos, "CDS", str(i), str(transcript_data['total_cds'])
        current_pos += length
    
    # Check 3'UTR
    for start, end in transcript_data['utr3']:
        length = end - start + 1
        if mrna_pos <= current_pos + length - 1:
            return start + (mrna_pos - current_pos), "3'UTR", ".", "."
        current_pos += length
    
    return None, "unknown", ".", "."

def process_single_file(args):
    """Process a single cleaveland file - designed for multiprocessing"""
    cleaveland_file, transcripts, mirna_targets, output_dir = args
    
    try:
        # Parse the cleaveland file
        results = parse_cleaveland_batch(cleaveland_file, mirna_targets)
        
        if not results:
            print(f"Warning: No results found in {cleaveland_file}")
            return []
        
        processed_results = []
        
        for result in results:
            transcript_id = result['transcript_id']
            
            if transcript_id not in transcripts:
                print(f"Warning: Transcript {transcript_id} not found in annotation")
                continue
            
            tx_data = transcripts[transcript_id]
            genomic_pos, region, cds_segment, total_cds = map_position(tx_data, result['slice_pos'])
            
            processed_result = {
                'miRNA': result['miRNA'],
                'transcript': transcript_id,
                'cleavage_pos': result['slice_pos'],
                'genome_coord': genomic_pos,
                'region': region,
                'cds_segment': cds_segment,
                'total_cds': total_cds,
                'binding_start': result['binding_start'] if result['binding_start'] else '.',
                'binding_end': result['binding_end'] if result['binding_end'] else '.',
                'strand': tx_data['strand'],
                'tx_start': tx_data['start'],
                'tx_end': tx_data['end'],
                'source_file': os.path.basename(cleaveland_file)
            }
            
            processed_results.append(processed_result)
        
        return processed_results
        
    except Exception as e:
        print(f"Error processing {cleaveland_file}: {e}")
        return []

def create_region_distribution_csv(all_results, output_file):
    """Create CSV file with region distribution summary"""
    region_counts = Counter(result['region'] for result in all_results)
    total_sites = len(all_results)
    
    with open(output_file, 'w') as f:
        f.write("Region,Count,Percentage\n")
        for region, count in sorted(region_counts.items()):
            percentage = (count / total_sites) * 100
            f.write(f"{region},{count},{percentage:.2f}\n")

def create_bar_chart(all_results, output_file):
    """Create bar chart of slice site distribution"""
    region_counts = Counter(result['region'] for result in all_results)
    
    plt.figure(figsize=(10, 6))
    regions = list(region_counts.keys())
    counts = list(region_counts.values())
    
    # Create bar chart with different colors
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    bars = plt.bar(regions, counts, color=colors[:len(regions)])
    
    # Add value labels on bars
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{count}', ha='center', va='bottom', fontsize=10)
    
    plt.title('miRNA Cleavage Site Distribution by Region', fontsize=14, fontweight='bold')
    plt.xlabel('Genomic Region', fontsize=12)
    plt.ylabel('Number of Cleavage Sites', fontsize=12)
    plt.xticks(rotation=45)
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def create_pie_chart(all_results, output_file):
    """Create pie chart of slice site distribution"""
    region_counts = Counter(result['region'] for result in all_results)
    
    plt.figure(figsize=(10, 8))
    regions = list(region_counts.keys())
    counts = list(region_counts.values())
    
    # Create pie chart with custom colors
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    # Calculate percentages
    total = sum(counts)
    
    # Create pie chart
    wedges, texts, autotexts = plt.pie(counts, labels=regions, autopct='%1.1f%%',
                                      colors=colors[:len(regions)], startangle=90)
    
    # Enhance text appearance
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(10)
    
    # Add count information to labels
    for i, (region, count) in enumerate(zip(regions, counts)):
        texts[i].set_text(f'{region}\n(n={count})')
        texts[i].set_fontsize(10)
    
    plt.title('miRNA Cleavage Site Distribution by Region', fontsize=14, fontweight='bold')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def get_cleaveland_files(input_path):
    """Get cleaveland files from either a directory or a single file"""
    cleaveland_files = []
    
    if os.path.isfile(input_path):
        # Single file provided
        cleaveland_files.append(input_path)
        print(f"Processing single file: {input_path}")
    elif os.path.isdir(input_path):
        # Directory provided - find all cleaveland files
        for ext in ['*.txt', '*.out', '*.cleaveland']:
            cleaveland_files.extend(glob.glob(os.path.join(input_path, ext)))
        print(f"Processing directory: {input_path}")
    else:
        print(f"Error: {input_path} is neither a file nor a directory")
        return []
    
    return cleaveland_files

def main():
    if len(sys.argv) != 5:
        print("Usage: ./batch_mirna_annotator.py <cleaveland_input> <targets_file> <annotation.gff3> <output_file>")
        print("  cleaveland_input: Directory containing cleaveland files OR single cleaveland file")
        print("  targets_file: miRNA to target mappings file")
        print("  annotation.gff3: GFF3 annotation file")
        print("  output_file: Output file for combined results")
        print("\nAdditional outputs will be created:")
        print("  - slice_region_distribution_summary.csv")
        print("  - slice_site_distribution_bar.png")
        print("  - slice_site_distribution_pie.png")
        sys.exit(1)
    
    cleaveland_input, targets_file, gff3_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    
    # Find all cleaveland files (handles both files and directories)
    cleaveland_files = get_cleaveland_files(cleaveland_input)
    
    if not cleaveland_files:
        print(f"No cleaveland files found in {cleaveland_input}")
        sys.exit(1)
    
    print(f"Found {len(cleaveland_files)} cleaveland file(s) to process")
    
    # Parse annotation and targets files once
    print("Parsing annotation file...")
    transcripts = parse_gff3(gff3_file)
    print(f"Loaded {len(transcripts)} transcripts")
    
    print("Parsing targets file...")
    mirna_targets = parse_mirna_targets(targets_file)
    print(f"Loaded {len(mirna_targets)} miRNA-target mappings")
    
    # Prepare arguments for multiprocessing
    process_args = [(f, transcripts, mirna_targets, os.path.dirname(cleaveland_input)) for f in cleaveland_files]
    
    # Process files in parallel
    num_cores = min(mp.cpu_count(), len(cleaveland_files))
    print(f"Processing files using {num_cores} cores...")
    
    all_results = []
    
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        for i, results in enumerate(executor.map(process_single_file, process_args)):
            all_results.extend(results)
            if (i + 1) % 10 == 0:  # Progress update every 10 files
                print(f"Processed {i + 1}/{len(cleaveland_files)} files")
    
    if not all_results:
        print("No results to process. Exiting.")
        sys.exit(1)
    
    # Write main output table (maintains original functionality)
    print(f"Writing {len(all_results)} results to {output_file}")
    
    with open(output_file, 'w') as f:
        # Write header
        f.write("miRNA,Transcript,Cleavage_Pos,Genome_Coord,Region,CDS_Segment,Total_CDS,"
                "Binding_Start,Binding_End,Strand,Tx_Start,Tx_End,Source_File\n")
        
        # Write results
        for result in all_results:
            f.write(f"{result['miRNA']},{result['transcript']},{result['cleavage_pos']},"
                   f"{result['genome_coord']},{result['region']},{result['cds_segment']},"
                   f"{result['total_cds']},{result['binding_start']},{result['binding_end']},"
                   f"{result['strand']},{result['tx_start']},{result['tx_end']},{result['source_file']}\n")
    
    # Create additional output files
    output_dir = os.path.dirname(output_file) if os.path.dirname(output_file) else '.'
    
    # 1. Region distribution CSV
    csv_file = os.path.join(output_dir, 'slice_region_distribution_summary.csv')
    print(f"Creating region distribution summary: {csv_file}")
    create_region_distribution_csv(all_results, csv_file)
    
    # 2. Bar chart
    bar_file = os.path.join(output_dir, 'slice_site_distribution_bar.png')
    print(f"Creating bar chart: {bar_file}")
    create_bar_chart(all_results, bar_file)
    
    # 3. Pie chart
    pie_file = os.path.join(output_dir, 'slice_site_distribution_pie.png')
    print(f"Creating pie chart: {pie_file}")
    create_pie_chart(all_results, pie_file)
    
    print(f"Analysis complete! Results written to {output_file}")
    
    # Print summary statistics
    regions = {}
    mirnas = set()
    transcripts_processed = set()
    
    for result in all_results:
        region = result['region']
        regions[region] = regions.get(region, 0) + 1
        mirnas.add(result['miRNA'])
        transcripts_processed.add(result['transcript'])
    
    print(f"\nSummary:")
    print(f"Total sites: {len(all_results)}")
    print(f"Unique miRNAs: {len(mirnas)}")
    print(f"Unique transcripts: {len(transcripts_processed)}")
    print(f"Region distribution:")
    for region, count in sorted(regions.items()):
        print(f"  {region}: {count}")
    
    print(f"\nAdditional files created:")
    print(f"  - {csv_file}")
    print(f"  - {bar_file}")
    print(f"  - {pie_file}")

if __name__ == "__main__":
    main()
