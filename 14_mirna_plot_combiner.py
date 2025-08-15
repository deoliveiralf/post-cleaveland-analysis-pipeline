#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import os
import argparse
from pathlib import Path

def combine_mirna_plots(table_file, mirna_dir, targets_dir, correlation_dir, output_dir="combined_plots", figsize=(15, 5)):
    """
    Combine three different plots for miRNA-target pairs from separate directories.
    
    Parameters:
    -----------
    table_file : str
        Path to the tab-separated file containing miRNA ID and Target ID columns
    mirna_dir : str
        Directory containing miRNA heatmap PNG files
    targets_dir : str
        Directory containing targets heatmap PNG files
    correlation_dir : str
        Directory containing correlation plot PNG files
    output_dir : str
        Directory to save the combined plots
    figsize : tuple
        Figure size (width, height) in inches
    """
    
    # Read the miRNA-target pairs table
    try:
        df = pd.read_csv(table_file, sep='\t')
        print(f"Loaded {len(df)} miRNA-target pairs from {table_file}")
    except Exception as e:
        print(f"Error reading table file: {e}")
        return
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Process each miRNA-target pair
    successful_combinations = 0
    missing_files = []
    
    for idx, row in df.iterrows():
        mirna_id = row['miRNA ID']
        target_id = row['Target ID']
        
        # Generate filenames based on the pattern
        mirna_heatmap = f"{mirna_id}_logFC_heatmap.png"
        targets_heatmap = f"{mirna_id}_logFC_heatmap.png"  # Same pattern as miRNA?
        correlation_plot = f"{mirna_id}_{target_id}_combined.png"
        
        # Full paths to different directories
        mirna_path = os.path.join(mirna_dir, mirna_heatmap)
        targets_path = os.path.join(targets_dir, targets_heatmap)
        correlation_path = os.path.join(correlation_dir, correlation_plot)
        
        # Check if all files exist
        files_exist = [
            os.path.exists(mirna_path),
            os.path.exists(targets_path),
            os.path.exists(correlation_path)
        ]
        
        if not all(files_exist):
            missing = []
            if not files_exist[0]: missing.append(mirna_heatmap)
            if not files_exist[1]: missing.append(targets_heatmap)
            if not files_exist[2]: missing.append(correlation_plot)
            missing_files.extend(missing)
            print(f"Skipping {mirna_id} - {target_id}: Missing files {missing}")
            continue
        
        try:
            # Load images
            img_mirna = mpimg.imread(mirna_path)
            img_targets = mpimg.imread(targets_path)
            img_correlation = mpimg.imread(correlation_path)
            
            # Create figure with custom layout
            fig = plt.figure(figsize=figsize)
            gs = GridSpec(1, 3, figure=fig, width_ratios=[1, 1, 1.2])
            
            # Plot miRNA heatmap
            ax1 = fig.add_subplot(gs[0, 0])
            ax1.imshow(img_mirna)
            ax1.set_title(f'{mirna_id}\nLog Fold Change', fontsize=10, pad=10)
            ax1.axis('off')
            
            # Plot targets heatmap
            ax2 = fig.add_subplot(gs[0, 1])
            ax2.imshow(img_targets)
            ax2.set_title(f'{mirna_id} Targets\nLog Fold Changes', fontsize=10, pad=10)
            ax2.axis('off')
            
            # Plot correlation
            ax3 = fig.add_subplot(gs[0, 2])
            ax3.imshow(img_correlation)
            ax3.set_title(f'Correlation Analysis\n{mirna_id} vs {target_id}', fontsize=10, pad=10)
            ax3.axis('off')
            
            # Adjust layout
            plt.tight_layout()
            
            # Save combined plot
            output_filename = f"{mirna_id}_{target_id}_combined_analysis.png"
            output_path = os.path.join(output_dir, output_filename)
            plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()
            
            successful_combinations += 1
            print(f"✓ Created: {output_filename}")
            
        except Exception as e:
            print(f"Error processing {mirna_id} - {target_id}: {e}")
            continue
    
    # Summary
    print(f"\n=== SUMMARY ===")
    print(f"Total pairs processed: {len(df)}")
    print(f"Successful combinations: {successful_combinations}")
    print(f"Failed combinations: {len(df) - successful_combinations}")
    
    if missing_files:
        print(f"\nMissing files ({len(set(missing_files))} unique):")
        for file in sorted(set(missing_files)):
            print(f"  - {file}")
    
    print(f"\nCombined plots saved in: {output_dir}/")

def create_sample_table(filename="mirna_targets_sample.txt"):
    """
    Create a sample miRNA-targets table for testing
    """
    sample_data = """miRNA ID	Target ID
novel404_Chr09_39038	Sevir.1G015900.1
novel404_Chr09_39038	Sevir.9G194400.1
novel404_Chr09_39038	Sevir.9G459700.1
miR160_Chr02_7942	Sevir.1G296300.1
miR160_Chr02_7942	Sevir.2G320400.1
miR160_Chr02_7942	Sevir.4G268700.5
miR160_Chr02_7942	Sevir.4G268700.6
miR160_Chr02_7942	Sevir.7G179200.2
miR160_Chr02_7942	Sevir.9G219100.1
miR160_Chr02_7942	Sevir.9G518100.2"""
    
    with open(filename, 'w') as f:
        f.write(sample_data)
    print(f"Sample table created: {filename}")

def main():
    """
    Main function to handle command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Combine miRNA-target plots from separate PNG directories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s mirna_targets.txt mirna_heatmap/ targets_heatmap/ correlation/
  %(prog)s mirna_targets.txt mirna_heatmap/ targets_heatmap/ correlation/ -o results/
  %(prog)s mirna_targets.txt mirna_heatmap/ targets_heatmap/ correlation/ --figsize 18 6
        """
    )
    
    parser.add_argument('table_file', 
                       help='Tab-separated file with miRNA ID and Target ID columns')
    
    parser.add_argument('mirna_dir', 
                       help='Directory containing miRNA heatmap PNG files')
    
    parser.add_argument('targets_dir', 
                       help='Directory containing targets heatmap PNG files')
    
    parser.add_argument('correlation_dir', 
                       help='Directory containing correlation plot PNG files')
    
    parser.add_argument('-o', '--output-dir', 
                       default='combined_plots', 
                       help='Output directory for combined plots (default: combined_plots)')
    
    parser.add_argument('--figsize', 
                       nargs=2, 
                       type=float, 
                       default=[15, 5], 
                       metavar=('WIDTH', 'HEIGHT'),
                       help='Figure size in inches (default: 15 5)')
    
    parser.add_argument('--create-sample', 
                       action='store_true',
                       help='Create a sample miRNA-targets table for testing')
    
    args = parser.parse_args()
    
    # Create sample table if requested
    if args.create_sample:
        create_sample_table("mirna_targets_sample.txt")
        print("Sample table created. You can now run:")
        print(f"  {parser.prog} mirna_targets_sample.txt mirna_heatmap/ targets_heatmap/ correlation/")
        return
    
    # Check if table file exists
    if not os.path.exists(args.table_file):
        print(f"Error: Table file '{args.table_file}' not found!")
        return
    
    # Check if all directories exist
    directories = {
        'miRNA heatmap': args.mirna_dir,
        'targets heatmap': args.targets_dir,
        'correlation': args.correlation_dir
    }
    
    for name, directory in directories.items():
        if not os.path.exists(directory):
            print(f"Error: {name} directory '{directory}' not found!")
            return
    
    print("Starting miRNA-target plot combination...")
    print(f"Table file: {args.table_file}")
    print(f"miRNA heatmap directory: {args.mirna_dir}")
    print(f"Targets heatmap directory: {args.targets_dir}")
    print(f"Correlation directory: {args.correlation_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Figure size: {args.figsize[0]} x {args.figsize[1]} inches")
    print("-" * 50)
    
    # Run the combination process
    combine_mirna_plots(
        table_file=args.table_file,
        mirna_dir=args.mirna_dir,
        targets_dir=args.targets_dir,
        correlation_dir=args.correlation_dir,
        output_dir=args.output_dir,
        figsize=tuple(args.figsize)
    )

if __name__ == "__main__":
    main()
