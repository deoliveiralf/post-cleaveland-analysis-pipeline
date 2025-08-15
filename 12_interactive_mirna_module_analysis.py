#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import argparse
from datetime import datetime
from matplotlib.colors import LinearSegmentedColormap

# Configuration
SAMPLE_PAIRS = [
    ('Third_Internode', 'Fifth_Internode'),
    ('Young_Leaves', 'Mature_Leaves'), 
    ('Young_Roots', 'Mature_Roots')
]
TOP_N_MODULES = 10
SIGNIFICANCE_LEVEL = 0.05

# Custom color maps
FC_CMAP = LinearSegmentedColormap.from_list('diverging', ['#2166ac', '#f7f7f7', '#b2182b'], N=256)
CORR_CMAP = LinearSegmentedColormap.from_list('correlation', ['#006837', '#a7d26d', '#ffffbf'], N=256)

class Logger:
    """Dual logger for console and file output."""
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, 'w')

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.terminal.flush()
        self.log.flush()

def setup_logging(log_file, output_dir):
    """Initialize logging to file and console."""
    os.makedirs(output_dir, exist_ok=True)
    sys.stdout = Logger(log_file)
    print(f"\n=== Analysis started at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ===\n")

def load_data(mirna_modules_file, mirnas_cpm_file, transcripts_tpm_file, cleavage_sites_file=None):
    """Load and preprocess input files."""
    print("Loading data files...")
    files = [mirna_modules_file, mirnas_cpm_file, transcripts_tpm_file]
    for f in files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Missing file: {f}")
    
    modules = pd.read_csv(mirna_modules_file, sep='\t')
    mirnas = pd.read_csv(mirnas_cpm_file, sep='\t')
    transcripts = pd.read_csv(transcripts_tpm_file, sep='\t')
    
    # Load cleavage sites if provided
    cleavage_sites = None
    if cleavage_sites_file and os.path.exists(cleavage_sites_file):
        cleavage_sites = pd.read_csv(cleavage_sites_file)
        print(f"Loaded cleavage sites for {len(cleavage_sites)} miRNA-target pairs with additional annotations")
    
    # Convert to numeric
    for df in [mirnas, transcripts]:
        df.iloc[:,1:] = df.iloc[:,1:].apply(pd.to_numeric, errors='coerce').fillna(0)
    
    print(f"Loaded {len(modules)} miRNA-target pairs across {len(mirnas)} miRNAs and {len(transcripts)} transcripts.")
    return modules, mirnas, transcripts, cleavage_sites

def calculate_log_fold_change(df, group1, group2):
    """Calculate log2 fold change between sample groups."""
    cols1 = [c for c in df.columns if group1 in c]
    cols2 = [c for c in df.columns if group2 in c]
    mean1 = df[cols1].mean(axis=1)
    mean2 = df[cols2].mean(axis=1)
    return np.log2((mean2 + 0.1) / (mean1 + 0.1))

def pearson_correlation(x, y):
    """Calculate Pearson correlation and p-value manually (without scipy)."""
    try:
        x = np.array(x, dtype=float)
        y = np.array(y, dtype=float)
        
        # Remove NaN values
        valid_mask = (~np.isnan(x)) & (~np.isnan(y)) & (~np.isinf(x)) & (~np.isinf(y))
        x_clean = x[valid_mask]
        y_clean = y[valid_mask]
        
        n = len(x_clean)
        if n < 3:  # Minimum 3 points for meaningful correlation
            return np.nan, np.nan
            
        # Calculate means
        mean_x = np.mean(x_clean)
        mean_y = np.mean(y_clean)
        
        # Calculate correlation coefficient
        numerator = np.sum((x_clean - mean_x) * (y_clean - mean_y))
        sum_sq_x = np.sum((x_clean - mean_x)**2)
        sum_sq_y = np.sum((y_clean - mean_y)**2)
        
        denominator = np.sqrt(sum_sq_x * sum_sq_y)
        
        if denominator == 0 or sum_sq_x == 0 or sum_sq_y == 0:
            return np.nan, np.nan
            
        r = numerator / denominator
        
        # Ensure r is within valid bounds
        r = np.clip(r, -1, 1)
        
        # Calculate p-value using t-distribution approximation
        if abs(r) >= 0.999:  # Very high correlation
            p_value = 0.001
        else:
            t_stat = r * np.sqrt((n - 2) / (1 - r**2))
            # Simplified p-value approximation (two-tailed test)
            df = n - 2
            t_abs = abs(t_stat)
            # Very simplified approximation - for better accuracy, implement full t-distribution
            p_value = 2 * np.exp(-0.717 * t_abs - 0.416 * t_abs**2) if t_abs < 3 else 0.001
            p_value = min(1.0, max(0.0, p_value))
                
        return r, p_value
    except Exception as e:
        print(f"Error in correlation calculation: {e}")
        return np.nan, np.nan

def linear_regression(x, y):
    """Manual linear regression for plotting trend lines."""
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)
    valid_mask = (~np.isnan(x)) & (~np.isnan(y))
    x = x[valid_mask]
    y = y[valid_mask]
    
    if len(x) < 2:
        return 0, 0
    
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]
    return m, c

def plot_correlation_scatter(mirna_values, target_values, ax, sample_pair, r_value, p_value):
    """Plot miRNA vs target expression with regression line and r value."""
    # Convert to numpy arrays and filter NaN values
    mirna_values = np.array(mirna_values, dtype=float)
    target_values = np.array(target_values, dtype=float)
    valid_mask = (~np.isnan(mirna_values)) & (~np.isnan(target_values))
    mirna_values = mirna_values[valid_mask]
    target_values = target_values[valid_mask]
    
    if len(mirna_values) == 0:
        ax.text(0.5, 0.5, "No valid data points", ha='center', va='center', fontsize=12)
        ax.set_title(f'{sample_pair}\n(No valid data)', fontsize=12)
        return
    
    # Calculate regression line
    m, c = linear_regression(mirna_values, target_values)
    
    # Plot scatter points with larger size
    ax.scatter(mirna_values, target_values, s=80, alpha=0.7, color='#1f77b4', edgecolors='black', linewidth=0.5)
    
    # Plot regression line
    if len(mirna_values) > 1:
        x_vals = np.array([min(mirna_values), max(mirna_values)])
        y_vals = c + m * x_vals
        ax.plot(x_vals, y_vals, color='#d62728', linewidth=3)
    
    # Add r value annotation with significance star if needed
    sig_star = '*' if (not np.isnan(p_value)) and (p_value < SIGNIFICANCE_LEVEL) else ''
    ax.text(0.05, 0.95, f"r = {r_value:.2f}{sig_star}", 
            transform=ax.transAxes, fontsize=12, fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.9, edgecolor='black'))
    
    ax.set_xlabel('miRNA CPM (log2)', fontsize=12)
    ax.set_ylabel('Target TPM (log2)', fontsize=12)
    ax.set_title(f'{sample_pair[:15]}...\nmiRNA vs Target' if len(sample_pair) > 15 else f'{sample_pair}\nmiRNA vs Target', fontsize=12)
    ax.tick_params(labelsize=10)
    ax.grid(True, alpha=0.3)

def create_module_table(mirna, data, output_dir):
    """Create detailed table for each module with expression values and correlations."""
    print(f"  Creating detailed table for {mirna}...")
    
    # Get basic info
    targets = data['targets']
    mirna_expr = data['mirna_expr']
    target_expr = data['target_expr']
    cleavage_info = data.get('cleavage_info', {})
    
    # Get all common samples
    all_samples = sorted(set(mirna_expr.columns[1:]) & set(target_expr.columns[1:]))
    
    # Create comprehensive table (without cleavage info)
    table_data = []
    
    # Add miRNA row
    mirna_row = {
        'Type': 'miRNA',
        'ID': mirna,
        'Description': f'miRNA module with {len(targets)} targets'
    }
    
    # Add miRNA expression values
    for sample in all_samples:
        mirna_row[sample] = mirna_expr[sample].values[0] if sample in mirna_expr.columns else np.nan
    
    # Add fold changes for miRNA
    for g1, g2 in SAMPLE_PAIRS:
        pair_name = f'{g2}/{g1}'
        mirna_row[f'FC_{pair_name}'] = data['fold_changes'][f'{pair_name}_miRNA']
    
    table_data.append(mirna_row)
    
    # Add target rows (without cleavage info)
    for i, target in enumerate(targets):
        target_row = {
            'Type': 'Target',
            'ID': target,
            'Description': f'Target #{i+1} of {mirna}'
        }
        
        # Add target expression values
        target_df = target_expr[target_expr['GeneID'] == target]
        if not target_df.empty:
            target_expr_row = target_df.iloc[0]
            for sample in all_samples:
                target_row[sample] = target_expr_row[sample] if sample in target_expr_row.index else np.nan
        else:
            # Target not found in expression data
            for sample in all_samples:
                target_row[sample] = np.nan
            print(f"    Warning: Target {target} not found in expression data")
        
        # Add fold changes for target
        for j, (g1, g2) in enumerate(SAMPLE_PAIRS):
            pair_name = f'{g2}/{g1}'
            if i < len(data['fold_changes'][f'{pair_name}_Targets']):
                target_row[f'FC_{pair_name}'] = data['fold_changes'][f'{pair_name}_Targets'][i]
            else:
                target_row[f'FC_{pair_name}'] = np.nan
        
        table_data.append(target_row)
    
    # Create DataFrame
    df = pd.DataFrame(table_data)
    
    # Create correlation summary table (without cleavage info)
    corr_data = []
    for i, target in enumerate(targets):
        corr_row = {
            'miRNA': mirna,
            'Target': target,
            'Target_Index': i
        }
        
        # Add correlations for each sample pair
        for g1, g2 in SAMPLE_PAIRS:
            pair_name = f'{g2}/{g1}'
            if i < len(data['correlations'][pair_name]):
                corr_row[f'Correlation_{pair_name}'] = data['correlations'][pair_name][i]
                corr_row[f'P_value_{pair_name}'] = data['p_values'][pair_name][i]
                corr_row[f'Significant_{pair_name}'] = (
                    data['p_values'][pair_name][i] < SIGNIFICANCE_LEVEL 
                    if not np.isnan(data['p_values'][pair_name][i]) else False
                )
            else:
                corr_row[f'Correlation_{pair_name}'] = np.nan
                corr_row[f'P_value_{pair_name}'] = np.nan
                corr_row[f'Significant_{pair_name}'] = False
        
        # Add summary statistics
        all_corrs = []
        for g1, g2 in SAMPLE_PAIRS:
            pair_name = f'{g2}/{g1}'
            if i < len(data['correlations'][pair_name]):
                all_corrs.append(data['correlations'][pair_name][i])
        
        valid_corrs = [c for c in all_corrs if not np.isnan(c)]
        
        if valid_corrs:
            corr_row['Mean_Correlation'] = np.mean(valid_corrs)
            corr_row['Mean_Abs_Correlation'] = np.mean(np.abs(valid_corrs))
            corr_row['Min_Correlation'] = np.min(valid_corrs)
            corr_row['Max_Correlation'] = np.max(valid_corrs)
        else:
            corr_row['Mean_Correlation'] = np.nan
            corr_row['Mean_Abs_Correlation'] = np.nan
            corr_row['Min_Correlation'] = np.nan
            corr_row['Max_Correlation'] = np.nan
        
        corr_data.append(corr_row)
    
    corr_df = pd.DataFrame(corr_data)
    
    # Create dedicated cleavage site table - include all targets, with NA for those without cleavage info
    cleavage_data = []
    for target in targets:
        if target in cleavage_info:
            # Target has cleavage info - add all sites with all available information
            for site_info in cleavage_info[target]:
                cleavage_data.append({
                    'miRNA': mirna,
                    'Target': target,
                    'Cleavage_Position': site_info['Cleavage_Pos'],
                    'Genome_Coordinate': site_info['Genome_Coord'],
                    'Region': site_info['Region'],
                    'CDS_Segment': site_info['CDS_Segment'],
                    'Total_CDS': site_info['Total_CDS'],
                    'Binding_Start': site_info['Binding_Start'],
                    'Binding_End': site_info['Binding_End'],
                    'Strand': site_info['Strand'],
                    'Tx_Start': site_info['Tx_Start'],
                    'Tx_End': site_info['Tx_End'],
                    'Source_File': site_info['Source_File']
                })
        else:
            # Target has no cleavage info - add one row with NA values
            cleavage_data.append({
                'miRNA': mirna,
                'Target': target,
                'Cleavage_Position': 'NA',
                'Genome_Coordinate': 'NA',
                'Region': 'NA',
                'CDS_Segment': 'NA',
                'Total_CDS': 'NA',
                'Binding_Start': 'NA',
                'Binding_End': 'NA',
                'Strand': 'NA',
                'Tx_Start': 'NA',
                'Tx_End': 'NA',
                'Source_File': 'NA'
            })
    
    cleavage_df = pd.DataFrame(cleavage_data)

    # Create module summary table
    summary_data = {
        'miRNA': mirna,
        'Target_Count': len(targets),
        'Mean_Abs_Correlation': data['mean_abs_corr'],
        'Mean_Neg_Correlation': data['mean_neg_corr'],
        'Available_Samples': len(all_samples),
        'Sample_Pairs_Analyzed': len(SAMPLE_PAIRS)
    }
    
    # Add fold change summaries
    for g1, g2 in SAMPLE_PAIRS:
        pair_name = f'{g2}/{g1}'
        summary_data[f'miRNA_FC_{pair_name}'] = data['fold_changes'][f'{pair_name}_miRNA']
        target_fcs = data['fold_changes'][f'{pair_name}_Targets']
        summary_data[f'Mean_Target_FC_{pair_name}'] = np.mean(target_fcs)
        summary_data[f'Median_Target_FC_{pair_name}'] = np.median(target_fcs)
    
    summary_df = pd.DataFrame([summary_data])
    
    # Save tables
    safe_mirna_name = mirna.replace('/', '_').replace('\\', '_')
    
    # Main expression table
    expr_file = f"{output_dir}/{safe_mirna_name}_expression_table.csv"
    df.to_csv(expr_file, index=False)
    
    # Correlation table
    corr_file = f"{output_dir}/{safe_mirna_name}_correlation_table.csv"
    corr_df.to_csv(corr_file, index=False)
    
    # Dedicated cleavage site table (always created)
    cleavage_file = f"{output_dir}/{safe_mirna_name}_cleavage_sites.csv"
    cleavage_df.to_csv(cleavage_file, index=False)
    print(f"    Saved detailed cleavage sites table: {cleavage_file}")
    
    # Summary table
    summary_file = f"{output_dir}/{safe_mirna_name}_summary_table.csv"
    summary_df.to_csv(summary_file, index=False)
    
    print(f"    Saved tables: {expr_file}, {corr_file}, {cleavage_file}, {summary_file}")
    
    return df, corr_df, summary_df, cleavage_df

def prepare_module_data(modules, mirnas, transcripts, cleavage_sites=None):
    """Organize data by miRNA-transcript modules and calculate metrics."""
    print("Preparing module data...")
    module_dict = {}
    
    for mirna, group in modules.groupby('miRNA_Chr ID'):
        targets = group['Target ID'].unique()
        mirna_expr = mirnas[mirnas['GeneID'] == mirna]
        target_expr = transcripts[transcripts['GeneID'].isin(targets)]
        
        if not mirna_expr.empty and not target_expr.empty:
            # Get cleavage information for this miRNA-target pairs with all available annotations
            cleavage_info = {}
            if cleavage_sites is not None:
                for target in targets:
                    sites = cleavage_sites[(cleavage_sites['miRNA'] == mirna) & 
                                         (cleavage_sites['Transcript'] == target)]
                    if not sites.empty:
                        # Store all site information as a list of dictionaries
                        cleavage_info[target] = sites.to_dict('records')
            
            fc_data, correlations, p_values = {}, {}, {}
            all_corrs = []
            all_neg_corrs = []
            
            for g1, g2 in SAMPLE_PAIRS:
                pair_name = f'{g2}/{g1}'
                mirna_fc = calculate_log_fold_change(mirna_expr, g1, g2).values[0]
                target_fcs, corrs, ps = [], [], []
                
                # Get all available samples for this miRNA and targets
                all_common_samples = sorted(set(mirna_expr.columns[1:]) & set(target_expr.columns[1:]))
                
                # Filter samples to only those belonging to the current sample pair (g1, g2)
                pair_samples = [c for c in all_common_samples if g1 in c or g2 in c]
                
                print(f"    Processing {pair_name}: Using samples {pair_samples}")
                
                for target_idx, (_, target_row) in enumerate(target_expr.iterrows()):
                    target = target_row['GeneID']
                    target_fc = calculate_log_fold_change(pd.DataFrame(target_row).T, g1, g2).values[0]
                    target_fcs.append(target_fc)
                    
                    # Extract values for correlation calculation
                    mirna_vals = mirna_expr[pair_samples].values.flatten()
                    target_vals = target_row[pair_samples].values
                    
                    r, p = pearson_correlation(mirna_vals, target_vals)
                    corrs.append(r)
                    ps.append(p)
                    if not np.isnan(r):
                        all_corrs.append(abs(r))
                        if r < 0:  # Collect negative correlations
                            all_neg_corrs.append(r)
                
                fc_data[f'{pair_name}_miRNA'] = mirna_fc
                fc_data[f'{pair_name}_Targets'] = target_fcs
                correlations[pair_name] = corrs
                p_values[pair_name] = ps
            
            module_dict[mirna] = {
                'targets': targets,
                'mirna_expr': mirna_expr,
                'target_expr': target_expr,
                'fold_changes': fc_data,
                'correlations': correlations,
                'p_values': p_values,
                'mean_abs_corr': np.nanmean(all_corrs) if all_corrs else 0,
                'mean_neg_corr': np.nanmean(all_neg_corrs) if all_neg_corrs else np.nan,
                'cleavage_info': cleavage_info
            }
    
    print(f"Prepared {len(module_dict)} valid miRNA-target modules.")
    return module_dict

def add_cell_borders(ax, data_matrix):
    """Add borders to heatmap cells."""
    for i in range(data_matrix.shape[0]+1):
        ax.axhline(i, color='black', linewidth=0.8)
    for j in range(data_matrix.shape[1]+1):
        ax.axvline(j, color='black', linewidth=0.8)

def plot_individual_heatmaps(module_data, output_dir):
    """Generate individual heatmaps and correlation plots for all modules."""
    print("\nGenerating individual plots:")
    
    for i, (mirna, data) in enumerate(module_data.items(), 1):
        print(f"  {i}. Processing {mirna}...")
        targets = data['targets']
        
        # Create module tables
        create_module_table(mirna, data, output_dir)
        
        # Dynamic figure sizing based on number of targets
        base_height = 10
        target_height = max(8, len(targets) * 0.6)  # Minimum 8 inches, scale with targets
        total_height = base_height + target_height
        
        # Create figure with improved layout and larger size
        fig = plt.figure(figsize=(24, total_height))
        gs = fig.add_gridspec(3, 4, width_ratios=[4, 4, 2, 2], 
                            height_ratios=[2, len(targets) * 0.5, 6],
                            wspace=0.25, hspace=0.4)
        
        # miRNA Fold Change Heatmap (with colorbar)
        ax1 = fig.add_subplot(gs[0, :2])
        mirna_fc = [data['fold_changes'][f'{g2}/{g1}_miRNA'] for g1, g2 in SAMPLE_PAIRS]
        sns.heatmap([mirna_fc], cmap=FC_CMAP, center=0, annot=True, 
                    fmt=".2f", cbar=True, cbar_kws={'location': 'right', 'shrink': 0.6},
                    annot_kws={'size': 14, 'weight': 'bold'}, 
                    xticklabels=[f'{g2}/{g1}' for g1, g2 in SAMPLE_PAIRS], 
                    yticklabels=['miRNA'], ax=ax1)
        ax1.set_title(f'miRNA Fold Change (Mean |r| = {data["mean_abs_corr"]:.2f})', 
                     fontsize=16, fontweight='bold', pad=20)
        ax1.tick_params(labelsize=12)
        add_cell_borders(ax1, np.array([mirna_fc]))
        
        # Targets Fold Change Heatmap
        ax2 = fig.add_subplot(gs[1, :2])
        target_fc_matrix = np.array([data['fold_changes'][f'{g2}/{g1}_Targets'] for g1, g2 in SAMPLE_PAIRS]).T
        sns.heatmap(target_fc_matrix, cmap=FC_CMAP, center=0, annot=True, 
                    fmt=".2f", cbar=False, annot_kws={'size': 10, 'weight': 'bold'},
                    xticklabels=[f'{g2}/{g1}' for g1, g2 in SAMPLE_PAIRS], 
                    yticklabels=targets, ax=ax2)
        ax2.set_title('Target Fold Changes', fontsize=16, fontweight='bold', pad=20)
        ax2.tick_params(labelsize=12)
        add_cell_borders(ax2, target_fc_matrix)
        
        # Correlation Heatmap with significance stars
        ax3 = fig.add_subplot(gs[2, :2])
        corr_matrix = np.array([data['correlations'][f'{g2}/{g1}'] for g1, g2 in SAMPLE_PAIRS]).T
        p_matrix = np.array([data['p_values'][f'{g2}/{g1}'] for g1, g2 in SAMPLE_PAIRS]).T
        
        # Create annotation matrix with stars for significant correlations
        annot_matrix = []
        for i in range(corr_matrix.shape[0]):
            row = []
            for j in range(corr_matrix.shape[1]):
                val = corr_matrix[i,j]
                p = p_matrix[i,j]
                star = '*' if (not np.isnan(p)) and (p < SIGNIFICANCE_LEVEL) else ''
                row.append(f"{val:.2f}{star}")
            annot_matrix.append(row)
        annot_matrix = np.array(annot_matrix)
        
        sns.heatmap(corr_matrix, cmap=CORR_CMAP, center=0, 
                    annot=annot_matrix, fmt="", vmin=-1, vmax=1, cbar=True,
                    cbar_kws={'location': 'right', 'shrink': 0.6},
                    annot_kws={'size': 10, 'weight': 'bold'},
                    xticklabels=[f'{g2}/{g1}' for g1, g2 in SAMPLE_PAIRS], 
                    yticklabels=targets, ax=ax3)
        ax3.set_title('Pearson Correlation Coefficients (* p<0.05)', fontsize=16, fontweight='bold', pad=20)
        ax3.tick_params(labelsize=12)
        add_cell_borders(ax3, corr_matrix)
        
        # Correlation Scatter Plots - distributed vertically on the right side
        for j, (g1, g2) in enumerate(SAMPLE_PAIRS):
            # Position scatter plots vertically on the right with better spacing
            if j == 0:
                ax = fig.add_subplot(gs[0, 2:])
            elif j == 1:
                ax = fig.add_subplot(gs[1, 2:])
            else:
                ax = fig.add_subplot(gs[2, 2:])
            
            pair_name = f'{g2}/{g1}'
            
            # Get all available samples and filter to only this sample pair
            all_common_samples = sorted(set(data['mirna_expr'].columns[1:]) & 
                                      set(data['target_expr'].columns[1:]))
            pair_samples = [c for c in all_common_samples if g1 in c or g2 in c]
            
            # Find the target with the highest absolute correlation for this sample pair
            best_target_idx = 0
            best_corr = 0
            correlations_for_pair = data['correlations'][pair_name]
            
            for idx, corr in enumerate(correlations_for_pair):
                if not np.isnan(corr) and abs(corr) > abs(best_corr):
                    best_corr = corr
                    best_target_idx = idx
            
            # Plot the best target for this sample pair using only samples from this pair
            target_values = data['target_expr'].iloc[best_target_idx][pair_samples].values
            mirna_values = data['mirna_expr'][pair_samples].values.flatten()
            
            # Get the actual correlation and p-value for this target and sample pair
            r_value = correlations_for_pair[best_target_idx] if best_target_idx < len(correlations_for_pair) else np.nan
            p_value = data['p_values'][pair_name][best_target_idx] if best_target_idx < len(data['p_values'][pair_name]) else np.nan
            
            # Add target name to the plot title
            target_name = data['targets'][best_target_idx] if best_target_idx < len(data['targets']) else 'Target'
            plot_correlation_scatter(mirna_values, target_values, ax, 
                                   f"{pair_name}\n({target_name})", r_value, p_value)
        
        # Improve overall layout and margins
        plt.subplots_adjust(left=0.08, right=0.95, top=0.93, bottom=0.07)
        
        # Save with higher DPI and better format
        plt.savefig(f"{output_dir}/{mirna.replace('/', '_')}_analysis.png", 
                    dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        plt.close()

def generate_top_modules_report(module_data, top_n, output_dir):
    """Generate summary of top modules by correlation strength."""
    print(f"\nGenerating top {top_n} modules report...")
    ranked_modules = sorted(
        module_data.items(),
        key=lambda x: x[1]['mean_abs_corr'], 
        reverse=True
    )[:top_n]
    
    # Filter modules with negative correlations for the second plot
    negative_modules = [(mirna, data) for mirna, data in module_data.items() 
                       if not np.isnan(data['mean_neg_corr'])]
    ranked_negative_modules = sorted(
        negative_modules,
        key=lambda x: abs(x[1]['mean_neg_corr']), 
        reverse=True
    )[:top_n]
    
    # Create figure with improved sizing and layout
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 14))
    
    # First plot: Top modules by mean absolute correlation
    bars1 = ax1.bar(
        range(len(ranked_modules)),
        [mod[1]['mean_abs_corr'] for mod in ranked_modules],
        color=plt.cm.viridis(np.linspace(0, 1, len(ranked_modules))),
        edgecolor='black',
        linewidth=1
    )
    
    ax1.set_title(f'Top {top_n} miRNA Modules by Mean Absolute Pearson Correlation', 
                  fontsize=18, fontweight='bold', pad=25)
    ax1.set_ylabel('Mean |r| across targets', fontsize=14, labelpad=15)
    ax1.set_xlabel('miRNA Module', fontsize=14, labelpad=15)
    ax1.set_xticks(range(len(ranked_modules)))
    ax1.set_xticklabels([mod[0] for mod in ranked_modules], rotation=45, ha='right', fontsize=12)
    ax1.tick_params(axis='y', labelsize=12)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Add correlation examples as annotations
    for i, (mod, bar) in enumerate(zip(ranked_modules, bars1)):
        best_corr = max(
            np.abs(np.concatenate(list(mod[1]['correlations'].values()))),
            default=0
        )
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height/2, 
                f"max |r|={best_corr:.2f}", 
                ha='center', va='center', color='white', fontweight='bold', fontsize=10)
    
    # Second plot: Top modules with negative correlations
    if ranked_negative_modules:
        bars2 = ax2.bar(
            range(len(ranked_negative_modules)),
            [mod[1]['mean_neg_corr'] for mod in ranked_negative_modules],
            color=plt.cm.Reds_r(np.linspace(0, 1, len(ranked_negative_modules))),
            edgecolor='black',
            linewidth=1
        )
        
        ax2.set_title(f'Top {top_n} miRNA Modules by Mean Negative Pearson Correlation', 
                      fontsize=18, fontweight='bold', pad=25)
        ax2.set_ylabel('Mean negative r across targets', fontsize=14, labelpad=15)
        ax2.set_xlabel('miRNA Module', fontsize=14, labelpad=15)
        ax2.set_xticks(range(len(ranked_negative_modules)))
        ax2.set_xticklabels([mod[0] for mod in ranked_negative_modules], rotation=45, ha='right', fontsize=12)
        ax2.tick_params(axis='y', labelsize=12)
        ax2.grid(True, alpha=0.3, axis='y')
        
        # Add correlation examples as annotations
        for i, (mod, bar) in enumerate(zip(ranked_negative_modules, bars2)):
            # Find the most negative correlation
            all_corrs = np.concatenate(list(mod[1]['correlations'].values()))
            most_negative = min(all_corrs[~np.isnan(all_corrs)], default=0)
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height/2, 
                    f"min r={most_negative:.2f}", 
                    ha='center', va='center', color='white', fontweight='bold', fontsize=10)
    else:
        ax2.text(0.5, 0.5, 'No modules with negative correlations found', 
                ha='center', va='center', transform=ax2.transAxes, fontsize=16, fontweight='bold')
        ax2.set_title('Top miRNA Modules by Negative Correlations - None Found', 
                      fontsize=18, fontweight='bold')
        ax2.set_xlim(0, 1)
        ax2.set_ylim(0, 1)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.1, right=0.95, top=0.93, bottom=0.2)
    
    # Save with higher DPI and better format
    plt.savefig(f"{output_dir}/top_{top_n}_modules_summary.png", 
                dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    
    # Save to CSV - including negative correlations data
    top_df = pd.DataFrame({
        'miRNA': [mod[0] for mod in ranked_modules],
        'Target_Count': [len(mod[1]['targets']) for mod in ranked_modules],
        'Mean_Abs_Correlation': [mod[1]['mean_abs_corr'] for mod in ranked_modules],
        'Mean_Neg_Correlation': [mod[1]['mean_neg_corr'] for mod in ranked_modules],
        'Max_Abs_Correlation': [
            max(np.abs(np.concatenate(list(mod[1]['correlations'].values()))), default=0)
            for mod in ranked_modules
        ]
    })
    top_df.to_csv(f"{output_dir}/top_{top_n}_modules.csv", index=False)
    
    # Save negative correlations data separately if any exist
    if ranked_negative_modules:
        neg_df = pd.DataFrame({
            'miRNA': [mod[0] for mod in ranked_negative_modules],
            'Target_Count': [len(mod[1]['targets']) for mod in ranked_negative_modules],
            'Mean_Neg_Correlation': [mod[1]['mean_neg_corr'] for mod in ranked_negative_modules],
            'Most_Negative_Correlation': [
                min(np.concatenate(list(mod[1]['correlations'].values()))[
                    ~np.isnan(np.concatenate(list(mod[1]['correlations'].values())))], default=0)
                for mod in ranked_negative_modules
            ]
        })
        neg_df.to_csv(f"{output_dir}/top_{top_n}_negative_modules.csv", index=False)
        print(f"Saved negative correlations summary to '{output_dir}/top_{top_n}_negative_modules.csv'")
    
    print(f"Saved top {top_n} modules summary to '{output_dir}/top_{top_n}_modules.*'")

def main():
    parser = argparse.ArgumentParser(
        description="Enhanced miRNA module analysis with tables and visualizations."
    )
    parser.add_argument(
        "inputs", nargs=3, metavar=("MODULES", "MIRNAS", "TRANSCRIPTS"),
        help="Input files: mirna_modules_renamed.txt mirnas_CPM.txt transcripts_TPM.txt"
    )
    parser.add_argument(
        "--cleavage-sites", "-c", 
        help="Optional file with cleavage site information"
    )
    parser.add_argument(
        "--output-dir", "-o", default="module_heatmaps",
        help="Output directory for results (default: module_heatmaps)"
    )
    parser.add_argument(
        "--top-n", type=int, default=TOP_N_MODULES,
        help=f"Number of top modules to report (default: {TOP_N_MODULES})"
    )
    parser.add_argument(
        "--log-file", default="miRNA_module_analysis.log",
        help="Log file path (default: miRNA_module_analysis.log, will be inside output dir)"
    )

    args = parser.parse_args()

    # Log file will be inside the output directory
    log_file_path = os.path.join(args.output_dir, args.log_file)
    setup_logging(log_file_path, args.output_dir)

    try:
        print("=== miRNA Module Analyzer with Tables ===")
        print(f"Configuration:\n- Sample pairs: {SAMPLE_PAIRS}\n- Top modules: {args.top_n}\n")
        print(f"Inputs:\n  Modules: {args.inputs[0]}\n  miRNAs: {args.inputs[1]}\n  Transcripts: {args.inputs[2]}")
        print(f"Output directory: {args.output_dir}\n")

        modules, mirnas, transcripts, cleavage_sites = load_data(*args.inputs, args.cleavage_sites)
        module_data = prepare_module_data(modules, mirnas, transcripts, cleavage_sites)

        if not module_data:
            print("Error: No valid miRNA-target modules found!")
            return

        plot_individual_heatmaps(module_data, args.output_dir)
        generate_top_modules_report(module_data, args.top_n, args.output_dir)

        print("\n=== Analysis completed successfully ===")
        print(f"Log file saved to: {log_file_path}")

    except Exception as e:
        print(f"\n!!! Error: {str(e)} !!!", file=sys.stderr)
        raise
    finally:
        if isinstance(sys.stdout, Logger):
            sys.stdout.log.close()
            sys.stdout = sys.stdout.terminal

if __name__ == "__main__":
    main()
