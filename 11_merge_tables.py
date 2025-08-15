#!/usr/bin/env python3
import pandas as pd
import sys
from datetime import datetime
import logging
import argparse

def remove_redundant_rows(df):
    """Remove rows where all columns have identical values to previous row"""
    duplicates = df.duplicated()
    cleaned_df = df[~duplicates]
    num_removed = len(df) - len(cleaned_df)
    if num_removed > 0:
        logging.info(f"Removed {num_removed} redundant row(s)")
    return cleaned_df

def load_cpm_filter_list(cpm_file, column_name=None):
    """Load list of genes/transcripts from CPM table first column"""
    try:
        # Try reading with headers first
        if cpm_file.endswith('.csv'):
            cpm_df = pd.read_csv(cpm_file)
        else:
            cpm_df = pd.read_csv(cpm_file, sep='\t')
        
        # Check if first row looks like data rather than headers
        first_col_first_row = str(cpm_df.iloc[0, 0])
        if (first_col_first_row.startswith(('novel', 'ENSG', 'AT', 'Os', 'Zm')) or 
            any(char.isdigit() for char in first_col_first_row)):
            logging.info("CPM table appears to have no headers, re-reading without headers...")
            # Re-read without headers
            if cpm_file.endswith('.csv'):
                cpm_df = pd.read_csv(cpm_file, header=None)
            else:
                cpm_df = pd.read_csv(cpm_file, sep='\t', header=None)
            # Use first column (index 0)
            column_index = 0
            column_name = f"Column_{column_index}"
        else:
            # Use headers normally
            if column_name is None:
                column_index = 0
                column_name = cpm_df.columns[0]
            else:
                column_index = list(cpm_df.columns).index(column_name)
        
        # Extract the filter list from the appropriate column
        if isinstance(column_name, str) and column_name.startswith('Column_'):
            # No headers case
            filter_list = set(cpm_df.iloc[:, column_index].dropna().astype(str).unique())
        else:
            # Headers case
            filter_list = set(cpm_df[column_name].dropna().astype(str).unique())
            
        logging.info(f"Loaded {len(filter_list)} unique entries from CPM file column '{column_name}'")
        return filter_list, column_name
    except Exception as e:
        logging.error(f"Error reading CPM file: {e}")
        raise

def apply_cpm_filter(df, filter_list, target_column='Target ID'):
    """Filter dataframe based on CPM list"""
    original_count = len(df)
    
    # Try exact match first
    filtered_df = df[df[target_column].isin(filter_list)]
    
    # If no matches, try base gene ID matching (remove .1, .2 etc.)
    if len(filtered_df) == 0:
        logging.info("No exact matches found, trying base gene ID matching...")
        df_copy = df.copy()
        df_copy['base_id'] = df_copy[target_column].str.replace(r'\.\d+$', '', regex=True)
        filtered_df = df_copy[df_copy['base_id'].isin(filter_list)]
        if len(filtered_df) > 0:
            filtered_df = filtered_df.drop(columns=['base_id'])
    
    # If still no matches, try reverse - remove extensions from filter list
    if len(filtered_df) == 0:
        logging.info("Trying reverse matching (removing extensions from filter list)...")
        filter_base_list = {gene.split('.')[0] for gene in filter_list}
        df_copy = df.copy()
        df_copy['base_id'] = df_copy[target_column].str.replace(r'\.\d+$', '', regex=True)
        filtered_df = df_copy[df_copy['base_id'].isin(filter_base_list)]
        if len(filtered_df) > 0:
            filtered_df = filtered_df.drop(columns=['base_id'])
    
    filtered_count = len(filtered_df)
    logging.info(f"CPM filtering: {original_count} -> {filtered_count} entries ({filtered_count/original_count*100:.1f}% retained)")
    
    return filtered_df

def try_to_excel(df, filename):
    """Attempt to export to Excel with graceful fallback"""
    try:
        import openpyxl  # Try importing first
        df.to_excel(filename, index=False)
        logging.info(f"Excel output saved to {filename}")
        return True
    except ImportError:
        logging.warning("openpyxl not installed - skipping Excel output")
        return False
    except Exception as e:
        logging.error(f"Failed to create Excel file: {e}")
        return False

def generate_html_with_filters(df, filename):
    """Generate enhanced HTML with entry count and column filters"""
    # Get column names for filter generation
    columns = df.columns.tolist()
    
    # Create filter dropdowns HTML
    filter_html = ""
    for i, col in enumerate(columns):
        unique_values = sorted(df[col].dropna().astype(str).unique())
        options = "".join([f'<option value="{val}">{val}</option>' for val in unique_values])
        filter_html += f'''
        <div class="filter-group">
            <label for="filter_{i}">{col}:</label>
            <select id="filter_{i}" class="column-filter" data-column="{i}">
                <option value="">All</option>
                {options}
            </select>
            <div class="nan-filter-individual">
                <input type="checkbox" id="hide-nan_{i}" class="nan-checkbox-individual" data-column="{i}">
                <label for="hide-nan_{i}">Hide NaN/Empty</label>
            </div>
        </div>'''
    
    # Convert DataFrame to HTML table with proper IDs
    table_html = df.to_html(index=False, classes='dataframe', border=0, table_id='data-table')
    
    html_content = f'''<!DOCTYPE html>
<html>
<head>
    <title>miRNA Target Analysis Results</title>
    <style>
        body {{ 
            font-family: Arial, sans-serif; 
            margin: 20px;
            background-color: #f8f9fa;
        }}
        
        .header {{
            background-color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 20px;
        }}
        
        h1 {{ 
            color: #2e6c80; 
            margin: 0 0 10px 0;
        }}
        
        .timestamp {{ 
            color: #666; 
            font-size: 0.9em; 
            margin-bottom: 10px;
        }}
        
        .entry-count {{
            background-color: #e8f4f8;
            color: #2e6c80;
            padding: 10px 15px;
            border-radius: 5px;
            font-weight: bold;
            display: inline-block;
            margin-bottom: 20px;
        }}
        
        .filters-container {{
            background-color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 20px;
        }}
        
        .filters-title {{
            color: #2e6c80;
            font-weight: bold;
            margin-bottom: 15px;
            font-size: 1.1em;
        }}
        
        .filters-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
        }}
        
        .filter-group {{
            display: flex;
            flex-direction: column;
        }}
        
        .filter-group label {{
            font-weight: bold;
            color: #333;
            margin-bottom: 5px;
            font-size: 0.9em;
        }}
        
        .nan-filter-individual {{
            display: flex;
            align-items: center;
            gap: 5px;
            margin-top: 8px;
            padding: 6px 8px;
            background-color: #f8f9fa;
            border-radius: 3px;
            border: 1px solid #e9ecef;
        }}
        
        .nan-filter-individual input[type="checkbox"] {{
            width: 14px;
            height: 14px;
            cursor: pointer;
        }}
        
        .nan-filter-individual label {{
            margin: 0;
            cursor: pointer;
            font-size: 0.8em;
            color: #6c757d;
            font-weight: normal;
        }}
        
        .column-filter {{
            padding: 8px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 0.9em;
            background-color: white;
        }}
        
        .column-filter:focus {{
            outline: none;
            border-color: #2e6c80;
            box-shadow: 0 0 5px rgba(46, 108, 128, 0.3);
        }}
        
        .filter-controls {{
            display: flex;
            align-items: center;
            gap: 15px;
            margin-top: 15px;
            flex-wrap: wrap;
        }}
        
        .clear-filters {{
            background-color: #2e6c80;
            color: white;
            border: none;
            padding: 10px 20px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 0.9em;
        }}
        
        .clear-filters:hover {{
            background-color: #245a6b;
        }}
        
        .table-container {{
            background-color: white;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        
        table {{ 
            border-collapse: collapse; 
            width: 100%; 
        }}
        
        th {{ 
            background-color: #2e6c80; 
            color: white; 
            text-align: left; 
            padding: 12px 8px;
            font-weight: bold;
            position: sticky;
            top: 0;
            z-index: 10;
        }}
        
        td {{ 
            padding: 10px 8px; 
            border-bottom: 1px solid #eee;
            vertical-align: top;
        }}
        
        tr:nth-child(even) {{ 
            background-color: #f8f9fa; 
        }}
        
        tr:hover {{ 
            background-color: #e6f7ff; 
        }}
        
        .hidden-row {{
            display: none !important;
        }}
        
        .results-summary {{
            background-color: #fff3cd;
            border: 1px solid #ffeaa7;
            padding: 10px 15px;
            border-radius: 4px;
            margin-bottom: 15px;
            font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>miRNA Target Analysis Results</h1>
        <div class="timestamp">Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
        <div class="entry-count">Total Entries: {len(df)}</div>
    </div>
    
    <div class="filters-container">
        <div class="filters-title">Filter Data</div>
        <div class="results-summary" id="filter-summary" style="display: none;">
            Showing <span id="visible-count">{len(df)}</span> of {len(df)} entries
        </div>
        <div class="filters-grid">
            {filter_html}
        </div>
        <div class="filter-controls">
            <button class="clear-filters" onclick="clearAllFilters()">Clear All Filters</button>
        </div>
    </div>
    
    <div class="table-container">
        {table_html}
    </div>

    <script>
        // Initialize filter functionality
        document.addEventListener('DOMContentLoaded', function() {{
            const filters = document.querySelectorAll('.column-filter');
            const nanCheckboxes = document.querySelectorAll('.nan-checkbox-individual');
            const table = document.getElementById('data-table');
            const rows = table.querySelectorAll('tbody tr');
            const filterSummary = document.getElementById('filter-summary');
            const visibleCount = document.getElementById('visible-count');
            
            filters.forEach(filter => {{
                filter.addEventListener('change', applyFilters);
            }});
            
            nanCheckboxes.forEach(checkbox => {{
                checkbox.addEventListener('change', applyFilters);
            }});
            
            function isEmptyOrNaN(value) {{
                if (!value) return true;
                const trimmed = value.trim();
                return trimmed === '' || 
                       trimmed.toLowerCase() === 'nan' || 
                       trimmed.toLowerCase() === 'null' ||
                       trimmed === '-' ||
                       trimmed === 'n/a' ||
                       trimmed === 'na';
            }}
            
            function applyFilters() {{
                const activeFilters = [];
                const activeNaNFilters = [];
                
                filters.forEach(filter => {{
                    if (filter.value) {{
                        activeFilters.push({{
                            column: parseInt(filter.dataset.column),
                            value: filter.value.toLowerCase()
                        }});
                    }}
                }});
                
                nanCheckboxes.forEach(checkbox => {{
                    if (checkbox.checked) {{
                        activeNaNFilters.push(parseInt(checkbox.dataset.column));
                    }}
                }});
                
                let visibleRowCount = 0;
                
                rows.forEach(row => {{
                    let showRow = true;
                    
                    // Apply column filters
                    activeFilters.forEach(filterObj => {{
                        const cell = row.cells[filterObj.column];
                        if (cell && !cell.textContent.toLowerCase().includes(filterObj.value)) {{
                            showRow = false;
                        }}
                    }});
                    
                    // Apply individual NaN filters
                    activeNaNFilters.forEach(columnIndex => {{
                        const cell = row.cells[columnIndex];
                        if (cell && isEmptyOrNaN(cell.textContent)) {{
                            showRow = false;
                        }}
                    }});
                    
                    if (showRow) {{
                        row.classList.remove('hidden-row');
                        visibleRowCount++;
                    }} else {{
                        row.classList.add('hidden-row');
                    }}
                }});
                
                // Update summary
                if (activeFilters.length > 0 || activeNaNFilters.length > 0) {{
                    filterSummary.style.display = 'block';
                    visibleCount.textContent = visibleRowCount;
                }} else {{
                    filterSummary.style.display = 'none';
                }}
            }}
            
            window.clearAllFilters = function() {{
                filters.forEach(filter => {{
                    filter.value = '';
                }});
                nanCheckboxes.forEach(checkbox => {{
                    checkbox.checked = false;
                }});
                rows.forEach(row => {{
                    row.classList.remove('hidden-row');
                }});
                filterSummary.style.display = 'none';
            }}
        }});
    </script>
</body>
</html>'''
    
    with open(filename, 'w') as f:
        f.write(html_content)

def main():
    parser = argparse.ArgumentParser(description='Merge miRNA tables and generate reports')
    parser.add_argument('mirna_table', help='Input miRNA table file')
    parser.add_argument('function_table', help='Input function table file')
    parser.add_argument('annotation_table', help='Input annotation table file')
    
    # CPM filtering options
    parser.add_argument('--cpm-filter', help='CPM table file to filter results')
    parser.add_argument('--cpm-column', help='Column name in CPM file to use for filtering (default: first column)')
    parser.add_argument('--target-column', default='miRNA_Chr ID', 
                      help='Column in merged table to match against CPM list (default: miRNA_Chr ID)')
    
    # Output options
    parser.add_argument('--format', choices=['all', 'html', 'tsv', 'xlsx', 'json'], 
                      nargs='+', default=['all'],
                      help='Output format(s) to generate')
    parser.add_argument('--output-prefix', default='miRNA_results',
                      help='Prefix for output files')
    parser.add_argument('--skip-excel', action='store_true',
                      help='Skip Excel output even if requested')
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('miRNA_analysis.log'),
            logging.StreamHandler()
        ]
    )
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    try:
        logging.info("Reading input files...")
        mirna_df = pd.read_csv(args.mirna_table, sep='\t')
        function_df = pd.read_csv(args.function_table, sep='\t')
        annotation_df = pd.read_csv(args.annotation_table, sep='\t', comment='#')
    except Exception as e:
        logging.error(f"Error reading files: {e}")
        sys.exit(1)

    # Process tables
    function_df = function_df[['Setaria viridis ID', 'Function']].rename(columns={'Setaria viridis ID': 'Gene'})
    annotation_df = annotation_df.iloc[:, [1, 11, 13, 15]]
    annotation_df.columns = ['Gene', 'Best-hit-arabi-defline', 'Best-hit-clamy-defline', 'Best-hit-rice-defline']

    # Merge tables
    logging.info("Merging tables...")
    mirna_df['Target_base'] = mirna_df['Target ID'].str.replace(r'\.\d+$', '', regex=True)
    merged_df = pd.merge(mirna_df, function_df, on='Gene', how='left')
    merged_df = pd.merge(merged_df, annotation_df, on='Gene', how='left')

    # Handle Target_base matches
    temp_function = function_df.rename(columns={'Gene': 'Target_base'})
    merged_df_target = pd.merge(mirna_df, temp_function, on='Target_base', how='left')
    
    temp_annotation = annotation_df.copy()
    temp_annotation['Gene'] = temp_annotation['Gene'].str.replace(r'\.\d+$', '', regex=True)
    temp_annotation = temp_annotation.rename(columns={'Gene': 'Target_base'})
    merged_df_target = pd.merge(merged_df_target, temp_annotation, on='Target_base', how='left')
    
    for col in ['Function', 'Best-hit-arabi-defline', 'Best-hit-clamy-defline', 'Best-hit-rice-defline']:
        merged_df[col] = merged_df[col].combine_first(merged_df_target[col])

    # Clean up
    merged_df = merged_df.drop(columns=['Target_base'])
    merged_df = remove_redundant_rows(merged_df)

    # Apply CPM filtering if provided
    if args.cpm_filter:
        logging.info(f"Applying CPM filtering from {args.cpm_filter}")
        filter_list, filter_column = load_cpm_filter_list(args.cpm_filter, args.cpm_column)
        merged_df = apply_cpm_filter(merged_df, filter_list, args.target_column)
        
        if len(merged_df) == 0:
            logging.warning("No entries remained after CPM filtering!")
            sys.exit(1)
        
        # Update output prefix to indicate filtering
        args.output_prefix = f"{args.output_prefix}_CPM_filtered"

    # Generate outputs
    base_filename = f"{args.output_prefix}_{timestamp}"
    
    if 'all' in args.format or 'tsv' in args.format:
        tsv_file = f"{base_filename}.tsv"
        merged_df.to_csv(tsv_file, sep='\t', index=False)
        logging.info(f"TSV output saved to {tsv_file}")

    if 'all' in args.format or 'html' in args.format:
        html_file = f"{base_filename}.html"
        generate_html_with_filters(merged_df, html_file)
        logging.info(f"Enhanced HTML output saved to {html_file}")

    if ('all' in args.format or 'xlsx' in args.format) and not args.skip_excel:
        xlsx_file = f"{base_filename}.xlsx"
        if not try_to_excel(merged_df, xlsx_file):
            logging.info("Note: Install openpyxl with 'pip install openpyxl' for Excel support")

    if 'all' in args.format or 'json' in args.format:
        json_file = f"{base_filename}.json"
        merged_df.to_json(json_file, orient="records")
        logging.info(f"JSON output saved to {json_file}")

    logging.info("Processing complete!")

if __name__ == "__main__":
    main()
