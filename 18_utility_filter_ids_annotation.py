#!/usr/bin/env python3
"""
Script to filter an annotation table based on a list of IDs,
ignoring variant numbers after the last dot. Includes pretty output options.
"""

##USAGE:
# Basic usage: ./filter_ids_annotation.py -i id_list.txt -a annotation_table.tsv
# Generate HTML report: ./filter_ids_annotation.py -i id_list.txt -a annotation_table.tsv --html --format fancy_grid
# Set maximum width for better readability with long text: ./filter_ids_annotation.py -i id_list.txt -a annotation_table.tsv --max_width 40
# if your IDs are in the third column of annotation table (index 2): ./filter_ids_annotation.py -i id_list.txt -a annotation_table.tsv -c 2

import pandas as pd
import re
import argparse
import sys
from tabulate import tabulate
import os
import textwrap

def parse_arguments():
    parser = argparse.ArgumentParser(description='Filter annotation table based on ID list, ignoring variant numbers.')
    parser.add_argument('-i', '--id_list', required=True, help='Path to file containing list of IDs')
    parser.add_argument('-a', '--annotation', required=True, help='Path to annotation table file')
    parser.add_argument('-o', '--output', default='filtered_annotation.tsv', help='Output file path (default: filtered_annotation.tsv)')
    parser.add_argument('-s', '--separator', default='\t', help='Field separator in annotation table (default: tab)')
    parser.add_argument('-c', '--column', default=1, type=int, help='0-based index of the column containing Setaria IDs (default: 1)')
    parser.add_argument('--html', action='store_true', help='Generate an HTML file for prettier viewing')
    parser.add_argument('--no_console', action='store_true', help='Disable console table output')
    parser.add_argument('--format', default='pretty', choices=['plain', 'simple', 'github', 'grid', 'fancy_grid', 'pipe', 'orgtbl', 'jira'],
                        help='Table format for console output (default: pretty)')
    parser.add_argument('--max_width', type=int, default=0, help='Maximum width for table columns (0 for no limit)')

    return parser.parse_args()

def wrap_text_in_columns(df, max_width):
    """Wrap text in columns to improve readability"""
    if max_width <= 0:
        return df

    # Create a copy to avoid modifying the original
    wrapped_df = df.copy()

    # Apply wrapping to string columns only
    for col in wrapped_df.select_dtypes(include=['object']).columns:
        wrapped_df[col] = wrapped_df[col].apply(
            lambda x: textwrap.fill(str(x), width=max_width) if pd.notnull(x) else x
        )

    return wrapped_df

def generate_html(df, output_file):
    """Generate a styled HTML file for the DataFrame"""
    html_file = os.path.splitext(output_file)[0] + '.html'

    # Define CSS styles for better visualization
    styles = [
        dict(selector="table", props=[
            ("border-collapse", "collapse"),
            ("font-family", "Arial, sans-serif"),
            ("width", "100%"),
            ("margin", "20px 0"),
        ]),
        dict(selector="th", props=[
            ("background-color", "#4CAF50"),
            ("color", "white"),
            ("font-weight", "bold"),
            ("text-align", "left"),
            ("padding", "10px"),
            ("border", "1px solid #ddd"),
        ]),
        dict(selector="td", props=[
            ("padding", "8px"),
            ("border", "1px solid #ddd"),
            ("text-align", "left"),
        ]),
        dict(selector="tr:nth-child(even)", props=[
            ("background-color", "#f2f2f2"),
        ]),
        dict(selector="tr:hover", props=[
            ("background-color", "#ddd"),
        ]),
        dict(selector="caption", props=[
            ("font-size", "1.2em"),
            ("font-weight", "bold"),
            ("margin-bottom", "10px"),
            ("text-align", "left"),
        ]),
    ]

    # Create styled DataFrame
    styled_df = df.style.set_table_styles(styles)

    # Generate HTML with caption
    html_content = styled_df.set_caption(
        f"Filtered Annotation Table - {len(df)} entries"
    ).to_html()

    # Add viewport meta tag for mobile responsiveness
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Filtered Annotation Results</title>
</head>
<body>
{html_content}
</body>
</html>
"""

    # Write to file
    with open(html_file, 'w') as f:
        f.write(html_content)

    return html_file

def main():
    args = parse_arguments()

    # Read ID list
    try:
        with open(args.id_list, 'r') as f:
            id_list = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        sys.stderr.write(f"Error: ID list file '{args.id_list}' not found.\n")
        sys.exit(1)

    # Extract base IDs (remove variant number after last dot)
    base_ids = [re.sub(r'\.\d+$', '', id_) for id_ in id_list]
    print(f"Loaded {len(id_list)} IDs ({len(set(base_ids))} unique base IDs)")

    # Read annotation table
    try:
        df = pd.read_csv(args.annotation, sep=args.separator)
    except FileNotFoundError:
        sys.stderr.write(f"Error: Annotation table file '{args.annotation}' not found.\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"Error reading annotation table: {e}\n")
        sys.exit(1)

    # Get the actual column name
    try:
        column_name = df.columns[args.column]
    except IndexError:
        sys.stderr.write(f"Error: Column index {args.column} is out of bounds. Table has {len(df.columns)} columns.\n")
        sys.exit(1)

    print(f"Using column '{column_name}' for ID matching")

    # Extract base IDs from the specified column in the table
    df['Base_ID'] = df[column_name].astype(str).replace(r'\.\d+$', '', regex=True)

    # Filter rows where the base ID matches any in our list
    filtered_df = df[df['Base_ID'].isin(base_ids)]

    # Drop the temporary column we created
    filtered_df = filtered_df.drop(columns=['Base_ID'])

    # Save the result to TSV
    filtered_df.to_csv(args.output, sep=args.separator, index=False)
    print(f"Found {len(filtered_df)} matching entries out of {len(df)} total rows.")
    print(f"Results saved to {args.output}")

    # Generate HTML output if requested
    if args.html:
        html_file = generate_html(filtered_df, args.output)
        print(f"HTML report generated: {html_file}")

    # Display table in console if not disabled
    if not args.no_console and not filtered_df.empty:
        # Prepare the DataFrame for display
        display_df = wrap_text_in_columns(filtered_df, args.max_width)

        # Print the table
        print("\nFiltered Results:")
        print(tabulate(display_df, headers="keys", tablefmt=args.format, showindex=False))

if __name__ == "__main__":
    main()
