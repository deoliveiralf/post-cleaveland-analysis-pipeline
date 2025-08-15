#!/usr/bin/env python3
import os
import re
import argparse
import logging
from datetime import datetime
from pathlib import Path
import base64

def setup_logging():
    """Setup logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('html_png_integration.log'),
            logging.StreamHandler()
        ]
    )

def find_png_files(png_directory):
    """Find all PNG files in the directory and create a mapping"""
    png_files = {}
    png_dir = Path(png_directory)
    
    if not png_dir.exists():
        logging.error(f"PNG directory not found: {png_directory}")
        return png_files
    
    for png_file in png_dir.glob("*.png"):
        # Extract the base name by removing _analysis.png suffix
        filename_base = png_file.stem
        if filename_base.endswith("_analysis"):
            filename_base = filename_base[:-9]  # Remove "_analysis"
        
        # Store the exact format as primary key
        png_files[filename_base] = png_file
        
        logging.debug(f"Found PNG: {filename_base} -> {png_file}")
    
    unique_files = len(list(png_dir.glob('*.png')))
    logging.info(f"Found {unique_files} PNG files in {png_directory}")
    logging.info(f"Created {len(png_files)} mapping entries")
    
    # Show actual PNG keys found
    png_keys = list(png_files.keys())[:10]
    logging.info(f"Sample PNG file keys: {png_keys}")
    return png_files

def encode_image_to_base64(image_path):
    """Convert image to base64 string for embedding"""
    try:
        with open(image_path, "rb") as img_file:
            return base64.b64encode(img_file.read()).decode('utf-8')
    except Exception as e:
        logging.error(f"Error encoding image {image_path}: {e}")
        return None

def extract_mirna_ids_from_html(html_content):
    """Extract miRNA_Chr IDs from the HTML table to understand the format"""
    mirna_ids = set()
    
    # Find all table rows
    row_pattern = r'<tr[^>]*>(.*?)</tr>'
    rows = re.findall(row_pattern, html_content, re.DOTALL)
    
    for row in rows:
        # Extract all cells from the row
        cell_pattern = r'<td[^>]*>(.*?)</td>'
        cells = re.findall(cell_pattern, row, re.DOTALL)
        
        # Look for miRNA_Chr ID format in any column that contains miR/MIR and Chr
        for cell in cells:
            cell_text = re.sub(r'<[^>]+>', '', cell).strip()
            if cell_text and ('miR' in cell_text or 'MIR' in cell_text) and 'Chr' in cell_text:
                mirna_ids.add(cell_text)
                logging.debug(f"Found miRNA_Chr ID: {cell_text}")
    
    return mirna_ids

def create_image_modal_html():
    """Create HTML/CSS/JS for image modal popup"""
    return '''
    <style>
        .view-module-btn {
            background-color: #2e6c80;
            color: white;
            border: none;
            padding: 6px 12px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 0.8em;
            margin: 2px;
            transition: background-color 0.3s;
        }
        
        .view-module-btn:hover {
            background-color: #245a6b;
        }
        
        .view-module-btn:disabled {
            background-color: #ccc;
            cursor: not-allowed;
        }
        
        .no-module-text {
            color: #666;
            font-style: italic;
            font-size: 0.8em;
        }
        
        /* Modal styles */
        .modal {
            display: none;
            position: fixed;
            z-index: 1000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(0,0,0,0.8);
            overflow: auto;
        }
        
        .modal-content {
            background-color: #fefefe;
            margin: 2% auto;
            padding: 20px;
            border-radius: 8px;
            width: 95%;
            max-width: 1200px;
            position: relative;
            max-height: 90vh;
            overflow-y: auto;
        }
        
        .close {
            color: #aaa;
            float: right;
            font-size: 28px;
            font-weight: bold;
            position: absolute;
            right: 15px;
            top: 10px;
            cursor: pointer;
        }
        
        .close:hover,
        .close:focus {
            color: #000;
            text-decoration: none;
        }
        
        .modal-image {
            width: 100%;
            height: auto;
            display: block;
            margin: 10px auto;
        }
        
        .modal-title {
            color: #2e6c80;
            margin-bottom: 15px;
            padding-right: 40px;
        }
        
        .download-btn {
            background-color: #28a745;
            color: white;
            border: none;
            padding: 10px 20px;
            border-radius: 4px;
            cursor: pointer;
            margin-top: 10px;
            font-size: 0.9em;
        }
        
        .download-btn:hover {
            background-color: #218838;
        }
        
        .modal-loading {
            text-align: center;
            padding: 50px;
            color: #666;
        }
        
        .spinner {
            border: 4px solid #f3f3f3;
            border-top: 4px solid #2e6c80;
            border-radius: 50%;
            width: 40px;
            height: 40px;
            animation: spin 1s linear infinite;
            margin: 20px auto;
        }
        
        @keyframes spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }
    </style>
    
    <!-- Modal HTML -->
    <div id="imageModal" class="modal">
        <div class="modal-content">
            <span class="close">&times;</span>
            <h2 class="modal-title" id="modalTitle">miRNA Module Analysis</h2>
            <div class="modal-loading" id="modalLoading">
                <div class="spinner"></div>
                <p>Loading module analysis...</p>
            </div>
            <div id="modalImageContainer" style="display: none;">
                <img class="modal-image" id="modalImage" src="" alt="Module Analysis">
                <button class="download-btn" id="downloadBtn">Download Image</button>
            </div>
        </div>
    </div>
    
    <script>
        // Modal functionality
        const modal = document.getElementById('imageModal');
        const modalImage = document.getElementById('modalImage');
        const modalTitle = document.getElementById('modalTitle');
        const modalLoading = document.getElementById('modalLoading');
        const modalImageContainer = document.getElementById('modalImageContainer');
        const downloadBtn = document.getElementById('downloadBtn');
        const closeBtn = document.querySelector('.close');
        
        // Store PNG file mappings (will be populated by Python)
        const pngMappings = {PNG_MAPPINGS_PLACEHOLDER};
        
        function findImageKey(mirnaId) {
            console.log(`Looking for image key for: ${mirnaId}`);
            console.log('Available PNG keys:', Object.keys(pngMappings));
            
            // Try exact match first
            if (pngMappings[mirnaId]) {
                console.log(`Found exact match: ${mirnaId}`);
                return mirnaId;
            }
            
            // Try case-insensitive match
            const lowerMirna = mirnaId.toLowerCase();
            for (const key in pngMappings) {
                if (key.toLowerCase() === lowerMirna) {
                    console.log(`Found case-insensitive match: ${key}`);
                    return key;
                }
            }
            
            console.log(`No match found for: ${mirnaId}`);
            return null;
        }
        
        function showModuleImage(mirnaId) {
            modal.style.display = 'block';
            modalTitle.textContent = `miRNA Module Analysis: ${mirnaId}`;
            modalLoading.style.display = 'block';
            modalImageContainer.style.display = 'none';
            
            // Find the correct image key
            const imageKey = findImageKey(mirnaId);
            
            if (imageKey && pngMappings[imageKey]) {
                // Load the image
                const img = new Image();
                img.onload = function() {
                    modalImage.src = this.src;
                    modalLoading.style.display = 'none';
                    modalImageContainer.style.display = 'block';
                    
                    // Set up download functionality
                    downloadBtn.onclick = function() {
                        const link = document.createElement('a');
                        link.download = `${mirnaId.replace(/[\/\\]/g, '_')}_analysis.png`;
                        link.href = modalImage.src;
                        link.click();
                    };
                };
                img.onerror = function() {
                    modalLoading.innerHTML = '<p style="color: red;">Error loading image</p>';
                };
                
                // Check if it's base64 data or file path
                if (typeof pngMappings[imageKey] === 'string' && pngMappings[imageKey].length > 100) {
                    // Assume it's base64 data
                    img.src = `data:image/png;base64,${pngMappings[imageKey]}`;
                } else {
                    // Assume it's a file path
                    img.src = pngMappings[imageKey];
                }
            } else {
                modalLoading.innerHTML = 
                    `<p style="color: #666;">No module analysis available for: ${mirnaId}</p>
                     <p style="font-size: 0.8em; margin-top: 10px;">Available PNG files: ${Object.keys(pngMappings).slice(0,5).join(', ')}${Object.keys(pngMappings).length > 5 ? '...' : ''}</p>
                     <p style="font-size: 0.8em;">Total PNG files: ${Object.keys(pngMappings).length}</p>`;
            }
        }
        
        // Close modal events
        closeBtn.onclick = function() {
            modal.style.display = 'none';
        };
        
        window.onclick = function(event) {
            if (event.target === modal) {
                modal.style.display = 'none';
            }
        };
        
        // ESC key to close modal
        document.addEventListener('keydown', function(event) {
            if (event.key === 'Escape' && modal.style.display === 'block') {
                modal.style.display = 'none';
            }
        });
    </script>'''

def integrate_png_links(html_file, png_directory, output_file, embed_images=False):
    """Integrate PNG links into the HTML table"""
    
    # Find PNG files
    png_files = find_png_files(png_directory)
    
    if not png_files:
        logging.warning("No PNG files found. Proceeding without module links.")
    
    # Read the HTML file
    try:
        with open(html_file, 'r', encoding='utf-8') as f:
            html_content = f.read()
    except Exception as e:
        logging.error(f"Error reading HTML file {html_file}: {e}")
        return False
    
    # Extract miRNA IDs from HTML for debugging
    mirna_ids_in_html = extract_mirna_ids_from_html(html_content)
    logging.info(f"Found {len(mirna_ids_in_html)} unique miRNA IDs in HTML")
    
    # Show samples for comparison
    html_sample = sorted(list(mirna_ids_in_html))[:10]
    png_sample = sorted(list(png_files.keys()))[:10]
    
    logging.info("=== FORMAT COMPARISON ===")
    logging.info(f"HTML miRNA IDs (sample): {html_sample}")
    logging.info(f"PNG keys (sample): {png_sample}")
    
    # Check for exact matches
    exact_matches = set(mirna_ids_in_html) & set(png_files.keys())
    logging.info(f"Exact matches found: {len(exact_matches)}")
    if exact_matches:
        exact_sample = sorted(list(exact_matches))[:5]
        logging.info(f"Exact match examples: {exact_sample}")
    
    # Show mismatches for debugging
    html_only = set(mirna_ids_in_html) - set(png_files.keys())
    png_only = set(png_files.keys()) - set(mirna_ids_in_html)
    
    if html_only:
        html_only_sample = sorted(list(html_only))[:5]
        logging.info(f"miRNA IDs in HTML but no PNG (sample): {html_only_sample}")
    
    if png_only:
        png_only_sample = sorted(list(png_only))[:5]
        logging.info(f"PNG files with no HTML match (sample): {png_only_sample}")
    
    # Prepare PNG mappings for JavaScript
    png_mappings = {}
    if embed_images:
        logging.info("Embedding images as base64...")
        for mirna_name, png_path in png_files.items():
            base64_data = encode_image_to_base64(png_path)
            if base64_data:
                png_mappings[mirna_name] = base64_data
                logging.debug(f"Embedded image for {mirna_name}")
    else:
        # For non-embedded version, use relative paths
        html_dir = Path(html_file).parent
        for mirna_name, png_path in png_files.items():
            try:
                relative_path = png_path.relative_to(html_dir)
                png_mappings[mirna_name] = str(relative_path)
            except ValueError:
                # If relative path fails, use absolute path
                png_mappings[mirna_name] = str(png_path)
    
    # Create the modal HTML
    modal_html = create_image_modal_html()
    
    # Replace the PNG mappings placeholder
    png_mappings_js = str(png_mappings).replace("'", '"')
    modal_html = modal_html.replace('{PNG_MAPPINGS_PLACEHOLDER}', png_mappings_js)
    
    # Find the table and add the module column to header
    # Look for table header row
    header_patterns = [
        r'(<tr[^>]*>.*?<th.*?</tr>)',  # Row with th tags
        r'(<tr[^>]*>.*?</tr>)'        # First row (fallback)
    ]
    
    header_match = None
    for pattern in header_patterns:
        header_match = re.search(pattern, html_content, re.DOTALL)
        if header_match and '<th' in header_match.group(1):
            break
    
    if not header_match:
        logging.error("Could not find table header row")
        return False
    
    header_row = header_match.group(1)
    
    # Add new header for module column
    new_header = '<th>Module Analysis</th>'
    new_header_row = header_row.replace('</tr>', new_header + '</tr>')
    html_content = html_content.replace(header_row, new_header_row)
    
    # Keep track of matched entries for statistics
    matched_count = 0
    total_data_rows = 0
    
    # Find and process all data rows
    def add_module_column(match):
        nonlocal matched_count, total_data_rows
        row_content = match.group(1)
        total_data_rows += 1
        
        # Extract all cells from the row
        cell_pattern = r'<td[^>]*>(.*?)</td>'
        cells = re.findall(cell_pattern, row_content, re.DOTALL)
        
        if len(cells) < 1:
            return match.group(0)  # Return original if no cells
        
        # Look for miRNA_Chr ID in any cell (since we don't know exact column position)
        mirna_chr_id = None
        for cell in cells:
            cell_text = re.sub(r'<[^>]+>', '', cell).strip()
            if cell_text and ('miR' in cell_text or 'MIR' in cell_text) and 'Chr' in cell_text:
                mirna_chr_id = cell_text
                break
        
        # Create the module column content
        if mirna_chr_id:
            # Check if we have a matching PNG file
            if mirna_chr_id in png_files:
                module_cell = f'<td><button class="view-module-btn" onclick="showModuleImage(\'{mirna_chr_id}\')">View Module</button></td>'
                matched_count += 1
                logging.debug(f"Added button for: {mirna_chr_id}")
            else:
                module_cell = '<td><span class="no-module-text">No module</span></td>'
                logging.debug(f"No PNG found for: {mirna_chr_id}")
        else:
            module_cell = '<td><span class="no-module-text">-</span></td>'
        
        # Add the new cell to the row
        return f'<tr>{row_content}{module_cell}</tr>'
    
    # Find and process all data rows (rows with <td> tags)
    data_row_pattern = r'<tr[^>]*>([^<]*(?:<td[^>]*>.*?</td>[^<]*)+)</tr>'
    html_content = re.sub(data_row_pattern, add_module_column, html_content, flags=re.DOTALL)
    
    # Insert the modal HTML before the closing body tag
    body_end = html_content.rfind('</body>')
    if body_end != -1:
        html_content = html_content[:body_end] + modal_html + html_content[body_end:]
    else:
        # If no body tag, append at the end
        html_content += modal_html
    
    # Update the title and add integration info
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    integration_info = f'''
        <div class="integration-info" style="background-color: #d4edda; border: 1px solid #c3e6cb; padding: 10px; border-radius: 4px; margin-bottom: 15px;">
            <strong>Integration Info:</strong> Module analysis images integrated on {timestamp}.<br>
            Found {len(png_files)} PNG files, {matched_count}/{total_data_rows} table entries matched.<br>
            Match rate: {(matched_count/total_data_rows)*100:.1f}% - Click "View Module" buttons to see detailed analysis.
        </div>
    '''
    
    # Insert integration info after the first opening tag
    body_start = html_content.find('<body')
    if body_start != -1:
        body_tag_end = html_content.find('>', body_start) + 1
        html_content = html_content[:body_tag_end] + integration_info + html_content[body_tag_end:]
    else:
        # Try to insert after any existing div
        first_tag = re.search(r'<[^/!][^>]*>', html_content)
        if first_tag:
            insert_pos = first_tag.end()
            html_content = html_content[:insert_pos] + integration_info + html_content[insert_pos:]
        else:
            html_content = integration_info + html_content
    
    # Write the modified HTML
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        logging.info(f"Integrated HTML saved to {output_file}")
        
        # Detailed summary
        logging.info("=== INTEGRATION SUMMARY ===")
        logging.info(f"Total PNG files found: {len(png_files)}")
        logging.info(f"Total table rows processed: {total_data_rows}")
        logging.info(f"Successfully matched entries: {matched_count}")
        logging.info(f"Match rate: {(matched_count/total_data_rows)*100:.1f}%")
        
        if matched_count == 0:
            logging.error("NO MATCHES FOUND! Check format differences above.")
        elif matched_count < total_data_rows * 0.5:
            logging.warning("Low match rate detected. Check format comparison above.")
        else:
            logging.info("Good match rate achieved!")
        
        return True
    except Exception as e:
        logging.error(f"Error writing output file {output_file}: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Integrate PNG module files with HTML table')
    parser.add_argument('html_file', help='Input HTML file from merge_tables.py')
    parser.add_argument('png_directory', help='Directory containing PNG module files')
    parser.add_argument('--output', '-o', help='Output HTML file name', 
                       default=None)
    parser.add_argument('--embed-images', action='store_true',
                       help='Embed PNG images as base64 in HTML (creates larger file but fully self-contained)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Determine output filename
    if args.output is None:
        html_path = Path(args.html_file)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.output = html_path.parent / f"{html_path.stem}_with_modules_{timestamp}.html"
    
    logging.info("Starting HTML-PNG integration...")
    logging.info(f"Input HTML: {args.html_file}")
    logging.info(f"PNG Directory: {args.png_directory}")
    logging.info(f"Output HTML: {args.output}")
    logging.info(f"Embed Images: {args.embed_images}")
    
    # Check if input files exist
    if not os.path.exists(args.html_file):
        logging.error(f"HTML file not found: {args.html_file}")
        return 1
    
    if not os.path.exists(args.png_directory):
        logging.error(f"PNG directory not found: {args.png_directory}")
        return 1
    
    # Perform integration
    success = integrate_png_links(args.html_file, args.png_directory, args.output, args.embed_images)
    
    if success:
        logging.info("Integration completed successfully!")
        logging.info(f"Open {args.output} in a web browser to view the integrated results.")
        return 0
    else:
        logging.error("Integration failed!")
        return 1

if __name__ == "__main__":
    exit(main())
