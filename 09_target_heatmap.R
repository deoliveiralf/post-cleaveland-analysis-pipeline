# =============================================================================
# TARGET GENE HEATMAP GENERATOR - Based on miRNA targets (FIXED VERSION)
# =============================================================================
# This script generates heatmaps for target genes grouped by their miRNA
# It calculates logFC values and creates both expression and logFC heatmaps

# Load required libraries
library(edgeR)
library(ComplexHeatmap)
library(RColorBrewer)
library(viridis)
library(grid)
library(circlize)

# =============================================================================
# CONFIGURATION AND DATA LOADING
# =============================================================================

setwd("C:/Users/Leandro/OneDrive - usp.br/LIGNINLAB/Projeto FAPESP CNPq - Setaria miRNAs/R analysis/mirna_and_target_plots_EdgeR")
getwd()
dir()


# Create output directories
output_dir <- "Target_Gene_Heatmaps_by_miRNA"
logfc_output_dir <- "Target_Gene_LogFC_Heatmaps_by_miRNA"

for (dir in c(output_dir, logfc_output_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat(paste("Created output directory:", dir, "\n"))
  }
}

# Function to safely load data files
safe_load_file <- function(filename, required = TRUE) {
  if (!file.exists(filename)) {
    if (required) {
      stop(paste("Required file not found:", filename))
    } else {
      cat(paste("Optional file not found:", filename, "\n"))
      return(NULL)
    }
  }
  
  tryCatch({
    data <- read.table(filename, header = TRUE, row.names = 1, sep = "\t", fill = TRUE)
    cat(paste("Successfully loaded", filename, "- Dimensions:", nrow(data), "x", ncol(data), "\n"))
    return(data)
  }, error = function(e) {
    if (required) {
      stop(paste("Error loading required file:", filename, "\nError:", e$message))
    } else {
      cat(paste("Error loading optional file:", filename, "\nError:", e$message, "\n"))
      return(NULL)
    }
  })
}

# Load target gene expression data
cat("Loading target gene expression data...\n")
target_expression_file <- "transcripts_TPM.txt"  # Adjust filename as needed
if (!file.exists(target_expression_file)) {
  # Try alternative filenames
  possible_files <- c("target_genes_CPM.txt", "targets_TPM.txt", "targets_CPM.txt", 
                      "gene_expression_TPM.txt", "gene_expression_CPM.txt")
  
  file_found <- FALSE
  for (filename in possible_files) {
    if (file.exists(filename)) {
      target_expression_file <- filename
      file_found <- TRUE
      break
    }
  }
  
  if (!file_found) {
    cat("Available files in directory:\n")
    print(list.files(pattern = "*.txt"))
    stop(paste("Target gene expression file not found. Tried:", 
               paste(c("transcripts_TPM.txt", possible_files), collapse = ", ")))
  }
}

cat(paste("Using target expression file:", target_expression_file, "\n"))
target_cpm_data <- safe_load_file(target_expression_file, required = TRUE)

# Convert to matrix to avoid warnings
target_cpm_data <- as.matrix(target_cpm_data)
cat("Converted expression data to matrix format\n")

# Load miRNA-target relationships
cat("Loading miRNA-target relationships...\n")
target_table <- tryCatch({
  read.table("renamed_mirna_target_modules_table.txt", header = TRUE, sep = "\t")
}, error = function(e) {
  stop(paste("Error loading miRNA-target table:", e$message))
})

# Debug: Check the structure of target_table
cat("=== DEBUGGING TARGET TABLE ===\n")
cat("Target table dimensions:", nrow(target_table), "x", ncol(target_table), "\n")
cat("Available columns in target table:\n")
print(colnames(target_table))
cat("\nTarget table structure:\n")
str(target_table)
cat("\nFirst few rows:\n")
print(head(target_table))

# FIXED: Specify the correct column names based on your data
mirna_col <- "miRNA.ID"
target_col <- "Target.ID"  # This should contain your target gene identifiers

cat(paste("Using miRNA column:", mirna_col, "\n"))
cat(paste("Using target gene column:", target_col, "\n"))

# Verify columns exist
if (!mirna_col %in% colnames(target_table)) {
  stop(paste("miRNA column not found:", mirna_col))
}
if (!target_col %in% colnames(target_table)) {
  stop(paste("Target column not found:", target_col))
}

# Convert to character vectors to avoid factor issues
target_table[[mirna_col]] <- as.character(target_table[[mirna_col]])
target_table[[target_col]] <- as.character(target_table[[target_col]])

# Remove any rows with NA values
complete_rows <- complete.cases(target_table[, c(mirna_col, target_col)])
target_table <- target_table[complete_rows, ]

cat(paste("After cleaning:", nrow(target_table), "miRNA-target relationships\n"))
cat("Sample of cleaned target table:\n")
print(head(target_table[, c(mirna_col, target_col)]))

# Check if target IDs match expression data row names
target_ids <- unique(target_table[[target_col]])
expression_ids <- rownames(target_cpm_data)
matching_ids <- intersect(target_ids, expression_ids)

cat(paste("Total unique targets in table:", length(target_ids), "\n"))
cat(paste("Total genes in expression data:", length(expression_ids), "\n"))
cat(paste("Matching targets found:", length(matching_ids), "\n"))

if (length(matching_ids) == 0) {
  cat("\nWARNING: No matching targets found between target table and expression data!\n")
  cat("Sample target IDs from table:\n")
  print(head(target_ids, 10))
  cat("\nSample expression row names:\n")
  print(head(expression_ids, 10))
  cat("\nPlease check if target IDs need to be converted or mapped.\n")
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Function to extract group information from sample names
extract_sample_groups <- function(sample_names) {
  groups <- gsub("[0-9]$", "", sample_names)
  return(groups)
}

# Function to get samples for a specific comparison
get_comparison_samples <- function(comparison_name, all_samples) {
  parts <- strsplit(comparison_name, "_vs_")[[1]]
  if (length(parts) != 2) {
    stop(paste("Invalid comparison name:", comparison_name))
  }
  
  group1 <- parts[1]
  group2 <- parts[2]
  
  sample_groups <- extract_sample_groups(all_samples)
  samples_group1 <- all_samples[sample_groups == group1]
  samples_group2 <- all_samples[sample_groups == group2]
  
  comparison_samples <- c(samples_group1, samples_group2)
  comparison_groups <- c(rep(group1, length(samples_group1)), 
                         rep(group2, length(samples_group2)))
  
  return(list(samples = comparison_samples, groups = comparison_groups))
}

# Function to calculate logFC between groups
calculate_logfc <- function(expression_data, group1_samples, group2_samples) {
  # Ensure we're working with numeric data
  if (is.matrix(expression_data) || is.data.frame(expression_data)) {
    group1_values <- as.numeric(expression_data[group1_samples])
    group2_values <- as.numeric(expression_data[group2_samples])
  } else {
    group1_values <- as.numeric(expression_data[group1_samples])
    group2_values <- as.numeric(expression_data[group2_samples])
  }
  
  # Remove NA values
  group1_values <- group1_values[!is.na(group1_values)]
  group2_values <- group2_values[!is.na(group2_values)]
  
  if (length(group1_values) == 0 || length(group2_values) == 0) {
    return(NA)
  }
  
  # Calculate group means (in log2 space)
  group1_mean <- mean(group1_values, na.rm = TRUE)
  group2_mean <- mean(group2_values, na.rm = TRUE)
  
  # LogFC = log2(group2/group1) = log2(group2) - log2(group1)
  # Since data is already in log2(TPM+1), we can subtract directly
  logfc <- group2_mean - group1_mean
  
  return(logfc)
}

# Function to perform t-test between groups
calculate_pvalue <- function(expression_data, group1_samples, group2_samples) {
  # Ensure we're working with numeric data
  if (is.matrix(expression_data) || is.data.frame(expression_data)) {
    group1_values <- as.numeric(expression_data[group1_samples])
    group2_values <- as.numeric(expression_data[group2_samples])
  } else {
    group1_values <- as.numeric(expression_data[group1_samples])
    group2_values <- as.numeric(expression_data[group2_samples])
  }
  
  # Remove NA values
  group1_values <- group1_values[!is.na(group1_values)]
  group2_values <- group2_values[!is.na(group2_values)]
  
  if (length(group1_values) < 2 || length(group2_values) < 2) {
    return(NA)
  }
  
  # Check if variances are equal (roughly)
  if (var(group1_values) == 0 && var(group2_values) == 0) {
    return(if(mean(group1_values) == mean(group2_values)) 1 else 0)
  }
  
  tryCatch({
    t_test <- t.test(group2_values, group1_values)
    return(t_test$p.value)
  }, error = function(e) {
    return(NA)
  })
}

# Enhanced color palette function
create_enhanced_color_palette <- function(categories, palette_type = "groups") {
  n_categories <- length(unique(categories))
  
  if (palette_type == "groups") {
    if (n_categories == 1) {
      colors <- c("#A8E6CF")
    } else if (n_categories == 2) {
      colors <- c("#A8E6CF", "#FFE0B2")
    } else if (n_categories == 3) {
      colors <- c("#A8E6CF", "#FFE0B2", "#D1C4E9")
    } else if (n_categories == 4) {
      colors <- c("#A8E6CF", "#FFE0B2", "#D1C4E9", "#FFCDD2")
    } else {
      colors <- c("#A8E6CF", "#FFE0B2", "#D1C4E9", "#FFCDD2", "#B3E5FC", 
                  "#DCEDC8", "#FFF9C4", "#F8BBD0", "#E1F5FE", "#F3E5F5")[1:n_categories]
    }
  }
  
  names(colors) <- unique(categories)
  return(colors)
}

# Expression color scheme
create_expression_color_scheme <- function(expression_data) {
  data_range <- range(expression_data, na.rm = TRUE)
  
  # Handle cases where all values are the same or invalid
  if (is.infinite(data_range[1]) || is.infinite(data_range[2]) || data_range[1] == data_range[2]) {
    data_range <- c(0, 1)
  }
  
  range_extension <- (data_range[2] - data_range[1]) * 0.1
  if (range_extension == 0) range_extension <- 0.5
  extended_range <- c(data_range[1] - range_extension, data_range[2] + range_extension)
  
  colors <- c("#0D47A1", "#1976D2", "#42A5F5", "#90CAF9", "#E3F2FD", 
              "#FFFFFF", "#FFEBEE", "#FFCDD2", "#EF5350", "#E53935", "#C62828")
  breaks <- seq(extended_range[1], extended_range[2], length.out = length(colors))
  
  return(circlize::colorRamp2(breaks, colors))
}

# LogFC color scheme
create_logfc_color_scheme <- function(logfc_values) {
  # Remove NA values for calculation
  valid_logfc <- logfc_values[!is.na(logfc_values)]
  
  if (length(valid_logfc) == 0) {
    return(circlize::colorRamp2(c(-1, 0, 1), c("#BDBDBD", "#F5F5F5", "#BDBDBD")))
  }
  
  # Get the maximum absolute value for symmetric scaling
  max_abs_logfc <- max(abs(valid_logfc))
  
  # Ensure minimum range for visualization
  if (max_abs_logfc < 0.5) {
    max_abs_logfc <- 0.5
  }
  
  # Create symmetric breaks
  breaks <- c(-max_abs_logfc, -max_abs_logfc/2, 0, max_abs_logfc/2, max_abs_logfc)
  
  # Blue to red color scheme for logFC
  colors <- c("#0D47A1", "#64B5F6", "#FFFFFF", "#FF8A65", "#D32F2F")
  
  return(circlize::colorRamp2(breaks, colors))
}

# =============================================================================
# MAIN FUNCTION TO GENERATE TARGET HEATMAPS
# =============================================================================

generate_target_heatmaps <- function(target_tpm_data, target_table, 
                                     output_dir, logfc_output_dir,
                                     pvalue_threshold = 0.05) {
  
  # Define comparisons based on your specific groups
  all_samples <- colnames(target_tpm_data)
  sample_groups <- extract_sample_groups(all_samples)
  unique_groups <- unique(sample_groups)
  
  cat("Available sample groups:\n")
  print(table(sample_groups))
  
  # Define your specific comparisons
  comparisons <- list(
    "Mature_Leaves_vs_Young_Leaves" = c("Mature_Leaves", "Young_Leaves"),
    "Third_Internode_vs_Fifth_Internode" = c("Third_Internode", "Fifth_Internode"),
    "Mature_Roots_vs_Young_Roots" = c("Mature_Roots", "Young_Roots")
  )
  
  # Filter comparisons to only include those with available samples
  valid_comparisons <- list()
  for (comp_name in names(comparisons)) {
    groups <- comparisons[[comp_name]]
    if (all(groups %in% unique_groups)) {
      valid_comparisons[[comp_name]] <- groups
    } else {
      cat(paste("Skipping comparison", comp_name, "- groups not found in data\n"))
    }
  }
  
  if (length(valid_comparisons) == 0) {
    stop("No valid comparisons found based on available sample groups")
  }
  
  cat("Valid comparisons:\n")
  print(names(valid_comparisons))
  
  # Convert to log2(TPM+1) if not already
  log_tpm_data <- log2(target_tpm_data + 1)
  
  # Group targets by miRNA using the correct column names
  mirna_targets <- split(target_table[[target_col]], target_table[[mirna_col]])
  
  cat(paste("Found", length(mirna_targets), "miRNAs with targets\n"))
  
  # Summary results
  summary_results <- data.frame()
  
  # Process each miRNA and its targets
  for (mirna_id in names(mirna_targets)) {
    target_genes <- mirna_targets[[mirna_id]]
    target_genes <- unique(target_genes)  # Remove duplicates
    
    # Filter targets that exist in expression data
    available_targets <- intersect(target_genes, rownames(log_tpm_data))
    
    if (length(available_targets) == 0) {
      cat(paste("No expression data found for targets of", mirna_id, "\n"))
      next
    }
    
    cat(paste("\nProcessing", mirna_id, "with", length(available_targets), "targets\n"))
    
    # Extract expression data for these targets
    target_expression <- log_tpm_data[available_targets, , drop = FALSE]
    
    # Generate expression heatmaps for each comparison
    for (comp_name in names(valid_comparisons)) {
      cat(paste("  Creating expression heatmap for", comp_name, "\n"))
      
      tryCatch({
        # Get samples for this comparison
        comp_info <- get_comparison_samples(comp_name, all_samples)
        comp_samples <- comp_info$samples
        comp_groups <- comp_info$groups
        
        if (length(comp_samples) == 0) {
          cat(paste("    No samples found for comparison:", comp_name, "\n"))
          next
        }
        
        # Extract expression for comparison samples
        comp_target_expression <- target_expression[, comp_samples, drop = FALSE]
        
        # Create annotation
        annotation_col <- data.frame(
          Group = comp_groups,
          row.names = comp_samples
        )
        
        # Create color palette
        group_colors <- create_enhanced_color_palette(comp_groups, "groups")
        
        # Create filename
        clean_mirna <- gsub("[^A-Za-z0-9_-]", "_", mirna_id)
        clean_comp <- gsub("[^A-Za-z0-9_]", "_", comp_name)
        filename <- file.path(output_dir, paste0(clean_mirna, "_targets_", clean_comp, "_expression_heatmap.png"))
        
        # Determine plot dimensions
        n_targets <- nrow(comp_target_expression)
        n_samples <- ncol(comp_target_expression)
        
        plot_width <- max(8, n_samples * 0.5 + 4)
        plot_height <- max(6, n_targets * 0.3 + 3)
        
        # Create heatmap
        png(filename, width = plot_width, height = plot_height, units = "in", res = 300)
        
        # Always create heatmaps, even for single targets
        color_mapping <- create_expression_color_scheme(as.matrix(comp_target_expression))
        
        col_ha <- HeatmapAnnotation(
          Group = annotation_col$Group,
          col = list(Group = group_colors),
          annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
        )
        
        ht <- Heatmap(
          as.matrix(comp_target_expression),  # Ensure matrix format
          name = "log2(TPM+1)",
          col = color_mapping,
          top_annotation = col_ha,
          show_row_names = TRUE,
          show_column_names = TRUE,      
          row_names_gp = gpar(fontsize = max(8, min(12, 100/n_targets))),
          column_names_gp = gpar(fontsize = 10),
          column_title = paste(mirna_id, "Targets -", gsub("_", " ", comp_name)),
          column_title_gp = gpar(fontsize = 12, fontface = "bold"),
          row_title = paste("Target Genes (n=", n_targets, ")", sep = ""),
          row_title_gp = gpar(fontsize = 11, fontface = "bold"),
          cluster_rows = if(n_targets > 1) TRUE else FALSE,  # Don't cluster single targets
          cluster_columns = FALSE,
          rect_gp = gpar(col = "white", lwd = 0.5),
          border = TRUE,
          width = unit(max(4, n_samples * 0.5), "cm"),
          height = unit(max(2, n_targets * 0.5), "cm")
        )
        
        draw(ht)
        
        dev.off()
        cat(paste("    Expression heatmap saved:", filename, "\n"))
        
      }, error = function(e) {
        cat(paste("    Error creating expression heatmap:", e$message, "\n"))
      })
    }
    
    # Generate logFC heatmap
    cat(paste("  Creating logFC heatmap for", mirna_id, "\n"))
    
    tryCatch({
      # Calculate logFC for each comparison and target
      logfc_matrix <- matrix(NA, nrow = length(available_targets), 
                             ncol = length(valid_comparisons),
                             dimnames = list(available_targets, names(valid_comparisons)))
      
      pvalue_matrix <- matrix(NA, nrow = length(available_targets), 
                              ncol = length(valid_comparisons),
                              dimnames = list(available_targets, names(valid_comparisons)))
      
      for (comp_name in names(valid_comparisons)) {
        groups <- valid_comparisons[[comp_name]]
        group1_samples <- all_samples[sample_groups == groups[1]]
        group2_samples <- all_samples[sample_groups == groups[2]]
        
        if (length(group1_samples) > 0 && length(group2_samples) > 0) {
          for (target_gene in available_targets) {
            # Extract expression for specific target gene
            target_expr <- log_tpm_data[target_gene, ]
            
            # Calculate logFC
            logfc_val <- calculate_logfc(target_expr, group1_samples, group2_samples)
            logfc_matrix[target_gene, comp_name] <- logfc_val
            
            # Calculate p-value
            pval <- calculate_pvalue(target_expr, group1_samples, group2_samples)
            pvalue_matrix[target_gene, comp_name] <- pval
          }
        }
      }
      
      # Create significance matrix
      sig_matrix <- ifelse(pvalue_matrix < pvalue_threshold & !is.na(pvalue_matrix), "*", "")
      
      # Create filename
      clean_mirna <- gsub("[^A-Za-z0-9_-]", "_", mirna_id)
      filename <- file.path(logfc_output_dir, paste0(clean_mirna, "_logFC_heatmap.png"))
      
      # Prepare matrices for display
      display_logfc <- logfc_matrix
      display_logfc[is.na(display_logfc)] <- 0
      
      # Create heatmap
      n_targets <- nrow(display_logfc)
      n_comparisons <- ncol(display_logfc)
      
      plot_width <- max(8, n_comparisons * 1.5 + 4)
      plot_height <- max(6, n_targets * 0.3 + 3)
      
      png(filename, width = plot_width, height = plot_height, units = "in", res = 300)
      
      # Always create heatmaps, even for single targets
      logfc_colors <- create_logfc_color_scheme(as.vector(logfc_matrix))
      
      # Cell annotation function
      cell_fun <- function(j, i, x, y, width, height, fill) {
        if (i <= nrow(logfc_matrix) && j <= ncol(logfc_matrix)) {
          if (is.na(logfc_matrix[i, j])) {
            grid.text("NA", x, y, gp = gpar(fontsize = 8, col = "gray50"))
          } else if (sig_matrix[i, j] == "*") {
            grid.text("*", x, y, gp = gpar(fontsize = 12, fontface = "bold", col = "gold"))
          }
        }
      }
      
      ht <- Heatmap(
        as.matrix(display_logfc),  # Ensure matrix format
        name = "logFC",
        col = logfc_colors,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = max(8, min(12, 100/n_targets))),
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 45,
        column_title = paste(mirna_id, "Targets - Log Fold Changes"),
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        row_title = paste("Target Genes (n=", n_targets, ")", sep = ""),
        row_title_gp = gpar(fontsize = 11, fontface = "bold"),
        cluster_rows = if(n_targets > 1) TRUE else FALSE,  # Don't cluster single targets
        cluster_columns = FALSE,
        cell_fun = cell_fun,
        rect_gp = gpar(col = "white", lwd = 0.5),
        border = TRUE,
        width = unit(max(4, n_comparisons * 1.5), "cm"),
        height = unit(max(2, n_targets * 0.5), "cm")
      )
      
      draw(ht)
      
      # Add legend
      grid.text(paste("* = p-value <", pvalue_threshold), 
                x = 0.5, y = 0.08, 
                gp = gpar(fontsize = 10, fontface = "bold"))
      
      dev.off()
      cat(paste("    LogFC heatmap saved:", filename, "\n"))
      
      # Add to summary
      for (comp_name in names(valid_comparisons)) {
        for (target_gene in available_targets) {
          summary_row <- data.frame(
            miRNA_ID = mirna_id,
            Target_Gene = target_gene,
            Comparison = comp_name,
            logFC = ifelse(is.na(logfc_matrix[target_gene, comp_name]), NA, 
                           round(logfc_matrix[target_gene, comp_name], 3)),
            PValue = ifelse(is.na(pvalue_matrix[target_gene, comp_name]), NA, 
                            round(pvalue_matrix[target_gene, comp_name], 4)),
            Significant = ifelse(is.na(pvalue_matrix[target_gene, comp_name]), "No", 
                                 ifelse(pvalue_matrix[target_gene, comp_name] < pvalue_threshold, "Yes", "No")),
            stringsAsFactors = FALSE
          )
          
          if (nrow(summary_results) == 0) {
            summary_results <- summary_row
          } else {
            summary_results <- rbind(summary_results, summary_row)
          }
        }
      }
      
    }, error = function(e) {
      cat(paste("    Error creating logFC heatmap:", e$message, "\n"))
    })
  }
  
  # Save summary
  summary_file <- file.path(output_dir, "Target_Heatmaps_Summary.csv")
  write.csv(summary_results, summary_file, row.names = FALSE)
  cat(paste("\nSummary saved to:", summary_file, "\n"))
  
  # Print statistics
  cat("\n=== SUMMARY STATISTICS ===\n")
  cat(paste("Total target-comparison pairs processed:", nrow(summary_results), "\n"))
  if (nrow(summary_results) > 0) {
    cat(paste("Significant results:", sum(summary_results$Significant == "Yes", na.rm = TRUE), "\n"))
    cat("miRNA distribution:\n")
    print(table(summary_results$miRNA_ID))
  }
  
  return(summary_results)
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

cat("Starting target gene heatmap generation...\n")
cat("==========================================\n")

# Generate target heatmaps
summary_results <- generate_target_heatmaps(
  target_tpm_data = target_cpm_data,
  target_table = target_table,
  output_dir = output_dir,
  logfc_output_dir = logfc_output_dir,
  pvalue_threshold = 0.05
)

cat("\n==========================================\n")
cat("Target gene heatmap generation completed!\n")
cat(paste("Expression heatmaps: Check the", output_dir, "directory\n"))
cat(paste("LogFC heatmaps: Check the", logfc_output_dir, "directory\n"))
cat("\nGenerated outputs:\n")
cat("1. Expression heatmaps for target genes grouped by miRNA\n")
cat("2. LogFC heatmaps for target genes grouped by miRNA\n")
cat("3. Summary table with statistics\n")
cat("4. All targets shown with * marking significant ones\n")