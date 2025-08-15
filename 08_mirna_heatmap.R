# =============================================================================
# Set the directory and check
#==============================================================================
setwd("C:/Users/Leandro/OneDrive - usp.br/LIGNINLAB/Projeto FAPESP CNPq - Setaria miRNAs/R analysis/mirna_and_target_plots_EdgeR")
getwd()
dir()

# =============================================================================
# INDIVIDUAL miRNA HEATMAP GENERATOR - FIXED VERSION
# =============================================================================
# This script generates individual heatmaps for each miRNA found in the 
# renamed_mirna_target file, organized by EdgeR comparison results
# Fixed version that shows all results and handles missing data properly

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

# Create output directories (no subfolders for individual miRNAs)
output_dir <- "Individual_miRNA_Heatmaps_Enhanced"
logfc_output_dir <- "Individual_miRNA_LogFC_Heatmaps"

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
    cat(paste("Successfully loaded", filename, "\n"))
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

# Load CPM data with error handling
cat("Loading CPM data...\n")
cpm_data <- safe_load_file("renamed_normalized_counts_CPM.txt", required = TRUE)

# Load target miRNAs with error handling
cat("Loading target miRNAs...\n")
target_mirnas <- tryCatch({
  read.table("renamed_mirna_target_modules_table.txt", header = TRUE, sep = "\t")
}, error = function(e) {
  stop(paste("Error loading target miRNAs file:", e$message))
})

if (!"miRNA.ID" %in% colnames(target_mirnas)) {
  stop("Target miRNAs file must contain a column named 'miRNA.ID'")
}
target_mirna_ids <- unique(target_mirnas$miRNA.ID)
cat(paste("Found", length(target_mirna_ids), "unique target miRNAs\n"))

# Load DE results with error handling
cat("Loading differential expression results...\n")
de_files <- c(
  "Mature_Leaves_vs_Young_Leaves" = "DE_results_Mature_Leaves_vs_Young_Leaves_PValue_0.05.txt",
  "Mature_Roots_vs_Young_Roots" = "DE_results_Mature_Roots_vs_Young_Roots_PValue_0.05.txt",
  "Third_Internode_vs_Fifth_Internode" = "DE_results_Third_Internode_vs_Fifth_Internode_PValue_0.05.txt"
)

de_results <- list()
for (comp_name in names(de_files)) {
  filename <- de_files[comp_name]
  if (file.exists(filename)) {
    tryCatch({
      de_results[[comp_name]] <- read.table(filename, header = TRUE, row.names = 1, sep = "\t", fill = TRUE)
      cat(paste("Loaded", nrow(de_results[[comp_name]]), "results for", comp_name, "\n"))
    }, error = function(e) {
      cat(paste("Error loading", filename, ":", e$message, "\n"))
      de_results[[comp_name]] <- NULL
    })
  } else {
    cat(paste("Warning: File not found:", filename, "\n"))
  }
}

# Check if we have any DE results loaded
if (length(de_results) == 0) {
  stop("No differential expression results were successfully loaded. Cannot proceed.")
}

# =============================================================================
# FIXED HELPER FUNCTIONS
# =============================================================================

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

# Enhanced expression color scheme
create_expression_color_scheme <- function(expression_data, expression_category) {
  data_range <- range(expression_data, na.rm = TRUE)
  
  # Handle cases where all values are the same or invalid
  if (is.infinite(data_range[1]) || is.infinite(data_range[2]) || data_range[1] == data_range[2]) {
    data_range <- c(0, 1)
  }
  
  range_extension <- (data_range[2] - data_range[1]) * 0.1
  if (range_extension == 0) range_extension <- 0.5
  extended_range <- c(data_range[1] - range_extension, data_range[2] + range_extension)
  
  if (expression_category == "No expression (all zeros)") {
    colors <- c("#F5F5F5", "#E0E0E0", "#BDBDBD")
    breaks <- seq(extended_range[1], extended_range[2], length.out = 3)
  } else if (expression_category == "Constant expression") {
    colors <- c("#E3F2FD", "#90CAF9", "#1976D2")
    breaks <- seq(extended_range[1], extended_range[2], length.out = 3)
  } else if (expression_category == "Very low expression") {
    colors <- c("#FFF8E1", "#FFE082", "#FF8F00")
    breaks <- seq(extended_range[1], extended_range[2], length.out = 3)
  } else {
    colors <- c("#0D47A1", "#1976D2", "#42A5F5", "#90CAF9", "#E3F2FD", 
                "#FFFFFF", "#FFEBEE", "#FFCDD2", "#EF5350", "#E53935", "#C62828")
    breaks <- seq(extended_range[1], extended_range[2], length.out = length(colors))
  }
  
  return(circlize::colorRamp2(breaks, colors))
}

# FIXED: LogFC color scheme that handles NA values properly
create_logfc_color_scheme <- function(logfc_values) {
  # Remove NA values for calculation
  valid_logfc <- logfc_values[!is.na(logfc_values)]
  
  if (length(valid_logfc) == 0) {
    # Gray for no data - create a simple color scheme
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

# Function to format p-values for display
format_pvalue <- function(pval) {
  if (is.na(pval)) return("NA")
  if (pval < 0.001) return("< 0.001")
  if (pval < 0.01) return(sprintf("%.3f", pval))
  return(sprintf("%.3f", pval))
}

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

# FIXED: Function to determine significance status (shows all, marks significant)
get_significance_status <- function(mirna_id, comparison_name, de_results, 
                                    pvalue_threshold = 0.05, logfc_threshold = 0, 
                                    use_fdr = FALSE) {
  if (!comparison_name %in% names(de_results)) {
    return("Not tested")
  }
  
  if (!mirna_id %in% rownames(de_results[[comparison_name]])) {
    return("Not found")
  }
  
  result <- de_results[[comparison_name]][mirna_id, ]
  
  if (is.na(result$FDR) || is.na(result$logFC)) {
    return("No data")
  }
  
  # Show regulation direction regardless of significance
  if (result$logFC > logfc_threshold) {
    base_status <- "Upregulated"
  } else if (result$logFC < -logfc_threshold) {
    base_status <- "Downregulated"
  } else {
    base_status <- "No change"
  }
  
  # Add significance marker based on p-value (not FDR)
  pval <- if ("PValue" %in% colnames(result)) {
    result$PValue
  } else if ("pvalue" %in% colnames(result)) {
    result$pvalue
  } else if ("P.Value" %in% colnames(result)) {
    result$P.Value
  } else {
    NA
  }
  
  if (!is.na(pval) && pval < pvalue_threshold) {
    return(paste0(base_status, " *"))
  }
  
  return(base_status)
}

# FIXED: Function to create logFC matrix (handles all miRNAs, not just significant)
create_logfc_matrix <- function(mirna_ids, de_results) {
  comparison_names <- names(de_results)
  
  # Initialize matrix with NA values
  logfc_matrix <- matrix(NA, 
                         nrow = length(mirna_ids), 
                         ncol = length(comparison_names),
                         dimnames = list(mirna_ids, comparison_names))
  
  # Fill matrix with logFC values (all values, not just significant)
  for (comp_name in comparison_names) {
    if (comp_name %in% names(de_results)) {
      for (mirna_id in mirna_ids) {
        if (mirna_id %in% rownames(de_results[[comp_name]])) {
          logfc_val <- de_results[[comp_name]][mirna_id, "logFC"]
          # Include all logFC values, even if not significant
          if (!is.na(logfc_val)) {
            logfc_matrix[mirna_id, comp_name] <- logfc_val
          }
        }
      }
    }
  }
  
  return(logfc_matrix)
}

# FIXED: Function to create significance marker matrix (using p-value)
create_significance_matrix <- function(mirna_ids, de_results, pvalue_threshold = 0.05, use_fdr = FALSE) {
  comparison_names <- names(de_results)
  
  # Initialize matrix with empty strings
  sig_matrix <- matrix("", 
                       nrow = length(mirna_ids), 
                       ncol = length(comparison_names),
                       dimnames = list(mirna_ids, comparison_names))
  
  # Fill matrix with significance markers
  for (comp_name in comparison_names) {
    if (comp_name %in% names(de_results)) {
      for (mirna_id in mirna_ids) {
        if (mirna_id %in% rownames(de_results[[comp_name]])) {
          result <- de_results[[comp_name]][mirna_id, ]
          
          # Get p-value from different possible column names
          pval <- if ("PValue" %in% colnames(result)) {
            result$PValue
          } else if ("pvalue" %in% colnames(result)) {
            result$pvalue
          } else if ("P.Value" %in% colnames(result)) {
            result$P.Value
          } else {
            NA
          }
          
          if (!is.na(pval)) {
            if (pval < pvalue_threshold) {
              sig_matrix[mirna_id, comp_name] <- "*"  # Significant
            }
            # If not significant, leave empty (will show value without *)
          } else {
            sig_matrix[mirna_id, comp_name] <- "ND"  # No Data
          }
        } else {
          sig_matrix[mirna_id, comp_name] <- "NF"  # Not Found
        }
      }
    }
  }
  
  return(sig_matrix)
}

# =============================================================================
# FIXED MAIN HEATMAP GENERATION FUNCTION
# =============================================================================

generate_enhanced_mirna_heatmaps <- function(cpm_data, target_mirna_ids, de_results, 
                                             output_dir, logfc_output_dir,
                                             pvalue_threshold = 0.05, 
                                             logfc_threshold = 0,
                                             use_fdr = FALSE) {
  
  # Debug: Check input data structure
  cat("=== INPUT DATA VALIDATION ===\n")
  cat(paste("CPM data dimensions:", nrow(cpm_data), "x", ncol(cpm_data), "\n"))
  cat("Sample names:\n")
  print(colnames(cpm_data))
  cat("Sample groups extracted:\n")
  sample_groups <- extract_sample_groups(colnames(cpm_data))
  print(table(sample_groups))
  cat("\n")
  
  # Filter CPM data to only include target miRNAs
  target_mirnas_in_data <- intersect(target_mirna_ids, rownames(cpm_data))
  cat(paste("Found", length(target_mirnas_in_data), "target miRNAs in CPM data\n"))
  
  if (length(target_mirnas_in_data) == 0) {
    stop("No target miRNAs found in CPM data!")
  }
  
  # Convert CPM to log-CPM (add small pseudocount to avoid log(0))
  log_cpm_data <- log2(cpm_data + 1)
  
  # Create logFC matrices for heatmap visualization
  cat("\nCreating logFC matrices...\n")
  logfc_matrix <- create_logfc_matrix(target_mirnas_in_data, de_results)
  significance_matrix <- create_significance_matrix(target_mirnas_in_data, de_results, 
                                                    pvalue_threshold, use_fdr)
  
  # Debug: Check logFC matrix
  cat("LogFC matrix summary:\n")
  print(summary(as.vector(logfc_matrix)))
  cat("Number of non-NA values in logFC matrix:", sum(!is.na(logfc_matrix)), "\n")
  
  # Create summary data frame
  summary_results <- data.frame()
  
  # Generate enhanced heatmap for each target miRNA (original expression-based)
  for (mirna_id in target_mirnas_in_data) {
    cat(paste("\nProcessing miRNA:", mirna_id, "\n"))
    
    # Extract expression data for this miRNA
    mirna_expression <- log_cpm_data[mirna_id, ]
    mirna_expression <- as.numeric(mirna_expression)
    names(mirna_expression) <- colnames(log_cpm_data)
    
    # Handle NA values
    if (any(is.na(mirna_expression))) {
      cat(paste("Warning: Found NA values in", mirna_id, "- replacing with 0\n"))
      mirna_expression[is.na(mirna_expression)] <- 0
    }
    
    # Check expression characteristics
    expr_var <- var(mirna_expression, na.rm = TRUE)
    expr_max <- max(mirna_expression, na.rm = TRUE)
    expr_mean <- mean(mirna_expression, na.rm = TRUE)
    
    # Determine expression category
    if (all(mirna_expression == 0, na.rm = TRUE)) {
      expr_category <- "No expression (all zeros)"
    } else if (is.na(expr_var) || expr_var == 0) {
      expr_category <- "Constant expression"
    } else if (expr_max < 1) {
      expr_category <- "Very low expression"
    } else {
      expr_category <- "Variable expression"
    }
    
    cat(paste("  Expression category:", expr_category, 
              "| Mean:", round(expr_mean, 3), 
              "| Max:", round(expr_max, 3), 
              "| Var:", round(expr_var, 3), "\n"))
    
    # Generate enhanced heatmap for each comparison (expression-based)
    for (comparison_name in names(de_results)) {
      cat(paste("  Creating enhanced heatmap for comparison:", comparison_name, "\n"))
      
      # Get samples for this comparison
      tryCatch({
        comp_info <- get_comparison_samples(comparison_name, colnames(cpm_data))
        comp_samples <- comp_info$samples
        comp_groups <- comp_info$groups
        
        if (length(comp_samples) == 0) {
          cat(paste("    No samples found for comparison:", comparison_name, "\n"))
          next
        }
        
        # Extract expression data for comparison samples
        comp_expression <- mirna_expression[comp_samples]
        comp_expression <- as.numeric(comp_expression)
        names(comp_expression) <- comp_samples
        
        # Create annotation data frame
        annotation_col <- data.frame(
          Group = comp_groups,
          row.names = comp_samples
        )
        
        # Create enhanced color palette for groups
        group_colors <- create_enhanced_color_palette(comp_groups, "groups")
        
        # Get significance status and statistics
        sig_status <- get_significance_status(mirna_id, comparison_name, de_results, 
                                              pvalue_threshold, logfc_threshold, use_fdr)
        
        # Get comprehensive statistics if available
        logfc_val <- NA
        fdr_val <- NA
        pval_val <- NA
        stats_text <- ""
        
        if (comparison_name %in% names(de_results) && 
            mirna_id %in% rownames(de_results[[comparison_name]])) {
          result <- de_results[[comparison_name]][mirna_id, ]
          logfc_val <- result$logFC
          fdr_val <- result$FDR
          
          # Try to get p-value from different possible column names
          if ("PValue" %in% colnames(result)) {
            pval_val <- result$PValue
          } else if ("pvalue" %in% colnames(result)) {
            pval_val <- result$pvalue
          } else if ("P.Value" %in% colnames(result)) {
            pval_val <- result$P.Value
          }
          
          stats_text <- sprintf("logFC: %.3f | FDR: %s | P-value: %s", 
                                logfc_val, 
                                format_pvalue(fdr_val),
                                format_pvalue(pval_val))
        }
        
        # Create filename (directly in output_dir, no subfolder)
        clean_mirna <- gsub("[^A-Za-z0-9_]", "_", mirna_id)
        clean_comp <- gsub("[^A-Za-z0-9_]", "_", comparison_name)
        filename <- file.path(output_dir, paste0(clean_mirna, "_", clean_comp, "_enhanced_heatmap.png"))
        
        # Create enhanced color scheme for expression
        color_mapping <- create_expression_color_scheme(comp_expression, expr_category)
        
        # Determine legend title based on expression category
        legend_title <- switch(expr_category,
                               "No expression (all zeros)" = "log2(CPM+1)\n(No expression)",
                               "Constant expression" = "log2(CPM+1)\n(Constant)",
                               "Very low expression" = "log2(CPM+1)\n(Low expression)",
                               "log2(CPM+1)")
        
        # Create enhanced visualization
        if (length(comp_samples) <= 2) {
          # Enhanced bar plot for very few samples
          png(filename, width = 10, height = 8, units = "in", res = 300)
          par(mar = c(12, 6, 8, 3), bg = "white")
          
          # Choose colors based on expression category and groups
          if (expr_category == "No expression (all zeros)") {
            bar_colors <- rep("#BDBDBD", length(comp_expression))
          } else {
            bar_colors <- group_colors[comp_groups]
          }
          
          # Create bar plot with enhanced styling
          bp <- barplot(comp_expression, 
                        main = "",
                        ylab = "log2(CPM + 1)", 
                        las = 2,
                        col = bar_colors,
                        border = "white",
                        cex.lab = 1.3,
                        cex.axis = 1.1,
                        ylim = c(0, max(comp_expression) * 1.15))
          
          # Add enhanced title with multiple lines
          title(main = paste(mirna_id, "\n", gsub("_", " ", comparison_name)), 
                cex.main = 1.4, font.main = 2, line = 6)
          title(main = paste("Status:", sig_status, " | ", expr_category), 
                cex.main = 1.0, font.main = 1, line = 4.5, col = "gray30")
          
          # Add statistics text with enhanced formatting
          if (stats_text != "") {
            mtext(stats_text, side = 1, line = 9, cex = 1.0, font = 2, col = "darkblue")
          }
          
          # Add expression summary
          expr_info <- sprintf("Mean: %.3f | Range: %.3f - %.3f | Samples: %d", 
                               mean(comp_expression), min(comp_expression), 
                               max(comp_expression), length(comp_expression))
          mtext(expr_info, side = 1, line = 10.5, cex = 0.9, col = "gray40")
          
          # Add enhanced legend
          if (expr_category != "No expression (all zeros)" && length(unique(comp_groups)) > 1) {
            legend("topright", legend = names(group_colors), fill = group_colors, 
                   cex = 1.0, bg = "white", box.col = "gray80", 
                   title = "Groups", title.col = "black")
          }
          
          dev.off()
          
        } else {
          # Create enhanced ComplexHeatmap
          heatmap_matrix <- matrix(comp_expression, 
                                   nrow = 1, 
                                   dimnames = list(mirna_id, names(comp_expression)))
          
          # Ensure matrix has proper numeric values
          if (any(is.na(heatmap_matrix))) {
            heatmap_matrix[is.na(heatmap_matrix)] <- 0
          }
          
          # Create enhanced heatmap
          png(filename, width = 10, height = 8, units = "in", res = 300)
          
          tryCatch({
            # Enhanced column annotation with better styling
            col_ha <- HeatmapAnnotation(
              Group = annotation_col$Group,
              col = list(Group = group_colors),
              annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
              annotation_legend_param = list(
                Group = list(
                  title_gp = gpar(fontsize = 11, fontface = "bold"), 
                  labels_gp = gpar(fontsize = 10),
                  grid_height = unit(4, "mm"),
                  grid_width = unit(4, "mm")
                )
              ),
              border = TRUE
            )
            
            # Create enhanced heatmap
            ht <- Heatmap(
              heatmap_matrix,
              name = legend_title,
              col = color_mapping,
              top_annotation = col_ha,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 13, fontface = "bold"),
              column_names_gp = gpar(fontsize = 10),
              column_title = paste(mirna_id, "Expression"),
              column_title_gp = gpar(fontsize = 14, fontface = "bold"),
              rect_gp = gpar(col = "white", lwd = 1.5),
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 11, fontface = "bold"),
                labels_gp = gpar(fontsize = 10),
                grid_height = unit(4, "mm"),
                grid_width = unit(4, "mm"),
                legend_height = unit(6, "cm")
              ),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              width = unit(max(6, ncol(heatmap_matrix) * 1.0), "cm"),
              height = unit(2, "cm"),
              border = TRUE
            )
            
            # Draw heatmap
            draw(ht)
            
            # Add enhanced annotations at the bottom
            grid.text(gsub("_", " ", comparison_name), x = 0.5, y = 0.22, 
                      gp = gpar(fontsize = 12, fontface = "bold"))
            
            # Status and category
            status_text <- paste("Status:", sig_status, " | ", expr_category)
            grid.text(status_text, x = 0.5, y = 0.18, 
                      gp = gpar(fontsize = 11, col = "gray30"))
            
            # Statistics with enhanced formatting
            if (stats_text != "") {
              grid.text(stats_text, x = 0.5, y = 0.14, 
                        gp = gpar(fontsize = 11, fontface = "bold", col = "darkblue"))
            }
            
            # Expression summary
            expr_summary <- sprintf("Mean: %.3f | Max: %.3f | Variance: %.3f | Samples: %d",
                                    mean(comp_expression, na.rm = TRUE), 
                                    max(comp_expression, na.rm = TRUE),
                                    var(comp_expression, na.rm = TRUE),
                                    length(comp_expression))
            grid.text(expr_summary, x = 0.5, y = 0.10, 
                      gp = gpar(fontsize = 10, col = "gray40"))
            
          }, error = function(e) {
            cat(paste("      Error creating enhanced ComplexHeatmap:", e$message, "\n"))
            # Enhanced fallback plot
            plot(1:length(comp_expression), comp_expression, 
                 type = "b", pch = 19, col = group_colors[comp_groups], lwd = 2,
                 main = paste(mirna_id, "\n", gsub("_", " ", comparison_name), 
                              "\n", sig_status, "\n(", expr_category, ")"),
                 xlab = "Samples", ylab = "log2(CPM + 1)",
                 xaxt = "n", cex = 1.2, cex.main = 1.1)
            axis(1, at = 1:length(comp_expression), labels = names(comp_expression), las = 2, cex.axis = 0.9)
            
            # Add grid for better readability
            grid(col = "lightgray", lty = "dotted")
            
            # Re-plot points on top of grid
            points(1:length(comp_expression), comp_expression, 
                   pch = 19, col = group_colors[comp_groups], cex = 1.2)
            lines(1:length(comp_expression), comp_expression, col = "gray60", lwd = 1)
            
            # Add legend for fallback plot
            if (length(unique(comp_groups)) > 1) {
              legend("topright", legend = names(group_colors), col = group_colors, 
                     pch = 19, cex = 0.9, bg = "white")
            }
            
            # Add statistics text
            if (stats_text != "") {
              mtext(stats_text, side = 1, line = 6, cex = 0.9, col = "darkblue")
            }
          })
          
          dev.off()
        }
        
        cat(paste("    Enhanced heatmap saved:", filename, "\n"))
        
        # Add to summary with enhanced information
        summary_row <- data.frame(
          miRNA_ID = mirna_id,
          Comparison = comparison_name,
          Expression_Category = expr_category,
          Significance_Status = sig_status,
          N_Samples = length(comp_samples),
          Groups = paste(unique(comp_groups), collapse = " vs "),
          Mean_Expression = round(mean(comp_expression, na.rm = TRUE), 3),
          Max_Expression = round(max(comp_expression, na.rm = TRUE), 3),
          Expression_Variance = round(var(comp_expression, na.rm = TRUE), 3),
          Expression_Range = paste(round(range(comp_expression, na.rm = TRUE), 3), collapse = " - "),
          logFC = ifelse(is.na(logfc_val), NA, round(logfc_val, 3)),
          FDR = ifelse(is.na(fdr_val), NA, round(fdr_val, 4)),
          PValue = ifelse(is.na(pval_val), NA, round(pval_val, 4)),
          Filename = basename(filename),
          stringsAsFactors = FALSE
        )
        
        # Add to results
        if (!exists("summary_results") || nrow(summary_results) == 0) {
          summary_results <- summary_row
        } else {
          # Ensure consistent columns
          missing_cols_in_summary <- setdiff(names(summary_row), names(summary_results))
          missing_cols_in_row <- setdiff(names(summary_results), names(summary_row))
          
          if (length(missing_cols_in_summary) > 0) {
            for (col in missing_cols_in_summary) {
              summary_results[[col]] <- NA
            }
          }
          if (length(missing_cols_in_row) > 0) {
            for (col in missing_cols_in_row) {
              summary_row[[col]] <- NA
            }
          }
          
          summary_row <- summary_row[names(summary_results)]
          summary_results <- rbind(summary_results, summary_row)
        }
        
      }, error = function(e) {
        cat(paste("Error processing", mirna_id, "for", comparison_name, ":", e$message, "\n"))
      })
    }
    
    # Generate logFC-based heatmap for this miRNA
    cat(paste("  Creating logFC heatmap for", mirna_id, "\n"))
    
    # Extract logFC values for this miRNA
    mirna_logfc <- logfc_matrix[mirna_id, ]
    mirna_sig <- significance_matrix[mirna_id, ]
    
    # Check if we have any valid logFC values
    valid_logfc <- !is.na(mirna_logfc)
    if (!any(valid_logfc)) {
      cat(paste("    No valid logFC data for", mirna_id, "- skipping logFC heatmap\n"))
      next
    }
    
    # Create filename (directly in logfc_output_dir, no subfolder)
    clean_mirna <- gsub("[^A-Za-z0-9_]", "_", mirna_id)
    filename <- file.path(logfc_output_dir, paste0(clean_mirna, "_logFC_heatmap.png"))
    
    # Prepare matrix for heatmap
    heatmap_matrix <- matrix(mirna_logfc, nrow = 1, 
                             dimnames = list(mirna_id, names(mirna_logfc)))
    
    # Create display matrix (replace NA with 0 for visualization)
    display_matrix <- heatmap_matrix
    display_matrix[is.na(display_matrix)] <- 0
    
    # Create logFC color scheme
    logfc_colors <- create_logfc_color_scheme(mirna_logfc)
    
    # Get statistics for annotation
    stats_text <- ""
    for (comp_name in names(mirna_logfc)) {
      if (!is.na(mirna_logfc[comp_name]) && comp_name %in% names(de_results)) {
        if (mirna_id %in% rownames(de_results[[comp_name]])) {
          result <- de_results[[comp_name]][mirna_id, ]
          pval_val <- if ("PValue" %in% colnames(result)) result$PValue else 
            if ("pvalue" %in% colnames(result)) result$pvalue else
              if ("P.Value" %in% colnames(result)) result$P.Value else NA
          logfc_val <- result$logFC
          
          if (!is.na(logfc_val)) {
            sig_symbol <- if (!is.na(pval_val) && pval_val < pvalue_threshold) "*" else ""
            stats_text <- paste(stats_text, 
                                sprintf("%s: %.2f%s", gsub("_", " ", comp_name), logfc_val, sig_symbol), 
                                sep = " | ")
          }
        }
      }
    }
    stats_text <- sub("^\\s*\\|\\s*", "", stats_text)
    
    # Calculate dynamic dimensions
    n_comparisons <- length(mirna_logfc)
    base_width <- 8  # inches
    base_height <- 6 # inches
    comparison_width <- 0.5 # inches per comparison
    
    # Create heatmap with improved dimensions
    png(filename, 
        width = max(base_width, n_comparisons * comparison_width), 
        height = base_height,
        units = "in", 
        res = 300)
    
    tryCatch({
      # Create cell annotation function
      cell_fun <- function(j, i, x, y, width, height, fill) {
        if (i <= nrow(heatmap_matrix) && j <= ncol(heatmap_matrix)) {
          original_val <- heatmap_matrix[i, j]
          sig_val <- mirna_sig[j]
          
          if (is.na(original_val)) {
            grid.text("NA", x, y, gp = gpar(fontsize = 10, fontface = "bold", col = "gray30"))
          } else {
            grid.text(sprintf("%.2f", original_val), x, y, 
                      gp = gpar(fontsize = 9, fontface = "bold", 
                                col = if(abs(original_val) > 1) "white" else "black"))
            
            if (!is.na(sig_val) && sig_val == "*") {
              grid.text("*", x, y + unit(0.15, "npc"), 
                        gp = gpar(fontsize = 12, fontface = "bold", col = "yellow"))
            }
          }
        }
      }
      
      # Create the heatmap with improved layout
      ht <- Heatmap(
        display_matrix,
        name = "logFC",
        col = logfc_colors,
        show_row_names = FALSE,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 10, fontface = "bold"),
        column_names_rot = 45,
        column_title = paste(mirna_id, "- Log Fold Change"),
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        rect_gp = gpar(col = "white", lwd = 1.5),
        cell_fun = cell_fun,
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 10, fontface = "bold"),
          labels_gp = gpar(fontsize = 8),
          grid_height = unit(4, "mm"),
          grid_width = unit(4, "mm"),
          legend_height = unit(3, "cm"),
          at = c(-2, -1, 0, 1, 2),
          labels = c("-2", "-1", "0", "1", "2")
        ),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        width = unit(ncol(display_matrix) * 0.5, "in"),
        height = unit(1, "in"),
        border = TRUE
      )
      
      # Draw heatmap
      draw(ht)
      
      # Add annotations with improved layout
      grid.text(paste("Log Fold Change values (* = p-value <", pvalue_threshold, ")"), 
                x = 0.5, y = 0.20, 
                gp = gpar(fontsize = 10, fontface = "bold"))
      
      if (stats_text != "") {
        # Split long text into multiple lines if needed
        if (nchar(stats_text) > 60) {
          stats_parts <- strsplit(stats_text, " \\| ")[[1]]
          for (i in seq_along(stats_parts)) {
            grid.text(stats_parts[i], x = 0.5, y = 0.15 - (i-1)*0.05, 
                      gp = gpar(fontsize = 9, col = "darkblue", fontface = "bold"))
          }
        } else {
          grid.text(stats_text, x = 0.5, y = 0.15, 
                    gp = gpar(fontsize = 9, col = "darkblue", fontface = "bold"))
        }
      }
      
      # Add color interpretation
      grid.text("Blue: Downregulated | White: No change | Red: Upregulated", 
                x = 0.5, y = 0.08, 
                gp = gpar(fontsize = 8, col = "gray40"))
      
      # Add summary
      n_significant <- sum(mirna_sig == "*", na.rm = TRUE)
      n_total <- sum(!is.na(mirna_logfc))
      summary_info <- sprintf("Significant comparisons: %d out of %d", n_significant, n_total)
      grid.text(summary_info, x = 0.5, y = 0.03, 
                gp = gpar(fontsize = 8, col = "gray50"))
      
    }, error = function(e) {
      cat(paste("Error creating ComplexHeatmap for", mirna_id, ":", e$message, "\n"))
      cat("Creating fallback heatmap using base R...\n")
      
      # FIXED: Enhanced fallback using base R - properly handle single row matrix
      par(mar = c(10, 6, 6, 4), bg = "white")
      
      # FIXED: Check if we have valid data for plotting
      if (length(mirna_logfc) == 0 || all(is.na(mirna_logfc))) {
        plot(1, 1, type = "n", 
             main = paste(mirna_id, "\nNo LogFC Data Available"),
             xlab = "", ylab = "logFC", axes = FALSE)
        text(1, 1, "No data", cex = 1.5, col = "gray50")
        dev.off()
        return()
      }
      
      # FIXED: Create a simple bar plot instead of image() for single miRNA
      # Convert to numeric and handle NA values
      plot_values <- as.numeric(mirna_logfc)
      plot_names <- names(mirna_logfc)
      
      # Replace NA with 0 for plotting, but track where they were
      na_positions <- is.na(plot_values)
      plot_values[na_positions] <- 0
      
      # Create color vector based on values
      colors_vec <- rep("lightgray", length(plot_values))
      colors_vec[plot_values > 0] <- "#FF6B6B"  # Red for positive
      colors_vec[plot_values < 0] <- "#4ECDC4"  # Blue for negative
      colors_vec[na_positions] <- "#CCCCCC"    # Gray for NA
      
      # Create the bar plot
      bp <- barplot(plot_values,
                    names.arg = plot_names,
                    main = paste(mirna_id, "\nLog Fold Change across Comparisons"),
                    ylab = "logFC",
                    col = colors_vec,
                    border = "white",
                    las = 2,
                    cex.names = 0.8,
                    cex.main = 1.2,
                    ylim = c(min(plot_values, na.rm = TRUE) * 1.2, 
                             max(plot_values, na.rm = TRUE) * 1.2))
      
      # Add horizontal line at y = 0
      abline(h = 0, col = "black", lwd = 1, lty = 2)
      
      # Add value labels on bars
      for (i in 1:length(plot_values)) {
        if (na_positions[i]) {
          text(bp[i], plot_values[i], "NA", 
               pos = 3, cex = 0.9, font = 2, col = "gray30")
        } else {
          text(bp[i], plot_values[i], sprintf("%.2f", plot_values[i]), 
               pos = if(plot_values[i] >= 0) 3 else 1, 
               cex = 0.8, font = 2)
        }
        
        # Add significance stars
        if (!na_positions[i] && !is.na(mirna_sig[i]) && mirna_sig[i] == "*") {
          text(bp[i], plot_values[i], "*", 
               pos = if(plot_values[i] >= 0) 3 else 1, 
               offset = 0.5, cex = 1.5, font = 2, col = "gold")
        }
      }
      
      # Add legend
      legend("topright", 
             legend = c("Upregulated", "Downregulated", "No data", "* Significant"),
             fill = c("#FF6B6B", "#4ECDC4", "#CCCCCC", "white"),
             border = c("white", "white", "white", "white"),
             cex = 0.8)
      
      # Add statistics at bottom
      if (stats_text != "") {
        mtext(stats_text, side = 1, line = 8, cex = 0.9, col = "darkblue")
      }
      
      # Add grid for better readability
      grid(col = "lightgray", lty = "dotted", lwd = 0.5)
      
      # Redraw bars on top of grid
      barplot(plot_values,
              names.arg = plot_names,
              col = colors_vec,
              border = "white",
              las = 2,
              add = TRUE,
              axes = FALSE)
    })
    
    dev.off()
    
    cat(paste("    LogFC heatmap saved:", filename, "\n"))
  }
  
  # Generate overview logFC heatmap for all miRNAs
  cat("\nCreating overview logFC heatmap for all target miRNAs...\n")
  
  # FIXED: Better validation of input matrices
  if (is.null(logfc_matrix) || nrow(logfc_matrix) == 0 || ncol(logfc_matrix) == 0) {
    cat("Invalid or empty logFC matrix for overview heatmap\n")
    return()
  }
  
  # Filter out miRNAs with all NA values (but keep those with some data)
  has_some_data <- apply(!is.na(logfc_matrix), 1, any)
  valid_mirnas <- rownames(logfc_matrix)[has_some_data]
  
  if (length(valid_mirnas) == 0) {
    cat("No valid logFC data found for overview heatmap\n")
    
    # Create a message file instead
    message_file <- file.path(logfc_output_dir, "No_LogFC_Data_Available.txt")
    writeLines(c("No logFC data available for overview heatmap",
                 paste("Checked", length(target_mirnas_in_data), "miRNAs"),
                 paste("All values were NA"),
                 "This may indicate that the DE analysis files are missing or incomplete."),
               message_file)
    return()
  }
  
  cat(paste("Creating overview heatmap for", length(valid_mirnas), "miRNAs with data\n"))
  
  # FIXED: Subset matrices with proper dimension handling
  overview_logfc <- logfc_matrix[valid_mirnas, , drop = FALSE]
  overview_sig <- significance_matrix[valid_mirnas, , drop = FALSE]
  
  # FIXED: Ensure we maintain matrix structure
  if (!is.matrix(overview_logfc)) {
    overview_logfc <- as.matrix(overview_logfc)
  }
  if (!is.matrix(overview_sig)) {
    overview_sig <- as.matrix(overview_sig)
  }
  
  # For display, replace NA with 0 (but keep track of original NAs)
  display_logfc <- overview_logfc
  display_logfc[is.na(display_logfc)] <- 0
  
  # Create filename (directly in logfc_output_dir)
  filename <- file.path(logfc_output_dir, "All_miRNAs_LogFC_Overview_Heatmap.png")
  
  # Create logFC color scheme
  logfc_colors <- create_logfc_color_scheme(as.vector(overview_logfc))
  
  # Calculate dimensions
  n_mirnas <- nrow(display_logfc)
  n_comparisons <- ncol(display_logfc)
  
  # FIXED: Better dimension validation
  if (n_mirnas == 0 || n_comparisons == 0) {
    cat("Invalid dimensions for overview heatmap\n")
    return()
  }
  
  # Calculate dynamic dimensions
  base_width <- 10  # inches
  base_height <- 8  # inches
  mirna_height <- 0.3 # inches per miRNA
  comparison_width <- 0.5 # inches per comparison
  
  # Create heatmap
  png(filename, 
      width = max(base_width, n_comparisons * comparison_width), 
      height = max(base_height, n_mirnas * mirna_height),
      units = "in", 
      res = 300)
  
  tryCatch({
    # FIXED: Create cell annotation for significance and NA values
    cell_fun <- function(j, i, x, y, width, height, fill) {
      # FIXED: Ensure we're within bounds and handle matrix indexing properly
      if (i <= nrow(overview_logfc) && j <= ncol(overview_logfc)) {
        if (is.na(overview_logfc[i, j])) {
          # Show "NA" for missing data
          grid.text("NA", x, y, gp = gpar(fontsize = 8, fontface = "bold", col = "gray50"))
        } else if (!is.na(overview_sig[i, j]) && overview_sig[i, j] == "*") {
          # Show * for significant values
          grid.text("*", x, y, gp = gpar(fontsize = 10, fontface = "bold", col = "white"))
        }
      }
    }
    
    # Create heatmap
    ht <- Heatmap(
      display_logfc,
      name = "logFC",
      col = logfc_colors,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = gpar(fontsize = max(6, min(12, 120/n_mirnas))),
      column_names_gp = gpar(fontsize = 11),
      column_names_rot = 45,
      column_title = "miRNA Log Fold Changes - Overview (All Results)",
      column_title_gp = gpar(fontsize = 16, fontface = "bold"),
      row_title = "miRNAs",
      row_title_gp = gpar(fontsize = 14, fontface = "bold"),
      rect_gp = gpar(col = "white", lwd = 0.5),
      cell_fun = cell_fun,
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 10),
        grid_height = unit(4, "mm"),
        grid_width = unit(4, "mm"),
        legend_height = unit(8, "cm")
      ),
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      clustering_distance_rows = "euclidean",
      clustering_method_rows = "ward.D2",
      width = unit(max(8, n_comparisons * 1.5), "cm"),
      height = unit(max(6, n_mirnas * 0.4), "cm"),
      border = TRUE
    )
    
    # Draw heatmap
    draw(ht)
    
    # Add annotations
    grid.text(paste("Overview of", n_mirnas, "miRNAs across", n_comparisons, "comparisons"), 
              x = 0.5, y = 0.08, 
              gp = gpar(fontsize = 12, fontface = "bold"))
    
    grid.text(paste("* = p-value <", pvalue_threshold, "| NA = No data | Blue: Down | Red: Up"), 
              x = 0.5, y = 0.04, 
              gp = gpar(fontsize = 10, col = "gray40"))
    
    # Add summary statistics
    n_significant <- sum(overview_sig == "*", na.rm = TRUE)
    n_upregulated <- sum(overview_logfc > 0 & overview_sig == "*", na.rm = TRUE)
    n_downregulated <- sum(overview_logfc < 0 & overview_sig == "*", na.rm = TRUE)
    n_missing <- sum(is.na(overview_logfc))
    
    summary_text <- sprintf("Significant: %d total (%d up, %d down) | Missing: %d", 
                            n_significant, n_upregulated, n_downregulated, n_missing)
    grid.text(summary_text, x = 0.5, y = 0.12, 
              gp = gpar(fontsize = 11, col = "darkblue", fontface = "bold"))
    
  }, error = function(e) {
    cat(paste("Error creating overview logFC heatmap:", e$message, "\n"))
    
    # FIXED: Fallback plot that handles the data properly
    par(mar = c(12, 10, 6, 3))
    
    # FIXED: Create a simple heatmap using base R with proper dimension handling
    if (nrow(display_logfc) > 0 && ncol(display_logfc) > 0) {
      # Create color palette
      col_palette <- colorRampPalette(c("#0D47A1", "white", "#D32F2F"))(50)
      
      # FIXED: Proper matrix transposition and dimension checking
      plot_matrix <- t(as.matrix(display_logfc))
      
      # Ensure we have the right dimensions
      if (nrow(plot_matrix) == ncol(display_logfc) && ncol(plot_matrix) == nrow(display_logfc)) {
        # Create the heatmap
        image(1:ncol(plot_matrix), 1:nrow(plot_matrix), 
              plot_matrix, 
              col = col_palette,
              main = "miRNA LogFC Overview (All Results)",
              xlab = "", ylab = "",
              axes = FALSE)
        
        # Add axes
        axis(1, at = 1:ncol(plot_matrix), 
             labels = colnames(plot_matrix), las = 2, cex.axis = 0.8)
        axis(2, at = 1:nrow(plot_matrix), 
             labels = rownames(plot_matrix), las = 2, cex.axis = 0.7)
        
        # Add grid lines
        abline(h = 1:nrow(plot_matrix) + 0.5, col = "white", lwd = 0.5)
        abline(v = 1:ncol(plot_matrix) + 0.5, col = "white", lwd = 0.5)
      } else {
        plot(1, 1, type = "n", main = "Matrix Dimension Error")
        text(1, 1, paste("Dimension mismatch:", 
                         nrow(plot_matrix), "x", ncol(plot_matrix)), 
             cex = 1.2, col = "red")
      }
    } else {
      plot(1, 1, type = "n", main = "No LogFC Data Available for Overview")
      text(1, 1, "No valid logFC data found", cex = 1.5, col = "gray50")
    }
  })
  
  dev.off()
  
  cat(paste("Overview logFC heatmap saved:", filename, "\n"))
  
  # Create summary table
  if (nrow(overview_logfc) > 0) {
    logfc_summary <- data.frame(
      miRNA_ID = rownames(overview_logfc),
      overview_logfc,
      stringsAsFactors = FALSE
    )
    
    summary_file <- file.path(logfc_output_dir, "LogFC_Summary_Table.csv")
    write.csv(logfc_summary, summary_file, row.names = FALSE)
    cat(paste("LogFC summary table saved:", summary_file, "\n"))
  }
  
  # Save enhanced summary
  summary_file <- file.path(output_dir, "Enhanced_miRNA_Heatmaps_Summary.csv")
  write.csv(summary_results, summary_file, row.names = FALSE)
  cat(paste("\nEnhanced summary saved to:", summary_file, "\n"))
  
  # Print summary statistics
  cat("\n=== ENHANCED SUMMARY STATISTICS ===\n")
  cat(paste("Total enhanced heatmaps generated:", nrow(summary_results), "\n"))
  
  if (nrow(summary_results) > 0) {
    cat("Expression category distribution:\n")
    print(table(summary_results$Expression_Category))
    
    cat("\nSignificance status distribution:\n")
    print(table(summary_results$Significance_Status))
    
    cat("\nComparison distribution:\n")
    print(table(summary_results$Comparison))
    
    # Statistics for significant results
    sig_results <- summary_results[grepl("\\*", summary_results$Significance_Status), ]
    if (nrow(sig_results) > 0) {
      cat(paste("\nSignificant results (marked with *):", nrow(sig_results), "\n"))
      cat("LogFC range for significant results:", 
          paste(round(range(sig_results$logFC, na.rm = TRUE), 3), collapse = " to "), "\n")
    }
  }
  
  return(summary_results)
}

# =============================================================================
# MAIN EXECUTION WITH FIXES
# =============================================================================

cat("Starting FIXED individual miRNA heatmap generation with LogFC analysis...\n")
cat("===========================================================================\n")

# Generate the enhanced heatmaps (both expression and logFC-based)
summary_results <- generate_enhanced_mirna_heatmaps(
  cpm_data = cpm_data,
  target_mirna_ids = target_mirna_ids,
  de_results = de_results,
  output_dir = output_dir,
  logfc_output_dir = logfc_output_dir,
  pvalue_threshold = 0.05,
  logfc_threshold = 0,  # Changed to 0 to show all results
  use_fdr = FALSE  # Using p-value instead of FDR for significance
)

cat("\n===========================================================================\n")
cat("FIXED individual miRNA heatmap generation completed!\n")
cat(paste("Expression-based heatmaps: Check the", output_dir, "directory\n"))
cat(paste("LogFC-based heatmaps: Check the", logfc_output_dir, "directory\n"))
cat("\nGenerated outputs:\n")
cat("1. Individual expression heatmaps for each miRNA and comparison\n")
cat("2. Individual logFC heatmaps for each miRNA across all comparisons\n")
cat("3. Overview logFC heatmap showing all miRNAs together\n")
cat("4. Summary tables with statistics and file information\n")
cat("5. All results shown (not filtered) with * marking significant ones (p-value < 0.05)\n")
