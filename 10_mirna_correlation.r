# miRNA-Target Correlation Analysis Script (Improved with Debugging)
# This script performs Pearson's correlation analysis for miRNA-target pairs
# with specific pairwise comparisons and generates tables and plots

setwd("C:/Users/Leandro/OneDrive - usp.br/LIGNINLAB/Projeto FAPESP CNPq - Setaria miRNAs/R analysis/mirna_correlation")
cat("Current working directory:", getwd(), "\n")


# Load required libraries
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(dplyr)) install.packages("dplyr")
if (!require(gridExtra)) install.packages("gridExtra")
if (!require(grid)) install.packages("grid")
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

# Create output directories with full paths
current_dir <- getwd()
tables_dir <- file.path(current_dir, "tables")
plots_dir <- file.path(current_dir, "plots")

dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

cat("Working directory:", current_dir, "\n")
cat("Tables will be saved to:", tables_dir, "\n")
cat("Plots will be saved to:", plots_dir, "\n\n")

# Read input files with error handling
cat("Reading input files...\n")

# Read miRNA-target pairs
if (!file.exists("renamed_mirna_target_modules_table.txt")) {
  stop("Error: renamed_mirna_target_modules_table.txt not found!")
}
mirna_targets <- read.table("renamed_mirna_target_modules_table.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Read", nrow(mirna_targets), "miRNA-target pairs\n")
print(head(mirna_targets))

# Read transcript expression data
if (!file.exists("transcripts_TPM.txt")) {
  stop("Error: transcripts_TPM.txt not found!")
}
transcripts <- read.table("transcripts_TPM.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
cat("Read transcript data for", nrow(transcripts), "genes with", ncol(transcripts), "samples\n")
cat("Sample names:", colnames(transcripts)[1:min(5, ncol(transcripts))], "...\n")

# Read miRNA expression data
if (!file.exists("renamed_mirnas_CPM.txt")) {
  stop("Error: renamed_mirnas_CPM.txt not found!")
}
mirnas <- read.table("renamed_mirnas_CPM.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
cat("Read miRNA data for", nrow(mirnas), "miRNAs with", ncol(mirnas), "samples\n")
cat("Sample names:", colnames(mirnas)[1:min(5, ncol(mirnas))], "...\n\n")

# Check for column name mismatches
cat("Checking sample name consistency...\n")
transcript_samples <- colnames(transcripts)
mirna_samples <- colnames(mirnas)

if (!identical(transcript_samples, mirna_samples)) {
  cat("Warning: Sample names differ between transcript and miRNA files\n")
  cat("Transcript samples:", paste(transcript_samples[1:min(5, length(transcript_samples))], collapse = ", "), "...\n")
  cat("miRNA samples:", paste(mirna_samples[1:min(5, length(mirna_samples))], collapse = ", "), "...\n")
  
  # Use common samples
  common_samples <- intersect(transcript_samples, mirna_samples)
  cat("Using", length(common_samples), "common samples\n")
  transcripts <- transcripts[, common_samples]
  mirnas <- mirnas[, common_samples]
}

# Define sample groupings for pairwise comparisons
sample_groups <- list(
  "Third_vs_Fifth_Internode" = list(
    group1 = c("Third_Internode1", "Third_Internode2", "Third_Internode3"),
    group2 = c("Fifth_Internode1", "Fifth_Internode2", "Fifth_Internode3"),
    name1 = "Third_Internode",
    name2 = "Fifth_Internode"
  ),
  "Mature_vs_Young_Leaves" = list(
    group1 = c("Mature_Leaves1", "Mature_Leaves2", "Mature_Leaves3"),
    group2 = c("Young_Leaves1", "Young_Leaves2", "Young_Leaves3"),
    name1 = "Mature_Leaves",
    name2 = "Young_Leaves"
  ),
  "Mature_vs_Young_Roots" = list(
    group1 = c("Mature_Roots1", "Mature_Roots2", "Mature_Roots3"),
    group2 = c("Young_Roots1", "Young_Roots2", "Young_Roots3"),
    name1 = "Mature_Roots",
    name2 = "Young_Roots"
  )
)

# Check which samples are available for each comparison
cat("Checking sample availability for comparisons...\n")
available_samples <- colnames(transcripts)
for (comp_name in names(sample_groups)) {
  comp <- sample_groups[[comp_name]]
  available_group1 <- intersect(comp$group1, available_samples)
  available_group2 <- intersect(comp$group2, available_samples)
  cat(comp_name, ":\n")
  cat("  Group1 (", comp$name1, "):", length(available_group1), "of", length(comp$group1), "samples available\n")
  cat("  Group2 (", comp$name2, "):", length(available_group2), "of", length(comp$group2), "samples available\n")
}
cat("\n")

# Function to perform correlation analysis for a miRNA-target pair
perform_correlation_analysis <- function(mirna_id, target_id) {
  
  # Extract expression data
  mirna_expr <- as.numeric(mirnas[mirna_id, ])
  target_expr <- as.numeric(transcripts[target_id, ])
  
  # Check if both genes have expression data
  if (all(is.na(mirna_expr)) || all(is.na(target_expr))) {
    warning(paste("Missing expression data for", mirna_id, "or", target_id))
    return(NULL)
  }
  
  # Create results list for all comparisons
  all_results <- list()
  
  # Perform correlation for each pairwise comparison
  for (comparison_name in names(sample_groups)) {
    comparison <- sample_groups[[comparison_name]]
    
    # Get indices for the samples in each group
    group1_indices <- which(colnames(mirnas) %in% comparison$group1)
    group2_indices <- which(colnames(mirnas) %in% comparison$group2)
    
    # Skip if not enough samples in either group
    if (length(group1_indices) == 0 || length(group2_indices) == 0) {
      cat("    Skipping", comparison_name, "- missing sample groups\n")
      next
    }
    
    # Combine indices for this comparison
    selected_indices <- c(group1_indices, group2_indices)
    
    # Extract expression values for selected samples
    mirna_selected <- mirna_expr[selected_indices]
    target_selected <- target_expr[selected_indices]
    
    # Remove any NA values
    valid_indices <- !is.na(mirna_selected) & !is.na(target_selected)
    mirna_clean <- mirna_selected[valid_indices]
    target_clean <- target_selected[valid_indices]
    
    # Skip if not enough valid data points
    if (length(mirna_clean) < 3) {
      cat("    Skipping", comparison_name, "- only", length(mirna_clean), "valid data points\n")
      next
    }
    
    # Check for zero variance (which causes correlation to fail)
    if (var(mirna_clean) == 0 || var(target_clean) == 0) {
      cat("    Skipping", comparison_name, "- zero variance in expression data\n")
      next
    }
    
    # Perform Pearson correlation with error handling
    cor_test <- tryCatch({
      cor.test(mirna_clean, target_clean, method = "pearson")
    }, error = function(e) {
      cat("    Error in correlation test for", comparison_name, ":", e$message, "\n")
      return(NULL)
    }, warning = function(w) {
      cat("    Warning in correlation test for", comparison_name, ":", w$message, "\n")
      # Return the result even with warnings (like zero standard deviation)
      cor.test(mirna_clean, target_clean, method = "pearson")
    })
    
    # Skip if correlation test failed
    if (is.null(cor_test)) {
      next
    }
    
    # Create sample labels
    sample_labels <- c(rep(comparison$name1, length(group1_indices)), 
                       rep(comparison$name2, length(group2_indices)))
    sample_labels <- sample_labels[valid_indices]
    
    # Store results
    comparison_result <- list(
      comparison = comparison_name,
      mirna_expr = mirna_clean,
      target_expr = target_clean,
      sample_labels = sample_labels,
      correlation = cor_test$estimate,
      p_value = cor_test$p.value,
      conf_int = cor_test$conf.int
    )
    
    all_results[[comparison_name]] <- comparison_result
    cat("    ", comparison_name, ": r =", round(cor_test$estimate, 4), ", p =", format(cor_test$p.value, digits = 3), "\n")
  }
  
  return(all_results)
}

# Function to create correlation table
create_correlation_table <- function(results, mirna_id, target_id) {
  
  table_data <- data.frame(
    Comparison = character(),
    miRNA_ID = character(),
    Target_ID = character(),
    Correlation_r = numeric(),
    P_value = numeric(),
    CI_lower = numeric(),
    CI_upper = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (comparison_name in names(results)) {
    result <- results[[comparison_name]]
    
    table_data <- rbind(table_data, data.frame(
      Comparison = comparison_name,
      miRNA_ID = mirna_id,
      Target_ID = target_id,
      Correlation_r = round(result$correlation, 4),
      P_value = result$p_value,
      CI_lower = round(result$conf_int[1], 4),
      CI_upper = round(result$conf_int[2], 4),
      stringsAsFactors = FALSE
    ))
  }
  
  return(table_data)
}

# Function to create individual scatter plots
create_scatter_plots <- function(results, mirna_id, target_id) {
  
  plots <- list()
  
  for (comparison_name in names(results)) {
    result <- results[[comparison_name]]
    
    # Create data frame for plotting
    plot_data <- data.frame(
      miRNA_Expression = result$mirna_expr,
      Target_Expression = result$target_expr,
      Sample_Group = result$sample_labels
    )
    
    # Create the plot
    p <- ggplot(plot_data, aes(x = miRNA_Expression, y = Target_Expression, color = Sample_Group)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
      labs(
        title = paste("Correlation:", mirna_id, "vs", target_id),
        subtitle = paste(comparison_name, 
                         "\nr =", round(result$correlation, 4), 
                         ", p =", format(result$p_value, scientific = TRUE, digits = 3)),
        x = paste("miRNA Expression (", mirna_id, ")", sep = ""),
        y = paste("Target Expression (", target_id, ")", sep = ""),
        color = "Sample Group"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "bottom"
      ) +
      scale_color_brewer(type = "qual", palette = "Set1")
    
    plots[[comparison_name]] <- p
  }
  
  return(plots)
}

# Function to create combined scatter plot with all comparisons
create_combined_scatter_plot <- function(results, mirna_id, target_id) {
  
  # Check if we have results for all comparisons
  if (length(results) == 0) {
    return(NULL)
  }
  
  # Create individual plots for combining
  individual_plots <- list()
  
  for (comparison_name in names(results)) {
    result <- results[[comparison_name]]
    
    # Create data frame for plotting
    plot_data <- data.frame(
      miRNA_Expression = result$mirna_expr,
      Target_Expression = result$target_expr,
      Sample_Group = result$sample_labels
    )
    
    # Create clean comparison name for title
    clean_comp_name <- gsub("_", " ", comparison_name)
    clean_comp_name <- gsub("vs", "vs", clean_comp_name)
    
    # Create individual plot
    p <- ggplot(plot_data, aes(x = miRNA_Expression, y = Target_Expression, color = Sample_Group)) +
      geom_point(size = 2.5, alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", size = 0.8) +
      labs(
        title = clean_comp_name,
        subtitle = paste("r =", round(result$correlation, 4), 
                         ", p =", format(result$p_value, scientific = TRUE, digits = 3)),
        x = "miRNA Expression",
        y = "Target Expression",
        color = "Sample Group"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.position = "bottom",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.margin = margin(t = 5, b = 5),
        panel.border = element_rect(colour = "grey80", fill = NA, size = 0.5)
      ) +
      scale_color_brewer(type = "qual", palette = "Set1")
    
    individual_plots[[comparison_name]] <- p
  }
  
  # Load required library for combining plots
  if (!require(gridExtra)) {
    install.packages("gridExtra")
    library(gridExtra)
  }
  if (!require(grid)) {
    install.packages("grid")
    library(grid)
  }
  
  # Create main title
  main_title <- paste("Correlation Analysis:", mirna_id, "vs", target_id)
  
  # Combine plots based on number of available comparisons
  if (length(individual_plots) == 1) {
    combined_plot <- individual_plots[[1]] +
      labs(title = main_title, subtitle = individual_plots[[1]]$labels$subtitle)
  } else if (length(individual_plots) == 2) {
    combined_plot <- grid.arrange(
      grobs = individual_plots,
      ncol = 2,
      top = textGrob(main_title, gp = gpar(fontsize = 14, fontface = "bold"))
    )
  } else {
    # 3 plots arranged in a row
    combined_plot <- grid.arrange(
      grobs = individual_plots,
      ncol = 3,
      top = textGrob(main_title, gp = gpar(fontsize = 14, fontface = "bold"))
    )
  }
  
  return(combined_plot)
}

# Main analysis loop
cat("Starting miRNA-target correlation analysis...\n")
cat("Number of miRNA-target pairs to analyze:", nrow(mirna_targets), "\n\n")

# Check if any target genes exist in the data
cat("Checking gene ID availability...\n")
mirna_ids_available <- sum(mirna_targets$miRNA.ID %in% rownames(mirnas))
target_ids_available <- sum(mirna_targets$Target.ID %in% rownames(transcripts))
cat("miRNA IDs found in expression data:", mirna_ids_available, "of", nrow(mirna_targets), "\n")
cat("Target IDs found in expression data:", target_ids_available, "of", nrow(mirna_targets), "\n\n")

# Show some examples of IDs
cat("Example miRNA IDs in pairs:", paste(head(mirna_targets$miRNA.ID, 3), collapse = ", "), "\n")
cat("Example miRNA IDs in data:", paste(head(rownames(mirnas), 3), collapse = ", "), "\n")
cat("Example target IDs in pairs:", paste(head(mirna_targets$Target.ID, 3), collapse = ", "), "\n")
cat("Example target IDs in data:", paste(head(rownames(transcripts), 3), collapse = ", "), "\n\n")

processed_pairs <- 0
successful_pairs <- 0

# Initialize comprehensive results table
all_correlations <- data.frame(
  miRNA_ID = character(),
  Target_ID = character(),
  Comparison = character(),
  Correlation_r = numeric(),
  P_value = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  Significance = character(),
  N_samples = integer(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(mirna_targets)) {
  mirna_id <- mirna_targets$miRNA.ID[i]
  target_id <- mirna_targets$Target.ID[i]
  
  cat("Processing pair", i, "of", nrow(mirna_targets), ":", mirna_id, "vs", target_id, "\n")
  
  # Check if both genes exist in the expression data
  if (!(mirna_id %in% rownames(mirnas))) {
    cat("  Warning: miRNA", mirna_id, "not found in expression data\n")
    next
  }
  
  if (!(target_id %in% rownames(transcripts))) {
    cat("  Warning: Target", target_id, "not found in expression data\n")
    next
  }
  
  processed_pairs <- processed_pairs + 1
  
  # Perform correlation analysis
  results <- perform_correlation_analysis(mirna_id, target_id)
  
  if (is.null(results) || length(results) == 0) {
    cat("  Skipping due to insufficient data\n")
    next
  }
  
  successful_pairs <- successful_pairs + 1
  
  # Create file names (clean special characters)
  clean_mirna <- gsub("[^A-Za-z0-9_.]", "_", mirna_id)
  clean_target <- gsub("[^A-Za-z0-9_.]", "_", target_id)
  base_filename <- paste(clean_mirna, clean_target, sep = "_")
  
  # Create and save correlation table
  correlation_table <- create_correlation_table(results, mirna_id, target_id)
  table_filename <- file.path(tables_dir, paste(base_filename, "_correlation.csv", sep = ""))
  
  tryCatch({
    write.csv(correlation_table, table_filename, row.names = FALSE)
    cat("  ??? Saved table:", table_filename, "\n")
  }, error = function(e) {
    cat("  ??? Error saving table:", e$message, "\n")
  })
  
  # Add results to comprehensive table
  for (comparison_name in names(results)) {
    result <- results[[comparison_name]]
    
    # Handle NA values in correlation and p-value
    cor_value <- result$correlation
    p_val <- result$p_value
    
    # Check for NA values and handle them
    if (is.na(cor_value)) {
      cor_value <- NA
    }
    
    if (is.na(p_val)) {
      p_val <- NA
      significance <- "NA"
    } else {
      # Determine significance level
      if (p_val < 0.001) {
        significance <- "***"
      } else if (p_val < 0.01) {
        significance <- "**"
      } else if (p_val < 0.05) {
        significance <- "*"
      } else {
        significance <- "ns"
      }
    }
    
    # Handle confidence intervals
    ci_lower <- ifelse(is.na(result$conf_int[1]), NA, result$conf_int[1])
    ci_upper <- ifelse(is.na(result$conf_int[2]), NA, result$conf_int[2])
    
    # Add row to comprehensive table
    new_row <- data.frame(
      miRNA_ID = mirna_id,
      Target_ID = target_id,
      Comparison = comparison_name,
      Correlation_r = round(cor_value, 4),
      P_value = p_val,
      CI_lower = round(ci_lower, 4),
      CI_upper = round(ci_upper, 4),
      Significance = significance,
      N_samples = length(result$mirna_expr),
      stringsAsFactors = FALSE
    )
    
    all_correlations <- rbind(all_correlations, new_row)
  }
  
  # Create and save individual scatter plots
  plots <- create_scatter_plots(results, mirna_id, target_id)
  
  for (comparison_name in names(plots)) {
    clean_comparison <- gsub("[^A-Za-z0-9_]", "_", comparison_name)
    plot_filename <- file.path(plots_dir, paste(base_filename, "_", clean_comparison, ".png", sep = ""))
    
    tryCatch({
      ggsave(plot_filename, plots[[comparison_name]], 
             width = 10, height = 8, dpi = 300, units = "in")
      cat("  ??? Saved individual plot:", plot_filename, "\n")
    }, error = function(e) {
      cat("  ??? Error saving individual plot:", e$message, "\n")
    })
  }
  
  # Create and save combined scatter plot
  combined_plot <- create_combined_scatter_plot(results, mirna_id, target_id)
  
  if (!is.null(combined_plot)) {
    combined_plot_filename <- file.path(plots_dir, paste(base_filename, "_combined.png", sep = ""))
    
    tryCatch({
      # Adjust width based on number of comparisons
      plot_width <- ifelse(length(results) == 3, 15, ifelse(length(results) == 2, 12, 8))
      plot_height <- 6
      
      ggsave(combined_plot_filename, combined_plot, 
             width = plot_width, height = plot_height, dpi = 300, units = "in")
      cat("  ??? Saved combined plot:", combined_plot_filename, "\n")
    }, error = function(e) {
      cat("  ??? Error saving combined plot:", e$message, "\n")
    })
  }
  
  cat("  Completed analysis for", mirna_id, "vs", target_id, "\n\n")
}

cat("Analysis completed!\n")
cat("Summary:\n")
cat("- Total miRNA-target pairs:", nrow(mirna_targets), "\n")
cat("- Pairs with matching IDs:", processed_pairs, "\n")
cat("- Pairs successfully analyzed:", successful_pairs, "\n")
cat("- Results saved in:\n")
cat("  - Tables:", tables_dir, "\n")
cat("  - Plots:", plots_dir, "\n")

# Create and save comprehensive correlation table
if (nrow(all_correlations) > 0) {
  comprehensive_table_filename <- file.path(tables_dir, "all_correlations_comprehensive.csv")
  tryCatch({
    write.csv(all_correlations, comprehensive_table_filename, row.names = FALSE)
    cat("??? Saved comprehensive correlation table:", comprehensive_table_filename, "\n")
  }, error = function(e) {
    cat("??? Error saving comprehensive table:", e$message, "\n")
  })
  
  # Create summary statistics table
  summary_stats <- all_correlations %>%
    group_by(Comparison) %>%
    summarise(
      N_pairs = n(),
      Mean_r = round(mean(Correlation_r, na.rm = TRUE), 4),
      Median_r = round(median(Correlation_r, na.rm = TRUE), 4),
      SD_r = round(sd(Correlation_r, na.rm = TRUE), 4),
      Min_r = round(min(Correlation_r, na.rm = TRUE), 4),
      Max_r = round(max(Correlation_r, na.rm = TRUE), 4),
      Significant_p0.05 = sum(P_value < 0.05, na.rm = TRUE),
      Significant_p0.01 = sum(P_value < 0.01, na.rm = TRUE),
      Significant_p0.001 = sum(P_value < 0.001, na.rm = TRUE),
      Positive_correlations = sum(Correlation_r > 0, na.rm = TRUE),
      Negative_correlations = sum(Correlation_r < 0, na.rm = TRUE),
      Strong_correlations_abs_gt_0.7 = sum(abs(Correlation_r) > 0.7, na.rm = TRUE),
      Moderate_correlations_abs_0.3_to_0.7 = sum(abs(Correlation_r) >= 0.3 & abs(Correlation_r) <= 0.7, na.rm = TRUE),
      Weak_correlations_abs_lt_0.3 = sum(abs(Correlation_r) < 0.3, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Add percentage columns
  summary_stats <- summary_stats %>%
    mutate(
      Percent_significant_p0.05 = round((Significant_p0.05 / N_pairs) * 100, 2),
      Percent_significant_p0.01 = round((Significant_p0.01 / N_pairs) * 100, 2),
      Percent_significant_p0.001 = round((Significant_p0.001 / N_pairs) * 100, 2),
      Percent_positive = round((Positive_correlations / N_pairs) * 100, 2),
      Percent_negative = round((Negative_correlations / N_pairs) * 100, 2),
      Percent_strong = round((Strong_correlations_abs_gt_0.7 / N_pairs) * 100, 2),
      Percent_moderate = round((Moderate_correlations_abs_0.3_to_0.7 / N_pairs) * 100, 2),
      Percent_weak = round((Weak_correlations_abs_lt_0.3 / N_pairs) * 100, 2)
    )
  
  # Save summary table
  summary_table_filename <- file.path(tables_dir, "correlation_summary_statistics.csv")
  tryCatch({
    write.csv(summary_stats, summary_table_filename, row.names = FALSE)
    cat("??? Saved summary statistics table:", summary_table_filename, "\n")
  }, error = function(e) {
    cat("??? Error saving summary table:", e$message, "\n")
  })
  
  # Display summary in console
  cat("\n=== CORRELATION ANALYSIS SUMMARY ===\n")
  cat("Total correlations analyzed:", nrow(all_correlations), "\n\n")
  
  for (i in 1:nrow(summary_stats)) {
    comp <- summary_stats$Comparison[i]
    cat("Comparison:", comp, "\n")
    cat("  - Number of pairs:", summary_stats$N_pairs[i], "\n")
    cat("  - Mean correlation:", summary_stats$Mean_r[i], "\n")
    cat("  - Significant correlations (p < 0.05):", summary_stats$Significant_p0.05[i], 
        "(", summary_stats$Percent_significant_p0.05[i], "%)\n")
    cat("  - Strong correlations (|r| > 0.7):", summary_stats$Strong_correlations_abs_gt_0.7[i], 
        "(", summary_stats$Percent_strong[i], "%)\n")
    cat("  - Positive correlations:", summary_stats$Positive_correlations[i], 
        "(", summary_stats$Percent_positive[i], "%)\n")
    cat("  - Negative correlations:", summary_stats$Negative_correlations[i], 
        "(", summary_stats$Percent_negative[i], "%)\n\n")
  }
  
  # Create overall summary across all comparisons
  overall_summary <- data.frame(
    Total_correlations = nrow(all_correlations),
    Mean_r_all = round(mean(all_correlations$Correlation_r, na.rm = TRUE), 4),
    Median_r_all = round(median(all_correlations$Correlation_r, na.rm = TRUE), 4),
    SD_r_all = round(sd(all_correlations$Correlation_r, na.rm = TRUE), 4),
    Total_significant_p0.05 = sum(all_correlations$P_value < 0.05, na.rm = TRUE),
    Total_significant_p0.01 = sum(all_correlations$P_value < 0.01, na.rm = TRUE),
    Total_significant_p0.001 = sum(all_correlations$P_value < 0.001, na.rm = TRUE),
    Total_positive = sum(all_correlations$Correlation_r > 0, na.rm = TRUE),
    Total_negative = sum(all_correlations$Correlation_r < 0, na.rm = TRUE),
    Total_strong = sum(abs(all_correlations$Correlation_r) > 0.7, na.rm = TRUE),
    Total_moderate = sum(abs(all_correlations$Correlation_r) >= 0.3 & abs(all_correlations$Correlation_r) <= 0.7, na.rm = TRUE),
    Total_weak = sum(abs(all_correlations$Correlation_r) < 0.3, na.rm = TRUE)
  )
  
  # Add percentages to overall summary
  overall_summary <- overall_summary %>%
    mutate(
      Percent_significant_p0.05 = round((Total_significant_p0.05 / Total_correlations) * 100, 2),
      Percent_significant_p0.01 = round((Total_significant_p0.01 / Total_correlations) * 100, 2),
      Percent_significant_p0.001 = round((Total_significant_p0.001 / Total_correlations) * 100, 2),
      Percent_positive = round((Total_positive / Total_correlations) * 100, 2),
      Percent_negative = round((Total_negative / Total_correlations) * 100, 2),
      Percent_strong = round((Total_strong / Total_correlations) * 100, 2),
      Percent_moderate = round((Total_moderate / Total_correlations) * 100, 2),
      Percent_weak = round((Total_weak / Total_correlations) * 100, 2)
    )
  
  # Save overall summary
  overall_summary_filename <- file.path(tables_dir, "correlation_overall_summary.csv")
  tryCatch({
    write.csv(overall_summary, overall_summary_filename, row.names = FALSE)
    cat("??? Saved overall summary table:", overall_summary_filename, "\n")
  }, error = function(e) {
    cat("??? Error saving overall summary:", e$message, "\n")
  })
  
  cat("=== OVERALL SUMMARY (ALL COMPARISONS) ===\n")
  cat("Total correlations:", overall_summary$Total_correlations, "\n")
  cat("Overall mean correlation:", overall_summary$Mean_r_all, "\n")
  cat("Significant correlations (p < 0.05):", overall_summary$Total_significant_p0.05, 
      "(", overall_summary$Percent_significant_p0.05, "%)\n")
  cat("Strong correlations (|r| > 0.7):", overall_summary$Total_strong, 
      "(", overall_summary$Percent_strong, "%)\n")
  cat("Positive vs Negative:", overall_summary$Total_positive, "vs", overall_summary$Total_negative, "\n")
  
} else {
  cat("No correlations were successfully calculated.\n")
}

# List created files
if (dir.exists(tables_dir)) {
  table_files <- list.files(tables_dir, pattern = "\\.csv$")
  cat("- Created", length(table_files), "table files\n")
}

if (dir.exists(plots_dir)) {
  plot_files <- list.files(plots_dir, pattern = "\\.png$")
  cat("- Created", length(plot_files), "plot files\n")
}