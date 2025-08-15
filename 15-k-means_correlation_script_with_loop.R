library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(readr)
library(stringr)

# Set working directory to folder containing your CSV files
setwd("C:/Users/Leandro/OneDrive - usp.br/LIGNINLAB/Projeto FAPESP CNPq - Setaria miRNAs/R analysis/k-means_correlation/table_expression_2_samples")

# Get list of all CSV files in directory
csv_files <- list.files(pattern = "\\.csv$")

# Define the processing function (updated with your working version)
process_mirna_target_data <- function(file_path) {
  # Read the data
  data <- read_csv(file_path)
  
  # Extract sample information
  sample_cols <- colnames(data)[4:21]  # Adjust if your column structure differs
  
  # Process expression data
  processed_data <- data %>%
    select(Type, ID, Description, all_of(sample_cols)) %>%
    pivot_longer(cols = all_of(sample_cols), names_to = "Sample", values_to = "Expression") %>%
    mutate(
      Tissue = case_when(
        grepl("Fifth_Internode", Sample) ~ "Fifth_Internode",
        grepl("Third_Internode", Sample) ~ "Third_Internode", 
        grepl("Mature_Leaves", Sample) ~ "Mature_Leaves",
        grepl("Young_Leaves", Sample) ~ "Young_Leaves",
        grepl("Mature_Roots", Sample) ~ "Mature_Roots",
        grepl("Young_Roots", Sample) ~ "Young_Roots"
      ),
      Replicate = as.numeric(str_extract(Sample, "\\d+$")),
      Tissue_ordered = factor(Tissue, levels = c(
        "Young_Leaves", "Mature_Leaves", 
        "Young_Roots", "Mature_Roots",
        "Third_Internode", "Fifth_Internode"
      )),
      Log_Expression = log2(Expression + 1),
      Gene_ID = ID
    ) %>%
    group_by(Gene_ID) %>%
    mutate(Z_score = scale(Log_Expression)[,1]) %>%
    ungroup()
  
  # Create summary statistics
  summary_data <- processed_data %>%
    group_by(Type, Gene_ID, Description, Tissue_ordered) %>%
    summarise(
      Mean_Z_score = mean(Z_score, na.rm = TRUE),
      SE_Z_score = sd(Z_score, na.rm = TRUE) / sqrt(n()),
      Lower_CI = Mean_Z_score - 1.96 * SE_Z_score,
      Upper_CI = Mean_Z_score + 1.96 * SE_Z_score,
      N_samples = n(),
      .groups = "drop"
    ) %>%
    mutate(Tissue_position = as.numeric(Tissue_ordered))
  
  return(list(processed = processed_data, summary = summary_data))
}

create_mirna_target_plot <- function(mirna_data, target_data, processed_data) {
  # Calculate dynamic layout parameters
  n_targets <- length(unique(target_data$Gene_ID))
  targets_per_row <- min(6, max(4, ceiling(sqrt(n_targets)))) # Between 4-6 columns
  n_rows <- ceiling(n_targets / targets_per_row)
  
  # Calculate dynamic dimensions (in inches)
  base_width <- 16
  base_height <- 9
  height_per_row <- 3.5
  plot_width <- base_width
  plot_height <- base_height + (n_rows * height_per_row)
  
  tissue_labels <- c("YL", "ML", "YR", "MR", "TI", "FI")
  
  paired_plots <- lapply(unique(target_data$Gene_ID), function(target_id) {
    target_subset <- target_data %>% filter(Gene_ID == target_id)
    
    # Correlation calculation - MUST come first
    cor_result <- tryCatch({
      cor_data <- processed_data %>%
        filter(Type == "miRNA" | (Type == "Target" & Gene_ID == target_id)) %>%
        select(Sample, Tissue, Type, Z_score) %>%
        pivot_wider(names_from = Type, values_from = Z_score) %>%
        filter(complete.cases(.))
      
      if (nrow(cor_data) < 3) {
        list(r = NA, p_value = "NA (n<3)")
      } else {
        ct <- cor.test(cor_data$miRNA, cor_data$Target)
        list(
          r = round(ct$estimate, 2),
          p_value = ifelse(ct$p.value < 0.001, "< 0.001", round(ct$p.value, 3))
        )
      }
    }, error = function(e) {
      list(r = NA, p_value = "Error")
    })
    
    # Define cor_label HERE after cor_result is created
    cor_label <- if(is.na(cor_result$r)) {
      "r = NA (no data)"
    } else {
      paste0("r = ", cor_result$r, ", p = ", cor_result$p_value)
    }
    
    # miRNA plot
    p_mirna <- ggplot(mirna_data, aes(x = Tissue_position)) +
      geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), 
                  alpha = 0.2, fill = "#E69F00") +
      geom_line(aes(y = Mean_Z_score), color = "#E69F00", linewidth = 1) +
      geom_point(aes(y = Mean_Z_score), color = "#E69F00", size = 1) +
      theme_minimal() +
      scale_x_continuous(breaks = 1:6, labels = tissue_labels) +
      labs(title = "miRNA", y = "Z-score") +
      ylim(-3, 3) +
      theme(
        axis.title.x = element_blank(),
        plot.title = element_text(color = "#E69F00", face = "bold", size = 12)
      )
    
    # Target plot - NOW cor_label is defined
    p_target <- ggplot(target_subset, aes(x = Tissue_position)) +
      geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), 
                  alpha = 0.2, fill = "#56B4E9") +
      geom_line(aes(y = Mean_Z_score), color = "#56B4E9", linewidth = 1) +
      geom_point(aes(y = Mean_Z_score), color = "#56B4E9", size = 1) +
      scale_x_continuous(breaks = 1:6, labels = tissue_labels) +
      labs(
        title = paste0(target_id, "\n", cor_label),  # Now this works
        x = "Tissue", 
        y = "Z-score",
        tag = "Target"
      ) +
      ylim(-3, 3) +
      theme_minimal() +
      theme(
        plot.title = element_text(color = "#56B4E9", face = "bold", size = 15, hjust = 0.5, margin = margin(b = 5)),
        plot.tag = element_text(color = "#56B4E9", face = "bold", size = 12, 
                                margin = margin(t = 5, l = 10)),
        plot.tag.position = c(0.12, 0.98)
      )
    
    (p_mirna / p_target) + plot_layout(heights = c(1, 1.3))
  })
  
  # Final plot
  final_plot <- wrap_plots(paired_plots, ncol = targets_per_row) +
    plot_annotation(
      title = paste("Expression Analysis:", unique(mirna_data$Gene_ID), "and its Target Genes"),
      subtitle = paste("Z-score normalized expression showing", n_targets, "targets with Pearson correlation"),
      caption = paste("miRNA (orange, top) and its ", length(unique(target_data$Gene_ID)), 
                      " target genes (blue, bottom panels) across 6 tissue types.\n",
                      "YL=Young Leaves, ML=Mature Leaves, YR=Young Roots, MR=Mature Roots, TI=Third Internode, FI=Fifth Internode.\n",
                      "Error bars represent 95% confidence intervals across biological replicates."),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        plot.caption = element_text(size = 10, hjust = 0.5)
      )
    )
  
  if (!dir.exists("output_plots")) dir.create("output_plots")
  output_file <- paste0("output_plots/", make.names(unique(mirna_data$Gene_ID)), 
                        "_", "corr_all_samples.png")
  
  ggsave(output_file, final_plot, 
         width = plot_width, 
         height = plot_height,
         dpi = 300, bg = "white")
  
  return(output_file)
}

# Process all files in the directory
results <- lapply(csv_files, function(file) {
  tryCatch({
    # Process data
    processed <- process_mirna_target_data(file)
    
    # Separate miRNA and target data
    mirna_data <- processed$summary %>% filter(Type == "miRNA")
    target_data <- processed$summary %>% filter(Type == "Target")
    
    # Skip if no targets found
    if (nrow(target_data) == 0) {
      message("No targets found in file: ", file)
      return(NULL)
    }
    
    # Create and save plot
    output_path <- create_mirna_target_plot(
      mirna_data, 
      target_data, 
      processed$processed
    )
    
    message("Successfully processed: ", file, " -> Saved to: ", output_path)
    return(output_path)
    
  }, error = function(e) {
    message("Error processing file: ", file, "\n", e$message)
    return(NULL)
  })
})

# Print summary
successful_files <- unlist(results[!sapply(results, is.null)])
cat("\n=== PROCESSING SUMMARY ===\n")
cat("Total files processed:", length(csv_files), "\n")
cat("Successfully processed:", length(successful_files), "\n")
cat("Failed files:", length(csv_files) - length(successful_files), "\n")
if (length(successful_files) > 0) {
  cat("Output files saved in directory: output_plots/\n")
  cat("Example output:", successful_files[1], "\n")
}