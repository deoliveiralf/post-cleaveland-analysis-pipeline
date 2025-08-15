# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(readr)
library(stringr)  # Added for string operations

# Read the data (assuming the file is in your working directory)
setwd("C:/Users/Leandro/OneDrive - usp.br/LIGNINLAB/Projeto FAPESP CNPq - Setaria miRNAs/R analysis/k-means_correlation/table_expression")
# Replace with the actual path to your CSV file
data <- read_csv("novel359_Chr02_5549_expression_table.csv")

# Extract sample information
sample_cols <- colnames(data)[4:21]  # Columns with expression data
print("Sample columns found:")
print(sample_cols)

# Create a function to process and reshape the data
process_expression_data <- function(data) {
  
  # Reshape data from wide to long format
  long_data <- data %>%
    select(Type, ID, Description, all_of(sample_cols)) %>%
    pivot_longer(cols = all_of(sample_cols), 
                 names_to = "Sample", 
                 values_to = "Expression") %>%
    mutate(
      # Extract tissue type and replicate information
      Tissue = case_when(
        grepl("Fifth_Internode", Sample) ~ "Fifth_Internode",
        grepl("Third_Internode", Sample) ~ "Third_Internode", 
        grepl("Mature_Leaves", Sample) ~ "Mature_Leaves",
        grepl("Young_Leaves", Sample) ~ "Young_Leaves",
        grepl("Mature_Roots", Sample) ~ "Mature_Roots",
        grepl("Young_Roots", Sample) ~ "Young_Roots"
      ),
      Replicate = as.numeric(stringr::str_extract(Sample, "\\d+$")),
      
      # Create ordered tissue factor for plotting
      Tissue_ordered = factor(Tissue, levels = c(
        "Young_Leaves", "Mature_Leaves", 
        "Young_Roots", "Mature_Roots",
        "Third_Internode", "Fifth_Internode"
      )),
      
      # Log transform expression (add 1 to handle zeros)
      Log_Expression = log2(Expression + 1),
      
      # Z-score standardization per gene
      Gene_ID = ID
    ) %>%
    group_by(Gene_ID) %>%
    mutate(
      Z_score = scale(Log_Expression)[,1]
    ) %>%
    ungroup()
  
  return(long_data)
}

# Process the data
processed_data <- process_expression_data(data)

# Create summary statistics for plotting
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
  mutate(
    # Create tissue position for x-axis
    Tissue_position = as.numeric(Tissue_ordered)
  )

# Separate miRNA and target data
mirna_data <- summary_data %>% filter(Type == "miRNA")
target_data <- summary_data %>% filter(Type == "Target")

create_mirna_target_plot <- function(mirna_data, target_data, processed_data, targets_per_row = 6) {
  
  # Create tissue labels
  tissue_labels <- c("YL", "ML", "YR", "MR", "TI", "FI")
  
  # Create paired plots
  paired_plots <- lapply(unique(target_data$Gene_ID), function(target_id) {
    
    # Filter data for this target
    target_subset <- target_data %>% filter(Gene_ID == target_id)
    
    # SAFELY Calculate correlation stats
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
    
    # miRNA plot with small dots (size = 1)
    p_mirna <- ggplot(mirna_data, aes(x = Tissue_position)) +
      geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), 
                  alpha = 0.2, fill = "#E69F00") +
      geom_line(aes(y = Mean_Z_score), color = "#E69F00", linewidth = 1) +
      geom_point(aes(y = Mean_Z_score), color = "#E69F00", size = 1) + # Small dots
      theme_minimal() +
      labs(title = "miRNA", y = "Z-score") +
      ylim(-3, 3) +
      theme(
        axis.title.x = element_blank(),
        plot.title = element_text(color = "#E69F00", face = "bold", size = 10)
      )
    
    # Create correlation label
    cor_label <- if(is.na(cor_result$r)) {
      "r = NA (no data)"
    } else {
      paste0("r = ", cor_result$r, ", p = ", cor_result$p_value)
    }
    
    # Target plot with small dots (size = 1)
    p_target <- ggplot(target_subset, aes(x = Tissue_position)) +
      geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), 
                  alpha = 0.2, fill = "#56B4E9") +
      geom_line(aes(y = Mean_Z_score), color = "#56B4E9", linewidth = 1) +
      geom_point(aes(y = Mean_Z_score), color = "#56B4E9", size = 1) + # Small dots
      scale_x_continuous(breaks = 1:6, labels = tissue_labels) +
      labs(
        title = paste0(target_id, "\n", cor_label),
        x = "Tissue", 
        y = "Z-score"
      ) +
      ylim(-3, 3) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 10, hjust = 0.5, margin = margin(b = 5))
      )
    
    # Combine plots
    (p_mirna / p_target) + 
      plot_layout(heights = c(1, 1.3))
  })
  
  # Arrange all plots
  wrap_plots(paired_plots, ncol = targets_per_row) +
    plot_annotation(
      title = "Expression Analysis: miRNA and Target Genes",
      subtitle = "Pearson correlation (r) and p-values shown for each target (when calculable)",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      )
    )
}

# Run with error handling
tryCatch({
  final_plot <- create_mirna_target_plot(mirna_data, target_data, processed_data, targets_per_row = 6)
  ggsave("mirna_target_plot_with_cor.png", final_plot, 
         width = 18, height = 30, dpi = 300, bg = "white")
}, error = function(e) {
  message("Error creating plot: ", e$message)
})


# Make sure to include processed_data when calling the function
final_plot <- create_mirna_target_plot(mirna_data, target_data, processed_data, targets_per_row = 6)

# Save with increased height and adjusted DPI
ggsave("mirna_target_expression_plot.png", final_plot, 
       width = 18, height = 36, dpi = 300, bg = "white")

# Add overall title and caption
final_plot <- final_plot + 
  plot_annotation(
    title = "Expression Analysis: miRNA novel242_Chr06_28920 and its Target Genes",
    subtitle = "Z-score normalized expression across different tissue types and developmental stages",
    caption = paste0("miRNA (orange, top) and its ", length(unique(target_data$Gene_ID)), 
                     " target genes (blue, bottom panels) across 6 tissue types.\n",
                     "YL=Young Leaves, ML=Mature Leaves, YR=Young Roots, MR=Mature Roots, TI=Third Internode, FI=Fifth Internode.\n",
                     "Error bars represent 95% confidence intervals across biological replicates."),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      plot.caption = element_text(size = 10, hjust = 0)
    )
  )

# Display the plot
print(final_plot)

# Save the plot (increased height to accommodate paired plots)
ggsave("mirna_target_expression_plot.png", final_plot, 
       width = 18, height = 32, dpi = 300, bg = "white")

# Create a correlation plot between miRNA and targets
create_correlation_plot <- function(processed_data) {
  
  # Get miRNA expression
  mirna_expr <- processed_data %>%
    filter(Type == "miRNA") %>%
    select(Sample, Tissue, miRNA_Expression = Expression)
  
  # Get target expression and merge with miRNA
  correlation_data <- processed_data %>%
    filter(Type == "Target") %>%
    left_join(mirna_expr, by = c("Sample", "Tissue")) %>%
    group_by(Gene_ID, Tissue) %>%
    summarise(
      Target_Mean = mean(Expression, na.rm = TRUE),
      miRNA_Mean = mean(miRNA_Expression, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(Target_Mean) & !is.na(miRNA_Mean))
  
  # Create correlation plot
  cor_plot <- correlation_data %>%
    ggplot(aes(x = log2(miRNA_Mean + 1), y = log2(Target_Mean + 1))) +
    geom_point(aes(color = Tissue), alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    scale_color_brewer(type = "qual", palette = "Set2") +
    theme_minimal() +
    labs(
      title = "miRNA vs Target Expression Correlation",
      x = "miRNA Expression (log2)",
      y = "Target Expression (log2)",
      color = "Tissue Type"
    ) +
    theme(
      legend.position = "bottom"
    )
  
  return(cor_plot)
}

# Create correlation plot
correlation_plot <- create_correlation_plot(processed_data)
print(correlation_plot)

# Save correlation plot
ggsave("mirna_target_correlation.png", correlation_plot, 
       width = 10, height = 8, dpi = 300, bg = "white")

# Print summary statistics
cat("\n=== DATA SUMMARY ===\n")
cat("miRNA ID:", unique(mirna_data$Gene_ID), "\n")
cat("Number of target genes:", length(unique(target_data$Gene_ID)), "\n")
cat("Target gene IDs (first 5):", paste(head(unique(target_data$Gene_ID), 5), collapse = ", "), "\n")
cat("Tissue types:", paste(levels(processed_data$Tissue_ordered), collapse = ", "), "\n")
cat("Total samples per tissue:", unique(summary_data$N_samples)[1], "\n")

# Print expression ranges
cat("\n=== EXPRESSION RANGES ===\n")
mirna_range <- processed_data %>% filter(Type == "miRNA") %>% pull(Expression) %>% range()
target_range <- processed_data %>% filter(Type == "Target") %>% pull(Expression) %>% range()
cat("miRNA expression range:", round(mirna_range, 2), "\n")
cat("Target expression range:", round(target_range, 2), "\n")