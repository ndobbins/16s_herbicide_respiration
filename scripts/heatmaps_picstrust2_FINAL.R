# ============================================================================
# PICRUsT2 Carbon Pathways Visualization
# Author: Nathan Dobbins
# Date: 2025-08-12
# Description: Heatmap visualization of carbon cycling pathways from PICRUsT2 output
# ============================================================================

# Set working directory
setwd("C:/Users/Nathan/Documents/Practicum/herbicide_soil_DOBBINS/")

# ============================================================================
# INSTALL AND LOAD PACKAGES
# ============================================================================

# Install packages if not already installed
if (!require("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


# Load required libraries
library(readr)
library(pheatmap)
library(tidyverse)
library(DESeq2)

# ============================================================================
# DEFINE CARBON PATHWAYS AND HELPER FUNCTIONS
# ============================================================================

# Define carbon pathways of interest
carbon_pathways <- c(
  "superpathway of C1 compounds oxidation to CO2",
  "formaldehyde oxidation I", 
  "methyl-coenzyme M oxidation to CO2",
  "glucose degradation (oxidative)",
  "octane oxidation"
)

# Function to filter pathways based on keywords
filter_carbon_pathways <- function(pathway_matrix, keywords) {
  pathway_names <- rownames(pathway_matrix)
  
  # Use exact matching instead of pattern matching
  matching_indices <- which(pathway_names %in% keywords)
  
  return(pathway_matrix[matching_indices, , drop = FALSE])
}

# Function to calculate p-values for all carbon pathways
# *** CODED USING CLAUDE SONNET 4***
calculate_pathway_pvalues <- function(picrust2_data, metadata, carbon_pathways, control_group = "Control") {
  
  # Filter for carbon pathways
  carbon_matrix <- filter_carbon_pathways(picrust2_data, carbon_pathways)
  
  # Get all treatment groups
  treatment_groups <- unique(metadata$Group)
  control_group_name <- paste(control_group, "0", sep = "_")
  
  # Initialize results
  results <- data.frame(
    Pathway = character(),
    Treatment = character(),
    Control_Mean = numeric(),
    Treatment_Mean = numeric(),
    Fold_Change = numeric(),
    P_Value = numeric(),
    Significance = character(),
    stringsAsFactors = FALSE
  )
  
  # Get control samples
  control_samples <- rownames(metadata)[metadata$Group == control_group_name]
  
  # Loop through each pathway and each treatment
  for(pathway in rownames(carbon_matrix)) {
    for(treatment in treatment_groups) {
      if(treatment != control_group_name) {
        
        # Get treatment samples
        treatment_samples <- rownames(metadata)[metadata$Group == treatment]
        
        # Get pathway values for control and treatment
        control_values <- carbon_matrix[pathway, control_samples]
        treatment_values <- carbon_matrix[pathway, treatment_samples]
        
        # Calculate means
        control_mean <- mean(control_values, na.rm = TRUE)
        treatment_mean <- mean(treatment_values, na.rm = TRUE)
        
        # Calculate fold change
        fold_change <- treatment_mean - control_mean  # Since data is log-transformed
        
        # Perform t-test
        if(length(control_values) > 1 && length(treatment_values) > 1) {
          t_test <- t.test(control_values, treatment_values)
          p_val <- t_test$p.value
        } else {
          p_val <- NA
        }
        
        # Assign significance levels
        if(!is.na(p_val)) {
          if(p_val < 0.001) significance <- "***"
          else if(p_val < 0.01) significance <- "**"
          else if(p_val < 0.05) significance <- "*"
          else significance <- ""
        } else {
          significance <- ""
        }
        
        # Store results
        results <- rbind(results, data.frame(
          Pathway = pathway,
          Treatment = treatment,
          Control_Mean = control_mean,
          Treatment_Mean = treatment_mean,
          Fold_Change = fold_change,
          P_Value = p_val,
          Significance = significance
        ))
      }
    }
  }
  
  return(results)
}

# Function to create significance matrix
create_significance_matrix <- function(pvalue_results, pathway_order, treatment_order) {
  # Initialize empty matrix
  sig_matrix <- matrix("", nrow = length(pathway_order), ncol = length(treatment_order))
  rownames(sig_matrix) <- pathway_order
  colnames(sig_matrix) <- treatment_order
  
  # Fill in significance values
  for(i in 1:nrow(pvalue_results)) {
    pathway <- pvalue_results$Pathway[i]
    treatment <- pvalue_results$Treatment[i]
    significance <- pvalue_results$Significance[i]
    
    if(pathway %in% pathway_order && treatment %in% treatment_order) {
      sig_matrix[pathway, treatment] <- significance
    }
  }
  
  return(sig_matrix)
}

# Save pheatmap function
# Source: https://rdrr.io/github/GrahamHamilton/pipelineTools/src/R/save_pheatmap.R
save_pheatmap_png <- function(plot = NULL,
                              filename = NULL,
                              path = NULL,
                              width = 1200,
                              height = 1000,
                              res = 150) {
  # Path to save file
  if (!is.null(path)){
    file_path <- file.path(path, filename)  # Remove the fsep parameter
  }else{
    file_path <- filename
  }
  
  png(file_path, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(plot$gtable)
  dev.off()
}

# ============================================================================
# LOAD AND PREPARE DATA
# ============================================================================

# Load PICRUsT2 data
picrust2_df <- read.delim(
  "C:/Users/Nathan/picrust2_out_pipeline_04/pathways_out/path_abun_unstrat_descrip.tsv.gz",
  header = TRUE,
  row.names = 2,
  check.names = FALSE
)
head(picrust2_df)

# Remove the first column (which is 'pathway' now, since it wasn't used as row names)
picrust2_df <- picrust2_df[, -1]

# Convert to matrix for pheatmap
picrust2_matrix <- as.matrix(picrust2_df)

head(picrust2_df)

# Load metadata
map <- read.csv("data/MapFileHerb(in).csv", row.names = 1, check.names = FALSE)

# Make sure all sample names are trimmed
rownames(map) <- trimws(rownames(map))
colnames(picrust2_df) <- trimws(colnames(picrust2_df))

# Create a 'Group' column from Herbicide and Concentration
map$Group <- paste(map$Herbicide, map$Concentration, sep = "_")

# ============================================================================
# HEATMAP FOR ALL TREATMENTS
# ============================================================================

# Define target herbicide
herbicide_of_interest <- "2-4 D"

# Get samples from that herbicide AND from Control
samples_to_plot <- rownames(map)[
  map$Herbicide == herbicide_of_interest | map$Herbicide == "Control" | map$Herbicide == "Glyphosate" 
]

# Display Selected Sample IDs and Count
cat("Number of samples selected:", length(samples_to_plot), "\n")
print(samples_to_plot)

# Subset abundance matrix and map
picrust2_subset <- picrust2_df[, samples_to_plot, drop = FALSE]
map_subset <- map[samples_to_plot, , drop = FALSE]

# Transpose and align
picrust2_t <- t(picrust2_subset)
picrust2_t <- picrust2_t[rownames(map_subset), , drop = FALSE]

# Group-averaging
group_levels <- unique(map_subset$Group)
averaged_list <- lapply(group_levels, function(g) {
  rows <- which(map_subset$Group == g)
  colMeans(picrust2_t[rows, , drop = FALSE])
})
names(averaged_list) <- group_levels
averaged_matrix <- do.call(cbind, averaged_list)

# Filter for carbon pathways instead of top 50 variable
carbon_matrix <- filter_carbon_pathways(averaged_matrix, carbon_pathways)

# Transform and normalize
ko_log <- log10(carbon_matrix + 1)
ko_scaled <- t(scale(t(ko_log)))
ko_scaled <- ko_scaled[apply(ko_scaled, 1, function(x) all(is.finite(x))), ]
ko_scaled <- ko_scaled[, apply(ko_scaled, 2, function(x) all(is.finite(x)))]

# Reorder data for heatmap columns
desired_order <- c("Control_0", "2-4 D_0.125", "2-4 D_0.25", "2-4 D_0.5", "2-4 D_1", "Glyphosate_0.125", "Glyphosate_0.25", "Glyphosate_0.5", "Glyphosate_1")

# Reorder columns
ko_scaled <- ko_scaled[, desired_order]

# Final heatmap
pheatmap(ko_scaled, 
         main = "Carbon Pathways Across Herbicide Treatments",
         fontsize = 8,
         cluster_cols = FALSE,
         angle_col = 45,
         fontsize_row = 7)

# ============================================================================
# HEATMAP FOR ALL TREATMENTS NORMALIZED TO CONTROL
# ============================================================================

# Find the Control column index
control_colname <- grep("^Control", colnames(averaged_matrix), value = TRUE)[1]

# Subtract Control column values from all columns (center on Control)
averaged_matrix_norm <- sweep(averaged_matrix, 1, averaged_matrix[, control_colname], FUN = "-")

# Filter for carbon pathways
carbon_matrix_norm <- filter_carbon_pathways(averaged_matrix_norm, carbon_pathways)

# Scale rows if you want (optional since centering done)
ko_scaled <- t(scale(t(carbon_matrix_norm)))
ko_scaled <- ko_scaled[apply(ko_scaled, 1, function(x) all(is.finite(x))), ]
ko_scaled <- ko_scaled[, apply(ko_scaled, 2, function(x) all(is.finite(x)))]

# Reorder data for heatmap columns
desired_order <- c("Control_0", "2-4 D_0.125", "2-4 D_0.25", "2-4 D_0.5", "2-4 D_1", "Glyphosate_0.125", "Glyphosate_0.25", "Glyphosate_0.5", "Glyphosate_1")

# Reorder columns
ko_scaled <- ko_scaled[, desired_order]

# Plot heatmap, Normalized to control
heatmap_all <- pheatmap(ko_scaled,
                        main = "All Treatments",
                        fontsize = 8,
                        cluster_cols = FALSE,
                        angle_col = 45,
                        fontsize_row = 6
)
heatmap_all

# Capture the row order more directly
pathway_order_for_all <- heatmap_all$tree_row$labels[heatmap_all$tree_row$order]

print("Captured pathway order:")
print(pathway_order_for_all)

save_pheatmap_png(heatmap_all, filename = "carbon_pathways_heatmap_all.png", path = "Figures",width = 1700, height = 1000, res = 300)

# ============================================================================
# HEATMAP FOR 2-4D TREATMENT
# ============================================================================

# Define target herbicide
herbicide_of_interest <- "2-4 D"
# Select samples that are either Control or 2-4 D treatments ONLY
samples_to_plot <- rownames(map)[map$Herbicide %in% c("Control", herbicide_of_interest)]
# Display Selected Sample IDs and Count
cat("Number of samples selected:", length(samples_to_plot), "\n")
print(samples_to_plot)
# Subset picrust2 data and map metadata for selected samples
picrust2_subset <- picrust2_df[, samples_to_plot, drop = FALSE]
map_subset <- map[samples_to_plot, , drop = FALSE]
# Transpose to have samples as rows, pathways as columns
picrust2_t <- t(picrust2_subset)
picrust2_t <- picrust2_t[rownames(map_subset), , drop = FALSE]

# CALCULATE P-VALUES BEFORE GROUP AVERAGING (using individual samples)
# Transpose so pathways are rows, samples are columns
picrust2_for_pvalues <- t(picrust2_t)
# Then call the function
pvalue_results_24d <- calculate_pathway_pvalues(
  picrust2_data = picrust2_for_pvalues,  # Transposed data
  metadata = map_subset,
  carbon_pathways = carbon_pathways
)

# Group-averaging by Group (Control and each 2-4 D concentration)
group_levels <- unique(map_subset$Group)
averaged_list <- lapply(group_levels, function(g) {
  rows <- which(map_subset$Group == g)
  colMeans(picrust2_t[rows, , drop = FALSE])
})
names(averaged_list) <- group_levels
averaged_matrix <- do.call(cbind, averaged_list)
# Normalize to Control group: subtract Control column from all columns
control_colname <- grep("^Control", colnames(averaged_matrix), value = TRUE)[1]
averaged_matrix_norm <- sweep(averaged_matrix, 1, averaged_matrix[, control_colname], FUN = "-")
# Filter for carbon pathways instead of top 50 variable
carbon_matrix_norm <- filter_carbon_pathways(averaged_matrix_norm, carbon_pathways)
# Scale rows
ko_scaled <- t(scale(t(carbon_matrix_norm)))
ko_scaled <- ko_scaled[apply(ko_scaled, 1, function(x) all(is.finite(x))), ]
ko_scaled <- ko_scaled[, apply(ko_scaled, 2, function(x) all(is.finite(x)))]
# Reorder data for heatmap columns
desired_order <- c("Control_0", "2-4 D_0.125", "2-4 D_0.25", "2-4 D_0.5", "2-4 D_1")
# Reorder columns
ko_scaled <- ko_scaled[, desired_order]
# Only keep pathways that exist in the current ko_scaled matrix
available_pathways <- rownames(ko_scaled)
ko_scaled_ordered <- ko_scaled

# CREATE SIGNIFICANCE MATRIX
sig_matrix_24d <- create_significance_matrix(
  pvalue_results = pvalue_results_24d,
  pathway_order = available_pathways,
  treatment_order = desired_order
)

# Change column names to show just percentages
colnames(ko_scaled_ordered) <- c("Control", "12.5%", "25%", "50%", "100%")

# Also update the significance matrix column names to match
colnames(sig_matrix_24d) <- c("Control", "12.5%", "25%", "50%", "100%")

# Create heatmap with ordered rows AND significance overlay
heatmap_24d <- pheatmap(ko_scaled_ordered,
                        main = "2,4-D",
                        fontsize = 28,
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        angle_col = 45,
                        fontsize_row = 26,
                        fontsize_col = 15,
                        display_numbers = sig_matrix_24d,
                        number_color = "black",
                        fontsize_number = 60,
                        border_color = "grey60")

heatmap_24d

save_pheatmap_png(heatmap_24d, 
                  filename = "carbon_pathways_heatmap_24d_sig.png", 
                  path = "Figures", 
                  width = 5000, 
                  height = 2700, 
                  res = 300)

# ============================================================================
# HEATMAP FOR GLYPHOSATE TREATMENT
# ============================================================================

# Define target herbicide
herbicide_of_interest <- "Glyphosate"
# Select samples that are either Control or Glyphosate treatments ONLY
samples_to_plot <- rownames(map)[map$Herbicide %in% c("Control", herbicide_of_interest)]
# Display Selected Sample IDs and Count
cat("Number of samples selected:", length(samples_to_plot), "\n")
print(samples_to_plot)
# Subset picrust2 data and map metadata for selected samples
picrust2_subset <- picrust2_df[, samples_to_plot, drop = FALSE]
map_subset <- map[samples_to_plot, , drop = FALSE]
# Transpose to have samples as rows, pathways as columns
picrust2_t <- t(picrust2_subset)
picrust2_t <- picrust2_t[rownames(map_subset), , drop = FALSE]

# CALCULATE P-VALUES BEFORE GROUP AVERAGING (using individual samples)
# Transpose so pathways are rows, samples are columns
picrust2_for_pvalues <- t(picrust2_t)
# Then call the function
pvalue_results_glyph <- calculate_pathway_pvalues(
  picrust2_data = picrust2_for_pvalues,  # Transposed data
  metadata = map_subset,
  carbon_pathways = carbon_pathways
)

# Group-averaging by Group (Control and each Glyphosate concentration)
group_levels <- unique(map_subset$Group)
averaged_list <- lapply(group_levels, function(g) {
  rows <- which(map_subset$Group == g)
  colMeans(picrust2_t[rows, , drop = FALSE])
})
names(averaged_list) <- group_levels
averaged_matrix <- do.call(cbind, averaged_list)
# Normalize to Control group: subtract Control column from all columns
control_colname <- grep("^Control", colnames(averaged_matrix), value = TRUE)[1]
averaged_matrix_norm <- sweep(averaged_matrix, 1, averaged_matrix[, control_colname], FUN = "-")
# Filter for carbon pathways instead of top 50 variable
carbon_matrix_norm <- filter_carbon_pathways(averaged_matrix_norm, carbon_pathways)
# Scale rows (optional)
ko_scaled <- t(scale(t(carbon_matrix_norm)))
ko_scaled <- ko_scaled[apply(ko_scaled, 1, function(x) all(is.finite(x))), ]
ko_scaled <- ko_scaled[, apply(ko_scaled, 2, function(x) all(is.finite(x)))]
# Reorder data for heatmap columns
desired_order <- c("Control_0", "Glyphosate_0.125", "Glyphosate_0.25", "Glyphosate_0.5", "Glyphosate_1")
# Reorder columns
ko_scaled <- ko_scaled[, desired_order]
# Only keep pathways that exist in the current ko_scaled matrix
available_pathways <- rownames(ko_scaled)
ko_scaled_ordered <- ko_scaled

# CREATE SIGNIFICANCE MATRIX
sig_matrix_glyph <- create_significance_matrix(
  pvalue_results = pvalue_results_glyph,
  pathway_order = available_pathways,
  treatment_order = desired_order
)

# Change column names to show just percentages
colnames(ko_scaled_ordered) <- c("Control", "12.5%", "25%", "50%", "100%")

# Also update the significance matrix column names to match
colnames(sig_matrix_glyph) <- c("Control", "12.5%", "25%", "50%", "100%")

# Create heatmap with ordered rows AND significance overlay
heatmap_glyph <- pheatmap(ko_scaled_ordered,
                          main = "Glyphosate",
                          fontsize = 28,
                          cluster_rows = FALSE,
                          cluster_cols = FALSE,
                          angle_col = 45,
                          fontsize_row = 26,
                          fontsize_col = 15,
                          display_numbers = sig_matrix_glyph,
                          number_color = "black",
                          fontsize_number = 60,
                          border_color = "grey60")

heatmap_glyph

save_pheatmap_png(heatmap_glyph, filename = "carbon_pathways_heatmap_glyph_sig.png", path = "Figures",width = 5000, 
                  height = 2700, 
                  res = 300)

