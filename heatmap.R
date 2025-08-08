# ===========================
# 6. CHROMATIN ACCESSIBILITY HEATMAP - UPDATED
# ===========================
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(SummarizedExperiment)
library(dplyr)
library(reshape2)

# Define your new comparisons
comparisons <- list(
  list(name = "RAS_vs_Proliferating", contrast = c("group", "RAS", "Proliferating")),
  list(name = "MafK45_RAS_vs_RAS", contrast = c("group", "MafK45_RAS", "RAS")),
  list(name = "Ras_MafK045_KD_vs_RAS", contrast = c("group", "Ras_MafK045_KD", "RAS")),
  list(name = "Ras_MafK138_KD_vs_RAS", contrast = c("group", "Ras_MafK138_KD", "RAS")),
  list(name = "MafK045_KD_Ras_vs_RAS", contrast = c("group", "MafK045_KD_Ras", "RAS")),
  list(name = "MafK138_KD_Ras_vs_RAS", contrast = c("group", "MafK138_KD_Ras", "RAS"))
)

# ===========================
# PERFORM DIFFERENTIAL ACCESSIBILITY ANALYSIS
# ===========================

# Initialize list to store results
da_results_list <- list()
DARs_list <- list()

# Loop through each comparison
for(i in 1:length(comparisons)) {
  comp <- comparisons[[i]]
  cat("Processing comparison:", comp$name, "\n")
  
  # Run DESeq2 analysis
  da_result <- results(dds, contrast = comp$contrast, alpha = 0.05)
  da_result <- as.data.frame(da_result)
  da_result$peak_id <- rownames(atac_count_data_matrix)
  
  # Store results
  da_results_list[[comp$name]] <- da_result
  
  # Extract significant DARs
  DARs <- subset(da_result, padj < 0.05 & abs(log2FoldChange) > 1)
  DARs$source <- comp$name
  DARs_list[[comp$name]] <- DARs
  
  cat("  - Total peaks tested:", nrow(da_result), "\n")
  cat("  - Significant DARs:", nrow(DARs), "\n")
  cat("  - Up-regulated:", sum(DARs$log2FoldChange > 1), "\n")
  cat("  - Down-regulated:", sum(DARs$log2FoldChange < -1), "\n\n")
}

# ===========================
# COMBINE ALL SIGNIFICANT DARs
# ===========================

# Combine all significant DARs
all_DARs_combined <- do.call(rbind, DARs_list)
cat("Total DARs from all comparisons:", nrow(all_DARs_combined), "\n")

# Keep only one row per unique peak_id (first occurrence)
all_DARs_unique <- all_DARs_combined[!duplicated(all_DARs_combined$peak_id), ]
cat("Unique DARs after removing duplicates:", nrow(all_DARs_unique), "\n")

# Clean up the data
sig_peaks <- dplyr::select(all_DARs_unique, -source)

# Extract genomic coordinates
sig_peaks$chr <- sapply(strsplit(sig_peaks$peak_id, "_"), "[", 1)
sig_peaks$start <- as.numeric(sapply(strsplit(sig_peaks$peak_id, "_"), "[", 2))
sig_peaks$end <- as.numeric(sapply(strsplit(sig_peaks$peak_id, "_"), "[", 3))

# ===========================
# PREPARE DATA FOR HEATMAP
# ===========================

# Extract count data for significant peaks only
logit_DARs <- row.names(atac_count_data) %in% sig_peaks$peak_id
table(logit_DARs)

# Variance Stabilizing Transformation
vsd_atac <- vst(dds, blind = FALSE)  
vsd_mat <- assay(vsd_atac)  
dar_signals <- vsd_mat[logit_DARs, ]

cat("Matrix dimensions for heatmap:", dim(dar_signals), "\n")

# Z-score normalization
z_score <- function(x) {(x - mean(x)) / sd(x)}
sig_peak_scaled <- t(apply(dar_signals, 1, z_score))

# Handle any NaN values (peaks with zero variance)
sig_peak_scaled[is.na(sig_peak_scaled)] <- 0

# ===========================
# CREATE ANNOTATIONS FOR HEATMAP
# ===========================

# Column annotation (samples) - update colors for your groups
unique_groups <- unique(sample_info$group)
n_groups <- length(unique_groups)

# Create color palette for your groups
group_colors <- RColorBrewer::brewer.pal(min(n_groups, 11), "Spectral")
names(group_colors) <- unique_groups

# If you have more than 11 groups, extend with additional colors
if(n_groups > 11) {
  additional_colors <- rainbow(n_groups - 11)
  group_colors <- c(group_colors, additional_colors)
  names(group_colors) <- unique_groups
}

col_anno <- HeatmapAnnotation(
  Group = sample_info$group,
  col = list(Group = group_colors),
  Rep_text = anno_text(
    sample_info$Replication,
    gp = gpar(fontsize = 8, col = "black",
              fill = c("rep1" = "gray", 
                       "rep2" = "orange",
                       "rep3" = "red")),
    location = 0.5,
    just = "center",
    rot = 0
  ),
  annotation_name_gp = gpar(fontsize = 10)
)

# Row annotation (peaks) - assuming you have annotated_peaks
if(exists("annotated_peaks")) {
  # Simplify annotation for visualization
  simplified_anno <- sapply(annotated_peaks$annotation, function(x) {
    if(grepl("Promoter", x)) return("Promoter")
    else if(grepl("5' UTR", x)) return("5' UTR") 
    else if(grepl("3' UTR", x)) return("3' UTR")
    else if(grepl("Exon", x)) return("Exon")
    else if(grepl("Intron", x)) return("Intron")
    else if(grepl("Downstream", x)) return("Downstream")
    else if(grepl("Distal", x)) return("Distal Intergenic")
    else return("Other")
  })
  
  # Create row annotation
  row_anno <- HeatmapAnnotation(
    Annotation = simplified_anno,
    Log2FC = annotated_peaks$log2FoldChange,
    Direction = ifelse(annotated_peaks$log2FoldChange > 0, "Up", "Down"),
    which = "row",
    col = list(
      Annotation = c("Promoter" = "#E31A1C", "5' UTR" = "#FF7F00", "3' UTR" = "#FFFF33",
                     "Exon" = "#1F78B4", "Intron" = "#33A02C", "Downstream" = "#6A3D9A",
                     "Distal Intergenic" = "#B15928", "Other" = "#999999"),
      Log2FC = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
      Direction = c("Up" = "red", "Down" = "blue")
    ),
    annotation_name_gp = gpar(fontsize = 10)
  )
} else {
  # Create minimal row annotation without genomic annotations
  row_anno <- HeatmapAnnotation(
    Log2FC = sig_peaks$log2FoldChange,
    Direction = ifelse(sig_peaks$log2FoldChange > 0, "Up", "Down"),
    which = "row",
    col = list(
      Log2FC = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
      Direction = c("Up" = "red", "Down" = "blue")
    ),
    annotation_name_gp = gpar(fontsize = 10)
  )
}

# ===========================
# CREATE MAIN HEATMAP
# ===========================

# Main heatmap
main_heatmap <- Heatmap(
  sig_peak_scaled,
  name = "Accessibility\n(Z-score)",
  
  # Color scheme
  col = colorRamp2(c(-2, 0, 2), c("#1f77b4", "white","#ff7f0e")),
  
  # Annotations
  top_annotation = col_anno,
  left_annotation = row_anno,
  
  # Clustering
  row_km = 4,
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  
  # Titles
  row_title = paste0("Significant ATAC-seq Peaks (n=", nrow(sig_peaks), ")"),
  column_title = "Samples",
  
  # Display options
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10),
  
  # Heatmap size
  width = unit(8, "cm"),
  height = unit(12, "cm")
)

# Draw the heatmap
pdf("chromatin_accessibility_heatmap_new_comparisons.pdf", width = 12, height = 10)
draw(main_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# Display in R
draw(main_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")

# ===========================
# STORE CLUSTER RESULTS
# ===========================

# Store cluster result
ht_obj <- draw(main_heatmap)

# Get row indices (ordered) grouped by clusters
cluster_assignments <- row_order(ht_obj)

# Check names and size of clusters
str(cluster_assignments)

# Extract cluster data and save
for(i in 1:length(cluster_assignments)) {
  cluster_data <- sig_peaks[cluster_assignments[[i]], ]
  filename <- paste0('cluster', i, '_data_new_comparisons.csv')
  write.csv(cluster_data, file = filename)
  cat("Cluster", i, ":", nrow(cluster_data), "peaks saved to", filename, "\n")
}

# ===========================
# METHOD 2: HEATMAP WITH PEAK CLUSTERING BY DIRECTION
# ===========================

# Split peaks by direction of change
up_peak_indices <- which(sig_peaks$log2FoldChange > 0)
down_peak_indices <- which(sig_peaks$log2FoldChange < 0)

# Create split annotation
split_anno <- rep("Down-regulated", nrow(sig_peak_scaled))
split_anno[up_peak_indices] <- "Up-regulated"

# Heatmap with split
split_heatmap <- Heatmap(
  sig_peak_scaled,
  name = "Accessibility\n(Z-score)",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  
  # Split by regulation direction
  row_split = split_anno,
  
  # Annotations
  top_annotation = col_anno,
  left_annotation = row_anno,
  
  # Clustering within splits
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  
  # Titles
  row_title = "ATAC-seq Peaks",
  column_title = "Samples",
  
  # Display
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10),
  row_title_gp = gpar(fontsize = 12),
  
  # Size
  width = unit(8, "cm"),
  height = unit(12, "cm")
)

# Draw split heatmap
pdf("chromatin_accessibility_split_heatmap_new_comparisons.pdf", width = 12, height = 10)
draw(split_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

draw(split_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")

# ===========================
# SAVE INDIVIDUAL COMPARISON RESULTS
# ===========================

# Save individual comparison results
for(i in 1:length(comparisons)) {
  comp_name <- comparisons[[i]]$name
  filename <- paste0("DA_results_", comp_name, ".csv")
  write.csv(da_results_list[[comp_name]], file = filename, row.names = TRUE)
  
  # Also save significant DARs for each comparison
  sig_filename <- paste0("Significant_DARs_", comp_name, ".csv")
  write.csv(DARs_list[[comp_name]], file = sig_filename, row.names = TRUE)
}

# ===========================
# SUMMARY STATISTICS
# ===========================

# Print summary information
cat("\n=== HEATMAP SUMMARY ===\n")
cat("Total comparisons analyzed:", length(comparisons), "\n")
cat("Comparison names:\n")
for(comp in comparisons) {
  cat("  -", comp$name, "\n")
}
cat("\nTotal unique significant peaks plotted:", nrow(sig_peaks), "\n")
cat("Up-regulated peaks:", sum(sig_peaks$log2FoldChange > 0), "\n")
cat("Down-regulated peaks:", sum(sig_peaks$log2FoldChange < 0), "\n")
cat("Samples included:", ncol(sig_peak_scaled), "\n")

if(exists("simplified_anno")) {
  cat("Annotation distribution:\n")
  print(table(simplified_anno))
}

# Save the processed data
heatmap_data_output <- data.frame(
  sig_peak_scaled, 
  log2FC = sig_peaks$log2FoldChange,
  padj = sig_peaks$padj,
  chr = sig_peaks$chr,
  start = sig_peaks$start,
  end = sig_peaks$end
)

if(exists("simplified_anno")) {
  heatmap_data_output$annotation <- simplified_anno
}

write.csv(heatmap_data_output, 
          "significant_peaks_heatmap_data_new_comparisons.csv", 
          row.names = TRUE)

cat("\n=== FILES SAVED ===\n")
cat("- chromatin_accessibility_heatmap_new_comparisons.pdf\n")
cat("- chromatin_accessibility_split_heatmap_new_comparisons.pdf\n")
cat("- significant_peaks_heatmap_data_new_comparisons.csv\n")
cat("- Individual DA results and cluster files\n")

print("Analysis completed successfully!")