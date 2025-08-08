#########################################################
# ATAC-seq Analysis Pipeline
# Author: Kehan Li
# Date: 2025-06-13
# Description: Pipeline for differential accessibility and motif analysis
#########################################################

# Load required libraries
library(data.table)
library(here)
library(DESeq2)
library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(tidyr)
library(corrplot)


# 1. DATA LOADING AND PREPARATION#####
library(here)
# Load ATAC-seq count data
setwd("~/Desktop/CZ/CZ_ras_ATAC")
atac_count_data <- read_csv("data/clean data/ATAC_data_CZ.csv")
atac_count_data <- as.data.frame(atac_count_data)
atac_count_data$peak_id <- paste(atac_count_data$Chr, 
                                 atac_count_data$Start, 
                                 atac_count_data$End, sep = "_")
peak_ids <- atac_count_data$peak_id
row.names(atac_count_data) <- atac_count_data$peak_id
atac_count_data_matrix <- atac_count_data %>%
  select(-c(1:3), -peak_id)
rownames(atac_count_data_matrix) <- rownames(atac_count_data)


# create sample info
sample_names <- colnames(atac_count_data)[4:15]  # Explicitly take columns 4-15

# Create the sample_info dataframe
sample_info <- data.frame(
  sample_ID = sample_names,
  Replication = c("Rep1", "Rep2", "Rep1", "Rep2", "Rep1", "Rep2", 
                  "Rep1", "Rep2", "Rep1", "Rep2", "Rep1", "Rep1"),
  group = c("RAS", "RAS", "MafK45_RAS", "MafK45_RAS", "Proliferating", "Proliferating",
            "Ras_MafK045_KD", "Ras_MafK045_KD", "Ras_MafK138_KD", "Ras_MafK138_KD",
            "MafK045_KD_Ras", "MafK138_KD_Ras")
)



# Ensure sample order matches between count matrix and metadata
sample_metadata <- sample_info[match(colnames(atac_count_data_matrix), sample_info$sample_ID), ]


# Create GenomicRanges object for peaks
peak_ranges <- GRanges(
  seqnames = atac_count_data$Chr,
  ranges = IRanges(start = atac_count_data$Start, end = atac_count_data$End),
  peak_id = atac_count_data$peak_id
)


# 2. QUALITY CONTROL AND FILTERING#####

# Calculate basic statistics
total_peaks <- nrow(atac_count_data_matrix)
total_reads_per_sample <- colSums(atac_count_data_matrix)
reads_per_peak <- rowSums(atac_count_data_matrix)

# Print basic statistics
cat("Total peaks:", total_peaks, "\n")
cat("Total reads per sample:\n")
print(total_reads_per_sample)

# Filter low-count peaks (keep peaks with at least 10 reads in at least 2 samples)
keep_peaks <- rowSums(atac_count_data_matrix >= 10) >= 2
filtered_count_matrix <- atac_count_data_matrix[keep_peaks, ]
filtered_peak_ranges <- peak_ranges[keep_peaks]

rownames(filtered_count_matrix)<- rownames(atac_count_data_matrix)

cat("Peaks after filtering:", nrow(filtered_count_matrix), "\n")

# Library size normalization factors
lib_sizes <- colSums(filtered_count_matrix)
norm_factors <- lib_sizes / mean(lib_sizes)


# 3. EXPLORATORY DATA ANALYSIS####
#rlog
library(ggrepel)
rld <- rlog(dds, blind = TRUE)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))

#vst transformation
library(SummarizedExperiment)
library(ggplot2)
vst.r <- vst(dds,blind = TRUE)

#plotPCA from DEseq###########
# Get PCA data with variance explained
data <- plotPCA(rld, intgroup = c("group"), returnData = TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

library(ggrepel)

# Create the PCA plot
p <- ggplot(data, aes(x = PC1, y = PC2, color = group, fill = group)) + 
  # Add points with larger size and transparency
  geom_point(size = 4, alpha = 0.8, stroke = 1.2, shape = 21) +
  
  # Add confidence ellipses (removed fill parameter to avoid warning)
  stat_ellipse(aes(color = group), type = "t", linetype = 2, 
               size = 0.8, alpha = 0.7) +
  
  # Add sample labels
  geom_text_repel(aes(label = name), size = 3, 
                  box.padding = 0.5, point.padding = 0.3) +
  
  # Enhanced color palette
  scale_color_brewer(type = "qual", palette = "Set2", name = "Group") +
  scale_fill_brewer(type = "qual", palette = "Set2", name = "Group") +
  
  # Add variance explained to axis labels
  labs(
    title = "Principal Component Analysis",
    subtitle = "ATAC-seq Data Visualization",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    caption = "Based on normalized ATAC-seq accessibility data"
  ) +
  
  # Enhanced theme
  theme_minimal(base_size = 12) +
  theme(
    # Plot appearance
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    plot.caption = element_text(size = 10, color = "gray50"),
    
    # Axes
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.line = element_line(color = "gray30", linewidth = 0.5),
    
    # Legend
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.5),
    legend.margin = margin(10, 10, 10, 10),
    
    # Panel
    panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.minor = element_line(color = "gray95", linewidth = 0.3),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    
    # Border
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.8)
  )

# Display the plot
print(p)

ggsave("PCA.pdf", p, width = 10, height = 8, device = "pdf")

# 4. DIFFERENTIAL ACCESSIBILITY ANALYSIS ####

# Create DESeq2 object
####customize the reference group####
sample_metadata$group <- factor(sample_metadata$group)  # convert to factor first
sample_metadata$group <- relevel(sample_metadata$group, ref = "RAS")  # set RAS as ref group

dds <- DESeqDataSetFromMatrix(
  countData = filtered_count_matrix,
  colData = sample_metadata,
  design = ~  group
)

# Run DESeq2 analysis
dds <- DESeq(dds)
resultsNames(dds)

# Get differential accessibility results
da_results <- results(dds, contrast = c("group", "ETOD7", "PROD2"), alpha= 0.05) # ETOD7 vs PROD2
da_results<- as.data.frame(da_results)
da_results$peak_id <- rownames(atac_count_data_matrix)

# Add peak coordinates
da_results$chr <- sapply(strsplit(da_results$peak_id, "_"), "[", 1)
da_results$start <- as.numeric(sapply(strsplit(da_results$peak_id, "_"), "[", 2))
da_results$end <- as.numeric(sapply(strsplit(da_results$peak_id, "_"), "[", 3))

# Filter significant peaks
sig_peaks <- subset(da_results, abs(log2FoldChange) > 1 & padj < 0.05)
up_peaks <- subset(sig_peaks, log2FoldChange > 0)
down_peaks <- subset(sig_peaks, log2FoldChange < 0)

cat("Significant peaks:", nrow(sig_peaks), "\n")
cat("Up-regulated peaks (ETOD7 vs PROD2):", nrow(up_peaks), "\n")
cat("Down-regulated peaks (ETOD7 vs PROD2):", nrow(down_peaks), "\n")

# Volcano plot
library(EnhancedVolcano)
volcano_plot <- EnhancedVolcano(da_results,
                                lab = rep("", nrow(da_results)),  
                                x = 'log2FoldChange',
                                y = 'padj',
                                title = 'Differential Accessibility Analysis (ETOD7 vs PROD2)',
                                subtitle = " ",
                                pCutoff = 0.05,
                                FCcutoff = 1.0,
                                pointSize = 1.0,
                                labSize = 0, 
                                col = c('grey30', 'forestgreen', 'royalblue', 'red2'), 
                                colAlpha = 0.7,
                                legendPosition = 'right',
                                legendLabSize = 12,
                                legendIconSize = 4.0,
                                drawConnectors = FALSE,  
                                xlim = c(min(da_results$log2FoldChange, na.rm = TRUE) * 1.1, 
                                         max(da_results$log2FoldChange, na.rm = TRUE) * 1.1),
                                ylim = c(0, max(-log10(da_results$padj), na.rm = TRUE) * 1.1),
                                gridlines.major = TRUE,
                                gridlines.minor = FALSE
)

print(volcano_plot)
ggsave("volcano_plot.pdf", 
       plot = volcano_plot,
       width = 12, height = 10, 
       device = "pdf",
       bg = "white")


# 5. PEAK ANNOTATION #####
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg38)

# Load transcript database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Create GRanges for significant peaks
sig_peak_ranges <- GRanges(
  seqnames = sig_peaks$chr,
  ranges = IRanges(start = sig_peaks$start, end = sig_peaks$end),
  peak_id = sig_peaks$peak_id,
  log2FoldChange = sig_peaks$log2FoldChange,
  padj = sig_peaks$padj
)


# Annotate peaks
peak_annotation <- annotatePeak(sig_peak_ranges, 
                                tssRegion = c(-3000, 3000),
                                TxDb = txdb,
                                annoDb = "org.Hs.eg.db")

# Plot annotation distribution
plotAnnoPie(peak_annotation)
plotAnnoBar(peak_annotation)
plotDistToTSS(peak_annotation, title = "Distribution of Peaks relative to TSS") #nearest transcription start site

# Get annotated data frame
annotated_peaks <- as.data.frame(peak_annotation)

#6. TF (TRANSCRIPTION FACTOR) BINDING SITE BIAS ANALYSIS####
library(GenomicRanges)
library(chromVAR)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)
library(devtools)
devtools::install_github("GreenleafLab/chromVARmotifs")
BiocManager::install("motifmatchr")
library(chromVARmotifs)
library(motifmatchr)

# Create GRanges object
peak_granges <- GRanges(
  seqnames = atac_count_data$Chr,
  ranges = IRanges(
    start = atac_count_data$Start,
    end = atac_count_data$End
  )
)

# Convert data.table to matrix 
count_matrix <- as.matrix(filtered_count_matrix)

# Create SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(counts = count_matrix),  # Use the converted matrix
  rowRanges = peak_granges,
  colData = sample_metadata
)

# Add GC Bias Correction
se <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)

# Match motifs to peaks
motif_ix <- matchMotifs(
  human_pwms_v2, 
  se,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

library(SummarizedExperiment)

# Compute TF Binding Deviations
dev <- computeDeviations(
  object = se,
  annotations = motif_ix
)

# Continue with your analysis
tf_z <- deviationScores(dev)
row_order <- order(apply(tf_z, 1, sd), decreasing = TRUE)
ordered_tfs <- tf_z[row_order, ]
head(ordered_tfs)

write.csv(ordered_tfs, file = 'TF_deviation_score_ordered_by_decreasing_sd.csv')

# Calculate variability across samples for each TF & Get top 30 variable TFs
top_tfs <- head(rownames(tf_z)[order(apply(tf_z, 1, sd), decreasing = TRUE)], 30)

###tf heatmap####
#Create heatmap
library(ComplexHeatmap)
library(circlize)
col_fun <- colorRamp2(
  c(-40, -20, 0, 20, 40), 
  c("#313695", "#74add1", "white", "#f46d43", "#a50026")
)

# Create heatmap with proper spacing for row names and legend
tf_heatmap <- Heatmap(
  tf_z[top_tfs, ],
  name = "Deviation score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45,
  row_names_max_width = unit(8, "cm"),  # Reserve space for long row names
  
  # Add grid lines by setting cell borders
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "black", lwd = 0.2, fill = NA))
  },
  
  # Configure legend with proper spacing
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    legend_height = unit(4, "cm"),
    legend_width = unit(1.2, "cm"),
    direction = "vertical",
    border = TRUE
  )
)

# Create PDF with proper dimensions and spacing
pdf("Heatmap_TF.pdf", width = 16, height = 10)

# Use proper layout with pushViewport to control spacing
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(0.75, 0.25), "npc"))))

# Draw heatmap in left viewport
pushViewport(viewport(layout.pos.col = 1))
draw(tf_heatmap, newpage = FALSE, show_heatmap_legend = FALSE)
popViewport()

# Draw legend in right viewport
pushViewport(viewport(layout.pos.col = 2))
lgd = Legend(
  col_fun = col_fun, 
  title = "Deviation score",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 9),
  legend_height = unit(4, "cm"),
  legend_width = unit(1.2, "cm"),
  direction = "vertical",
  border = TRUE
)
draw(lgd, x = unit(0.2, "npc"), y = unit(0.5, "npc"))
popViewport()

popViewport()
dev.off()

##tf motif enrichment plot.R##### run this code to get plots

# 7. FUNCTIONAL ENRICHMENT ANALYSIS####

# Get genes associated with significant peaks
sig_genes <- unique(annotated_peaks$SYMBOL[!is.na(annotated_peaks$SYMBOL)])

gene_ids <- bitr(sig_genes, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)

entrez_genes <- gene_ids$ENTREZID

# 1. GO enrichment analysis
go_enrichment <- enrichGO(gene = sig_genes,
                          OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)
dotplot(go_enrichment, showCategory = 20)
barplot(go_enrichment, showCategory = 15)


dotplot(go_enrichment, showCategory = 20) +
  geom_text(aes(label = Count), 
            color = "white", 
            size = 3, 
            fontface = "bold") +
  ggtitle("Top 20 GO Enrichment Terms- ETOD7 vs PROD2")

# 2. KEGG pathway analysis
kegg_enrichment <- enrichKEGG(gene = entrez_genes,
                              organism = "hsa",
                              pvalueCutoff = 0.05)
dotplot(kegg_enrichment, showCategory = 20)
barplot(kegg_enrichment, showCategory = 15)

# 3. Reactome pathway analysis
reactome_enrichment <- enrichPathway(gene = entrez_genes,
                                     pvalueCutoff = 0.05,
                                     readable = TRUE)
dotplot(reactome_enrichment, showCategory = 15)
barplot(reactome_enrichment, showCategory = 15)


#combine the results
pdf("GREAT_Results.pdf", width = 12, height = 8)
top_15 <- head(go_table[order(go_table$Binom_Adjp_BH), ], 15)

ggplot(top_15, aes(x = -log10(Binom_Adjp_BH), y = reorder(name, -log10(Binom_Adjp_BH)))) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  labs(title = "Top 15 GO Biological Process Terms",
       x = "-log10(Adjusted P-value)",
       y = "GO Terms") +
  theme_minimal()
plotTopTerms(great_result, ontology = "GO Biological Process", n_terms = 15, title = "Top GO Biological Process Terms")
plotRegionGeneAssociationGraphs(great_result)
dev.off()

# 8. SAVE RESULTS ####

# Save differential accessibility results
write.csv(da_results, "differential_accessibility_results.csv", row.names = FALSE)
write.csv(sig_peaks, "significant_peaks.csv", row.names = FALSE)
write.csv(annotated_peaks, "annotated_significant_peaks.csv", row.names = FALSE)

# Save normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, "normalized_counts.csv")

# Save plots
ggsave("pca_plot.pdf", pca_plot, width = 8, height = 6)
ggsave("volcano_plot.pdf", volcano_plot, width = 8, height = 6)


# 9. SESSION INFO & SUMMARY REPORT ####

sessionInfo()

# SUMMARY REPORT
cat("\n=== ATAC-seq Analysis Summary ===\n")
cat("Total peaks analyzed:", nrow(filtered_count_matrix), "\n")
cat("Significant peaks (|log2FC| > 1, padj < 0.05):", nrow(sig_peaks), "\n")
cat("Up-regulated peaks:", nrow(up_peaks), "\n")
cat("Down-regulated peaks:", nrow(down_peaks), "\n")
cat("Genes associated with significant peaks:", length(sig_genes), "\n")
cat("GO terms enriched:", nrow(go_enrichment@result[go_enrichment@result$qvalue < 0.05,]), "\n")
cat("KEGG pathways enriched:", nrow(kegg_enrichment@result[kegg_enrichment@result$qvalue < 0.05,]), "\n")
