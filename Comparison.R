# Comprehensive ATAC-seq Analysis for Three Comparisons
# Libraries
library(DESeq2)
library(EnhancedVolcano)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(igraph)
library(ggraph)
library(dplyr)
library(gprofiler2)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(JASPAR2020)
library(motifmatchr)


# Initial DESeq2 setup
dds <- DESeqDataSetFromMatrix(
  countData = filtered_count_matrix,
  colData = sample_metadata,
  design = ~ group
)

dds <- DESeq(dds)

# Define comparisons
comparisons <- list(
  list(name = "RAS_vs_Proliferating", contrast = c("group", "RAS", "Proliferating")),
  list(name = "MafK45_RAS_vs_RAS", contrast = c("group", "MafK45_RAS", "RAS")),
  list(name = "Ras_MafK045_KD_vs_RAS", contrast = c("group", "Ras_MafK045_KD", "RAS")),
  list(name = "Ras_MafK138_KD_vs_RAS", contrast = c("group", "Ras_MafK138_KD", "RAS")),
  list(name = "MafK045_KD_Ras_vs_RAS", contrast = c("group", "MafK045_KD_Ras", "RAS")),
  list(name = "MafK138_KD_Ras_vs_RAS", contrast = c("group", "MafK138_KD_Ras", "RAS"))
)


# Initialize summary table
summary_table <- data.frame(
  Comparison = character(),
  Significant_Peaks = numeric(),
  Up_regulated_Peaks = numeric(),
  Down_regulated_Peaks = numeric(),
  Associated_Genes = numeric(),
  GO_Terms_Enriched = numeric(),
  KEGG_Pathways_Enriched = numeric(),
  Reactome_Pathways_Enriched = numeric(),
  TF_Targets_Enriched = numeric(),
  stringsAsFactors = FALSE
)

# Load transcript database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Function to create comprehensive analysis for each comparison
analyze_comparison <- function(dds, comparison, atac_count_data_matrix) {
  comp_name <- comparison$name
  cat(paste("\n=== Analyzing", comp_name, "===\n"))
  
  # Create directories
  dir.create(comp_name, showWarnings = FALSE)
  
  # 1. GET DIFFERENTIAL ACCESSIBILITY RESULTS
  da_results <- results(dds, contrast = comparison$contrast, alpha = 0.05)
  da_results <- as.data.frame(da_results)
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
  cat("Up-regulated peaks:", nrow(up_peaks), "\n")
  cat("Down-regulated peaks:", nrow(down_peaks), "\n")
  
  # 2. CREATE VOLCANO PLOT
  volcano_plot <- EnhancedVolcano(da_results,
                                  lab = rep("", nrow(da_results)),
                                  x = 'log2FoldChange',
                                  y = 'padj',
                                  title = paste('Differential Accessibility Analysis', comp_name),
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
                                  gridlines.minor = FALSE)
  
  # 3. PEAK ANNOTATION
  if(nrow(sig_peaks) > 0) {
    sig_peak_ranges <- GRanges(
      seqnames = sig_peaks$chr,
      ranges = IRanges(start = sig_peaks$start, end = sig_peaks$end),
      peak_id = sig_peaks$peak_id,
      log2FoldChange = sig_peaks$log2FoldChange,
      padj = sig_peaks$padj
    )
    
    peak_annotation <- annotatePeak(sig_peak_ranges,
                                    tssRegion = c(-3000, 3000),
                                    TxDb = txdb,
                                    annoDb = "org.Hs.eg.db")
    
    annotated_peaks <- as.data.frame(peak_annotation)
    
    # Save annotation plots separately
    # Annotation pie chart
    pdf(file.path(comp_name, paste0(comp_name, "_annotation_pie.pdf")), width = 8, height = 6)
    plotAnnoPie(peak_annotation)
    dev.off()
    
    # Annotation bar chart
    pdf(file.path(comp_name, paste0(comp_name, "_annotation_bar.pdf")), width = 10, height = 6)
    plotAnnoBar(peak_annotation)
    dev.off()
    
    # Distance to TSS plot
    pdf(file.path(comp_name, paste0(comp_name, "_distance_to_TSS.pdf")), width = 10, height = 6)
    plotDistToTSS(peak_annotation, title = paste("Distribution of Peaks relative to TSS -", comp_name))
    dev.off()
    
    # 4. FUNCTIONAL ENRICHMENT ANALYSIS
    sig_genes <- unique(annotated_peaks$SYMBOL[!is.na(annotated_peaks$SYMBOL)])
    
    if(length(sig_genes) > 5) {
      # GO enrichment
      go_enrichment <- enrichGO(gene = sig_genes,
                                OrgDb = org.Hs.eg.db,
                                keyType = "SYMBOL",
                                ont = "BP",
                                pAdjustMethod = "BH",
                                qvalueCutoff = 0.05)
      
      # KEGG enrichment
      gene_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      if(nrow(gene_ids) > 5) {
        entrez_genes <- gene_ids$ENTREZID
        
        kegg_enrichment <- enrichKEGG(gene = entrez_genes,
                                      organism = "hsa",
                                      pvalueCutoff = 0.05)
        
        # Reactome enrichment
        reactome_enrichment <- enrichPathway(gene = entrez_genes,
                                             pvalueCutoff = 0.05,
                                             readable = TRUE)
      }
      
      # 5. TRANSCRIPTION FACTOR ENRICHMENT ANALYSIS
      
      # Method 1: Using gProfiler2 for TF enrichment
      gost_result <- NULL
      tf_results <- data.frame()
      
      tryCatch({
        gost_result <- gost(query = sig_genes,
                            organism = "hsapiens",
                            sources = c("TF", "REAC", "GO:BP"),
                            correction_method = "gSCS")
        
        if(!is.null(gost_result) && !is.null(gost_result$result)) {
          tf_results <- gost_result$result[gost_result$result$source == "TF", ]
          
          # Flatten list columns for CSV export
          if(nrow(tf_results) > 0) {
            # Convert list columns to character strings
            if("intersection" %in% colnames(tf_results)) {
              tf_results$intersection <- sapply(tf_results$intersection, function(x) {
                tryCatch({
                  if(is.null(x)) return("")
                  if(is.list(x)) return(paste(unlist(x), collapse = ","))
                  if(is.character(x)) return(paste(x, collapse = ","))
                  return(as.character(x))
                }, error = function(e) return(""))
              })
            }
            if("parents" %in% colnames(tf_results)) {
              tf_results$parents <- sapply(tf_results$parents, function(x) {
                tryCatch({
                  if(is.null(x)) return("")
                  if(is.list(x)) return(paste(unlist(x), collapse = ","))
                  if(is.character(x)) return(paste(x, collapse = ","))
                  return(as.character(x))
                }, error = function(e) return(""))
              })
            }
          }
        }
      }, error = function(e) {
        cat("Warning: TF enrichment analysis failed:", e$message, "\n")
        tf_results <<- data.frame()
      })
      
      # Method 2: Motif enrichment using chromVAR (if you have motif data)
      # This requires additional setup with motif databases
      
      # 6. CREATE GENE NETWORK PLOTS
      
      # Function to create gene network from enrichment results
      create_gene_network <- function(enrichment_result, title_prefix, max_terms = 10) {
        if(is.null(enrichment_result) || nrow(enrichment_result@result) == 0) {
          return(NULL)
        }
        
        # Get top terms
        top_terms <- head(enrichment_result@result[enrichment_result@result$qvalue < 0.05, ], max_terms)
        
        if(nrow(top_terms) == 0) {
          return(NULL)
        }
        
        # Create network data
        edges <- data.frame()
        nodes <- data.frame()
        
        for(i in 1:nrow(top_terms)) {
          term_genes <- unlist(strsplit(top_terms$geneID[i], "/"))
          term_name <- top_terms$Description[i]
          
          # For KEGG pathways, convert ENTREZ IDs to gene symbols
          if(title_prefix == "KEGG Pathway") {
            # Convert ENTREZ IDs to gene symbols
            tryCatch({
              gene_symbols <- mapIds(org.Hs.eg.db, keys = term_genes, 
                                     column = "SYMBOL", keytype = "ENTREZID", 
                                     multiVals = "first")
              term_genes <- gene_symbols[!is.na(gene_symbols)]
            }, error = function(e) {
              cat("Warning: Could not convert ENTREZ IDs to symbols for", title_prefix, "\n")
            })
          }
          
          # Add term node
          if(!term_name %in% nodes$name) {
            nodes <- rbind(nodes, data.frame(name = term_name, type = "term", 
                                             size = -log10(top_terms$qvalue[i])))
          }
          
          # Add gene nodes and edges
          for(gene in term_genes) {
            if(!is.na(gene) && gene != "" && !gene %in% nodes$name) {
              nodes <- rbind(nodes, data.frame(name = gene, type = "gene", size = 1))
            }
            if(!is.na(gene) && gene != "") {
              edges <- rbind(edges, data.frame(from = term_name, to = gene))
            }
          }
        }
        
        # Only create network if we have valid edges
        if(nrow(edges) == 0) {
          return(NULL)
        }
        
        # Create network plot
        g <- graph_from_data_frame(edges, vertices = nodes)
        
        network_plot <- ggraph(g, layout = "fr") +
          geom_edge_link(alpha = 0.3, color = "grey70") +
          geom_node_point(aes(color = type, size = size)) +
          geom_node_text(aes(label = ifelse(type == "gene", name, "")), 
                         size = 2, repel = TRUE, max.overlaps = 20) +
          scale_color_manual(values = c("term" = "red", "gene" = "blue")) +
          scale_size_continuous(range = c(1, 6)) +
          theme_void() +
          labs(title = paste(title_prefix, "Gene Network -", comp_name)) +
          theme(legend.position = "bottom")
        
        return(network_plot)
      }
      
      # Create network plots
      go_network <- create_gene_network(go_enrichment, "GO Biological Process")
      kegg_network <- if(exists("kegg_enrichment")) create_gene_network(kegg_enrichment, "KEGG Pathway") else NULL
      reactome_network <- if(exists("reactome_enrichment")) create_gene_network(reactome_enrichment, "Reactome Pathway") else NULL
      
      # 7. TF ENRICHMENT PLOT (consistent with other enrichment plots)
      tf_plot <- NULL
      if(nrow(tf_results) > 0) {
        tryCatch({
          top_tf <- head(tf_results[order(tf_results$p_value), ], 15)
          
          # Safely calculate gene counts
          gene_counts <- sapply(top_tf$intersection, function(x) {
            tryCatch({
              if(is.character(x) && x != "" && !is.na(x)) {
                length(unlist(strsplit(as.character(x), ",")))
              } else {
                1  # Default count if parsing fails
              }
            }, error = function(e) 1)
          })
          
          # Create a data frame similar to enrichment results for consistency
          tf_plot_data <- data.frame(
            Description = top_tf$term_name,
            Count = gene_counts,
            pvalue = top_tf$p_value,
            p.adjust = top_tf$p_value,
            qvalue = top_tf$p_value,
            GeneRatio = paste0(gene_counts, "/", length(sig_genes))
          )
          
          # Create dotplot style visualization for TF enrichment
          tf_plot <- ggplot(tf_plot_data, aes(x = Count, y = reorder(Description, Count))) +
            geom_point(aes(size = Count, color = -log10(pvalue))) +
            geom_text(aes(label = Count), color = "white", size = 3, fontface = "bold") +
            scale_color_gradient(low = "blue", high = "red", name = "-log10(pvalue)") +
            scale_size_continuous(range = c(3, 10), name = "Count") +
            labs(title = paste("Top 15 Transcription Factor Targets -", comp_name),
                 x = "Gene Count",
                 y = "TF Targets") +
            theme_minimal() +
            theme(axis.text.y = element_text(size = 10))
        }, error = function(e) {
          cat("Warning: TF plot creation failed:", e$message, "\n")
          tf_plot <- NULL
        })
      }
      
      # 8. COMBINE MAIN PLOTS INTO ONE PDF (excluding annotation plots)
      pdf(file.path(comp_name, paste0(comp_name, "_comprehensive_analysis.pdf")), 
          width = 16, height = 12)
      
      # Page 1: Volcano plot
      print(volcano_plot)
      
      # Page 2-4: Enrichment plots with gene counts
      if(!is.null(go_enrichment) && nrow(go_enrichment@result) > 0) {
        tryCatch({
          go_dot <- dotplot(go_enrichment, showCategory = 20) +
            geom_text(aes(label = Count), color = "white", size = 3, fontface = "bold") +
            ggtitle(paste("GO Enrichment -", comp_name))
          print(go_dot)
        }, error = function(e) {
          cat("Error creating GO dotplot:", e$message, "\n")
        })
      }
      
      if(exists("kegg_enrichment") && !is.null(kegg_enrichment) && nrow(kegg_enrichment@result) > 0) {
        tryCatch({
          kegg_dot <- dotplot(kegg_enrichment, showCategory = 15) +
            geom_text(aes(label = Count), color = "white", size = 3, fontface = "bold") +
            ggtitle(paste("KEGG Enrichment -", comp_name))
          print(kegg_dot)
        }, error = function(e) {
          cat("Error creating KEGG dotplot:", e$message, "\n")
        })
      }
      
      if(exists("reactome_enrichment") && !is.null(reactome_enrichment) && nrow(reactome_enrichment@result) > 0) {
        tryCatch({
          reactome_dot <- dotplot(reactome_enrichment, showCategory = 15) +
            geom_text(aes(label = Count), color = "white", size = 3, fontface = "bold") +
            ggtitle(paste("Reactome Enrichment -", comp_name))
          print(reactome_dot)
        }, error = function(e) {
          cat("Error creating Reactome dotplot:", e$message, "\n")
        })
      }
      
      # Page 5: TF enrichment
      if(!is.null(tf_plot)) {
        tryCatch({
          print(tf_plot)
        }, error = function(e) {
          cat("Error creating TF plot:", e$message, "\n")
        })
      }
      
      # Page 6-8: Network plots
      if(!is.null(go_network)) {
        tryCatch({
          print(go_network)
        }, error = function(e) {
          cat("Error creating GO network:", e$message, "\n")
        })
      }
      
      if(!is.null(kegg_network)) {
        tryCatch({
          print(kegg_network)
        }, error = function(e) {
          cat("Error creating KEGG network:", e$message, "\n")
        })
      }
      
      if(!is.null(reactome_network)) {
        tryCatch({
          print(reactome_network)
        }, error = function(e) {
          cat("Error creating Reactome network:", e$message, "\n")
        })
      }
      
      dev.off()
      
      # 9. SAVE RESULTS TO FILES
      write.csv(da_results, file.path(comp_name, paste0(comp_name, "_differential_accessibility_results.csv")), row.names = FALSE)
      write.csv(sig_peaks, file.path(comp_name, paste0(comp_name, "_significant_peaks.csv")), row.names = FALSE)
      write.csv(annotated_peaks, file.path(comp_name, paste0(comp_name, "_annotated_significant_peaks.csv")), row.names = FALSE)
      
      # Save enrichment results with error handling
      tryCatch({
        if(!is.null(go_enrichment) && nrow(go_enrichment@result) > 0) {
          write.csv(go_enrichment@result, file.path(comp_name, paste0(comp_name, "_GO_enrichment.csv")), row.names = FALSE)
        }
      }, error = function(e) cat("Error saving GO results:", e$message, "\n"))
      
      tryCatch({
        if(exists("kegg_enrichment") && !is.null(kegg_enrichment) && nrow(kegg_enrichment@result) > 0) {
          write.csv(kegg_enrichment@result, file.path(comp_name, paste0(comp_name, "_KEGG_enrichment.csv")), row.names = FALSE)
        }
      }, error = function(e) cat("Error saving KEGG results:", e$message, "\n"))
      
      tryCatch({
        if(exists("reactome_enrichment") && !is.null(reactome_enrichment) && nrow(reactome_enrichment@result) > 0) {
          write.csv(reactome_enrichment@result, file.path(comp_name, paste0(comp_name, "_Reactome_enrichment.csv")), row.names = FALSE)
        }
      }, error = function(e) cat("Error saving Reactome results:", e$message, "\n"))
      
      tryCatch({
        if(nrow(tf_results) > 0) {
          write.csv(tf_results, file.path(comp_name, paste0(comp_name, "_TF_enrichment.csv")), row.names = FALSE)
        }
      }, error = function(e) cat("Error saving TF results:", e$message, "\n"))
      
      # Summary report with data collection
      n_go_enriched <- if(!is.null(go_enrichment) && nrow(go_enrichment@result) > 0) {
        nrow(go_enrichment@result[go_enrichment@result$qvalue < 0.05,])
      } else { 0 }
      
      n_kegg_enriched <- if(exists("kegg_enrichment") && !is.null(kegg_enrichment) && nrow(kegg_enrichment@result) > 0) {
        nrow(kegg_enrichment@result[kegg_enrichment@result$qvalue < 0.05,])
      } else { 0 }
      
      n_reactome_enriched <- if(exists("reactome_enrichment") && !is.null(reactome_enrichment) && nrow(reactome_enrichment@result) > 0) {
        nrow(reactome_enrichment@result[reactome_enrichment@result$qvalue < 0.05,])
      } else { 0 }
      
      n_tf_enriched <- if(nrow(tf_results) > 0) {
        nrow(tf_results[tf_results$p_value < 0.05,])
      } else { 0 }
      
      # Add to summary table
      summary_table <<- rbind(summary_table, data.frame(
        Comparison = comp_name,
        Significant_Peaks = nrow(sig_peaks),
        Up_regulated_Peaks = nrow(up_peaks),
        Down_regulated_Peaks = nrow(down_peaks),
        Associated_Genes = length(sig_genes),
        GO_Terms_Enriched = n_go_enriched,
        KEGG_Pathways_Enriched = n_kegg_enriched,
        Reactome_Pathways_Enriched = n_reactome_enriched,
        TF_Targets_Enriched = n_tf_enriched,
        stringsAsFactors = FALSE
      ))
      
      cat(paste("\n=== Summary for", comp_name, "===\n"))
      cat("Significant peaks:", nrow(sig_peaks), "\n")
      cat("Up-regulated peaks:", nrow(up_peaks), "\n")
      cat("Down-regulated peaks:", nrow(down_peaks), "\n")
      cat("Associated genes:", length(sig_genes), "\n")
      cat("GO terms enriched:", n_go_enriched, "\n")
      cat("KEGG pathways enriched:", n_kegg_enriched, "\n")
      cat("Reactome pathways enriched:", n_reactome_enriched, "\n")
      cat("TF targets enriched:", n_tf_enriched, "\n")
      
    } else {
      cat("Not enough genes for enrichment analysis\n")
    }
    
  } else {
    cat("No significant peaks found\n")
  }
}

# Run analysis for all comparisons
for(comparison in comparisons) {
  analyze_comparison(dds, comparison, atac_count_data_matrix)
}



# Create and save comprehensive summary table
cat("\n=== COMPREHENSIVE SUMMARY TABLE ===\n")
print(summary_table)

# Save summary table
write.csv(summary_table, "ATAC_seq_analysis_summary_table.csv", row.names = FALSE)

cat("\n=== COMPLETE ATAC-seq ANALYSIS FINISHED ===\n")
cat("Results saved in separate directories for each comparison:\n")
for(comparison in comparisons) {
  cat(paste("-", comparison$name, "\n"))
}