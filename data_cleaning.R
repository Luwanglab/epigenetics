#ATAC seq data cleanning
library(readr)
data<- fread("~/Documents/CZ/CZ_ras_rna/data/raw data/atac_peak_counts.txt") #n =359662
data<- as.data.frame(data)
data_names<- read_csv("Desktop/sample_info.csv")

#keep expected columns
cols <- colnames(data)
keep_cols <- cols[grepl("Chr|Start|End|CZ|NEX", cols)]
exclude <- c("NEX_1", "NEX_2", "NEX_9", "NEX_10", "NEX_14")
file_names <- basename(keep_cols)
exclude_pattern <- paste0("^(", paste(exclude, collapse = "|"), ")_")
keep_final_cols <- keep_cols[!grepl(exclude_pattern, file_names)]

data_filtered <- data[, keep_final_cols]

# create new column names
new_names <- basename(colnames(data_filtered))               
new_names <- sub("_S\\d+.*$", "", new_names)                  
colnames(data_filtered) <- new_names

#change to meaningful names

data_names$Sample
data_names$Indicator

old_names <- colnames(data_filtered)
name_map <- setNames(data_names$Indicator, data_names$Sample)
new_names <- ifelse(old_names %in% names(name_map), name_map[old_names], old_names)
colnames(data_filtered) <- new_names

colnames(data_filtered) <- gsub("^IMR90_", "", colnames(data_filtered))
colnames(data_filtered) <- gsub("knockdown", "KD", colnames(data_filtered))

###Data cleaning####
# use edgeR filterByExpr function
library(edgeR)
non_gene_data <- data_filtered[, -(1:3)]
# create group info
group <- factor(c(
  rep("RAS", 2),
  rep("MafK45_RAS", 2),
  rep("Proliferating", 2),
  rep("Ras_MafK045_KD", 2),
  rep("Ras_MafK138_KD", 2),
  "MafK045_KD_Ras",
  "MafK138_KD_Ras"
))


# Use filterByExpr to filter low-signal chromatin accessibility regions
# For ATAC-seq data, this function will:
# 1. Identify peak regions with weak signals across all conditions  
# 2. Retain regions with sufficient read coverage in at least one condition
# Parameters:
# - non_gene_data: peak region count matrix (rows = peaks, columns = samples)
# - group: sample grouping to determine which peaks are active under specific conditions
# - min.count: minimum read count threshold (50), peaks must reach this to be considered "accessible"
keep <- filterByExpr(non_gene_data, group = group, min.count = 50)
data_filter <- data_filtered[keep, , drop = FALSE]

setwd("~/Desktop/CZ")
write_csv(data_filter, file= 'ATAC_data_CZ.csv') # n= 151670

