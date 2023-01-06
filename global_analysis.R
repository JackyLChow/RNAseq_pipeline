###---load packages-------------------------------------------------------------
# global
## data manipulation
library(readr)
library(dplyr)
library(tidyr)
library(spatstat)
library(matrixStats)
library(stringr)
## plotting
library(ggplot2)
library(cowplot)

# bioniformatics
## counts normalization, batch correction, and differential gene expression
library(limma)
library(sva)
## gene expression visualization
library(ComplexHeatmap)
library(circlize)
## weighted genetic correlation network analysis
library(WGCNA)
## gene set enrichment
library(org.Hs.eg.db)
library(ReactomePA)
library(clusterProfiler)

###---data import---------------------------------------------------------------
# parameters file
parameters <- data.frame(read_csv(""))
rownames(parameters) <- parameters[, "description"]

# checkpoints
checkpoint_file <- parameters["checkpoint_file", "value"]
if("checkpoint.csv" %in% list.files(parameters["output_folder", "value"])){
  checkpoint <- data.frame(read_csv(checkpoint_file))
} else {
  checkpoint <- data.frame(matrix(ncol = 2, nrow = 0))
  names(checkpoint) <- c("file", "date_created")
  checkpoint[1, ] <- c("checkpoint.csv", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  write.csv(checkpoint, checkpoint_file, row.names = F)
}
#write.table(data.frame("lol", "lmao"), file = checkpoint_file, sep = ",", append = T, col.names = F, row.names = F)

# scripts
normalization_and_batch_correction <- "~/Documents/BFX_proj/RNAseq_pipeline/R_code/data_load_clean_v1.R" # counts pre-processing
counts_distribution <- "~/Documents/BFX_proj/RNAseq_pipeline/R_code/counts_distribution_v1.R" # QC counts by sample
PCA_calculation <- "~/Documents/BFX_proj/RNAseq_pipeline/R_code/pca_v1.R" # PCA pipeline
heatmap_generation <- "~/Documents/BFX_proj/RNAseq_pipeline/R_code/heatmap_v1.R" # heatmap pipeline
WGCNA_calculation <- "~/Documents/BFX_proj/RNAseq_pipeline/R_code/wgcna_v1.R" # heatmap pipeline

###---raw counts processing-----------------------------------------------------
if(any(grepl("^expression", checkpoint$file))){
  cat("\nnormalized counts already exist\n")
  } else {
    cat("\nnormalizing counts\n")
    source(normalization_and_batch_correction)
    if(is.na(parameters["metadata_batch_column", "value"])){
      cat("\nno batch correction applied\n")
    } else {
      cat(paste0("\nbatch correction in ComBat by ", parameters["metadata_batch_column", "value"],"\n"))
    }
    cat("... complete ...")
  }

###---counts_distribution-------------------------------------------------------
if(any(grepl("sample_counts_distribution", checkpoint$file))){
  cat("\ncounts distribution results already exist\n")
} else {
  cat("plotting counts distribution")
  source(counts_distribution)
  cat("... complete ...")
}

###---PCA-----------------------------------------------------------------------
if(all(c("PCA_results.RDS") %in% checkpoint$file)){
  cat("\nPCA results already exist\n")
} else {
  cat(paste0("\ncalculating PCA and plotting by ",
             parameters["metadata_comparison_column", "value"], "\n"))
  source(PCA_calculation)
  cat("... complete ...")
}

###---heatmap-------------------------------------------------------------------
if(all(c("heatmap_results.RDS") %in% checkpoint$file)){
  cat("\nheatmaps already exist\n")
} else {
  cat(paste0("\ngenerating heatmap and annotating by:\ncomparison = ",
             parameters["metadata_comparison_column", "value"], "\n"))
  source(heatmap_generation)
  cat("... complete ...")
}

###---WGCNA---------------------------------------------------------------------
if(any(grepl("^WGCNA", checkpoint$file))){
  cat("\nWGCNA results already exist\n")
} else {
  cat("\ncalculating WGCNA\n")
  source(WGCNA_calculation)
  cat("... complete ...")
}










