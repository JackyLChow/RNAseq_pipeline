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
library(ggrepel)

# bioniformatics
## differential gene expression
library(limma)
## gene expression visualization
library(ComplexHeatmap)
library(circlize)
## gene set enrichment
library(org.Hs.eg.db)
library(msigdbr)
library(ReactomePA)
library(clusterProfiler)

###---data import---------------------------------------------------------------
# parameters file
parameters <- data.frame(read_csv("~/Documents/BFX_proj/RNAseq_pipeline/_input/Chow_PNAS_2020.csv"))
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

# scripts
DGE_calculation <- "~/Documents/BFX_proj/RNAseq_pipeline/R_code/limma_dge_v1.R" # differential gene expression and gsea

###---DGE-----------------------------------------------------------------------
if(any(grepl("^differential gene expression", checkpoint$file))){
  cat("\nDGE results already exist\n")
} else {
  cat("\ncalculating DGE and GSEA\n")
  source(DGE_calculation)
  cat("... complete ...")
}
