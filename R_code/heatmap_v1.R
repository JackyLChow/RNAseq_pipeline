################################################################################
#
# Heatmap generation
#
################################################################################

### data management ------------------------------------------------------------
# load data
counts_files <- list.files(parameters["output_folder", "value"])
counts_files <- counts_files[grepl("^expression", counts_files)] # identify processed expression data
counts_names <- str_remove(counts_files, ".RDS")
counts_data <- list()
for(i in 1:length(counts_names)){
  counts_data[[counts_names[i]]] <- readRDS(paste0("~/Documents/BFX_proj/RNAseq_pipeline/_output/Chow_PNAS_2020/", counts_files[i]))
  rm(i)
}
metadata <- readRDS("~/Documents/BFX_proj/RNAseq_pipeline/_output/Chow_PNAS_2020/metadata.RDS")
# clean up
rm(counts_files, counts_names) 

# output data
heatmap_out <- list()

### heatmap parameters ---------------------------------------------------------
heatmap_expression_colors_ <- colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow")) # heatmap colors

### heatmap pipeline -----------------------------------------------------------
for(counts_ in names(counts_data)){
  # counts mutation
  counts_matrix_ <- counts_data[[counts_]] # call counts
  counts_matrix_ <- counts_matrix_[order(rowVars(as.matrix(counts_matrix_)), decreasing = T)[1:500], ] # filter top 500 variant genes
  counts_matrix_ <- t(scale(t(counts_matrix_))) # z-score each gene across samples
  
  # annotation
  if(is.na(parameters["metadata_batch_column", "value"])){
    top_anno_ <- HeatmapAnnotation(comparison = metadata[, parameters["metadata_comparison_column", "value"]])
  } else {
    top_anno_ <- HeatmapAnnotation(comparison = metadata[, parameters["metadata_comparison_column", "value"]],
                                   batch = metadata[, parameters["metadata_batch_column", "value"]])
  }
  
  # generate clustered heatmap
  set.seed(415); heatmap_ <- Heatmap(counts_matrix_,
                                     column_title = paste0("Heatmap on ", counts_),
                                     col = heatmap_expression_colors_,
                                     name = "z-score",
                                     #height = unit(2.5, "mm") * nrow(h_mtx),
                                     #width = unit(2.5, "mm") * ncol(h_mtx),
                                     show_row_dend = T, show_column_dend = T,
                                     show_row_names = F, show_column_names = F,
                                     row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
                                     top_annotation = top_anno_,
                                     show_heatmap_legend = T)
  
  jpeg(paste0(parameters["output_folder", "value"], "heatmap_on_", counts_, ".jpg"), width = 800, height = 800)
  print(heatmap_)
  dev.off()
  checkpoint <- rbind(checkpoint, c(paste0("heatmap_on_", counts_, ".jpg"), format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))); write.csv(checkpoint, checkpoint_file, row.names = F)
  
  heatmap_out[[paste0(counts_, "_heatmap_output")]] <- draw(heatmap_) # export heatmap parameters
}

rm(list = ls()[grepl("_$", ls())])

### export results -------------------------------------------------------------
saveRDS(heatmap_out, paste0(parameters["output_folder", "value"], "heatmap_results.RDS"))
checkpoint <- rbind(checkpoint, c("heatmap_results.RDS", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))); write.csv(checkpoint, checkpoint_file, row.names = F)

rm(heatmap_out, metadata, counts_data) # clean up



                             


