################################################################################
#
# principle component analysis
#
################################################################################

### data management ------------------------------------------------------------
# load data
counts_files <- list.files(parameters["output_folder", "value"])
counts_files <- counts_files[grepl("^expression", counts_files)] # identify processed expression data
counts_names <- str_remove(counts_files, ".RDS")
counts_data <- list()
for(i in 1:length(counts_names)){
  counts_data[[counts_names[i]]] <- readRDS(paste0(parameters["output_folder", "value"], counts_files[i]))
  rm(i)
}
metadata <- readRDS(paste0(parameters["output_folder", "value"], "metadata.RDS"))
# clean up
rm(counts_files, counts_names) 

# output data
pca_out <- list()

### PCA pipeline ---------------------------------------------------------------
for(counts_ in names(counts_data)){
  counts_matrix_ <- counts_data[[counts_]] # call counts
  # run PCA
  pca_ <- prcomp(t(counts_matrix_))
  pca_out[[paste0(counts_, "_plot")]] <- data.frame(pca_$x, metadata) # export data in plottable format
  pca_out[[paste0(counts_, "_prcomp_output")]] <- pca_ # export prcomp output
  # plot PCA
  pca_var_ <- round(pca_$sdev^2/sum(pca_$sdev^2), 3) * 100; names(pca_var_) <- colnames(pca_$x) # calculate percent variance explained
  plot_df_ <- pca_out[[paste0(counts_, "_plot")]] # retrieve data frame for plot
  pc_1_2_ <- ggplot(plot_df_, aes(PC1, PC2, color = .data[[parameters["metadata_comparison_column", "value"]]])) +
    xlab(paste0("PC1 ", pca_var_[1], "%")) +
    ylab(paste0("PC2 ", pca_var_[2], "%")) +
    geom_point(alpha = 0.75, size = 2) +
    theme(aspect.ratio = 1)
  pc_3_4_ <- ggplot(plot_df_, aes(PC3, PC4, color = .data[[parameters["metadata_comparison_column", "value"]]])) +
    xlab(paste0("PC3 ", pca_var_[3], "%")) +
    ylab(paste0("PC4 ", pca_var_[4], "%")) +
    geom_point(alpha = 0.75, size = 2) +
    theme(aspect.ratio = 1)
  pc_5_6_ <- ggplot(plot_df_, aes(PC5, PC6, color = .data[[parameters["metadata_comparison_column", "value"]]])) +
    xlab(paste0("PC5 ", pca_var_[5], "%")) +
    ylab(paste0("PC6 ", pca_var_[6], "%")) +
    geom_point(alpha = 0.75, size = 2) +
    theme(aspect.ratio = 1)
  pc_plot_ <- plot_grid(pc_1_2_ + theme(legend.position = "none"),
                        pc_3_4_ + theme(legend.position = "none"),
                        pc_5_6_ + theme(legend.position = "none"), nrow = 1)
  legend_ <- get_legend(pc_1_2_ + theme(legend.box.margin = margin(0, 12, 0, 12)))
  pc_plot_ <- plot_grid(pc_plot_, legend_, rel_widths = c(3, .4))
  title_ <- ggdraw() +
    draw_label(paste0("PCA on ", counts_), fontface = 'bold', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  jpeg(paste0(parameters["output_folder", "value"], "PCA_on_", counts_, ".jpg"), width = 1200, height = 400)
  print(plot_grid(title_, pc_plot_, ncol = 1, rel_heights = c(0.1, 1)))
  dev.off()
  checkpoint <- rbind(checkpoint, c(paste0("PCA_on_", counts_, ".jpg"), format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))); write.csv(checkpoint, checkpoint_file, row.names = F)
  rm(list = ls()[grepl("_$", ls())])
}

### export results -------------------------------------------------------------
saveRDS(pca_out, paste0(parameters["output_folder", "value"], "PCA_results.RDS"))
checkpoint <- rbind(checkpoint, c("PCA_results.RDS", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))); write.csv(checkpoint, checkpoint_file, row.names = F)

rm(pca_out, counts_data, metadata)
