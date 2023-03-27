################################################################################
#
# counts distribution QC
#
################################################################################

### data management ------------------------------------------------------------
# load data
counts_matrix_ <- readRDS(parameters["counts_file", "value"])
counts_ <- parameters["counts_type", "value"]

### counts distribution --------------------------------------------------------
counts_dataframe_ <- data.frame()
for(i in colnames(counts_matrix_)){
  sample_data_ <- data.frame(sample = rep(i, nrow(counts_matrix_)),
                             counts = counts_matrix_[, i],
                             row.names = NULL)
  counts_dataframe_ <- rbind(counts_dataframe_, sample_data_)
  rm(i)
}
# plot counts boxplot
plot_counts_ <- ggplot(counts_dataframe_, aes(sample, counts)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

title_ <- ggdraw() +
  draw_label(paste0("Sample counts distribution ", counts_), fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
png(paste0(parameters["output_folder", "value"], "sample_counts_distribution_on_", counts_, ".png"), width = ncol(counts_matrix_)*10, height = 400)
print(plot_grid(title_, plot_counts_, ncol = 1, rel_heights = c(0.1, 1)))
dev.off()
checkpoint <- rbind(checkpoint, c(paste0("sample_counts_distribution_on_", counts_, ".png"), format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))); write.csv(checkpoint, checkpoint_file, row.names = F)
rm(list = ls()[grepl("_$", ls())])

### export results -------------------------------------------------------------
checkpoint <- rbind(checkpoint, c("sample_counts_distribution", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))); write.csv(checkpoint, checkpoint_file, row.names = F)