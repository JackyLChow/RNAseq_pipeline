# "All purpose" tool for initial heatmap generation
## something quick, rote, but also tidy and visually appealing

library(matrixStats)

################################################################################
#
# prep data for heatmap
#
################################################################################

### counts ---
h_mtx <- c_mtx[order(rowVars(as.matrix(c_mtx)), decreasing = T)[1:500], ]
h_mtx <- t(scale(t(h_mtx)))

### heatmap annotation ---
# heatmap colors
col_htmp <- colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow2"))

################################################################################
#
# bake the heatmap
#
################################################################################

set.seed(415); hm <- Heatmap(h_mtx,
                             col = col_htmp,
                             name = "z-score",
                             #height = unit(2.5, "mm") * nrow(h_mtx), width = unit(2.5, "mm") * ncol(h_mtx),
                             show_row_dend = T, show_column_dend = T,
                             show_row_names = F, show_column_names = F,
                             row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
                             top_annotation = anno_col,
                             show_heatmap_legend = T)


