# Rote analysis of bulkRNAseq data
library(readr)
library(stringr)
library(circlize)
library(pals)
library(ComplexHeatmap)
library(ggplot2)

################################################################################
#
# load and shape input data
#
################################################################################

### directories ---
in_dir <- "~/Documents/BFX_proj/RNAseq_pipeline/_input/Chow_PNAS_2020/"
out_dir <- "~/Documents/BFX_proj/RNAseq_pipeline/_output/Chow_PNAS_2020/"
### normalized counts ---
n_count_ <- data.frame(read_csv(paste0(in_dir, "Chow_PNAS_normcounts.csv")))
c_mtx <- n_count_[, 2:ncol(n_count_)]
rownames(c_mtx) <- n_count_$gene
### metadata ---
meta_ <- data.frame(read_csv(paste0(in_dir, "Chow_PNAS_meta_med.csv")))
rownames(meta_) <- str_replace(meta_$sample, "-", ".")
meta_ <- meta_[colnames(c_mtx), ]

### parameters for analysis ---
# WGCNA
dum <- c("age", "F", "M", "control", "SBRT") # dummified columns for WGCNA
meta_dum <- meta_[, dum]
# PCA, heatmaps
srt <- c("age", "sex", "treatment", "path_T_stage") # "smart" columns for plots
meta_srt <- meta_[, c("age", "sex", "treatment", "path_T_stage")]
# color legend
col_pts_ <- glasbey(length(unique(meta_srt$path_T_stage))); names(col_pts_) <- unique(meta_srt$path_T_stage)

col_leg <- list(
  age = colorRamp2(c(min(meta_srt$age), max(meta_srt$age)), c("green1", "green4")),
  sex = c("F" = "red3", "M" = "blue3"),
  treatment = c("SBRT" = "hotpink", "control" = "dodgerblue"),
  path_T_stage = col_pts_)

anno_col <- HeatmapAnnotation(age = meta_srt$age,
                              sex = meta_srt$sex,
                              treatment = meta_srt$treatment,
                              path_T_stage = meta_srt$path_T_stage,
                              col = col_leg,
                              show_legend = F,
                              annotation_name_side = "left",
                              annotation_name_gp= gpar(fontsize = 8),
                              simple_anno_size = unit(0.25, "cm"))

rm(list = ls()[grepl("_$", ls())]) # clean up

################################################################################
#
# summarize data
#
################################################################################

################################################################################
# color legend
################################################################################

source("~/Documents/BFX_proj/Miniport/R_code/color_legend.R")

png(paste0(out_dir, "fig/legend.png"), height = 60 * (length(lgd) + 1), width = 30 * max(lengths(col_leg)), res = 100)
draw(packLegend(list = lgd, direction = "vertical"))
dev.off()

################################################################################
# pca
# 
# output: pca_df, data.frame to plot PCA and srt; pca_var, percent variance
################################################################################

source("~/Documents/BFX_proj/RNAseq_pipeline/R_code/pca.R")

for(m_s in srt){
  png(paste0(out_dir, "fig/pca_", m_s, ".png"), height = 400, width = 400, res = 100)
  if(is.numeric(pca_df[, m_s])){
    print(ggplot(pca_df, aes_string("PC1", "PC2", color = m_s)) +
            geom_point() +
            xlab(paste0("PC1 ", pca_var["PC1"], "%")) + ylab(paste0("PC2 ", pca_var["PC2"], "%")) +
            ggtitle(m_s) +
            scale_color_gradient(low = attr(col_leg[[m_s]], "colors")[1], high = attr(col_leg[[m_s]], "colors")[2]) +
            theme_dark() +
            theme(aspect.ratio = 1,
                  legend.position = "none"))
  } else {
    print(ggplot(pca_df, aes_string("PC1", "PC2", color = m_s)) +
            geom_point() +
            xlab(paste0("PC1 ", pca_var["PC1"], "%")) + ylab(paste0("PC2 ", pca_var["PC2"], "%")) +
            ggtitle(m_s) +
            scale_color_manual(values = col_leg[[m_s]]) +
            theme_dark() +
            theme(aspect.ratio = 1,
                  legend.position = "none"))
  }
  dev.off()
}

# export PCA coordinates and percent variance
write.csv(pca_df, paste0(out_dir, "df/pca.csv"), row.names = T)
saveRDS(pca_var, paste0(out_dir, "df/pca_var.rds"))

################################################################################
# heatmap
# 
# output: hm, Heatmap
################################################################################

source("~/Documents/BFX_proj/RNAseq_pipeline/R_code/ap_heatmap.R")

png(paste0(out_dir, "fig/heatmap.png"), height = 500, width = 1000, res = 100)
draw(hm)
dev.off()

# export HeatmapList for additional downstream analysis
saveRDS(draw(hm), paste0(out_dir, "df/heatmap_list.rds"))

################################################################################
# wgcna
# 
# output: w_hm, correlation plot; mod_genes_long, genes and mod_gsea_long, pathways associated with gene modules
################################################################################

source("~/Documents/BFX_proj/RNAseq_pipeline/R_code/wgcna.R")

png(paste0(out_dir, "fig/wgcna_corplot.png"), height = 700, width = 500, res = 100)
draw(w_hm)
dev.off()

write.csv(mod_genes_long, paste0(out_dir, "df/mod_genes.csv"))
write.csv(mod_gsea_long, paste0(out_dir, "df/mod_gsea.csv"))
