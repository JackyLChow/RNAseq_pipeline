# Generalized WGCNA workflow
## input: normalized counts data, curated metadata for WGCNA
## output: images of WGCNA heatmaps, table of WGCNA gene modules, table of WGCNA GSEA results

# CRAN
library(readr)
library(dplyr)
library(WGCNA)
library(scico)
library(circlize)
library(matrixStats)
library(stringr)

# Bioconductor
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(ReactomePA)
library(clusterProfiler)

################################################################################
#
# parameters
#
################################################################################

### trim count matrix for speed, if NULL then no trim ---
trim_c_mtx <- 5000

### gene synonym reference ---
hs <- org.Hs.eg.db
hs <- AnnotationDbi::select(hs,
                            keys = rownames(c_mtx),
                            columns = c("ENTREZID"),
                            keytype = "SYMBOL")
hs <- hs[!duplicated(hs$SYMBOL), ]

################################################################################
#
# prepare for WGCNA
#
################################################################################

tnc_mtx <- t(c_mtx) # transpose counts for WGCNA

### WGCNA qc for genes and samples ---
gsg_w <- goodSamplesGenes(tnc_mtx)
tnc_mtx <- tnc_mtx[gsg_w$goodSamples, gsg_w$goodGenes]

# clean up
rm(gsg_w)

################################################################################
#
# run gently modified WGCNA workflow
#
################################################################################

### filter top variant genes for speed ---
if(is.numeric(trim_c_mtx)){
  tnc_mtx <- tnc_mtx[, order(colVars(as.matrix(tnc_mtx)), decreasing = T)[1:trim_c_mtx]]
}

rm(trim_c_mtx)

### assay possible thresholds ---
sft_w <- pickSoftThreshold(tnc_mtx, powerVector = c(c(1:10), seq(from = 12, to = 20, by = 2)))

### calculate gene similarity ---
adj_w <- adjacency(tnc_mtx, power = sft_w$powerEstimate) # calculate adjacency matrix; how connected each gene is to other genes
tom_w <- TOMsimilarity(adj_w) # calculate topological overlap
geneTree_w <- hclust(as.dist(1 - tom_w), method = "average") # calculate gene tree by dissimilarity

### assign genes to gene module; produces vector of module assignments for each gene ---
mods <- cutreeDynamic(geneTree_w, # clustering
                      distM = 1 - tom_w, # distance matrix
                      deepSplit = 2,
                      pamRespectsDendro = F,
                      minClusterSize = 30) # minimum cluster size is 30 genes
mods <- labels2colors(mods)

### correlate final modules to traits ---
MEs_w <- moduleEigengenes(tnc_mtx, colors = mods)$eigengenes # re-calculate eigengenes with merged ME colors
MEs_w <- orderMEs(MEs_w) # cluster MEs by similarity; put grey (unassigned) at end
moduleTraitCor_w <- cor(MEs_w, meta_dum, use = "p") # correlation of sample eigenvalue vs sample meta value; pairwiase complete observations
moduleTraitPvalue_w <- corPvalueStudent(moduleTraitCor_w, nrow(tnc_mtx)) # calculate P value of correlations
moduleTraitFDR_w <- matrix(p.adjust(moduleTraitPvalue_w, method = "BH"),
                           nrow = nrow(moduleTraitPvalue_w),
                           ncol = ncol(moduleTraitPvalue_w),
                           dimnames = list(rownames(moduleTraitPvalue_w),
                                           colnames(moduleTraitPvalue_w))) # adjust by FDR

### correlate genes to final modules ---
geneModuleMembership_w <- cor(tnc_mtx, MEs_w, use = "p")
MMPvalue_w <- corPvalueStudent(as.matrix(geneModuleMembership_w), nrow(tnc_mtx))
MMPFDR_w <- matrix(p.adjust(MMPvalue_w, method = "BH"),
                   nrow = nrow(geneModuleMembership_w),
                   dimnames = list(rownames(MMPvalue_w),
                                   colnames(MMPvalue_w))) # adjust by FDR
### clean up ---
mod_trait <- list(cor = moduleTraitCor_w, pval = moduleTraitPvalue_w, fdr = moduleTraitFDR_w)
gene_mod <- list(gMM = geneModuleMembership_w, pval = MMPvalue_w, fdr = MMPFDR_w)

rm(list = ls()[grepl("_w$", ls())])

################################################################################
#
# analyze gene modules by GSEA
#
################################################################################

mod_genes <- list()
mod_gsea <- list()

### for loop to pull genes and run three levels of gsea
for(mod_ in unique(mods)){
  genes_ <- colnames(tnc_mtx)[mods == mod_] # get symbols of module genes
  
  entrez_ <- hs$ENTREZID[hs$SYMBOL %in% genes_] # convert symbols to entrezid
  entrez_ <- entrez_[!is.na(entrez_)] # filter missing entrezid
  
  reactome_ <- data.frame(enrichPathway( # run reactome enrichment
    entrez_, hs$ENTREZID,
    organism = "human",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.5,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE)); if(nrow(reactome_) > 0){reactome_$source <- "Reactome"}
  
  go_ <- data.frame(enrichGO( # run GO enrichment
    entrez_,
    org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "MF",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.5,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
    pool = FALSE)); if(nrow(go_) > 0){go_$source <- "GO"}
  
  gsea_ <- rbind(reactome_, go_)
  
  # add values to list
  mod_genes[[paste0("ME",mod_)]] <- genes_
  mod_gsea[[paste0("ME",mod_)]] <- gsea_
  
  rm(list = ls()[grepl("_$", ls())])
}

### export genes and gsea --
mod_genes_long <- data.frame(row.names = names(mod_genes),
                             mod_name = names(mod_genes))
mod_genes_long$mod_genes <- NULL
for(mod_ in mod_genes_long$mod_name){
  mod_genes_long[mod_, "mod_genes"] <- paste(mod_genes[[mod_]], collapse = ", ")
  rm(mod_) # clean up
}

mod_gsea_long <- bind_rows(mod_gsea, .id = "column_label")

rm(mod_genes, mod_gsea, mods, hs)

################################################################################
#
# visualize gene module correlation to traits
#
################################################################################

### heatmap parameters ---
# heatmap fill scale
col_max <- round(max(abs(mod_trait$cor)), digits = 2)

col_fun <- circlize::colorRamp2(c(-col_max, -col_max/2, 0, col_max/2, col_max),
                                c(scico(5, palette = "vanimo")[1], scico(5, palette = "vanimo")[2],
                                  scico(5, palette = "vanimo")[3],
                                  scico(5, palette = "vanimo")[4], scico(5, palette = "vanimo")[5]))

### plot module eigenvalue vs trait ---
# plot heatmap1; correlation only
anno_row <- HeatmapAnnotation(ME = rownames(mod_trait$cor),
                              col = list(ME = setNames(str_remove(rownames(mod_trait$cor), "ME"), rownames(mod_trait$cor))),
                              which = "row",
                              show_legend = F)

set.seed(415); w_hm <- Heatmap(mod_trait$cor,
                               name = "Correlation",
                               col = col_fun,
                               right_annotation = anno_row,
                               show_row_names = T,
                               row_names_side = "right",
                               width = unit(ncol(mod_trait$cor) * 8, "mm"),
                               height = unit(nrow(mod_trait$cor) * 8, "mm"),
                               cluster_columns = F)

rm(col_max, col_fun, anno_row, mod_trait, gene_mod, tnc_mtx)






































