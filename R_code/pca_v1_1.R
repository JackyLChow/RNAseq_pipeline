################################################################################
#
# principle component analysis
#
################################################################################

### data management ------------------------------------------------------------
# load data
counts_matrix_ <- readRDS(parameters["counts_file", "value"])
counts_ <- parameters["counts_type", "value"]
metadata <- readRDS(parameters["metadata_file", "value"])
comparisons_ <- unlist(str_split(parameters["metadata_comparison_columns", "value"], ","))

# output data
pca_out <- list()

### functions ------------------------------------------------------------------
# function to replace ENTREZID with SYMBOL in GSEA outputs
entrezid_to_symbol <- function(gsea_results){
  gsea_out_ <- gsea_results
  gsea_out_[, "core_enrichment_symbol"] <- NULL
  if(nrow(gsea_out_) == 0){
    cat("\nno enriched pathways\n")
    return(gsea_out_)
  } else {
    for(i in rownames(gsea_out_)){
      entrezids_ <- unlist(str_split(gsea_out_[i, "core_enrichment"], "/"))
      symbols_ <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = entrezids_,
                                        columns = "SYMBOL",
                                        keytype = "ENTREZID")[, "SYMBOL"]
      gsea_out_[i, "core_enrichment_symbol"] <- paste(symbols_, collapse = "/")
    }
    return(gsea_out_)
  }
}

#-------------------------------------------------------------------------------
# PCA pipeline 
#-------------------------------------------------------------------------------

### run PCA --------------------------------------------------------------------
pca_ <- prcomp(t(counts_matrix_))
pca_out[[paste0(counts_, "_plot")]] <- data.frame(pca_$x, metadata) # export data in plottable format
pca_out[[paste0(counts_, "_prcomp_output")]] <- pca_ # export prcomp output

### PCA eigenvector GSEA -------------------------------------------------------
# extract eigenvectors for first three PCs
pca_eigenvectors_ <- pca_out[[paste0(counts_, "_prcomp_output")]][["rotation"]][, paste0("PC", 1:3)]

for(i in paste0("PC", 1:3)){
  # make ranked gene list
  ## extract ordered gene list by logFC
  ordered_genes_ <- pca_eigenvectors_[, i]
  ordered_genes_ <- data.frame(Symbol = names(ordered_genes_), eigen_value = ordered_genes_)
  ordered_genes_ <- ordered_genes_[order(ordered_genes_$eigen_value, decreasing = T), ]
  
  ## generate ENTREZID key from synonym table
  key_ <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = ordered_genes_$Symbol,
                                columns = c("ENTREZID"),
                                keytype = "SYMBOL")
  key_ <- key_[!duplicated(key_$SYMBOL), ] # remove duplicated symbols
  
  ## load hallmark pathways
  msig_h_ <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene) %>%
    dplyr::rename(ont = gs_name, gene = entrez_gene)
  
  ## assign ENTREZID to gene
  ordered_genes_ <- left_join(ordered_genes_, key_, by = c("Symbol" = "SYMBOL"))
  
  ## filter genes with duplicated and missing ENTREZID
  ordered_genes_ <- ordered_genes_[!is.na(ordered_genes_$ENTREZID), ]
  ordered_genes_ <- ordered_genes_[!duplicated(ordered_genes_$ENTREZID), ]
  
  ## make ranked list
  ranked_list_ <- ordered_genes_[, "eigen_value"]
  names(ranked_list_) <- as.character(ordered_genes_[, "ENTREZID"])
  
  # run pathway enrichments
  ## reactome pathway enrichment
  set.seed(415); reactomeGSEA_ <- gsePathway(ranked_list_,
                                             maxGSSize = 500,
                                             pvalueCutoff = 0.1) %>% data.frame()
  pca_out[[paste0(i, "_Reactome_GSEA")]] <- entrezid_to_symbol(reactomeGSEA_)
  
  ## KEGG pathway enrichment
  set.seed(415); keggGSEA_ <- gseKEGG(ranked_list_,
                                      organism = "hsa",
                                      pvalueCutoff = 0.1) %>% data.frame()
  pca_out[[paste0(i, "_KEGG_GSEA")]] <- entrezid_to_symbol(keggGSEA_)
  
  ## GO pathway enrichment
  set.seed(415); goGSEA_ <- gseGO(ranked_list_,
                                  OrgDb = org.Hs.eg.db,
                                  pvalueCutoff = 0.1) %>% data.frame()
  pca_out[[paste0(i, "_GO_GSEA")]] <- entrezid_to_symbol(goGSEA_)
  
  ## Hallmark pathway enrichment
  set.seed(415); hallmarkGSEA_ <- GSEA(ranked_list_, TERM2GENE = msig_h_, scoreType = "pos") %>% data.frame()
  pca_out[[paste0(i, "_Hallmark_GSEA")]] <- entrezid_to_symbol(hallmarkGSEA_)
}

### plot PCA -------------------------------------------------------------------
pca_var_ <- round(pca_$sdev^2/sum(pca_$sdev^2), 3) * 100; names(pca_var_) <- colnames(pca_$x) # calculate percent variance explained
plot_df_ <- pca_out[[paste0(counts_, "_plot")]] # retrieve data frame for plot
point_size_ <- 10/log2(nrow(plot_df_))

for(i in comparisons_){
  pc_1_2_ <- ggplot(plot_df_, aes(PC1, PC2, color = .data[[i]])) +
    xlab(paste0("PC1 ", pca_var_[1], "%")) +
    ylab(paste0("PC2 ", pca_var_[2], "%")) +
    geom_point(alpha = 0.75, size = point_size_) +
    theme(aspect.ratio = 1)
  pc_3_4_ <- ggplot(plot_df_, aes(PC3, PC4, color = .data[[i]])) +
    xlab(paste0("PC3 ", pca_var_[3], "%")) +
    ylab(paste0("PC4 ", pca_var_[4], "%")) +
    geom_point(alpha = 0.75, size = point_size_) +
    theme(aspect.ratio = 1)
  pc_5_6_ <- ggplot(plot_df_, aes(PC5, PC6, color = .data[[i]])) +
    xlab(paste0("PC5 ", pca_var_[5], "%")) +
    ylab(paste0("PC6 ", pca_var_[6], "%")) +
    geom_point(alpha = 0.75, size = point_size_) +
    theme(aspect.ratio = 1)
  pc_plot_ <- plot_grid(pc_1_2_ + theme(legend.position = "none"),
                        pc_3_4_ + theme(legend.position = "none"),
                        pc_5_6_ + theme(legend.position = "none"), nrow = 1)
  legend_ <- get_legend(pc_1_2_ + theme(legend.box.margin = margin(0, 12, 0, 12)))
  pc_plot_ <- plot_grid(pc_plot_, legend_, rel_widths = c(3, .4))
  title_ <- ggdraw() +
    draw_label(paste0("PCA on ", counts_), fontface = 'bold', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  png(paste0(parameters["output_folder", "value"], "PCA_on_", counts_, "_", i, ".png"), width = 1200, height = 400)
  print(plot_grid(title_, pc_plot_, ncol = 1, rel_heights = c(0.1, 1)))
  dev.off()
  checkpoint <- rbind(checkpoint, c(paste0("PCA_on_", counts_, "_", i, ".png"), format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))); write.csv(checkpoint, checkpoint_file, row.names = F)
  rm(i)
}

### export results -------------------------------------------------------------
saveRDS(pca_out, paste0(parameters["output_folder", "value"], "PCA_results.RDS"))
checkpoint <- rbind(checkpoint, c("PCA_results.RDS", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))); write.csv(checkpoint, checkpoint_file, row.names = F)

rm(list = ls()[grepl("_$", ls())])
rm(pca_out, metadata)
