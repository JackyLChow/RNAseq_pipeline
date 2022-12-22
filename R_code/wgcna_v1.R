################################################################################
#
# weighted genetic correlation network analysis
#
################################################################################

# input: normalized counts data, "smart" metadata
# output: images of WGCNA heatmaps, table of WGCNA gene modules

###---data management ----------------------------------------------------------
# load data
counts_files <- list.files(parameters["output_folder", "value"])
counts_files <- counts_files[grepl("^expression", counts_files)] # identify processed expression data
counts_names <- str_remove(counts_files, ".RDS")
counts_data <- list()
for(i in 1:length(counts_names)){
  # transpose counts data for WGCNA
  counts_data[[counts_names[i]]] <- t(readRDS(paste0("~/Documents/BFX_proj/RNAseq_pipeline/_output/Chow_PNAS_2020/", counts_files[i])))
  rm(i)
}
metadata <- readRDS("~/Documents/BFX_proj/RNAseq_pipeline/_output/Chow_PNAS_2020/metadata.RDS")
# clean up
rm(counts_files, counts_names) 

# dummify metadata for WGCNA
dummified_columns <- unlist(strsplit(parameters["metadata_wgcna_columns", "value"], ", "))
dummified_metadata <- list()
for(col_ in dummified_columns){
  dummified_metadata[[col_]] <- dummify(metadata[, col_]) # dummify column
  colnames(dummified_metadata[[col_]]) <- paste0(col_, "_", colnames(dummified_metadata[[col_]]))
  rm(col_)
}
rm(dummified_columns)
dummified_metadata <- as.data.frame(do.call(cbind.data.frame, dummified_metadata))
rownames(dummified_metadata) <- rownames(metadata)

# gene synonym table
gene_synonyms <- org.Hs.eg.db
gene_synonyms <- AnnotationDbi::select(gene_synonyms,
                                       keys = colnames(counts_data[[1]]),
                                       columns = c("ENTREZID"),
                                       keytype = "SYMBOL")
gene_synonyms <- gene_synonyms[!duplicated(gene_synonyms$SYMBOL), ]

# trim counts data for speed
trim_input_counts <- 5000

###---WGCNA pipeline------------------------------------------------------------

for(counts_ in names(counts_data)){
  ####---load data------------------------------------------------------------
  input_counts_ <- counts_data[[counts_]] # select counts data
  gsg_ <- goodSamplesGenes(input_counts_) # WGCNA gene and sample QC
  input_counts_ <- input_counts_[gsg_$goodSamples, gsg_$goodGenes] # filter bad genes and samples
  dummified_metadata_ <- dummified_metadata[rownames(input_counts_), ] # filter and match metadata rows
  rm(gsg_) # clean up
  
  # filter top variant genes for speed
  if(is.numeric(trim_input_counts)){
    input_counts_ <- input_counts_[, order(colVars(as.matrix(input_counts_)), decreasing = T)[1:trim_input_counts]]
  }
  
  ####---output file------------------------------------------------------------
  wgcna_out_ <- list()
  
  ####---assign genes to modules------------------------------------------------
  # assay possible thresholds
  soft_threshold_ <- pickSoftThreshold(input_counts_, powerVector = c(c(1:10), seq(from = 12, to = 20, by = 2)))
  
  # calculate gene similarity
  adjacency_ <- adjacency(input_counts_, power = soft_threshold_$powerEstimate) # calculate adjacency matrix; how connected each gene is to other genes
  topological_overlap_matrix_ <- TOMsimilarity(adjacency_) # calculate topological overlap
  gene_tree_ <- hclust(as.dist(1 - topological_overlap_matrix_), method = "average") # calculate gene tree by dissimilarity
  
  # assign genes to gene module; produces vector of module assignments for each gene
  gene_module_numeric_assignment_ <- cutreeDynamic(gene_tree_, # clustering
                                                   distM = 1 - topological_overlap_matrix_, # distance matrix
                                                   deepSplit = 2,
                                                   pamRespectsDendro = F,
                                                   minClusterSize = 30) # minimum cluster size is 30 genes
  gene_module_assignment_ <- labels2colors(gene_module_numeric_assignment_) # convert numeric to color

  # output
  wgcna_out_[["gene_module_assignment"]] <- data.frame("SYMBOL" = colnames(input_counts_),
                                                       "module_numeric_assignment" = gene_module_numeric_assignment_,
                                                       "module_assignment" = gene_module_assignment_)
  
  ####---calculate gene module significance-------------------------------------
  # correlate final modules to traits
  module_eigenvalues_ <- moduleEigengenes(input_counts_, colors = gene_module_assignment_)$eigengenes # re-calculate eigengenes with merged ME colors
  module_eigenvalues_ <- orderMEs(module_eigenvalues_) # cluster MEs by similarity; put grey (unassigned) at end
  module_trait_correlation_ <- cor(module_eigenvalues_, dummified_metadata_, use = "p") # correlation of sample eigenvalue vs sample meta value; pairwiase complete observations
  module_trait_pvalue_ <- corPvalueStudent(module_trait_correlation_, nrow(input_counts_)) # calculate P value of correlations
  module_trait_fdr_ <- matrix(p.adjust(module_trait_pvalue_, method = "BH"),
                              nrow = nrow(module_trait_pvalue_),
                              ncol = ncol(module_trait_pvalue_),
                              dimnames = list(rownames(module_trait_pvalue_),
                                              colnames(module_trait_pvalue_))) # adjust by FDR
  # output
  wgcna_out_[["module_trait_correlation"]] <- module_trait_correlation_
  wgcna_out_[["module_trait_pvalue"]] <- module_trait_pvalue_
  wgcna_out_[["module_trait_fdr"]] <- module_trait_fdr_
  
  # correlate genes to final modules
  gene_module_correlation_ <- cor(input_counts_, module_eigenvalues_, use = "p") # correlate gene expression to module
  gene_module_pvalue_ <- corPvalueStudent(as.matrix(gene_module_correlation_), nrow(input_counts_))
  gene_module_fdr_ <- matrix(p.adjust(gene_module_pvalue_, method = "BH"),
                             nrow = nrow(gene_module_pvalue_),
                             dimnames = list(rownames(gene_module_pvalue_),
                                             colnames(gene_module_pvalue_))) # adjust by FDR
  # output
  wgcna_out_[["gene_module_correlation"]] <- gene_module_correlation_
  wgcna_out_[["gene_module_pvalue"]] <- gene_module_pvalue_
  wgcna_out_[["gene_module_fdr"]] <- gene_module_fdr_
  
  ####---calculate module gene set enrichments----------------------------------
  module_genes_ <- list()
  module_gsea_ <- list()
  
  # for loop to pull genes and run gsea on multiple databases
  for(mod_ in unique(wgcna_out_[["gene_module_assignment"]][, "module_assignment"])){
    gene_module_assignment_ <- wgcna_out_[["gene_module_assignment"]] # call gene module table
    genes_ <- gene_module_assignment_$SYMBOL[which(gene_module_assignment_$module_assignment == mod_)] # filter module genes
    gene_entrez_ <- gene_synonyms$ENTREZID[gene_synonyms$SYMBOL %in% genes_] # convert symbols to entrezid
    gene_entrez_ <- gene_entrez_[!is.na(gene_entrez_)] # filter missing entrezid
    
    # universe is only genes from subset
    #universe_entrez_ <- gene_synonyms$ENTREZID[gene_synonyms$SYMBOL %in% gene_module_assignment_$SYMBOL]
    #universe_entrez_ <- universe_entrez_[!is.na(universe_entrez_)] # filter missing entrezid
    
    universe_entrez_ <- gene_synonyms$ENTREZID
    
    reactome_ <- data.frame(enrichPathway( # run reactome enrichment
      gene = gene_entrez_, universe = universe_entrez_,
      organism = "human",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.5,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      readable = FALSE)); if(nrow(reactome_) > 0){reactome_$source <- "Reactome"}
    
    go_ <- data.frame(enrichGO( # run GO enrichment
      gene = gene_entrez_, universe = universe_entrez_,
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
    
    kegg_ <- data.frame(enrichKEGG( # run KEGG enrichment
      gene = gene_entrez_, universe = universe_entrez_,
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.5,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500)); if(nrow(kegg_) > 0){kegg_$source <- "KEGG"}
    
    gsea_ <- rbind(reactome_, go_, kegg_)
    
    # add values to list
    module_genes_[[paste0("ME",mod_)]] <- genes_
    module_gsea_[[paste0("ME",mod_)]] <- gsea_
  }
  
  # export genes and gsea
  module_genes_long_ <- data.frame(row.names = names(module_genes_), module_name = names(module_genes_))
  module_genes_long_$module_genes <- NULL
  for(mod_ in module_genes_long_$module_name){
    module_genes_long_[mod_, "mod_genes"] <- paste(module_genes_[[mod_]], collapse = ", ")
  }
  wgcna_out_[["module_genes"]] <- module_genes_long_
  wgcna_out_[["module_gsea"]] <- bind_rows(module_gsea_, .id = "module_name")
  
  ####---plot correlation results-----------------------------------------------
  
  # call necessary data sets
  module_genes_ <- wgcna_out_$module_genes
  module_gsea_ <- wgcna_out_$module_gsea
  module_trait_correlation_ <- wgcna_out_$module_trait_correlation
  module_trait_fdr_ <- wgcna_out_$module_trait_fdr
  gene_module_correlation_ <- wgcna_out_$gene_module_correlation
  
  ### heatmap parameters ---
  # heatmap fill scale
  col_max_ <- round(max(abs(module_trait_correlation_)), digits = 2)
  
  col_fun_ <- circlize::colorRamp2(c(-col_max_, -col_max_/2, 0, col_max_/2, col_max_),
                                   c("purple1", "purple4", "black", "yellow4", "yellow1"))
  
  # heatmap module annotation with shortened GSEA and gene lists
  module_data_ <- data.frame(row.names = module_genes_$module_name, mod = module_genes_$module_name)
  module_data_$paths <- NULL
  module_data_$genes <- NULL
  
  for(mod_ in rownames(module_data_)){
    if(sum(module_gsea_$module_name == mod_) > 0){
      gsea_ <- module_gsea_[module_gsea_$module_name == mod_, ]
      gsea_$s_d <- paste0(gsea_$source, ": ", substr(gsea_$Description, 1, 40), "...")
      gsea_ <- gsea_[order(gsea_$pvalue), ]
      
      paths_ <- paste0(" ", gsea_$s_d[1], "\n ", gsea_$s_d[2])
      module_data_[mod_, "paths"] <- paths_
    }
    
    gene_ <- module_genes_[module_genes_$module_name == mod_, "mod_genes"] # trim module membership table to genes assigned to module
    gene_ <- unlist(strsplit(gene_, ", "))
    gene_ <- gene_module_correlation_[gene_, mod_]
    gene_ <- gene_[order(gene_, decreasing = T)] # order by membership
    gene_ <- paste0(" pos: ", paste(names(head(gene_, 5)), collapse = ", "), "\n",
                    " neg: ", paste(names(tail(gene_, 5)), collapse = ", "))
    module_data_[mod_, "genes"] <- gene_
  }
  
  module_data_ <- module_data_[rownames(module_trait_correlation_), ] # rearrange to match correlation matrix
  
  # plot heatmap1; correlation only
  anno_row_l_ <- HeatmapAnnotation(ME = rownames(module_trait_correlation_),
                                   col = list(ME = setNames(str_remove(rownames(module_trait_correlation_), "ME"), rownames(module_trait_correlation_))),
                                   which = "row",
                                   show_legend = F,
                                   show_annotation_name = c(ME = F))
  
  anno_row_r_ <- HeatmapAnnotation(GSEA = anno_text(module_data_$paths, gp = gpar(fontsize = 10, border = "black")),
                                  genes = anno_text(module_data_$genes, gp = gpar(fontsize = 10, border = "black")),
                                  which = "row")
  
  png(paste0(parameters["output_folder", "value"], "WGCNA_correlation_heatmap_", counts_, ".png"), height = nrow(module_trait_correlation_) * 37.5 + 500, width = ncol(module_trait_correlation_) * 37.5 + 1500, res = 100)
  set.seed(415); hm_ <- Heatmap(module_trait_correlation_,
                                name = "Correlation",
                                col = col_fun_,
                                right_annotation = anno_row_r_,
                                left_annotation = anno_row_l_,
                                show_row_names = T,
                                row_names_side = "left",
                                width = unit(ncol(module_trait_correlation_) * 8, "mm"),
                                height = unit(nrow(module_trait_correlation_) * 8, "mm"),
                                rect_gp = gpar(type = "none"),
                                cell_fun = function(j, i, x, y, width, height, fill) {
                                  if(module_trait_fdr_[i, j] < 0.01){ # greyed out does not pass FDR threshold
                                    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey90", fill = "white"))
                                  } else if(module_trait_fdr_[i, j] < 0.05) {
                                    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey90", fill = "grey75"))
                                  } else if(module_trait_fdr_[i, j] < 0.5) {
                                    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey90", fill = "grey25"))
                                  } else {
                                    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey90", fill = "black"))
                                  }
                                  grid.circle(x = x, y = y, r = unit(3, "mm"),
                                              gp = gpar(fill = col_fun_(module_trait_correlation_[i, j]), col = "grey80")) # color is correlation
                                },
                                cluster_columns = F)
  
  print(hm_)
  dev.off()
  checkpoint <- rbind(checkpoint, c(paste0("WGCNA_correlation_heatmap_", counts_, ".png"), format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))); write.csv(checkpoint, checkpoint_file, row.names = F)
  
  ###---save WGCNA outputs------------------------------------------------------
  saveRDS(wgcna_out_, paste0(parameters["output_folder", "value"], "WGCNA_results_", counts_, ".RDS"))
  checkpoint <- rbind(checkpoint, c(paste0("WGCNA_results_", counts_, ".RDS"), format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))); write.csv(checkpoint, checkpoint_file, row.names = F)
  
  rm(list = ls()[grepl("_$", ls())])
}

rm(dummified_metadata, metadata, counts_data, trim_input_counts)


