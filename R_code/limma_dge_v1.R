################################################################################
#
# differential gene expression
#
################################################################################

### data management ------------------------------------------------------------
# load data
counts_matrix_ <- readRDS(parameters["counts_file", "value"])
counts_ <- parameters["counts_type", "value"]
metadata <- readRDS(parameters["metadata_file", "value"])

# output data
dge_out <- list()

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

### voom normalization ---------------------------------------------------------
# use voom if raw counts are available
if(parameters["has_raw_counts", "value"] == "yes"){
  counts_matrix_ <- voom(counts_matrix_)
}

### dge limma ------------------------------------------------------------------
# set up design matrix
design_main_ <- factor(metadata[, parameters["metadata_comparison_column", "value"]])

# limma pipeline
design_ <- model.matrix(~ design_main_)
l_ <- lmFit(counts_matrix_, design_)
e_ <- eBayes(l_)
deg_ <- topTable(e_, n = Inf)
dge_out[["differential_gene_expression_resutls"]] <- deg_ # archive deg_ results

# volcano plot
## reformat deg_ table for next steps
deg_$Symbol <- rownames(deg_)
## pick genes to show
show_ <- c(rownames(deg_)[order(deg_[, "logFC"])[1:5]], # bottom 5 logFC
           rownames(deg_)[order(deg_[, "logFC"], decreasing = T)[1:5]], # top 5 logFC
           rownames(deg_)[order(deg_[, "adj.P.Val"])[1:5]]) # 5 most significant
deg_$show <- ifelse(deg_[, "Symbol"] %in% show_, "show", "no_show")

## make volcano plot
vp_ <- ggplot(deg_, aes(logFC, -log10(adj.P.Val))) +
  geom_point(data = deg_[deg_[, "show"] == "show", ], color = "grey60") +
  geom_point(data = deg_[deg_[, "show"] == "no_show", ], size = 0.1, color = "grey60") +
  geom_label_repel(data = deg_[deg_[, "show"] == "show", ], aes(label = Symbol), min.segment.length = 0) +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  xlab(paste0(levels(design_main_)[1], "    <-    logFC    ->    ", levels(design_main_)[2])) +
  theme(aspect.ratio = 1)
dge_out[["volcano_plot"]] <- vp_

png(paste0(parameters["output_folder", "value"], "volcano_plot.png"), width = 300, height = 300, res = 75)
print(vp_)
dev.off()
checkpoint <- rbind(checkpoint, c("volcano_plot.png", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))); write.csv(checkpoint, checkpoint_file, row.names = F)

### gsea -----------------------------------------------------------------------
# make ranked gene list
## extract ordered gene list by logFC
ordered_genes_ <- deg_[order(deg_$logFC, decreasing = T), ]

## generate ENTREZID key from synonym table
key_ <- AnnotationDbi::select(org.Hs.eg.db,
                              keys = ordered_genes_$Symbol,
                              columns = c("ENTREZID"),
                              keytype = "SYMBOL")
key_ <- key_[!duplicated(key_$SYMBOL), ] # remove duplicated symbols

## assign ENTREZID to gene
ordered_genes_ <- left_join(ordered_genes_, key_, by = c("Symbol" = "SYMBOL"))

## filter genes with duplicated and missing ENTREZID
ordered_genes_ <- ordered_genes_[!is.na(ordered_genes_$ENTREZID), ]
ordered_genes_ <- ordered_genes_[!duplicated(ordered_genes_$ENTREZID), ]

## make ranked list
ranked_list_ <- ordered_genes_[, "logFC"]
names(ranked_list_) <- as.character(ordered_genes_[, "ENTREZID"])

# run pathway enrichments
## reactome pathway enrichment
set.seed(415); reactomeGSEA_ <- gsePathway(ranked_list_,
                                           maxGSSize = 500,
                                           pvalueCutoff = 0.1) %>% data.frame()
dge_out[["Reactome_GSEA"]] <- entrezid_to_symbol(reactomeGSEA_)

## KEGG pathway enrichment
set.seed(415); keggGSEA_ <- gseKEGG(ranked_list_,
                                    organism = "hsa",
                                    pvalueCutoff = 0.1) %>% data.frame()
dge_out[["KEGG_GSEA"]] <- entrezid_to_symbol(keggGSEA_)

## GO pathway enrichment
set.seed(415); goGSEA_ <- gseGO(ranked_list_,
                                OrgDb = org.Hs.eg.db,
                                pvalueCutoff = 0.1) %>% data.frame()
dge_out[["GO_GSEA"]] <- entrezid_to_symbol(goGSEA_)


### export results -------------------------------------------------------------
saveRDS(dge_out, paste0(parameters["output_folder", "value"], "differential_gene_expression_results.RDS"))
checkpoint <- rbind(checkpoint, c("differential_gene_expression_results.RDS", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))); write.csv(checkpoint, checkpoint_file, row.names = F)

rm(list = ls()[grepl("_$", ls())])
rm(dge_out, metadata)
