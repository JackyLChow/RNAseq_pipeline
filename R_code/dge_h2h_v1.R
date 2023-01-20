library(DESeq2)

################################################################################
#
# compare DESeq2 to limma outputs, thought experiment
#
################################################################################

### data management ------------------------------------------------------------
# load data
counts_matrix_ <- readRDS(parameters["counts_file", "value"])
metadata <- readRDS(parameters["metadata_file", "value"])

# output data
dge_h2h_out <- list()

### dge limma ------------------------------------------------------------------
# set up design matrix
design_main_ <- factor(metadata[, parameters["metadata_comparison_column", "value"]])

# limma pipeline
design_ <- model.matrix(~ design_main_)
l_ <- lmFit(counts_matrix_, design_)
e_ <- eBayes(l_)
deg_ <- topTable(e_, n = Inf)

dge_h2h_out[["limma"]] <- deg_ # archive deg_ results

### dge DESeq2 -----------------------------------------------------------------
dds_ <- DESeqDataSetFromMatrix(round(counts_matrix_), # need to round to be accepted into DESeq2
                               metadata, 
                               design = formula(paste("~", parameters["metadata_comparison_column", "value"])))
dds_ <- DESeq(dds_)
deg_ <- data.frame(results(dds_))

dge_h2h_out[["DESeq2"]] <- deg_ # archive deg_ results

### plot l2fc estimates
lr_ <- dge_h2h_out[["limma"]]
lr_ <- data.frame(Symbol = rownames(lr_), limma_l2fc = lr_$logFC, limma_pval = lr_$P.Value, rownames = rownames(lr_))

dr_ <- dge_h2h_out[["DESeq2"]]
dr_ <- data.frame(Symbol = rownames(dr_), ds2_l2fc = dr_$log2FoldChange, ds2_pval = dr_$pvalue, rownames = rownames(dr_))

h2h_ <- inner_join(lr_, dr_, by = "Symbol")

ggplot(h2h_, aes(limma_l2fc, ds2_l2fc)) +
  geom_point(alpha = 0.05, size = 0.1)

ggplot(h2h_, aes(-log10(limma_pval), -log10(ds2_pval))) +
  geom_point(alpha = 0.05, size = 0.1)

### export results -------------------------------------------------------------
saveRDS(dge_h2h_out, paste0(parameters["output_folder", "value"], "dge_h2h_results.RDS"))

rm(list = ls()[grepl("_$", ls())])
rm(dge_out, metadata)
