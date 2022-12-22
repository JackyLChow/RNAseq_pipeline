################################################################################
#
# load and clean data
#
################################################################################

# load counts table
input_counts <- data.frame(read.csv(parameters["counts_file", "value"]))
counts_type <- parameters["counts_type", "value"]
duplicated_rows <- duplicated(input_counts[, "gene"]) # identify duplicate gene entries
input_counts <- input_counts[!duplicated_rows, ] # remove duplicated genes
input_counts <- data.frame(input_counts[, 2:ncol(input_counts)], row.names = input_counts[, "gene"]) # trim raw counts
rm(duplicated_rows)

# load metadata
metadata <- data.frame(read.csv(parameters["metadata_file", "value"], row.names = 1))
rownames(metadata) <- sub("-", ".", rownames(metadata))
saveRDS(metadata, paste0(parameters["output_folder", "value"], "metadata.RDS")) # export cleaned metadata
checkpoint <- rbind(checkpoint, c("metadata.RDS", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))); write.csv(checkpoint, checkpoint_file, row.names = F)

if(counts_type == "raw"){
  # voom normalization
  voom_counts <- voom(as.matrix(input_counts))
  saveRDS(voom_counts$E, paste0(parameters["output_folder", "value"], "expression_voom_counts.RDS")) # export normalized counts
  checkpoint <- rbind(checkpoint, c("expression_voom_counts.RDS", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))); write.csv(checkpoint, checkpoint_file, row.names = F)
  counts_data <<- voom_counts$E
  rm(voom_counts)
} else {
  counts_data <<- log2(as.matrix(input_counts) + 1)
  saveRDS(counts_data, paste0(parameters["output_folder", "value"], "expression_log2_", counts_type, "_plus_1_counts.RDS")) # export normalized counts
  checkpoint <- rbind(checkpoint, c(paste0("expression_log2_", counts_type, "_plus_1_counts.RDS"), format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))); write.csv(checkpoint, checkpoint_file, row.names = F)
}

# batch correction
if(!is.na(parameters["metadata_batch_column", "value"])){
  combat_counts <- ComBat(counts_data, # normalized counts
                          batch = metadata[, parameters["metadata_batch_column", "value"]], # batch is sex
                          mod = model.matrix(~1, data = metadata)) # model is blank for batch correction
  saveRDS(combat_counts, paste0(parameters["output_folder", "value"], "expression_combat_counts.RDS")) # export batch corrected counts
  checkpoint <- rbind(checkpoint, c("expression_combat_counts.RDS", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))); write.csv(checkpoint, checkpoint_file, row.names = F)
  rm(combat_counts)
}

rm(input_counts, metadata, counts_data)
