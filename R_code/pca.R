# pca

pca_ <- prcomp(t(c_mtx)) # run PCA

pca_var <- round(pca_$sdev^2/sum(pca_$sdev^2), 3) * 100; names(pca_var) <- colnames(pca_$x) # calculate percent variance explained

pca_df <- data.frame(pca_$x, meta_srt) # make data.frame for ggplot

rm(list = ls()[grepl("_$", ls())])