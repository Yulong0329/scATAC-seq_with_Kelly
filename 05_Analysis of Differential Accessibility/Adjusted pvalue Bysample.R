cluster_ids <- c("10", "6", "9", "7", "15", "8", "11", "3", "16", "12")
base_path <- "/home1/yulongqi/Feature_Matrix/Clusters_vs_the_rest_Bysample/FindMarker_bysample_0506/"

adjp_list <- list()


for (cl in cluster_ids) {
  file_path <- paste0(base_path, "markers_cluster", cl, "_vs_rest.rds")
  marker_df <- readRDS(file_path)
  df <- data.frame(
    peak_id = rownames(marker_df),
    adj_pval = marker_df$p_val_adj
  )
  colnames(df)[2] <- paste0("cluster", cl)
  
  adjp_list[[cl]] <- df
}

#Merge by peak id and retain all peaks
adjp_df <- Reduce(function(x, y) merge(x, y, by = "peak_id", all = TRUE), adjp_list)

head(adjp_df)

saveRDS(adjp_df, file = "/home1/yulongqi/Feature_Matrix/Clusters_vs_the_rest_Bysample/FindMarker_bysample_0506/all_clusters_adj_pvals_matrix.rds")
write.csv(adjp_df, file = "/home1/yulongqi/Feature_Matrix/Clusters_vs_the_rest_Bysample/FindMarker_bysample_0506/all_clusters_adj_pvals_matrix.csv", row.names = FALSE)


