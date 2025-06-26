library(ggplot2)

bysample <- read.csv("/home1/yulongqi/Feature_Matrix/Clusters_vs_the_rest_Bysample/FindMarker_bysample_0506/allclusters_sig.csv")
allsample <- read.csv("/home1/yulongqi/Feature_Matrix/Clusters_vs_the_rest_allsample/FindMarker_allsample_0506/allclusters_sig.csv")


# Step 1: Identify all cluster columns (excluding non-p-valued columns)）
cluster_cols_bysample <- grep("^cluster", colnames(bysample))

# Step 2: Force these columns to be converted to numeric types
bysample[cluster_cols_bysample] <- lapply(bysample[cluster_cols_bysample], as.numeric)

# Step 3: Make a significance judgment for each row (p.adj < 0.05 in any cluster)
bysample$allclusters <- apply(bysample[, cluster_cols_bysample], 1, function(x) {
  if (all(is.na(x))) {
    FALSE
  } else {
    any(x < 0.05, na.rm = TRUE)
  }
})

# Duplicate previous steps
cluster_cols_allsample <- grep("^cluster", colnames(allsample))
allsample[cluster_cols_allsample] <- lapply(allsample[cluster_cols_allsample], as.numeric)
allsample$allclusters <- apply(allsample[, cluster_cols_allsample], 1, function(x) {
  if (all(is.na(x))) {
    FALSE  # 全部是 NA，我们保守地认为不显著
  } else {
    any(x < 0.05, na.rm = TRUE)
  }
})




# Calculate the total number and the significant number
bysample_total <- nrow(bysample)         # 40138
bysample_sig   <- sum(bysample$allclusters)
table(bysample$allclusters)

allsample_total <- nrow(allsample)       # 11287
allsample_sig   <- sum(allsample$allclusters)



plot_data <- data.frame(
  method = c("by-sample", "all-samples"),
  total = c(bysample_total, allsample_total),
  significant = c(bysample_sig, allsample_sig)
)

plot_data$fill_group <- ifelse(plot_data$method == "All Samples",
                               "Significant peaks (All Samples)",
                               "Significant peaks (By Sample)")

fill_colors <- c("Significant peaks (All Samples)" = "skyblue",
                 "Significant peaks (By Sample)" = "salmon")
bg_data <- plot_data
bg_data$fill_group <- "Other"  

# Merge data for stacked plotting (significant + non-significant)
stacked_data <- rbind(
  data.frame(method = plot_data$method,
             yval = plot_data$significant,
             fill_group = plot_data$fill_group,
             label = paste0(round(100 * plot_data$significant / plot_data$total, 1), "%")),
  data.frame(method = plot_data$method,
             yval = plot_data$total - plot_data$significant,
             fill_group = "Other",
             label = NA)
)

stacked_colors <- c("Significant peaks (All Samples)" = "skyblue",
                    "Significant peaks (By Sample)" = "salmon",
                    "Other" = "white")

# Creating Plot
p<-ggplot(stacked_data, aes(x = method, y = yval, fill = fill_group)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(data = stacked_data[!is.na(stacked_data$label), ],
            aes(label = label),
            position = position_stack(vjust = 0.5),
            color = "white", size = 5, fontface = "bold") +
  scale_fill_manual(
    values = stacked_colors,
    name = NULL,
    labels = c(
      "Other Peaks",
      "Significant peaks (All Samples)",
      "Significant peaks (By Sample)"
    )
  ) +
  labs(y = "# peaks", x = NULL, title = "Proportion of Differentially Accessible Peaks Across Peak Calling Strategies") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "right")
# Save
ggsave("/home1/yulongqi/Feature_Matrix/significant_peak_proportions_0603.png", plot = p, width = 6, height = 4, dpi = 300)

