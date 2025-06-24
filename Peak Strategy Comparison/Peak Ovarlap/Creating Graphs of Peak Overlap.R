library(ggplot2)

###Load data
bysample <- readRDS("/Users/yulongqiu/Desktop/Biostats/Master_Thesis/peak_df_by_with_overlaps.rds")
allsample <- readRDS("/Users/yulongqiu/Desktop/Biostats/Master_Thesis/peak_df_all_with_overlaps.rds")


# By-sample 
bysample_table <- as.data.frame(table(bysample$overlaps))
colnames(bysample_table) <- c("overlaps", "count")
bysample_table$proportion <- bysample_table$count / sum(bysample_table$count)
bysample_table$overlaps <- factor(bysample_table$overlaps, levels = bysample_table$overlaps)

# All-sample 
allsample_table <- as.data.frame(table(allsample$overlaps))
colnames(allsample_table) <- c("overlaps", "count")
allsample_table$proportion <- allsample_table$count / sum(allsample_table$count)
allsample_table$overlaps <- factor(allsample_table$overlaps, levels = allsample_table$overlaps)

# Draw separately
p1 <- ggplot(bysample_table, aes(x = overlaps, y = count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  geom_text(aes(label = paste0(round(proportion * 100, 2), "%")),
            vjust = -0.5, size = 3.5) +
  labs(title = "By-sample Peak Overlap Counts",
       x = "Number of Overlaps in Query Set", y = "Count") +
  theme_minimal()

p2 <- ggplot(allsample_table, aes(x = overlaps, y = count)) +
  geom_bar(stat = "identity", fill = "salmon", color = "black") +
  geom_text(aes(label = paste0(round(proportion * 100, 2), "%")),
            vjust = -0.5, size = 3.5) +
  labs(title = "All-sample Peak Overlap Counts",
       x = "Number of Overlaps in Query Set", y = "Count") +
  theme_minimal()

print(p1)
print(p2)

# Merge the data for comparative plotting
bysample_table$Reference_Set <- "By-sample method"
bysample_table$overlaps <- as.character(bysample_table$overlaps)

allsample_table$Reference_Set <- "All-samples method"
allsample_table$overlaps <- as.character(allsample_table$overlaps)

plot_df <- rbind(bysample_table, allsample_table)
plot_df$overlaps <- factor(plot_df$overlaps,
                           levels = sort(unique(as.numeric(plot_df$overlaps))))

# the merged comparison chart
p <- ggplot(plot_df, aes(x = overlaps, y = count, fill = Reference_Set)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
  geom_text(aes(label = paste0(round(proportion * 100, 2), "%")),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3.2) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = c("By-sample method" = "salmon", "All-samples method" = "skyblue")) +
  labs(title = "Comparison of Peak Overlap Counts",
       x = "Number of Overlaps in Query Set",
       y = "Count",
       fill = "Reference Set") +
  theme_minimal()
print(p)

# Save graph
ggsave("/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Output/peak_overlap_comparison_barplot.png",
       plot = p, width = 8, height = 5, dpi = 300)

