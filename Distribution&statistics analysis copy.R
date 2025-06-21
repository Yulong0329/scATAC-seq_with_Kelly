library(GenomicRanges)
library(dplyr)
library(ggplot2)

# Load peak data
load("/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Peaks_allsample/All_sample_Narrowpeaks/combined_peaks2.RData")
load("/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Peaks_bysample/Combined_and_reduced_Granges/combined_peaks.RData")

# Create data frame with length and strategy labels
all_sample_df <- data.frame(length = width(combined.peaks2), strategy = "All Sample")
by_sample_df  <- data.frame(length = width(combined.peaks), strategy = "By Sample")
peak_data <- rbind(all_sample_df, by_sample_df)

# Optional: remove very long peaks
max_len <- 3000
peak_data <- peak_data %>% filter(length <= max_len)

# Set bin width and assign bins
bin_width <- 100
peak_data <- peak_data %>%
  mutate(bin_start = floor(length / bin_width) * bin_width)

# Compute percentage for each bin
peak_summary <- peak_data %>%
  group_by(strategy, bin_start) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(strategy) %>%
  mutate(percent = count / sum(count) * 100)

# Plot using facet_wrap to get separate panels
ggplot(peak_summary, aes(x = bin_start, y = percent)) +
  geom_col(fill = "steelblue", alpha = 0.8, width = bin_width * 0.9) +
  facet_wrap(~strategy, ncol = 1) +
  scale_x_continuous(
    name = "Peak Length (bp)",
    breaks = seq(0, 3000, by = 1000)
  ) +
  ylab("Percentage of Peaks (%)") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 0, size = 10)
  ) +
  ggtitle("Peak Length Distribution (Percentage)")








library(ggplot2)
library(dplyr)

# 你前面已经完成 peak_summary 的生成
# 确保 strategy 是因子类型，便于颜色映射
peak_summary$strategy <- factor(peak_summary$strategy, levels = c("All Sample", "By Sample"))

# 绘图
ggplot(peak_summary, aes(x = bin_start, y = percent, fill = strategy)) +
  geom_col(alpha = 0.8, width = bin_width * 0.9, show.legend = FALSE) +
  facet_grid(. ~ strategy) +  # 横向排列
  scale_fill_manual(values = c("All Sample" = "skyblue", "By Sample" = "salmon")) +
  scale_x_continuous(
    name = "Peak Length (bp)",
    breaks = seq(0, 3000, by = 1000)
  ) +
  ylab("Percentage of Peaks (%)") +
  ggtitle("Peak Length Distribution by Strategy (Percentage)") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )




#两张图合并在一张图中，能够看到占比对比
library(ggplot2)
library(dplyr)

# 假设你已经有 peak_summary，包括 bin_start, percent, strategy 列

# 把 strategy 转为 factor，设置显示顺序
peak_summary$strategy <- factor(peak_summary$strategy, levels = c("All Sample", "By Sample"))

# 修改绘图
ggplot(peak_summary, aes(x = bin_start, y = percent, fill = strategy)) +
  geom_col(alpha = 0.8, width = bin_width * 0.9, position = "dodge") +
  scale_fill_manual(
    name = "Strategy",
    values = c("All Sample" = "skyblue", "By Sample" = "salmon"),
    labels = c("All-samples method", "By-sample method")
  ) +
  scale_x_continuous(
    name = "Peak Length (bp)",
    breaks = seq(0, 3000, by = 1000)
  ) +
  ylab("Percentage of Peaks (%)") +
  ggtitle("Peak Length Distribution by Strategy (Percentage)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )



#制作了有legend的图
library(ggplot2)
library(dplyr)

# 假设 peak_summary 已含 bin_start, percent, strategy 列
peak_summary$strategy <- factor(peak_summary$strategy, levels = c("All Sample", "By Sample"))

p_combined<-ggplot(peak_summary, aes(x = bin_start, y = percent, fill = strategy)) +
  geom_col(alpha = 0.8, width = bin_width * 0.9, show.legend = TRUE) +
  facet_grid(. ~ strategy) +  # 横向排列分组图
  scale_fill_manual(
    name = "Strategy",
    values = c("All Sample" = "skyblue", "By Sample" = "salmon"),
    labels = c("All-samples method", "By-sample method")
  ) +
  scale_x_continuous(
    name = "Peak Length (bp)",
    breaks = seq(0, 3000, by = 1000)
  ) +
  ylab("Percentage of Peaks (%)") +
  ggtitle("Peak Length Distribution by Strategy (Percentage)") +
  theme_minimal() +
  theme(
    strip.text = element_blank(),  # 取消 facet 上方标签
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )
ggsave(filename = file.path(output_folder, "peak_length_distribution_Percentage_combined2.png"), 
       plot = p_combined, width = 8, height = 6, dpi = 300)








#Define the output folder for saving plots
output_folder <- "/Users/yulongqiu/Desktop/Biostats/Master_Thesis/Output"

#Making Histogram
library(ggplot2)

p <- ggplot(peak_data, aes(x=length, fill=strategy)) +
  geom_histogram(alpha=0.5, position="identity", bins=50) +
  theme_minimal() +
  labs(title="Comparison of Peak Length Distribution1", x="Peak Length", y="Count") +
  scale_fill_manual(values=c("skyblue", "salmon"))
ggsave(filename = file.path(output_folder, "peak_length_distribution.png"), 
       plot = p, width = 8, height = 6, dpi = 300)

#To show percentages in Y axis
peak_data <- peak_data %>%
       group_by(strategy) %>%
       mutate(percent = length / sum(length))
p1 <- ggplot(peak_data, aes(x=length, fill=strategy, y=..count../sum(..count..))) +
       geom_histogram(alpha=0.5, position="identity", bins=50) +
     theme_minimal() +
      labs(title="Comparison of Peak Length Distribution (Percentage)", 
                       x="Peak Length", y="Percentage") +
       scale_y_continuous(labels=scales::percent) +
     scale_fill_manual(values=c("skyblue", "salmon")) +
       facet_wrap(~strategy, scales = "free_y")
ggsave(filename = file.path(output_folder, "peak_length_distribution_Percentage.png"), 
       plot = p1, width = 8, height = 6, dpi = 300)


#合并的percentages两张图的结果
p_combined <- ggplot(peak_data, aes(x = length, fill = strategy, y = after_stat(count / sum(count)))) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 50) +
  theme_minimal() +
  labs(title = "Comparison of Peak Length Distribution (Percentage)", 
       x = "Peak Length", y = "Percentage") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("skyblue", "salmon"))
ggsave(filename = file.path(output_folder, "peak_length_distribution_Percentage_combined.png"), 
       plot = p_combined, width = 8, height = 6, dpi = 300)




























