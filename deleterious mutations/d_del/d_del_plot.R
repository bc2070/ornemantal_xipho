#!/usr/bin/env Rscript


setwd("/public4/group_crf/home/b20gengjlin3/Xiphophorus/NEW_DN_CDS/gradient/dN")



required_packages <- c("ggplot2", "ggpubr", "dplyr", "tidyr", "purrr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}


file_list <- list.files(pattern = "_cal_filtered.bed$")

if (length(file_list) < 1) {
  stop("no found  *_cal_filtered.bed ")
}

cat("reading ...\n")

all_data <- file_list %>%
  map_df(function(f) {

    df_temp <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    d_h <- as.numeric(df_temp[[5]])
    d_m <- as.numeric(df_temp[[6]])
    

    data.frame(
      Value = d_h - d_m,
      Group = gsub("_cal_filtered.bed", "", f)     )
  }) %>%
  drop_na(Value)

cat("Carry out statistical tests and select the significant items...\n")


stat_test <- compare_means(Value ~ Group, data = all_data, method = "wilcox.test")


sig_comparisons_df <- stat_test %>% filter(p < 0.05)


sig_list <- purrr::map2(sig_comparisons_df$group1, sig_comparisons_df$group2, ~c(.x, .y))


cat("chart...\n")

p <- ggplot(all_data, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.5, outlier.alpha = 0.3, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.6) +
  theme_bw() +
  labs(
    title = "Significant Genomic Differences (dH - dM)",
    subtitle = "Positive: closer to maculatus | Negative: closer to hellerii\n(Only showing p < 0.05 comparisons)",
    x = "Genomic Regions",
    y = "Index (d_hellerii - d_maculatus)"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

if (length(sig_list) > 0) {
  p <- p + stat_compare_means(
    comparisons = sig_list,
    method = "wilcox.test",
    label = "p.signif",           
    step_increase = 0.12,         
    tip_length = 0.01
  )
} else {
  cat("Note: There were no significant differences between any groups (p < 0.05), so no horizontal lines will be displayed.\n")
}


output_name <- "wilcox_sig_plot"
#ggsave(paste0(output_name, ".png"), p, width = 10, height = 7, dpi = 300)
ggsave(paste0(output_name, ".pdf"), p, width = 10, height = 7)


print(sig_comparisons_df[, c("group1", "group2", "p", "p.signif")])















required_packages <- c("ggplot2", "ggpubr", "dplyr", "tidyr", "purrr", "stringr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}


file_list <- list.files(pattern = "_cal_filtered.bed$")
if (length(file_list) < 1) stop("no found  *_cal_filtered.bed")

cat("reading ...\n")


all_data <- file_list %>%
  map_df(function(f) {
    df_temp <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    

    group_name <- gsub("_cal_filtered.bed", "", f)
    

    midpoint <- NA
    if (grepl("minor_", group_name)) {
      nums <- as.numeric(unlist(str_extract_all(group_name, "\\d+")))
      if (length(nums) >= 2) {
        midpoint <- mean(nums)
      }
    }
    

    data.frame(
      Value = as.numeric(df_temp[[5]]) - as.numeric(df_temp[[6]]),
      Group = group_name,
      Midpoint = midpoint,
      IsMinor = !is.na(midpoint)
    )
  }) %>%
  drop_na(Value)


non_minor_groups <- unique(all_data$Group[!all_data$IsMinor])
minor_groups_ordered <- all_data %>%
  filter(IsMinor) %>%
  distinct(Group, Midpoint) %>%
  arrange(Midpoint) %>%
  pull(Group)

all_data$Group <- factor(all_data$Group, levels = c(non_minor_groups, minor_groups_ordered))





stat_test <- compare_means(Value ~ Group, data = all_data, method = "wilcox.test")
sig_list <- stat_test %>% 
  filter(p < 0.05) %>% 
  purrr::map2(.$group1, .$group2, ~c(.x, .y))


minor_data <- all_data %>% filter(IsMinor)
if (nrow(minor_data) > 0) {

  cor_res <- cor.test(minor_data$Midpoint, minor_data$Value)
  lm_res <- lm(Value ~ Midpoint, data = minor_data)
  
  anno_text <- sprintf("Trend (Minor groups):\nPearson R = %.3f\np-value = %.3e", 
                       cor_res$estimate, cor_res$p.value)
  cat(anno_text, "\n")
}



p <- ggplot(all_data, aes(x = Group, y = Value, fill = IsMinor)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.5, outlier.alpha = 0.3, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  scale_fill_manual(values = c("TRUE" = "#69b3a2", "FALSE" = "#404080")) +   theme_bw() +
  labs(
    title = "Genomic Differences (dH - dM) with Trend Analysis",
    subtitle = "Minor series: 0-10 to 90-100 (grouped by Midpoint)",
    x = "Genomic Regions",
    y = "Index (d_hellerii - d_maculatus)"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )


if (nrow(minor_data) > 0) {

  minor_indices <- which(levels(all_data$Group) %in% minor_groups_ordered)
  

  p <- p + geom_smooth(
    data = subset(all_data, IsMinor),
    aes(x = as.numeric(Group), y = Value),
    method = "lm", color = "red", se = TRUE, inherit.aes = FALSE
  )
  

  p <- p + annotate("text", x = 1, y = max(all_data$Value)*0.9, 
                    label = anno_text, hjust = 0, color = "red", fontface = "italic")
}


print(p)

# ggsave("genomic_trend_analysis.pdf", width = 12, height = 7)
















required_packages <- c("ggplot2", "dplyr", "tidyr", "purrr", "stringr", "stats")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}


# setwd("/public4/group_crf/home/b20gengjlin3/Xiphophorus/NEW_DN_CDS/gradient/dN")


target_order <- c(
  "minor_0","minor_0-10", "minor_10-20", "minor_20-30", "minor_30-40", "minor_40-50", 
  "minor_50-60", "minor_60-70", "minor_70-80", "minor_80-90", "minor_90-100",
  "non_HF", "HF", "long_HF"
)


file_list <- list.files(pattern = "_cal_filtered.bed$")
if (length(file_list) < 1) stop("no found *_cal_filtered.bed ")

)
all_data <- file_list %>%
  map_df(function(f) {
    df_temp <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    group_name <- gsub("_cal_filtered.bed", "", f)
    data.frame(
      Value = as.numeric(df_temp[[5]]) - as.numeric(df_temp[[6]]),
      Group = group_name,
      IsMinor = grepl("minor_", group_name)
    )
  }) %>%
  drop_na(Value)


existing_groups <- unique(all_data$Group)
final_levels <- target_order[target_order %in% existing_groups]
others <- setdiff(existing_groups, final_levels)
all_data$Group <- factor(all_data$Group, levels = c(final_levels, others))

summary_stats <- all_data %>%
  group_by(Group) %>%
  summarise(
    Count = n(),
    Mean = mean(Value),
    Median = median(Value),
    SD = sd(Value)
  )
print(as.data.frame(summary_stats))

pairwise_results <- pairwise.wilcox.test(all_data$Value, all_data$Group, p.adjust.method = "fdr")
print(pairwise_results)

major_data <- all_data %>% filter(Group == "non_HF") %>% pull(Value)
minor_heavy_data <- all_data %>% 
  filter(Group %in% c("minor_50-60", "minor_60-70", "minor_70-80", "minor_80-90", "minor_90-100")) %>% 
  pull(Value)

if(length(major_data) > 0 && length(minor_heavy_data) > 0) {
  wt <- wilcox.test(major_data, minor_heavy_data)
  cat("Major (non_HF) Mean:", mean(major_data), "\n")
  cat("Minor (>50%) Mean:", mean(minor_heavy_data), "\n")
  cat("Mann-Whitney U statistic (W):", wt$statistic, "\n")
  cat("P-value:", wt$p.value, "\n")
} else {
  cat("data lost\n")
}


cat("\n chart...\n")
p <- ggplot(all_data, aes(x = Group, y = Value, fill = IsMinor)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.5, outlier.alpha = 0.3, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  scale_fill_manual(values = c("TRUE" = "#69b3a2", "FALSE" = "#404080")) + 
  theme_bw() +
  labs(title = "Genomic Differences (dH - dM)", x = "Genomic Regions", y = "Index (d_hellerii - d_maculatus)") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())

print(p)


 ggsave("genomic_boxplot_custom_order.pdf", width = 10, height = 6)
 