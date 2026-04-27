
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr") 


library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

setwd("/public4/group_crf/home/b20gengjlin3/Xiphophorus/NEW_DN_CDS/gradient/cds")


read_bed_as_df <- function(filepath) {
  df <- read.table(filepath, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  if (ncol(df) >= 3) {
    colnames(df)[1:3] <- c("Chrom", "Start", "End")
    df$Chrom <- as.character(df$Chrom) # 强制Chrom列为字符型
  } else {
    stop(paste("BED file", filepath, "does not have enough columns (expected at least 3)."))
  }
  return(df)
}



merge_overlapping_regions_df <- function(df) {
  if (nrow(df) == 0) {
    return(data.frame(Chrom = character(), Start = integer(), End = integer(), stringsAsFactors = FALSE))
  }
  
  df <- df %>% arrange(Chrom, Start) 
  
  merged_regions <- list()
  current_chrom <- df$Chrom[1]
  current_start <- df$Start[1]
  current_end <- df$End[1]
  
  for (i in 2:nrow(df)) {
    if (df$Chrom[i] == current_chrom && df$Start[i] <= current_end) {

      current_end <- max(current_end, df$End[i])
    } else {

      merged_regions <- c(merged_regions, list(c(current_chrom, current_start, current_end)))
      current_chrom <- df$Chrom[i]
      current_start <- df$Start[i]
      current_end <- df$End[i]
    }
  }

  merged_regions <- c(merged_regions, list(c(current_chrom, current_start, current_end)))
  

  merged_df <- as.data.frame(do.call(rbind, merged_regions), stringsAsFactors = FALSE)
  colnames(merged_df) <- c("Chrom", "Start", "End")
  merged_df$Start <- as.integer(merged_df$Start)
  merged_df$End <- as.integer(merged_df$End)
  
  return(merged_df)
}

calculate_total_overlap_df <- function(df1_regions, df2_regions) {
  total_overlap <- 0
  
  if (nrow(df1_regions) == 0 || nrow(df2_regions) == 0) {
    return(0)
  }
  

  df1_regions <- df1_regions %>% arrange(Chrom, Start)
  df2_regions <- df2_regions %>% arrange(Chrom, Start)
  

  chromosomes <- unique(c(df1_regions$Chrom, df2_regions$Chrom))
  
  for (chrom in chromosomes) {
    chrom_df1 <- df1_regions %>% filter(Chrom == chrom)
    chrom_df2 <- df2_regions %>% filter(Chrom == chrom)
    
    if (nrow(chrom_df1) == 0 || nrow(chrom_df2) == 0) {
      next 
    }
    

    idx1 <- 1
    idx2 <- 1
    while (idx1 <= nrow(chrom_df1) && idx2 <= nrow(chrom_df2)) {
      interval1_start <- chrom_df1$Start[idx1]
      interval1_end <- chrom_df1$End[idx1]
      interval2_start <- chrom_df2$Start[idx2]
      interval2_end <- chrom_df2$End[idx2]
      

      overlap_start <- max(interval1_start, interval2_start)
      overlap_end <- min(interval1_end, interval2_end)
      

      if (overlap_start < overlap_end) {
        total_overlap <- total_overlap + (overlap_end - overlap_start)
      }
      
      if (interval1_end < interval2_end) {
        idx1 <- idx1 + 1
      } else {
        idx2 <- idx2 + 1
      }
    }
  }
  return(total_overlap)
}


calculate_target_region_non_excluded_length <- function(target_bed_filepath, merged_exclusion_regions_df, desired_chroms) {

  target_regions_raw <- read_bed_as_df(target_bed_filepath)

  target_regions_filtered <- target_regions_raw %>%
    filter(Chrom %in% desired_chroms)
  

  if (nrow(target_regions_filtered) == 0) {
    return(0)
  }
  

  merged_target_regions <- merge_overlapping_regions_df(target_regions_filtered)
  

  total_target_length <- sum(merged_target_regions$End - merged_target_regions$Start)
  

  overlap_with_exclusion <- calculate_total_overlap_df(merged_target_regions, merged_exclusion_regions_df)
  

  non_excluded_length <- total_target_length - overlap_with_exclusion
  
  return(non_excluded_length)
}



output_csv_file <- "output2.csv" 
chrom_txt_file <- "chrom.txt"   
target_regions_dir <- "target_region" 
inversion_bed <- "inversion_pos.bed"
N_bed <- "N_pos.bed" 
low_posterior_bed <- "low_posterior.bed"

cds_len_file <- "cds_len.txt"     
gene_pos_file <- "gene_pos.gff"   


desired_chroms <- as.character(1:24)


num_resamples <- 1000 


target_region_files_fullpaths <- list.files(target_regions_dir, pattern = "\\.bed$", full.names = TRUE)
if (length(target_region_files_fullpaths) == 0) {
  stop(paste0("No .bed files found in the target regions directory: ", target_regions_dir, ". Please check the path and file existence."))
}
target_file_prefixes <- tools::file_path_sans_ext(basename(target_region_files_fullpaths))
cat("Automatically detected target region prefixes:\n")
print(target_file_prefixes)





genes_data <- read.csv(output_csv_file, header = TRUE, stringsAsFactors = FALSE)
genes_data$CDS_Length <- as.numeric(genes_data$CDS_Length) 


genes_long <- genes_data %>%
  separate_rows(Labels, sep = ",") %>%
  filter(!is.na(Labels) & Labels != "")


chrom_info <- read.table(chrom_txt_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                         colClasses = c("character", "integer", "integer", "character"))
chrom_info$Length <- chrom_info$End - chrom_info$Start + 1
original_chrom_info_rows <- nrow(chrom_info)
chrom_info <- chrom_info %>% filter(Chrom %in% desired_chroms)
cat(sprintf("Filtered out %d rows from chrom.txt (retained %d rows for desired chroms).\n",
            original_chrom_info_rows - nrow(chrom_info), nrow(chrom_info)))
if(nrow(chrom_info) == 0) stop("No desired chromosomes found in chrom.txt after filtering.")



print("Reading and merging exclusion regions...")
all_exclusion_regions_raw <- bind_rows(
  read_bed_as_df(inversion_bed),
  read_bed_as_df(N_bed),
  read_bed_as_df(low_posterior_bed)
)
original_exclusion_rows <- nrow(all_exclusion_regions_raw)
filtered_exclusion_regions <- all_exclusion_regions_raw %>%
  filter(Chrom %in% desired_chroms)
cat(sprintf("Filtered out %d rows from exclusion regions (retained %d rows for desired chroms).\n",
            original_exclusion_rows - nrow(filtered_exclusion_regions), nrow(filtered_exclusion_regions)))
merged_exclusion_regions_df <- merge_overlapping_regions_df(filtered_exclusion_regions)
print("Merged exclusion regions (data.frame) created:")
print(head(merged_exclusion_regions_df))


exclusion_length_by_chrom <- merged_exclusion_regions_df %>%
  group_by(Chrom) %>%
  summarise(ExcludedLength = sum(End - Start)) %>%
  ungroup()


chrom_lengths_non_excluded <- chrom_info %>%
  left_join(exclusion_length_by_chrom, by = "Chrom") %>%
  mutate(ExcludedLength = ifelse(is.na(ExcludedLength), 0, ExcludedLength)) %>%
  mutate(NonExcludedLength = Length - ExcludedLength) %>%
  select(Chrom, TotalChromLength = Length, ExcludedLength, NonExcludedLength)
print("Chromosome lengths after excluding masked regions:")
print(chrom_lengths_non_excluded)



print(paste0("\nCalculating non-excluded lengths for target regions in '", target_regions_dir, "'..."))

if (length(target_region_files_fullpaths) == 0) {
  warning(paste0("No .bed files found in the target regions directory: ", target_regions_dir))
} else {
  target_regions_non_excluded_lengths <- data.frame(
    RegionName = character(),
    NonExcludedLength = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (file_path in target_region_files_fullpaths) {
    region_name <- tools::file_path_sans_ext(basename(file_path))
    cat(paste0("Processing target region: ", region_name, "... "))
    
    non_excluded_len <- calculate_target_region_non_excluded_length(file_path, merged_exclusion_regions_df, desired_chroms)
    
    target_regions_non_excluded_lengths <- bind_rows(
      target_regions_non_excluded_lengths,
      data.frame(RegionName = region_name, NonExcludedLength = non_excluded_len, stringsAsFactors = FALSE)
    )
    cat(paste0("Non-excluded length: ", non_excluded_len, "\n"))
  }
  
  print("\nSummary of Non-Excluded Lengths for Target Regions:")
  print(target_regions_non_excluded_lengths)
}




print("\nCalculating CDS densities for each category with resampling...")

category_densities_summary <- data.frame(Category = character(),
                                         Mean_Resampled_CDS_Length_Sum = numeric(),
                                         Non_Excluded_Region_Length = numeric(),
                                         Mean_Density = numeric(),
                                         stringsAsFactors = FALSE)
category_cds_densities_distributions <- list()
category_raw_cds_lengths <- list()

for (prefix in target_file_prefixes) {
  cat(sprintf("Processing category: %s\n", prefix))
  
  current_category_genes <- genes_long %>% filter(Labels == prefix)
  cds_lengths_in_category <- current_category_genes$CDS_Length[!is.na(current_category_genes$CDS_Length)]
  
  if (length(cds_lengths_in_category) == 0) {
    cat(sprintf("  No valid CDS lengths found for category %s. Skipping density calculation.\n", prefix))
    next
  }
  
  category_raw_cds_lengths[[prefix]] <- cds_lengths_in_category
  
  target_file_path <- file.path(target_regions_dir, paste0(prefix, ".bed"))
  non_excluded_region_length <- NA
  
  if (!file.exists(target_file_path)) {
    cat(sprintf("  Warning: Target region file %s not found. Non-Excluded_Region_Length will be NA.\n", target_file_path))
  } else {
    target_regions_df_raw <- read_bed_as_df(target_file_path)
    original_target_rows <- nrow(target_regions_df_raw)
    target_regions_df <- target_regions_df_raw %>%
      filter(Chrom %in% desired_chroms)
    
    if (original_target_rows > nrow(target_regions_df)) {
      cat(sprintf("  Filtered out %d rows from target regions for %s (retained %d rows for desired chroms).\n",
                  original_target_rows - nrow(target_regions_df), prefix, nrow(target_regions_df)))
    }
    
    if(nrow(target_regions_df) == 0) {
      cat(sprintf("  No target regions left for desired chroms in category %s. Non-Excluded_Region_Length will be NA.\n", prefix))
    } else {
      merged_target_regions_df <- merge_overlapping_regions_df(target_regions_df)
      total_target_length <- sum(merged_target_regions_df$End - merged_target_regions_df$Start)
      excluded_from_target_length <- calculate_total_overlap_df(merged_target_regions_df, merged_exclusion_regions_df)
      non_excluded_region_length <- total_target_length - excluded_from_target_length
    }
  }
  
  if (is.na(non_excluded_region_length) || non_excluded_region_length <= 0) {
    cat(sprintf("  Non-Excluded_Region_Length for category %s is invalid or zero. Skipping density distribution calculation.\n", prefix))
    next
  }
  
  resampled_densities <- numeric(num_resamples)
  resampled_sums <- numeric(num_resamples)
  
  for (k in 1:num_resamples) {
    resampled_cds_lengths <- sample(cds_lengths_in_category, size = length(cds_lengths_in_category), replace = TRUE)
    resampled_sums[k] <- sum(resampled_cds_lengths)
    resampled_densities[k] <- resampled_sums[k] / non_excluded_region_length
  }
  
  mean_resampled_cds_sum <- mean(resampled_sums)
  mean_density <- mean(resampled_densities)
  
  category_densities_summary <- rbind(category_densities_summary, data.frame(
    Category = prefix,
    Mean_Resampled_CDS_Length_Sum = mean_resampled_cds_sum,
    Non_Excluded_Region_Length = non_excluded_region_length,
    Mean_Density = mean_density,
    stringsAsFactors = FALSE
  ))
  
  category_cds_densities_distributions[[prefix]] <- resampled_densities
}

print("\nCategory CDS Densities (Averaged from 1000 Resamples):")
print(category_densities_summary)




all_resampled_densities_for_plot <- data.frame(
  Category = character(),
  Density = numeric(),
  stringsAsFactors = FALSE
)

for (category_name in names(category_cds_densities_distributions)) {
  all_resampled_densities_for_plot <- rbind(all_resampled_densities_for_plot, data.frame(
    Category = category_name,
    Density = category_cds_densities_distributions[[category_name]],
    stringsAsFactors = FALSE
  ))
}
custom_order <- c("minor_0-10", "minor_10-20", "minor_20-30", "minor_30-40","minor_40-50","minor_50-60","minor_60-70","minor_70-80","minor_80-90","minor_90-100","HF","long_HF") 
all_resampled_densities_for_plot$Category <- factor(all_resampled_densities_for_plot$Category,
                                                    levels = custom_order)

# if (nrow(all_resampled_densities_for_plot) > 0) {
#   all_resampled_densities_for_plot$Category <- factor(all_resampled_densities_for_plot$Category,
#                                                       levels = names(sort(tapply(all_resampled_densities_for_plot$Density,
#                                                                                  all_resampled_densities_for_plot$Category,
#                                                                                  median), decreasing = FALSE)))
# } else {
#   cat("\nNo data available for plotting CDS density boxplot.\n")
# }

if (nrow(all_resampled_densities_for_plot) > 0) {
  p_boxplot <- ggplot(all_resampled_densities_for_plot, aes(x = Category, y = Density, fill = Category)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(color = "black", size = 0.5, alpha = 0.3, width = 0.2) +
    geom_hline(yintercept = total_cds_density, linetype = "dashed", color = "red", size = 1) +
    annotate("text", x = Inf, y = total_cds_density, label = paste0("Global CDS Density: ", sprintf("%.2e", total_cds_density)),
             vjust = -0.5, hjust = 1.1, color = "red", size = 4) +
    labs(title = paste0("CDS Density Distribution Across Categories (", num_resamples, " Resamples)"),
         x = "Gene Category",
         y = "CDS Density") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  
  print(p_boxplot)
  
  # ggsave("cds_density_boxplot.png", plot = p_boxplot, width = 12, height = 8, dpi = 300)
  # cat(paste0("\nCDS density boxplot saved as cds_density_boxplot.png (", num_resamples, " resamples)\n"))
}




perform_density_distribution_difference_test <- function(category1_name, category2_name, all_density_distributions_list, plot_results = TRUE) {
  
  data1 <- all_density_distributions_list[[category1_name]]
  data2 <- all_density_distributions_list[[category2_name]]
  
  if (length(data1) < 2 || length(data2) < 1) {
    warning(sprintf("Not enough density samples for %s (n=%d) or %s (n=%d). Skipping comparison.\n",
                    category1_name, length(data1), category2_name, length(data2)))
    return(list(
      category1 = category1_name,
      category2 = category2_name,
      mean_dist1 = NA,
      sd_dist1 = NA,
      median_dist2 = NA,
      p_value_raw = NA,
      p_value_two_tailed = NA,
      message = "Insufficient density samples"
    ))
  }
  
  mu1 <- mean(data1)
  sd1 <- sd(data1)
  
  if (sd1 == 0) {
    warning(paste("Standard deviation for density distribution of", category1_name, "is zero. Skipping comparison."))
    return(list(
      category1 = category1_name,
      category2 = category2_name,
      mean_dist1 = mu1,
      sd_dist1 = sd1,
      median_dist2 = median(data2),
      p_value_raw = NA,
      p_value_two_tailed = NA,
      message = "SD of density distribution 1 is zero"
    ))
  }
  
  median_data2 <- median(data2)
  p_value_raw <- pnorm(median_data2, mean = mu1, sd = sd1)
  p_value_two_tailed <- min(p_value_raw, 1 - p_value_raw) * 2
  
  if (plot_results) {
    cat(sprintf("\n--- Density Distribution Difference Test: %s (Reference) vs %s (Test) ---\n", category1_name, category2_name))
    cat(sprintf("Reference Category 1 (%s) Density Dist - Mean: %.4e, SD: %.4e\n", category1_name, mu1, sd1))
    cat(sprintf("Test Category 2 (%s) Density Dist - Median: %.4e\n", category2_name, median_data2))
    cat(sprintf("P-value (median of %s density dist in normal dist of %s - two-tailed): %.4f\n", category2_name, category1_name, p_value_two_tailed))
    
    plot_df <- data.frame(Value = c(data1, data2),
                          Category = c(rep(category1_name, length(data1)), rep(category2_name, length(data2))))
    
    p <- ggplot(plot_df, aes(x = Value, fill = Category, color = Category)) +
      geom_density(alpha = 0.4) +
      geom_vline(xintercept = median_data2, linetype = "dashed", color = "blue", size = 0.8) +
      annotate("text", x = median_data2, y = max(density(plot_df$Value)$y) * 0.9,
               label = paste("Median of", category2_name, ":", sprintf("%.2e", median_data2)),
               color = "blue", hjust = -0.1) +
      stat_function(fun = dnorm, args = list(mean = mu1, sd = sd1),
                    color = "red", linetype = "dotted", size = 1.2) +
      annotate("text", x = mu1, y = max(density(plot_df$Value)$y) * 0.8,
               label = paste("Normal fit for", category1_name),
               color = "red", hjust = -0.1) +
      labs(title = paste("Density Distribution Comparison:", category1_name, "vs", category2_name, "(P-value:", round(p_value_two_tailed, 4), ")"),
           x = "Resampled Density", y = "Density") +
      theme_minimal() +
      theme(legend.position = "bottom")
    print(p)
  }
  
  return(list(
    category1 = category1_name,
    category2 = category2_name,
    mean_dist1 = mu1,
    sd_dist1 = sd1,
    median_dist2 = median_data2,
    p_value_raw = p_value_raw,
    p_value_two_tailed = p_value_two_tailed
  ))
}

cat("\n--- Performing All-Pairs Density Distribution Difference Tests ---\n")

valid_categories_for_density_dist <- names(category_cds_densities_distributions)[sapply(category_cds_densities_distributions, function(x) length(x) >= 2)]

if (length(valid_categories_for_density_dist) < 2) {
  stop("Not enough valid categories (at least 2 with sufficient density samples) to perform all-pairs comparisons.")
}

all_combinations <- expand.grid(cat1_name = valid_categories_for_density_dist, cat2_name = valid_categories_for_density_dist, stringsAsFactors = FALSE) %>%
  filter(cat1_name != cat2_name) %>%
  rowwise() %>%
  filter(cat1_name < cat2_name) %>%
  ungroup()

all_test_results <- list()
pb <- txtProgressBar(min = 0, max = nrow(all_combinations), style = 3)

for (i in 1:nrow(all_combinations)) {
  cat1 <- all_combinations$cat1_name[i]
  cat2 <- all_combinations$cat2_name[i]
  
  result <- perform_density_distribution_difference_test(cat1, cat2, category_cds_densities_distributions, plot_results = FALSE)
  all_test_results[[i]] <- result
  setTxtProgressBar(pb, i)
}
close(pb)

results_df <- do.call(rbind, lapply(all_test_results, function(res) {
  as.data.frame(res, stringsAsFactors = FALSE)
})) %>%
  filter(!is.na(p_value_two_tailed))

results_df$p_adjusted <- p.adjust(results_df$p_value_two_tailed, method = "fdr")

alpha <- 0.05

significant_results <- results_df %>%
  filter(p_adjusted < alpha) %>%
  arrange(p_adjusted)

cat(sprintf("\n\n--- All-Pairs Density Distribution Comparison Results (%d valid comparisons) ---\n", nrow(results_df)))
if (nrow(significant_results) > 0) {
  cat(sprintf("\nSignificant differences found (adjusted p-value < %.2f):\n", alpha))
  print(significant_results)
  
  cat("\n--- Plotting significant density distribution differences (optional) ---\n")
  max_plots_to_show <- 5
  if (nrow(significant_results) > max_plots_to_show) {
    cat(sprintf("  (Showing top %d significant plots only)\n", max_plots_to_show))
    significant_results_to_plot <- significant_results %>% head(max_plots_to_show)
  } else {
    significant_results_to_plot <- significant_results
  }
  
  for (i in 1:nrow(significant_results_to_plot)) {
    cat1_sig <- significant_results_to_plot$category1[i]
    cat2_sig <- significant_results_to_plot$category2[i]
    perform_density_distribution_difference_test(cat1_sig, cat2_sig, category_cds_densities_distributions, plot_results = TRUE)
  }
  
} else {
  cat(sprintf("\nNo significant density distribution differences found at adjusted p-value < %.2f.\n", alpha))
}

cat("\n--- All Density Distribution Comparison Results (including non-significant) ---\n")
print(results_df %>% arrange(p_adjusted))


cat("\n\n--- Manual Comparison of Specific Categories (e.g., for detailed plot) ---\n")
cat("Available categories for comparison:\n")
print(names(category_cds_densities_distributions))

cat1_manual <- "HF"
cat2_manual <- "non_HF"

# cat1_manual <- readline("Enter the name of the reference category for manual comparison: ")
# cat2_manual <- readline("Enter the name of the test category for manual comparison: ")

if (cat1_manual %in% names(category_cds_densities_distributions) && cat2_manual %in% names(category_cds_densities_distributions)) {
  perform_density_distribution_difference_test(cat1_manual, cat2_manual, category_cds_densities_distributions, plot_results = TRUE)
} else {
  cat(sprintf("\nError: One or both categories '%s' and '%s' are not valid or do not have enough density samples for manual comparison.\n", cat1_manual, cat2_manual))
}