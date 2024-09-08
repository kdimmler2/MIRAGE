setwd('/scratch.global/dimml002/PaintedPatches/results')

library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(purrr)

######################################################################################
#   Reading and cleaning the data   ##################################################
######################################################################################

targets <- read.table('/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Sam_034_targets_liftover.txt')
#targets <- targets |> filter(V1 == chromosome)

df_discrete <- read.table('0.3_set1_AF_across_chrom_discrete.tsv', header=TRUE)
df_discrete <- df_discrete |> 
	filter(Masked == '-9') |>
	rename(AF = MAF) |>
	mutate(MAF = pmin(AF, 1 - AF))

#df_discrete <- df_discrete %>%
#  mutate(
#    target_site = ifelse(
#      ID %in% targets$V3, 'true_target',
#      ifelse(POS %in% targets$V2, 'target_position', 'flanking')
#    )
#  )

df_discrete <- df_discrete %>%
  mutate(
    target_site = ifelse(POS %in% targets$V2, 'target_position', 'flanking')
    )

df_discrete$Original <- as.numeric(df_discrete$Original)
df_discrete$Imputed <- as.numeric(df_discrete$Imputed)
df_discrete <- na.omit(df_discrete)

######################################################################################
#   Imputation Quality Score vs. MAF   ###############################################
######################################################################################

#Remove variants where not all genotypes are present
df_for_IQS <- df_discrete |>
  group_by(POS) |>
  filter(n_distinct(Original) >= 3, n_distinct(Imputed) >= 3) |>
  ungroup()

# Define the function to calculate IQS for each POS group with debug prints
calculate_IQS <- function(df) {
  df |>
    group_by(POS) |>
    group_split() |>
    map_dfr(~ {
      group_data <- .x
      
      # Print debug information for each group
      cat("Processing POS group:", unique(group_data$POS), "\n")
      print(group_data)
      
      true_genotypes <- group_data$Original
      imputed_genotypes <- group_data$Imputed
      
      # Ensure group_data has valid POS
      pos_value <- unique(group_data$POS)
      if (length(pos_value) == 0) {
        cat("Invalid POS for group, skipping.\n")
        return(NULL)
      }
      
      # Check for groups with insufficient unique values
      if (length(unique(true_genotypes)) == 1 || length(unique(imputed_genotypes)) == 1) {
        return(data.frame(POS = pos_value, IQS = NA))
      }
      
      # Create the 3x3 table using table
      contingency_table <- table(true_genotypes, imputed_genotypes)
      
      # Check if the contingency table is valid
      if (sum(contingency_table) == 0) {
        return(data.frame(POS = pos_value, IQS = NA))
      }
      
      # Convert table to matrix and ensure it's 3x3
      n <- matrix(0, nrow = 3, ncol = 3)
      n[1:dim(contingency_table)[1], 1:dim(contingency_table)[2]] <- contingency_table
      
      # Calculate observed proportion of agreement (Po)
      Po <- sum(diag(n)) / sum(n)
      
      # Calculate marginal frequencies
      true_marginals <- colSums(n)
      imputed_marginals <- rowSums(n)
      
      # Calculate chance agreement (Pc)
      Pc <- sum(true_marginals * imputed_marginals) / sum(n)^2
      
      # Calculate IQS
      IQS <- if (Pc == 1) 1 else (Po - Pc) / (1 - Pc)
      
      # Return the result for this group
      return(data.frame(POS = pos_value, IQS = IQS, MAF = unique(group_data$MAF), target_site = unique(group_data$target_site)))
    })
}


results_df <- calculate_IQS(df_for_IQS)

# Remove rows with NA values
cleaned_results_df <- results_df |>
  na.omit()

write.table(cleaned_results_df, '0.3_set1_IQS_across_chrom.tsv', sep='\t', row.names = FALSE, quote = FALSE)

#true_targets_df <- cleaned_results_df |> filter(target_site == "true_target")
flanking_df <- cleaned_results_df |> filter(target_site == "flanking")
targets_df <- cleaned_results_df |> filter(target_site == "target_position")

#print(true_targets_df)

#bins <- c(0, 0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)
#labels <- c("[0–0.005]", "[0.005–0.01]", "[0.01–0.02]", "[0.02–0.05]", "[0.05–0.10]", "[0.10–0.15]", "[0.15–0.20]", "[0.20–0.25]", "[0.25–0.30]", "[0.30–0.35]", "[0.35–0.40]", "[0.40–0.45]", "[0.45–0.50]")

# Define bins with size 0.005
bins <- seq(0, 0.5, by = 0.025)

# Create labels for the bins
labels <- paste0("[", bins[-length(bins)], "–", bins[-1], "]")

#true_targets_binned_maf_df <- true_targets_df |>
#  mutate(MAF_bin = cut(MAF, breaks = bins, labels = labels, include.lowest = TRUE))
#
## Calculate the average IQS values for each MAF bin
#true_targets_avg_IQS <- true_targets_binned_maf_df |>
#  group_by(MAF_bin) |>
#  summarise(!!set := mean(IQS, na.rm = TRUE))
#
#true_targets_avg_IQS <- true_targets_avg_IQS |> mutate('target_site' = 'true_target')

flanking_binned_maf_df <- flanking_df |>
  mutate(MAF_bin = cut(MAF, breaks = bins, labels = labels, include.lowest = TRUE))

# Calculate the average IQS values for each MAF bin
flanking_avg_IQS <- flanking_binned_maf_df |>
  group_by(MAF_bin) |>
  summarise(set1 = mean(IQS, na.rm = TRUE))

flanking_avg_IQS <- flanking_avg_IQS |> mutate('target_site' = 'flanking')

targets_binned_maf_df <- targets_df |>
  mutate(MAF_bin = cut(MAF, breaks = bins, labels = labels, include.lowest = TRUE))

# Calculate the average IQS values for each MAF bin
targets_avg_IQS <- targets_binned_maf_df |>
  group_by(MAF_bin) |>
  summarise(set1 = mean(IQS, na.rm = TRUE))

targets_avg_IQS <- targets_avg_IQS |> mutate('target_site' = 'target_position')

merged_df1 <- rbind(flanking_avg_IQS, targets_avg_IQS)

merged_df_0.2 <- merged_df_0.2 %>%
  left_join(merged_df1 %>% select(MAF_bin, target_site, set3), by = c("MAF_bin", "target_site"))

merged_df_0.2 <- merged_df_0.2 |> select(MAF_bin, target_site, set1, set2, set3)

#merged_df_0.05 |> merged_df_0.05 |> mutate(set3 = merged_df$set3[match(MAF_bin)])


###################################
merged_df_0.2 <- merged_df1
###################################


write.table(merged_df_0.2, '0.2_across_chrom_IQS.tsv', sep='\t', row.names = FALSE, quote = FALSE)


#df_list <- list(merged_df1, merged_df12)

#combined_df <- reduce(df_list, full_join, by = c('MAF_bin', 'target_site'))

write.table(merged_df, output_file, sep='\t', row.names = FALSE, quote = FALSE)
