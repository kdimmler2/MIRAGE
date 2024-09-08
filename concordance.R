library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(purrr)
library(patchwork)

setwd("/scratch.global/dimml002/PaintedPatches")

######################################################################################
#   Reading and cleaning the data set1   #############################################
######################################################################################

targets <- read.table('/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Sam_034_targets_liftover.txt', header=FALSE)
#targets <- targets |> filter(V1 %in% c('chr12', 'chr23', 'chr27'))

##########################
####   Masking 0.05   ####
##########################

df_dosage_0.05 <- read.table('0.05_set1_AF_across_chrom_dosage.tsv', header=TRUE)
R2s_0.05 <- read.table('0.05_set1_DR2_AR2_across_chrom_dosage.tsv', header=FALSE)
R2s_0.05 <- R2s_0.05 |> rename(ID = V3,
                             AR2 = V4,
                             DR2 = V5)
df_flipped_0.05 <- df_dosage_0.05 |> filter(MAF >= 0.5)
df_dosage_0.05 <- df_dosage_0.05 |> filter(MAF <= 0.5)
df_dosage_0.05 <- df_dosage_0.05 |> filter(Masked == '-9')
# df <- df |>
#   rename(AF = MAF) |>
#   mutate(MAF = pmin(AF, 1 - AF))

df_0.05 <- df_dosage_0.05 |>
  left_join(R2s_0.05, by = "ID")

df_0.05 <- df_0.05 |> mutate(target_site = ifelse(POS %in% targets$V2, 'target', 'flanking'))

df_0.05 <- df_0.05 |> mutate(Discordance = abs(Imputed - Original),
                              Concordance = 1 - (abs(Imputed - Original)),
                              genotype = case_when(
                              Original == 2 ~ 'alt',
                              Original == 1 ~ 'het',
                              Original == 0 ~ 'ref',))

##########################
####   Masking 0.1   #####
##########################

df_dosage_0.1 <- read.table('0.1_set3_AF_across_chrom_dosage.tsv', header=TRUE)
R2s_0.1 <- read.table('0.1_set3_DR2_AR2_across_chrom_dosage.tsv', header=FALSE)
R2s_0.1 <- R2s_0.1 |> rename(ID = V3,
                     AR2 = V4,
                     DR2 = V5)
df_flipped_0.1 <- df_dosage_0.1 |> filter(MAF >= 0.5)
df_dosage_0.1 <- df_dosage_0.1 |> filter(MAF <= 0.5)
df_dosage_0.1 <- df_dosage_0.1 |> filter(Masked == '-9')
# df <- df |>
#   rename(AF = MAF) |>
#   mutate(MAF = pmin(AF, 1 - AF))

df_0.1 <- df_dosage_0.1 |>
  left_join(R2s_0.1, by = "ID")

df_0.1 <- df_0.1 |> mutate(target_site = ifelse(POS %in% targets$V2, 'target', 'flanking'))

df_0.1 <- df_0.1 |> mutate(Discordance = abs(Imputed - Original),
                             Concordance = 1 - (abs(Imputed - Original)),
                             genotype = case_when(
                               Original == 2 ~ 'alt',
                               Original == 1 ~ 'het',
                               Original == 0 ~ 'ref',))

##########################
####   Masking 0.2   #####
##########################

df_dosage_0.2 <- read.table('0.2_set3_AF_across_chrom_dosage.tsv', header=TRUE)
R2s_0.2 <- read.table('0.2_set3_DR2_AR2_across_chrom_dosage.tsv', header=FALSE)
R2s_0.2 <- R2s_0.2 |> rename(ID = V3,
                             AR2 = V4,
                             DR2 = V5)
df_flipped_0.2 <- df_dosage_0.2 |> filter(MAF >= 0.5)
df_dosage_0.2 <- df_dosage_0.2 |> filter(MAF <= 0.5)
df_dosage_0.2 <- df_dosage_0.2 |> filter(Masked == '-9')
# df <- df |>
#   rename(AF = MAF) |>
#   mutate(MAF = pmin(AF, 1 - AF))

df_0.2 <- df_dosage_0.2 |>
  left_join(R2s_0.2, by = "ID")

df_0.2 <- df_0.2 |> mutate(target_site = ifelse(POS %in% targets$V2, 'target', 'flanking'))

df_0.2 <- df_0.2 |> mutate(Discordance = abs(Imputed - Original),
                             Concordance = 1 - (abs(Imputed - Original)),
                             genotype = case_when(
                               Original == 2 ~ 'alt',
                               Original == 1 ~ 'het',
                               Original == 0 ~ 'ref',))

##########################
####   Masking 0.3   #####
##########################

df_dosage_0.3 <- read.table('0.3_set3_AF_across_chrom_dosage.tsv', header=TRUE)
R2s_0.3 <- read.table('0.3_set3_DR2_AR2_across_chrom_dosage.tsv', header=FALSE)
R2s_0.3 <- R2s_0.3 |> rename(ID = V3,
                             AR2 = V4,
                             DR2 = V5)
df_flipped_0.3 <- df_dosage_0.3 |> filter(MAF >= 0.5)
df_dosage_0.3 <- df_dosage_0.3 |> filter(MAF <= 0.5)
df_dosage_0.3 <- df_dosage_0.3 |> filter(Masked == '-9')
# df <- df |>
#   rename(AF = MAF) |>
#   mutate(MAF = pmin(AF, 1 - AF))

df_0.3 <- df_dosage_0.3 |>
  left_join(R2s_0.3, by = "ID")

df_0.3 <- df_0.3 |> mutate(target_site = ifelse(POS %in% targets$V2, 'target', 'flanking'))

df_0.3 <- df_0.3 |> mutate(Discordance = abs(Imputed - Original),
                             Concordance = 1 - (abs(Imputed - Original)),
                             genotype = case_when(
                               Original == 2 ~ 'alt',
                               Original == 1 ~ 'het',
                               Original == 0 ~ 'ref',))

##########################
####   Masking 0.5   #####
##########################

df_dosage_0.5 <- read.table('0.5_set3_AF_across_chrom_dosage.tsv', header=TRUE)
R2s_0.5 <- read.table('0.5_set3_DR2_AR2_across_chrom_dosage.tsv', header=FALSE)
R2s_0.5 <- R2s_0.5 |> rename(ID = V3,
                             AR2 = V4,
                             DR2 = V5)
df_flipped_0.5 <- df_dosage_0.5 |> filter(MAF >= 0.5)
df_dosage_0.5 <- df_dosage_0.5 |> filter(MAF <= 0.5)
df_dosage_0.5 <- df_dosage_0.5 |> filter(Masked == '-9')
# df <- df |>
#   rename(AF = MAF) |>
#   mutate(MAF = pmin(AF, 1 - AF))

df_0.5 <- df_dosage_0.5 |>
  left_join(R2s_0.5, by = "ID")

df_0.5 <- df_0.5 |> mutate(target_site = ifelse(POS %in% targets$V2, 'target', 'flanking'))

df_0.5 <- df_0.5 |> mutate(Discordance = abs(Imputed - Original),
                           Concordance = 1 - (abs(Imputed - Original)),
                           genotype = case_when(
                             Original == 2 ~ 'alt',
                             Original == 1 ~ 'het',
                             Original == 0 ~ 'ref',))


######################################################################################
#   AF Binning   #####################################################################
######################################################################################

# bins <- c(0, 0.05, 0.15, 0.25, 0.35, 0.50)
# labels <- c("[0–0.05]", "[0.05–0.15]", "[0.15–0.25]", "[0.25–0.35]", "[0.35–0.50]")

# Define bins with size 0.005
bins <- seq(0, 0.5, by = 0.025)

# Create labels for the bins
labels <- paste0("[", bins[-length(bins)], "–", bins[-1], "]")

bin_df_0.3 <- df_0.3 %>%
  mutate(MAF_bin = cut(MAF, breaks = bins, labels = labels, include.lowest = TRUE))

######################################################################################
#   AR2 and DR2 Stuff set1   #########################################################
######################################################################################

##########################
####   Masking 0.05   ####
##########################

DR2_df_0.05 <- bin_df_0.05 |>
  group_by(MAF_bin) |>
  summarise(mean_DR2 = mean(DR2, na.rm = TRUE))

DR2_df_0.05_set3 <- DR2_df_0.05 |> rename(set1 = mean_DR2)

AR2_df_0.05 <- bin_df_0.05 |>
  group_by(MAF_bin) |>
  summarise(mean_AR2 = mean(AR2, na.rm = TRUE))

AR2_df_0.05_set3 <- AR2_df_0.05 |> rename(set1 = mean_AR2)

# Combine the data frames
df_list <- list(DR2_df_0.05_set1, DR2_df_0.05_set2, DR2_df_0.05_set3)
# Combine data frames by MAF_bin and target_site
combined_DR2_0.05 <- reduce(df_list, full_join, by = "MAF_bin")

write.table(combined_DR2_0.05, file = "combined_DR2_0.05.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

# Combine the data frames
df_list <- list(AR2_df_0.05_set1, AR2_df_0.05_set2, AR2_df_0.05_set3)
# Combine data frames by MAF_bin and target_site
combined_AR2_0.05 <- reduce(df_list, full_join, by = "MAF_bin")

write.table(combined_AR2_0.05, file = "combined_AR2_0.05.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

DR2_df_targets_0.05 <- bin_df |>
  group_by(MAF_bin, target_site) |>
  summarise(mean_DR2 = mean(DR2, na.rm = TRUE))

##########################
####   Masking 0.1   #####
##########################

# AF binning

DR2_df_0.1 <- bin_df_0.1 |>
  group_by(MAF_bin) |>
  summarise(mean_DR2 = mean(DR2, na.rm = TRUE))

DR2_df_0.1_set3 <- DR2_df_0.1 |> rename(set3 = mean_DR2)

AR2_df_0.1 <- bin_df_0.1 |>
  group_by(MAF_bin) |>
  summarise(mean_AR2 = mean(AR2, na.rm = TRUE))

AR2_df_0.1_set3 <- AR2_df_0.1 |> rename(set3 = mean_AR2)

# Combine the data frames
df_list <- list(DR2_df_0.1_set1, DR2_df_0.1_set2, DR2_df_0.1_set3)
# Combine data frames by MAF_bin and target_site
combined_DR2_0.1 <- reduce(df_list, full_join, by = "MAF_bin")

write.table(combined_DR2_0.1, file = "combined_DR2_0.1.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

# Combine the data frames
df_list <- list(AR2_df_0.1_set1, AR2_df_0.1_set2, AR2_df_0.1_set3)
# Combine data frames by MAF_bin and target_site
combined_AR2_0.1 <- reduce(df_list, full_join, by = "MAF_bin")

write.table(combined_AR2_0.1, file = "combined_AR2_0.1.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)


DR2_df_targets_0.1 <- bin_df |>
  group_by(MAF_bin, target_site) |>
  summarise(mean_DR2 = mean(DR2, na.rm = TRUE))

##########################
####   Masking 0.2   #####
##########################

# AF binning

DR2_df_0.2 <- bin_df_0.2 |>
  group_by(MAF_bin) |>
  summarise(mean_DR2 = mean(DR2, na.rm = TRUE))

DR2_df_0.2_set3 <- DR2_df_0.2 |> rename(set3 = mean_DR2)

AR2_df_0.2 <- bin_df_0.2 |>
  group_by(MAF_bin) |>
  summarise(mean_AR2 = mean(AR2, na.rm = TRUE))

AR2_df_0.2_set3 <- AR2_df_0.2 |> rename(set3 = mean_AR2)

# Combine the data frames
df_list <- list(DR2_df_0.2_set1, DR2_df_0.2_set2, DR2_df_0.2_set3)
# Combine data frames by MAF_bin and target_site
combined_DR2_0.2 <- reduce(df_list, full_join, by = "MAF_bin")

write.table(combined_DR2_0.2, file = "combined_DR2_0.2.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

# Combine the data frames
df_list <- list(AR2_df_0.2_set1, AR2_df_0.2_set2, AR2_df_0.2_set3)
# Combine data frames by MAF_bin and target_site
combined_AR2_0.2 <- reduce(df_list, full_join, by = "MAF_bin")

write.table(combined_AR2_0.2, file = "combined_AR2_0.2.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

DR2_df_targets_0.2 <- bin_df |>
  group_by(MAF_bin, target_site) |>
  summarise(mean_DR2 = mean(DR2, na.rm = TRUE))

##########################
####   Masking 0.3   #####
##########################

DR2_df_0.3 <- bin_df_0.3 |>
  group_by(MAF_bin) |>
  summarise(mean_DR2 = mean(DR2, na.rm = TRUE))

DR2_df_0.3_set3 <- DR2_df_0.3 |> rename(set3 = mean_DR2)

AR2_df_0.3 <- bin_df_0.3 |>
  group_by(MAF_bin) |>
  summarise(mean_AR2 = mean(AR2, na.rm = TRUE))

AR2_df_0.3_set3 <- AR2_df_0.3 |> rename(set3 = mean_AR2)

# Combine the data frames
df_list <- list(DR2_df_0.3_set1, DR2_df_0.3_set2, DR2_df_0.3_set3)
# Combine data frames by MAF_bin and target_site
combined_DR2_0.3 <- reduce(df_list, full_join, by = "MAF_bin")

write.table(combined_DR2_0.3, file = "combined_DR2_0.3.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

# Combine the data frames
df_list <- list(AR2_df_0.3_set1, AR2_df_0.3_set2, AR2_df_0.3_set3)
# Combine data frames by MAF_bin and target_site
combined_AR2_0.3 <- reduce(df_list, full_join, by = "MAF_bin")

write.table(combined_AR2_0.3, file = "combined_AR2_0.3.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

DR2_df_targets_0.3 <- bin_df |>
  group_by(MAF_bin, target_site) |>
  summarise(mean_DR2 = mean(DR2, na.rm = TRUE))

##########################
####   Masking 0.5   #####
##########################

DR2_df_0.5 <- bin_df_0.5 |>
  group_by(MAF_bin) |>
  summarise(mean_DR2 = mean(DR2, na.rm = TRUE))

DR2_df_0.5_set3 <- DR2_df_0.5 |> rename(set3 = mean_DR2)

AR2_df_0.5 <- bin_df_0.5 |>
  group_by(MAF_bin) |>
  summarise(mean_AR2 = mean(AR2, na.rm = TRUE))

AR2_df_0.5_set3 <- AR2_df_0.5 |> rename(set3 = mean_AR2)

# Combine the data frames
df_list <- list(DR2_df_0.5_set1, DR2_df_0.5_set2, DR2_df_0.5_set3)
# Combine data frames by MAF_bin and target_site
combined_DR2_0.5 <- reduce(df_list, full_join, by = "MAF_bin")

write.table(combined_DR2_0.5, file = "combined_DR2_0.5.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

# Combine the data frames
df_list <- list(AR2_df_0.5_set1, AR2_df_0.5_set2, AR2_df_0.5_set3)
# Combine data frames by MAF_bin and target_site
combined_AR2_0.5 <- reduce(df_list, full_join, by = "MAF_bin")

write.table(combined_AR2_0.5, file = "combined_AR2_0.5.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

DR2_df_targets_0.5 <- bin_df |>
  group_by(MAF_bin, target_site) |>
  summarise(mean_DR2 = mean(DR2, na.rm = TRUE))



##########################
####   AR2, DR2 Plots   ##
##########################

plot <- ggplot(DR2_df_0.1, aes(x = MAF_bin, y = mean_DR2, group=1)) +
  geom_point(color = "skyblue3", shape=1, size=2) +
  geom_line(color = "skyblue3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #ylim(0.9, 1.0) +
  geom_point(data = DR2_df_0.3, aes(x = MAF_bin, y = mean_DR2, group=2), color="darkorange", shape=2, size=2) +
  geom_line(data = DR2_df_0.3, aes(x = MAF_bin, y = mean_DR2, group=2), color="darkorange") +
  geom_point(data = DR2_df_0.2, aes(x = MAF_bin, y = mean_DR2, group=2), color="darkred", shape=5, size=2) +
  geom_line(data = DR2_df_0.2, aes(x = MAF_bin, y = mean_DR2, group=2), color="darkred") +
  geom_point(data = DR2_df_0.05, aes(x = MAF_bin, y = mean_DR2, group=2), color="darkgreen", shape=4, size=2) +
  geom_line(data = DR2_df_0.05, aes(x = MAF_bin, y = mean_DR2, group=2), color="darkgreen") +
  geom_point(data = DR2_df_0.5, aes(x = MAF_bin, y = mean_DR2, group=2), color="purple", shape=5, size=2) +
  geom_line(data = DR2_df_0.5, aes(x = MAF_bin, y = mean_DR2, group=2), color="purple") +
  labs(title = "Dosage R2 at different masked proportions",
       x = "Allele Frequency bin",
       y = "mean DR2")

plot

ggsave("DR2_masked.png", plot = plot, height = 10, width = 10, units = "in")

plot <- ggplot(AR2_df_0.1, aes(x = MAF_bin, y = mean_AR2, group=1)) +
  geom_point(color = "skyblue3", shape=1, size=2) +
  geom_line(color = "skyblue3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #ylim(0.9, 1.0) +
  geom_point(data = AR2_df_0.3, aes(x = MAF_bin, y = mean_AR2, group=2), color="darkorange", shape=2, size=2) +
  geom_line(data = AR2_df_0.3, aes(x = MAF_bin, y = mean_AR2, group=2), color="darkorange") +
  geom_point(data = AR2_df_0.2, aes(x = MAF_bin, y = mean_AR2, group=2), color="darkred", shape=5, size=2) +
  geom_line(data = AR2_df_0.2, aes(x = MAF_bin, y = mean_AR2, group=2), color="darkred") +
  geom_point(data = AR2_df_0.05, aes(x = MAF_bin, y = mean_AR2, group=2), color="darkgreen", shape=4, size=2) +
  geom_line(data = AR2_df_0.05, aes(x = MAF_bin, y = mean_AR2, group=2), color="darkgreen") +
  labs(title = "Allelic R2 at different masked proportions",
       x = "Allele Frequency bin",
       y = "mean AR2")

plot

ggsave("AR2_masked.png", plot = plot, height = 10, width = 10, units = "in")


######################################################################################
#   Accuracy   #######################################################################
######################################################################################

df1 <- read.table('results/0.5_set1_AF_across_chrom_discrete.tsv', header=TRUE)

df_flipped <- df1 |> filter(MAF >= 0.5)
df <- df1 |> filter(MAF <= 0.5)
df <- df |> filter(Masked == '-9')

# Define bins with size 0.005
bins <- seq(0, 0.5, by = 0.025)

# Create labels for the bins
labels <- paste0("[", bins[-length(bins)], "–", bins[-1], "]")

bin_df_0.5 <- df %>%
  mutate(MAF_bin = cut(MAF, breaks = bins, labels = labels, include.lowest = TRUE))

##########################
####   Masking 0.05   ####
##########################

proportion_0.05 <- bin_df_0.05 %>%
  group_by(MAF_bin) %>%
  summarize(proportion_match = mean(Imputed == Original))

proportion_0.05_split <- bin_df_0.05 %>%
  group_by(MAF_bin, genotype) %>%
  summarize(proportion_match = mean(Imputed == Original))

##########################
####   Masking 0.1   #####
##########################

proportion_0.1 <- bin_df_0.1 %>%
  group_by(MAF_bin) %>%
  summarize(proportion_match = mean(Imputed == Original))

proportion_0.1 <- bin_df_0.1 %>%
  group_by(MAF_bin, genotype) %>%
  summarize(proportion_match = mean(Imputed == Original))

##########################
####   Masking 0.2   #####
##########################

proportion_0.2 <- bin_df_0.2 %>%
  group_by(MAF_bin) %>%
  summarize(proportion_match = mean(Imputed == Original))

proportion_0.2 <- bin_df_0.2 %>%
  group_by(MAF_bin, genotype) %>%
  summarize(proportion_match = mean(Imputed == Original))

##########################
####   Masking 0.3   #####
##########################

proportion_0.3 <- bin_df_0.3 %>%
  group_by(MAF_bin) %>%
  summarize(proportion_match = mean(Imputed == Original))

proportion_0.3_split <- bin_df_0.3 %>%
  group_by(MAF_bin, genotype) %>%
  summarize(proportion_match = mean(Imputed == Original))


##########################
####   Masking 0.5   #####
##########################

proportion_0.5 <- bin_df_0.5 %>%
  group_by(MAF_bin) %>%
  summarize(proportion_match = mean(Imputed == Original))

proportion_0.3_split <- bin_df_0.3 %>%
  group_by(MAF_bin, genotype) %>%
  summarize(proportion_match = mean(Imputed == Original))


##########################
####   STITCH   ##########
##########################

proportion_stitch <- bin_df_stitch %>%
  group_by(MAF_bin) %>%
  summarize(proportion_match = mean(Imputed == Original))

write.table(proportion_stitch, 'stitch_accuracy_MAF.tsv', sep='\t', row.names = FALSE, quote = FALSE)

proportion_0.3_split <- bin_df_0.3 %>%
  group_by(MAF_bin, genotype) %>%
  summarize(proportion_match = mean(Imputed == Original))

##########################
####   Accuracy Plots   ##
##########################

plot <- ggplot(proportion_0.1, aes(x = MAF_bin, y = proportion_match, group=1)) +
  geom_point(color = "skyblue3", shape=1, size=2) +
  geom_line(color = "skyblue3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0.8, 1.0) +
  geom_point(data = proportion_0.3, aes(x = MAF_bin, y = proportion_match, group=2), color="darkorange", shape=2, size=2) +
  geom_line(data = proportion_0.3, aes(x = MAF_bin, y = proportion_match, group=2), color="darkorange") +
  geom_point(data = proportion_0.2, aes(x = MAF_bin, y = proportion_match, group=2), color="darkred", shape=5, size=2) +
  geom_line(data = proportion_0.2, aes(x = MAF_bin, y = proportion_match, group=2), color="darkred") +
  geom_point(data = proportion_0.05, aes(x = MAF_bin, y = proportion_match, group=2), color="darkgreen", shape=4, size=2) +
  geom_line(data = proportion_0.05, aes(x = MAF_bin, y = proportion_match, group=2), color="darkgreen") +
  geom_point(data = proportion_stitch, aes(x = MAF_bin, y = proportion_match, group=2), color="black") +
  geom_line(data = proportion_stitch, aes(x = MAF_bin, y = proportion_match, group=2), color="black") +
  labs(title = "Accuracy",
       x = "Allele Frequency bin",
       y = "mean Accuracy")

plot

ggsave("Accuracy_all_masked_stitch.png", plot = plot, height = 6, width = 10, units = "in")

plot <- ggplot(proportion_0.05, aes(x = MAF_bin, y = proportion_match, group=1)) +
  geom_point(color = "darkgreen", shape=4, size=2) +
  geom_line(color = "darkgreen", linetype="dashed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0.5, 1.0) +
  geom_point(data = proportion_0.05_split, aes(x = MAF_bin, y = proportion_match, group=genotype, color=genotype)) +
  geom_line(data = proportion_0.05_split, aes(x = MAF_bin, y = proportion_match, group=genotype, color=genotype)) +
  labs(title = "Accuracy for each genotype at 5% masking",
       x = "Allele Frequency bin",
       y = "mean Accuracy")

plot

ggsave("Accuracy_0.05_split_masked.png", plot = plot, height = 6, width = 10, units = "in")

plot <- ggplot(proportion_0.3, aes(x = MAF_bin, y = proportion_match, group=1)) +
  geom_point(color = "darkorange", shape=2, size=2) +
  geom_line(color = "darkorange", linetype="dashed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0.5, 1.0) +
  geom_point(data = proportion_0.3_split, aes(x = MAF_bin, y = proportion_match, group=genotype, color=genotype)) +
  geom_line(data = proportion_0.3_split, aes(x = MAF_bin, y = proportion_match, group=genotype, color=genotype)) +
  labs(title = "Accuracy for each genotype at 30% masking",
       x = "Allele Frequency bin",
       y = "mean Accuracy")

plot

ggsave("Accuracy_0.3_split_masked.png", plot = plot, height = 6, width = 10, units = "in")



######################################################################################
#   Reading and cleaning the data STITCH   ###########################################
######################################################################################

targets <- read.table('/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Sam_034_targets_liftover.txt', header=FALSE)
targets <- targets |> filter(V1 =='chr27')

df <- read.table('../EGGL_STITCH/K10_evaluation_results/dataframes/chr27/dosage.tsv', header=TRUE)
df_flipped <- df |> filter(MAF >= 0.5)
df <- df |> filter(MAF <= 0.5)

# df <- df |>
#   rename(AF = MAF) |>
#   mutate(MAF = pmin(AF, 1 - AF))

df <- df |> mutate(target_site = ifelse(POS %in% targets$V2, 'target', 'flanking'))

df <- df |> mutate(Discordance = abs(Imputed - Original),
                             Concordance = 1 - (abs(Imputed - Original)),
                             genotype = case_when(
                               Original == 2 ~ 'alt',
                               Original == 1 ~ 'het',
                               Original == 0 ~ 'ref',))


######################################################################################
#   AF Binning   #####################################################################
######################################################################################

# bins <- c(0, 0.05, 0.15, 0.25, 0.35, 0.50)
# labels <- c("[0–0.05]", "[0.05–0.15]", "[0.15–0.25]", "[0.25–0.35]", "[0.35–0.50]")

# Define bins with size 0.005
bins <- seq(0, 0.5, by = 0.025)

# Create labels for the bins
labels <- paste0("[", bins[-length(bins)], "–", bins[-1], "]")

bin_df_0.3 <- df_0.3 %>%
  mutate(MAF_bin = cut(MAF, breaks = bins, labels = labels, include.lowest = TRUE))

######################################################################################
#   AR2 and DR2 Stuff   ##############################################################
######################################################################################



##########################
####   AR2, DR2 Plots   ##
##########################

plot <- ggplot(DR2_df_0.1, aes(x = MAF_bin, y = mean_DR2, group=1)) +
  geom_point(color = "skyblue3", shape=1, size=2) +
  geom_line(color = "skyblue3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #ylim(0.9, 1.0) +
  geom_point(data = DR2_df_0.3, aes(x = MAF_bin, y = mean_DR2, group=2), color="darkorange", shape=2, size=2) +
  geom_line(data = DR2_df_0.3, aes(x = MAF_bin, y = mean_DR2, group=2), color="darkorange") +
  geom_point(data = DR2_df_0.2, aes(x = MAF_bin, y = mean_DR2, group=2), color="darkred", shape=5, size=2) +
  geom_line(data = DR2_df_0.2, aes(x = MAF_bin, y = mean_DR2, group=2), color="darkred") +
  geom_point(data = DR2_df_0.05, aes(x = MAF_bin, y = mean_DR2, group=2), color="darkgreen", shape=4, size=2) +
  geom_line(data = DR2_df_0.05, aes(x = MAF_bin, y = mean_DR2, group=2), color="darkgreen") +
  labs(title = "Dosage R2 at different masked proportions",
       x = "Allele Frequency bin",
       y = "mean DR2")

plot

ggsave("DR2_masked.png", plot = plot, height = 10, width = 10, units = "in")

plot <- ggplot(AR2_df_0.1, aes(x = MAF_bin, y = mean_AR2, group=1)) +
  geom_point(color = "skyblue3", shape=1, size=2) +
  geom_line(color = "skyblue3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #ylim(0.9, 1.0) +
  geom_point(data = AR2_df_0.3, aes(x = MAF_bin, y = mean_AR2, group=2), color="darkorange", shape=2, size=2) +
  geom_line(data = AR2_df_0.3, aes(x = MAF_bin, y = mean_AR2, group=2), color="darkorange") +
  geom_point(data = AR2_df_0.2, aes(x = MAF_bin, y = mean_AR2, group=2), color="darkred", shape=5, size=2) +
  geom_line(data = AR2_df_0.2, aes(x = MAF_bin, y = mean_AR2, group=2), color="darkred") +
  geom_point(data = AR2_df_0.05, aes(x = MAF_bin, y = mean_AR2, group=2), color="darkgreen", shape=4, size=2) +
  geom_line(data = AR2_df_0.05, aes(x = MAF_bin, y = mean_AR2, group=2), color="darkgreen") +
  labs(title = "Allelic R2 at different masked proportions",
       x = "Allele Frequency bin",
       y = "mean AR2")

plot

ggsave("AR2_masked.png", plot = plot, height = 10, width = 10, units = "in")


######################################################################################
#   Accuracy   #######################################################################
######################################################################################



##########################
####   Masking 0.05   ####
##########################

proportion_0.05 <- bin_df_0.05 %>%
  group_by(MAF_bin) %>%
  summarize(proportion_match = mean(Imputed == Original))

proportion_0.05_split <- bin_df_0.05 %>%
  group_by(MAF_bin, genotype) %>%
  summarize(proportion_match = mean(Imputed == Original))

##########################
####   Masking 0.1   #####
##########################

proportion_0.1 <- bin_df_0.1 %>%
  group_by(MAF_bin) %>%
  summarize(proportion_match = mean(Imputed == Original))

proportion_0.1 <- bin_df_0.1 %>%
  group_by(MAF_bin, genotype) %>%
  summarize(proportion_match = mean(Imputed == Original))

##########################
####   Masking 0.2   #####
##########################

proportion_0.2 <- bin_df_0.2 %>%
  group_by(MAF_bin) %>%
  summarize(proportion_match = mean(Imputed == Original))

proportion_0.2 <- bin_df_0.2 %>%
  group_by(MAF_bin, genotype) %>%
  summarize(proportion_match = mean(Imputed == Original))

##########################
####   Masking 0.3   #####
##########################

proportion_0.3 <- bin_df_0.3 %>%
  group_by(MAF_bin) %>%
  summarize(proportion_match = mean(Imputed == Original))

proportion_0.3_split <- bin_df_0.3 %>%
  group_by(MAF_bin, genotype) %>%
  summarize(proportion_match = mean(Imputed == Original))

##########################
####   Accuracy Plots   ##
##########################

plot <- ggplot(proportion_0.1, aes(x = MAF_bin, y = proportion_match, group=1)) +
  geom_point(color = "skyblue3", shape=1, size=2) +
  geom_line(color = "skyblue3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0.5, 1.0) +
  geom_point(data = proportion_0.3, aes(x = MAF_bin, y = proportion_match, group=2), color="darkorange", shape=2, size=2) +
  geom_line(data = proportion_0.3, aes(x = MAF_bin, y = proportion_match, group=2), color="darkorange") +
  geom_point(data = proportion_0.2, aes(x = MAF_bin, y = proportion_match, group=2), color="darkred", shape=5, size=2) +
  geom_line(data = proportion_0.2, aes(x = MAF_bin, y = proportion_match, group=2), color="darkred") +
  geom_point(data = proportion_0.05, aes(x = MAF_bin, y = proportion_match, group=2), color="darkgreen", shape=4, size=2) +
  geom_line(data = proportion_0.05, aes(x = MAF_bin, y = proportion_match, group=2), color="darkgreen") +
  labs(title = "Accuracy at different masked proportions",
       x = "Allele Frequency bin",
       y = "mean Accuracy")

plot

ggsave("Accuracy_all_masked.png", plot = plot, height = 6, width = 10, units = "in")

plot <- ggplot(proportion_0.05, aes(x = MAF_bin, y = proportion_match, group=1)) +
  geom_point(color = "darkgreen", shape=4, size=2) +
  geom_line(color = "darkgreen", linetype="dashed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0.5, 1.0) +
  geom_point(data = proportion_0.05_split, aes(x = MAF_bin, y = proportion_match, group=genotype, color=genotype)) +
  geom_line(data = proportion_0.05_split, aes(x = MAF_bin, y = proportion_match, group=genotype, color=genotype)) +
  labs(title = "Accuracy for each genotype at 5% masking",
       x = "Allele Frequency bin",
       y = "mean Accuracy")

plot

ggsave("Accuracy_0.05_split_masked.png", plot = plot, height = 6, width = 10, units = "in")

plot <- ggplot(proportion_0.3, aes(x = MAF_bin, y = proportion_match, group=1)) +
  geom_point(color = "darkorange", shape=2, size=2) +
  geom_line(color = "darkorange", linetype="dashed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0.5, 1.0) +
  geom_point(data = proportion_0.3_split, aes(x = MAF_bin, y = proportion_match, group=genotype, color=genotype)) +
  geom_line(data = proportion_0.3_split, aes(x = MAF_bin, y = proportion_match, group=genotype, color=genotype)) +
  labs(title = "Accuracy for each genotype at 30% masking",
       x = "Allele Frequency bin",
       y = "mean Accuracy")

plot

ggsave("Accuracy_0.3_split_masked.png", plot = plot, height = 6, width = 10, units = "in")

















######################################################################################
#   Correlation   ####################################################################
######################################################################################

#  correlation between true and imputed genotypes (Cor), which were coded as 0, 1, and 2 for genotypes AA, AB, and BB, respectively

df_discrete <- read.table('AF_across_chrom_discrete.tsv', header=TRUE)
df_discrete <- df_discrete |> rename(Discrete = Imputed)

df_for_cor <- df |>
  left_join(df_discrete |> select(ID, Sample, Discrete), by = c("ID", "Sample"))

df_cor_all <- bin_df |>
  group_by(MAF_bin) |>
  summarize(
    r2 = if(n() > 1 && sd(Discrete, na.rm = TRUE) != 0 && sd(Imputed, na.rm = TRUE) != 0) {
      (cor(Discrete, Imputed, use = "complete.obs", method = "pearson")^2)
    } else {
      NA
    },
    .groups = 'drop'
  )

plot <- ggplot(df_cor_all, aes(x = MAF_bin, y = r2)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0.5, 1.0)

plot

df_cor_targets <- bin_df |>
  group_by(MAF_bin, target_site) |>
  summarize(
    r2 = if(n() > 1 && sd(Discrete, na.rm = TRUE) != 0 && sd(Imputed, na.rm = TRUE) != 0) {
      (cor(Discrete, Imputed, use = "complete.obs", method = "pearson")^2)
    } else {
      NA
    },
    .groups = 'drop'
  )

plot <- ggplot(df_cor_targets, aes(x = MAF_bin, y = r2, color=target_site, group=target_site)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0.5, 1.0)

plot







df$Original <- as.numeric(df$Original)
df$Imputed <- as.numeric(df$Imputed)

df <- df |> mutate(Discordance = abs(Imputed - Original),
                   Concordance = 1 - (abs(Imputed - Original)),
                                      genotype = case_when(
                                      Original == 2 ~ 'alt',
                                      Original == 1 ~ 'het',
                                      Original == 0 ~ 'ref',))


df <- na.omit(df)

df_discrete <- df_discrete |> rename(Discrete = Imputed)

df_for_cor <- df %>%
  left_join(df_discrete %>% select(ID, Sample, Discrete), by = c("ID", "Sample"))



df_cor <- na.omit(df_cor)

plot <- ggplot(df_cor, aes(x = MAF, y = r2)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot

mean_discordance_df <- df %>%
  group_by(MAF, genotype) %>%
  summarise(mean_Discordance = mean(Discordance, na.rm = TRUE))

bins <- c(0, 0.05, 0.15, 0.25, 0.35, 0.50)
labels <- c("[0–0.05]", "[0.05–0.15]", "[0.15–0.25]", "[0.25–0.35]", "[0.35–0.50]")

# Define bins with size 0.005
bins <- seq(0, 0.5, by = 0.01)

# Create labels for the bins
labels <- paste0("[", bins[-length(bins)], "–", bins[-1], "]")

bin_df <- mean_discordance_df %>%
  mutate(MAF_bin = cut(MAF, breaks = bins, labels = labels, include.lowest = TRUE))

df_cor_all <- bin_df |>
  group_by(MAF_bin, target_site) |>
  summarize(
    r2 = if(n() > 1 && sd(Discrete, na.rm = TRUE) != 0 && sd(Imputed, na.rm = TRUE) != 0) {
      (cor(Discrete, Imputed, use = "complete.obs", method = "pearson")^2)
    } else {
      NA
    },
    .groups = 'drop'
  )

plot <- ggplot(df_cor, aes(x = MAF_bin, y = r2, color = genotype, group = genotype)) +
  geom_point() +
  geom_line() +
  geom_point(data = df_cor_all, aes(x = MAF_bin, y = r2), color = "blue", inherit.aes = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot <- ggplot(df_cor_all, aes(x = MAF_bin, y = r2, color = target_site, group = target_site)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot

plot <- ggplot(df_cor, aes(x = MAF_bin, y = r2, color=genotype, group=genotype)) +
  geom_point(df_cor_all, aes(x = MAF_bin, y = r2)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot

# Calculate the average IQS values for each MAF bin
split_avg_con <- all_df %>%
  group_by(MAF_bin, genotype) %>%
  summarise(mean_Concordance = mean(Concordance, na.rm = TRUE))

all_avg_con <- all_df %>%
  group_by(MAF_bin) %>%
  summarise(mean_Concordance = mean(Concordance, na.rm = TRUE))

plot <- ggplot(split_avg_con, aes(x = MAF_bin, y = mean_Concordance, color=genotype)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot




# Calculate the average IQS values for each MAF bin
split_avg_dis <- bin_df %>%
  group_by(MAF_bin, genotype) %>%
  summarise(mean_Discordance = mean(Discordance, na.rm = TRUE))

all_avg_con <- all_df %>%
  group_by(MAF_bin) %>%
  summarise(mean_Concordance = mean(Concordance, na.rm = TRUE))

plot <- ggplot(bin_df, aes(x = MAF_bin, y = mean_Discordance, color=genotype)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot
