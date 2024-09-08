library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbreak)

setwd('/scratch.global/dimml002/PaintedPatches/results')

df_0.05 <- read.table('combined_DR2_0.05.tsv', header=TRUE)
df_0.1 <- read.table('combined_DR2_0.1.tsv', header=TRUE)
df_0.2 <- read.table('combined_DR2_0.2.tsv', header=TRUE)
df_0.3 <- read.table('combined_DR2_0.3.tsv', header=TRUE)
df_0.5 <- read.table('combined_DR2_0.5.tsv', header=TRUE)

# Function to calculate mean and confidence interval
mean_ci <- function(x, conf = 0.95) {
  mean_x <- mean(x, na.rm = TRUE)
  stderr <- sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
  ci <- qt((1 + conf) / 2, df = length(na.omit(x)) - 1) * stderr
  return(data.frame(mean = mean_x, ci_lower = mean_x - ci, ci_upper = mean_x + ci))
}

# Calculate mean and CI for each row
ci_0.3 <- df_0.3 %>%
  rowwise() %>%
  mutate(mean_r2 = mean(c_across(set1:set3), na.rm = TRUE),
         ci = list(mean_ci(c_across(set1:set3)))) %>%
  unnest_wider(ci)

stitch_df <- read.table('/scratch.global/dimml002/EGGL_STITCH/three.chrom_dosage.tsv', header=TRUE)

stitch_df$Imputed <- as.numeric(stitch_df$Imputed)
stitch_df$Original <- as.numeric((stitch_df$Original))

stitch_df2 <- stitch_df |> filter(MAF < 0.5)

# Define bins with size 0.005
bins <- seq(0, 0.5, by = 0.025)

# Create labels for the bins
labels <- paste0("[", bins[-length(bins)], "–", bins[-1], "]")

bin_df_stitch <- stitch_df2 %>%
  mutate(MAF_bin = cut(MAF, breaks = bins, labels = labels, include.lowest = TRUE))

unique_imputed_count <- bin_df_stitch %>%
  distinct(MAF_bin, .keep_all = TRUE) %>%
  summarize(unique_imputed = n_distinct(Imputed))

stitch_DR2 <- bin_df_stitch %>%
  group_by(ID, MAF_bin) %>%
  summarize(DR2 = cor(Imputed, Original, use = "complete.obs")^2)

stitch_DR2_mean <- stitch_DR2 %>%
  group_by(MAF_bin) %>%
  summarize(mean_DR2 = mean(DR2))

write.table(stitch_DR2_mean, 'stitch_DR2_MAF_mean.tsv', sep='\t', row.names = FALSE, quote = FALSE)
stitch_DR2_mean <- read.table('stitch_DR2_MAF_mean.tsv', header=TRUE)

stitch_DR2_mean <- na.omit(stitch_DR2_mean)

na_count <- stitch_DR2 %>%
  filter(rowSums(is.na(.)) > 0) %>%
  nrow()

plot <- ggplot(ci_0.1, aes(x = MAF_bin, y = mean, group=1)) +
  geom_point(color = "skyblue3", shape=1, size=2) +
  geom_line(color = "skyblue3") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), color = "skyblue3", width = 0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #ylim(0.9, 1.0) +
  geom_point(data = ci_0.3, aes(x = MAF_bin, y = mean, group=2), color="darkorange", shape=2, size=2) +
  geom_line(data = ci_0.3, aes(x = MAF_bin, y = mean, group=2), color="darkorange") +
  geom_errorbar(data = ci_0.3, aes(x = MAF_bin, y = mean, ymin = ci_lower, ymax = ci_upper), color = "darkorange", width = 0.2) +
  geom_point(data = ci_0.2, aes(x = MAF_bin, y = mean, group=2), color="darkred", shape=5, size=2) +
  geom_line(data = ci_0.2, aes(x = MAF_bin, y = mean, group=2), color="darkred") +
  geom_errorbar(data = ci_0.2, aes(x = MAF_bin, y = mean, ymin = ci_lower, ymax = ci_upper), color = "darkred", width = 0.2) +
  geom_point(data = ci_0.05, aes(x = MAF_bin, y = mean, group=2), color="darkgreen", shape=4, size=2) +
  geom_line(data = ci_0.05, aes(x = MAF_bin, y = mean, group=2), color="darkgreen") +
  geom_errorbar(data = ci_0.05, aes(x = MAF_bin, y = mean, ymin = ci_lower, ymax = ci_upper), color = "darkgreen", width = 0.2) +
  geom_point(data = ci_0.5, aes(x = MAF_bin, y = mean, group=2), color="purple", shape=5, size=2) +
  geom_line(data = ci_0.5, aes(x = MAF_bin, y = mean, group=2), color="purple") +
  geom_errorbar(data = ci_0.5, aes(x = MAF_bin, y = mean, ymin = ci_lower, ymax = ci_upper), color = "purple", width = 0.2) +
  geom_point(data = stitch_DR2_mean, aes(x = MAF_bin, y = mean_DR2, group=2), color="black", size=2) +
  geom_line(data = stitch_DR2_mean, aes(x = MAF_bin, y = mean_DR2, group=2), color="black") +
  labs(title = "Dosage R2",
       x = "Allele Frequency bin",
       y = "mean DR2") +
  scale_y_continuous(
    breaks = seq(0,1, by = 0.05)  # Set y-axis labels at increments of 0.2
  )

plot

plot2 <- plot + 
  scale_y_break(c(0.575, 0.925), space=0.1)

plot2

ggsave("DR2_masked_stitch_version2.png", plot = plot2, height = 6, width = 10, units = "in")

df_0.05 <- read.table('combined_AR2_0.05.tsv', header=TRUE)
df_0.1 <- read.table('combined_AR2_0.1.tsv', header=TRUE)
df_0.2 <- read.table('combined_AR2_0.2.tsv', header=TRUE)
df_0.3 <- read.table('combined_AR2_0.3.tsv', header=TRUE)
df_0.5 <- read.table('combined_AR2_0.5.tsv', header=TRUE)

# Function to calculate mean and confidence interval
mean_ci <- function(x, conf = 0.95) {
  mean_x <- mean(x, na.rm = TRUE)
  stderr <- sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
  ci <- qt((1 + conf) / 2, df = length(na.omit(x)) - 1) * stderr
  return(data.frame(mean = mean_x, ci_lower = mean_x - ci, ci_upper = mean_x + ci))
}

# Calculate mean and CI for each row
ci_0.5 <- df_0.5 %>%
  rowwise() %>%
  mutate(mean_r2 = mean(c_across(set1:set3), na.rm = TRUE),
         ci = list(mean_ci(c_across(set1:set3)))) %>%
  unnest_wider(ci)

stitch_df <- read.table('/scratch.global/dimml002/EGGL_STITCH/three.chrom_discrete.tsv', header=TRUE)

stitch_df$Imputed <- as.numeric(stitch_df$Imputed)
stitch_df$Original <- as.numeric((stitch_df$Original))

stitch_df2 <- stitch_df |> filter(MAF < 0.5)

# Define bins with size 0.005
bins <- seq(0, 0.5, by = 0.025)

# Create labels for the bins
labels <- paste0("[", bins[-length(bins)], "–", bins[-1], "]")

bin_df_stitch <- stitch_df2 %>%
  mutate(MAF_bin = cut(MAF, breaks = bins, labels = labels, include.lowest = TRUE))

bin_df_stitch <- na.omit(bin_df_stitch)

# Filter out groups with insufficient complete data pairs
bin_df_stitch <- bin_df_stitch %>%
  group_by(ID, MAF_bin) %>%
  filter(!any(is.na(Imputed)) & !any(is.na(Original))) %>%
  summarize(DR2 = cor(Imputed, Original, use = "complete.obs")^2)

stitch_AR2 <- bin_df_stitch %>%
  group_by(ID, MAF_bin) %>%
  summarize(AR2 = cor(Imputed, Original, use = "complete.obs")^2)

stitch_AR2_mean <- stitch_AR2 %>%
  group_by(MAF_bin) %>%
  summarize(mean_AR2 = mean(AR2))

write.table(stitch_AR2_mean, 'stitch_AR2_MAF_mean.tsv', sep='\t', row.names = FALSE, quote = FALSE)

plot <- ggplot(ci_0.1, aes(x = MAF_bin, y = mean, group=1)) +
  geom_point(color = "skyblue3", shape=1, size=2) +
  geom_line(color = "skyblue3") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), color = "skyblue3", width = 0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #ylim(0.9, 1.0) +
  geom_point(data = ci_0.3, aes(x = MAF_bin, y = mean, group=2), color="darkorange", shape=2, size=2) +
  geom_line(data = ci_0.3, aes(x = MAF_bin, y = mean, group=2), color="darkorange") +
  geom_errorbar(data = ci_0.3, aes(x = MAF_bin, y = mean, ymin = ci_lower, ymax = ci_upper), color = "darkorange", width = 0.2) +
  geom_point(data = ci_0.2, aes(x = MAF_bin, y = mean, group=2), color="darkred", shape=5, size=2) +
  geom_line(data = ci_0.2, aes(x = MAF_bin, y = mean, group=2), color="darkred") +
  geom_errorbar(data = ci_0.2, aes(x = MAF_bin, y = mean, ymin = ci_lower, ymax = ci_upper), color = "darkred", width = 0.2) +
  geom_point(data = ci_0.05, aes(x = MAF_bin, y = mean, group=2), color="darkgreen", shape=4, size=2) +
  geom_line(data = ci_0.05, aes(x = MAF_bin, y = mean, group=2), color="darkgreen") +
  geom_errorbar(data = ci_0.05, aes(x = MAF_bin, y = mean, ymin = ci_lower, ymax = ci_upper), color = "darkgreen", width = 0.2) +
  geom_point(data = ci_0.5, aes(x = MAF_bin, y = mean, group=2), color="purple", shape=5, size=2) +
  geom_line(data = ci_0.5, aes(x = MAF_bin, y = mean, group=2), color="purple") +
  geom_errorbar(data = ci_0.5, aes(x = MAF_bin, y = mean, ymin = ci_lower, ymax = ci_upper), color = "purple", width = 0.2) +
  labs(title = "Allelic R2 at different masked proportions",
       x = "Allele Frequency bin",
       y = "mean AR2") +
  geom_point(data = stitch_AR2_mean, aes(x = MAF_bin, y = mean_AR2, group=2), color="black", size=2) +
  geom_line(data = stitch_AR2_mean, aes(x = MAF_bin, y = mean_AR2, group=2), color="black") +
  labs(title = "Dosage R2 at different masked proportions",
       x = "Allele Frequency bin",
       y = "mean DR2") +
  scale_y_continuous(
    breaks = seq(0,1, by = 0.05))  # Set y-axis labels at increments of 0.2

plot

ggsave("AR2_masked_stitch_notwork.png", plot = plot, height = 12, width = 20, units = "in")
