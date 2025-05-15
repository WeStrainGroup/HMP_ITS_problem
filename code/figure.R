library(Biostrings)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(PCAtools)
library(vegan)
library(dplyr)
library(tidyr)


# --- Figure 1 --- 
# Set color scheme for PE250 and PE300
col <- c(PE250 = "#e8998d", PE300 = "#6c9a8b")

# Read length distribution of raw reads
fq_raw_check <- read.delim("/data/wangxinyu/HMP_ITS/data/PRJNA356769_ITS2/1_dada2/check/PRJNA356769_ITS2_fq_stats_raw.tsv",
                           sep = "\t",
                           header = T)

# Create histogram of read lengths
his_reads <- ggplot(fq_raw_check, aes(x = avg_len)) +
    geom_histogram(aes(y = after_stat(count)/sum(after_stat(count))),
                   fill = "grey", color = "black", alpha = 0.75) +
    labs(x = "Average reads length of each sample",
         y = "Proportion") +
    theme_bw() + 
    theme(panel.border = element_rect(color = "black", linewidth = 1))
  
his_reads

# Compare the sequeincing depth between PE250 and PE300 sample
# drop sample with reads lower than 1000
fq_raw_check_filter <- fq_raw_check %>% filter(!grepl("\\.2\\.fastq\\.gz$", file))
fq_raw_check_filter <- fq_raw_check_filter %>% filter(num_seqs >= 1000)
fq_raw_check_filter$avg_len <- paste0("PE", fq_raw_check_filter$avg_len)
fq_raw_check_filter <- fq_raw_check_filter %>% select(file, num_seqs, avg_len)
head(fq_raw_check_filter)


depth <- ggplot(fq_raw_check_filter, aes(x = avg_len, y = num_seqs, fill = avg_len)) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
        geom_jitter(shape = 16, width = 0.25, size = 1.5, alpha = 0.3) +
        labs(x = "Sample type", y = "Sequencing depth") +
        stat_compare_means(
          method = "wilcox.test",  
          label = "p.format",          
          comparisons = list(c("PE250", "PE300")),
          tip.length = 0) +
        scale_fill_manual(values = col) + 
        coord_cartesian(ylim = c(0, 3.5e5)) +
        theme_bw() +
        theme(panel.border = element_rect(color = "black", linewidth = 1), 
              legend.position = "none",
              legend.title = element_blank()) 


ggsave(file.path("/data/wangxinyu/HMP_ITS/figure", "fig_1a.pdf"), depth, width = 4.5, height = 4.5)


# Read ASV sequences and analyze length distribution
PE250_fasta <- readDNAStringSet("/data/wangxinyu/HMP_ITS/data/PRJNA356769a_ITS2/1_dada2/itsx/PRJNA356769a_ITS2_asv_raw.fasta", format = "fasta")
PE300_fasta <- readDNAStringSet("/data/wangxinyu/HMP_ITS/data/PRJNA356769b_ITS2/1_dada2/itsx/PRJNA356769b_ITS2_asv_raw.fasta", format = "fasta")

PE250_fasta
PE300_fasta

# Calculate sequence lengths
PE250_fasta_len <- width(PE250_fasta)
PE300_fasta_len <- width(PE300_fasta)

table(PE250_fasta_len)
max(PE250_fasta_len) # 443

table(PE300_fasta_len)
max(PE300_fasta_len) # 497

# Combine length data for plotting
df <- rbind(data.frame(Length = PE250_fasta_len, Type = "PE250"), data.frame(Length = PE300_fasta_len, Type = "PE300"))

# Create density plot of ASV lengths
his_asv <- ggplot(df, aes(x = Length, y = ..density.., fill = Type)) +
      geom_density(linewidth = 0.25, alpha = 0.5, adjust = 1) +  
      scale_fill_manual(values = col) + 
      labs(
        x = "ASV length",
        y = "Density"
      ) +
      theme_bw() +
      geom_vline(xintercept = 400, linetype = "dashed", color = "black", linewidth = 0.5) + 
      theme(panel.border = element_rect(color = "black", linewidth = 1),
            legend.position = c(0.2, 0.8))

ggsave(file.path("/data/wangxinyu/HMP_ITS/figure", "fig_1b.pdf"), his_asv, width = 4.5, height = 4.5)


# Alpha diversity analysis
diversity <- read.csv("/data/wangxinyu/HMP_ITS/data/PRJNA356769_ITS2/3_phyloseq/fungi/ASV/PRJNA356769_ITS2_diversity_fungi_filterdp.csv", row.names = 1)
diversity <- diversity %>% select("Observed", "Shannon", "Simpson", "Fisher")
meta <- read.csv("/data/wangxinyu/HMP_ITS/data/meta.csv", row.names = 1)

# Merge diversity data with metadata
diversity_merged <- merge(diversity, meta, by = "row.names")
colnames(diversity_merged)[1] <- "sample"
head(diversity_merged)

# Reshape data for plotting
diversity_long <- diversity_merged %>%
  pivot_longer(cols = c("Observed", "Shannon", "Simpson"),
               names_to = "index",
               values_to = "value")

# Create boxplots for diversity indices
div_a <- ggplot(diversity_long, aes(x = type, y = value, fill = type)) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
        geom_jitter(shape = 16, width = 0.25, size = 1.5, alpha = 0.3) +
        facet_wrap(~ index, scales = "free_y", nrow = 1) +
        labs(x = "Sample type", y = "Diversity index") +
        scale_fill_manual(values = col) + 
        stat_compare_means(
          method = "wilcox.test",      
          label = "p.signif",          
          comparisons = list(c("PE250", "PE300")),
          tip.length = 0) +
        theme_bw() +
        theme(panel.border = element_rect(color = "black", linewidth = 1), 
              legend.position = "none",
              axis.text.x = element_blank(),
              legend.title = element_blank()) 

ggsave(file.path("/data/wangxinyu/HMP_ITS/figure", "fig_1c.pdf"), div_a, width = 4.5, height = 4.5)


# Beta diversity analysis
otu <- read.csv("/data/wangxinyu/HMP_ITS/data/PRJNA356769_ITS2/3_phyloseq/fungi/genus/PRJNA356769_ITS2_otu_fungi_gen_p0_clr_filterdp.csv", row.names = 1)
colnames(otu)
rownames(meta)

# Perform PCA analysis
set.seed(6311)
pca <- PCAtools::pca(otu, metadata = meta)

# Create PCA biplot
div_b0 <- biplot(pca, x = "PC1", y = "PC5",
                lab = NULL, colby = "type", colkey = col, pointSize = 1.5, ellipse = TRUE, ellipseLevel = 0.9,
                borderWidth = 0.5, axisLabSize = 11) +
  theme(legend.position = "none",
        legend.title = element_blank()) 
  

# Generate additional PCA plots
screeplot(pca)
pairsplot(pca, lab = NULL, colby = "type")
plotloadings(pca, labSize = 3)

# Permutation test for beta diversity
dist_matrix <- vegdist(t(otu), method = "euclidean")
dist_matrix

set.seed(6311)
adonis2(dist_matrix ~ type, data = meta, permutations = 999) # 0.002 **

# Add p-value annotation to PCA plot
div_b <- div_b0 + annotate("text", x = 5, y = 13, label = "P = 0.002 (perm = 999)", 
                 color = "black", size = 4)

ggsave(file.path("/data/wangxinyu/HMP_ITS/figure", "fig_1d.pdf"), div_b, width = 4.5, height = 4.5)

# Combine and save Figure 1
fig1 <- (depth + his_asv) / (div_a + div_b)
ggsave(file.path("/data/wangxinyu/HMP_ITS/figure", "fig_1.pdf"), fig1, width = 9, height = 9)



# --- Figure 2 --- 
# Read raw fastq statistics
fq_stats_raw <- read.delim("/data/wangxinyu/HMP_ITS/data/PRJNA356769_ITS2/1_dada2/check/PRJNA356769_ITS2_fq_stats_raw.tsv",
                           sep = "\t",
                           header = T)
fq_stats_raw <- fq_stats_raw %>% select("file", "num_seqs")
colnames(fq_stats_raw)[2] <- paste0(colnames(fq_stats_raw)[2], "_raw")
head(fq_stats_raw)

# Process R1 reads
fq_stats_F_cutITS86F <- read.delim("/data/wangxinyu/HMP_ITS/cutadapt/fq_stats_F_cutITS86F_relax.tsv",
                             sep = "\t",
                             header = T)

fq_stats_F_cutITS86F <- fq_stats_F_cutITS86F %>% select("file", "num_seqs")
colnames(fq_stats_F_cutITS86F)[2] <- paste0(colnames(fq_stats_F_cutITS86F)[2], "_cut")
head(fq_stats_F_cutITS86F)

# Merge R1 statistics
fq_stats_F_merged <- merge(fq_stats_raw, fq_stats_F_cutITS86F, by = "file")
head(fq_stats_F_merged)

# Filter samples with raw reads >= 1000
fq_stats_F_merged_filter <- fq_stats_F_merged %>% filter(num_seqs_raw >= 1000)

nrow(fq_stats_F_merged) # 324
nrow(fq_stats_F_merged_filter) # 310

# Calculate cut and uncut ratios for R1
fq_stats_F_merged_filter <- fq_stats_F_merged_filter %>%
  mutate(cut_ratio = num_seqs_cut / num_seqs_raw,
         uncut_ratio = 1 - cut_ratio,
         file = reorder(file, -cut_ratio))

# Reshape data for plotting
df_long <- fq_stats_F_merged_filter %>%
  pivot_longer(cols = c(cut_ratio, uncut_ratio),
               names_to = "type",
               values_to = "ratio") %>%
  mutate(type = factor(type, levels = c("uncut_ratio", "cut_ratio"),
         labels = c("ITS86F-", "ITS86F+")))

# Create stacked bar plot for R1
dis1 <- ggplot(df_long, aes(x = file, y = ratio, fill = type)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  scale_y_continuous(labels = scales::percent) +  
  labs(x = "R1 fastq files (num of reads > 1000, n = 310) ",
       y = "Percentage of reads",
       fill = "Reads type") +
  scale_fill_manual(values = c("ITS86F+" = "#e7ecef", "ITS86F-" = "#274c77")) +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "black") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.2, 0.2))


ggsave(file.path("/data/wangxinyu/HMP_ITS/figure", "fig_2b.pdf"), dis1, width = 4.5, height = 4.5)

# 统计 ratio > 0.75 的样本比例
high_ratio_file_count <- df_long %>%
  filter(type == "ITS86F+") %>%
  filter(ratio >= 0.75) %>%
  distinct(file) %>%
  nrow()
high_ratio_file_count / nrow(fq_stats_F_merged_filter) # 0.813 for default set and  0.903 for relaxed

# Process R2 reads
fq_stats_R_cutLR0B <- read.delim("/data/wangxinyu/HMP_ITS/cutadapt/fq_stats_R_cutLR0B_relax.tsv",
                             sep = "\t",
                             header = T)
fq_stats_R_cutLR0B <- fq_stats_R_cutLR0B %>% select("file", "num_seqs")
colnames(fq_stats_R_cutLR0B)[2] <- paste0(colnames(fq_stats_R_cutLR0B)[2], "_cut")
head(fq_stats_R_cutLR0B)

# Merge R2 statistics
fq_stats_R_merged <- merge(fq_stats_raw, fq_stats_R_cutLR0B, by = "file", all.y = TRUE)
head(fq_stats_R_merged)

# Filter samples with raw reads >= 1000
fq_stats_R_merged_filter <- fq_stats_R_merged %>% filter(num_seqs_raw >= 1000)

nrow(fq_stats_R_merged)
nrow(fq_stats_R_merged_filter) # 310

# Calculate cut and uncut ratios for R2
fq_stats_R_merged_filter <- fq_stats_R_merged_filter %>%
  mutate(cut_ratio = num_seqs_cut / num_seqs_raw,
         uncut_ratio = 1 - cut_ratio,
         file = reorder(file, -cut_ratio))

# Reshape data for plotting
df_long <- fq_stats_R_merged_filter %>%
  pivot_longer(cols = c(cut_ratio, uncut_ratio),
               names_to = "type",
               values_to = "ratio") %>%
  mutate(type = factor(type, levels = c("uncut_ratio", "cut_ratio"),
         labels = c("LR0B-", "LR0B+")))


# Create stacked bar plot for R2
dis2 <- ggplot(df_long, aes(x = file, y = ratio, fill = type)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  scale_y_continuous(labels = scales::percent) +  
  labs(x = "R2 fastq files (num of reads > 1000, n = 310) ",
       y = "Percentage of reads",
       fill = "Reads type") +
  scale_fill_manual(values = c("LR0B+" = "#e7ecef", "LR0B-" = "#274c77")) +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "black") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.2, 0.2))

ggsave(file.path("/data/wangxinyu/HMP_ITS/figure", "fig_2c.pdf"), dis2, width = 4.5, height = 4.5)

# 统计 ratio > 0.75 的样本比例
high_ratio_file_count <- df_long %>%
  filter(type == "LR0B+") %>%
  filter(ratio >= 0.75) %>%
  distinct(file) %>%
  nrow()
high_ratio_file_count / nrow(fq_stats_R_merged_filter) # 0.216 for default and 0.477 for relax

# Combine and save Figure 2
fig2 <- dis1 + dis2
ggsave(file.path("/data/wangxinyu/HMP_ITS/figure", "fig_2bc.pdf"), fig2, width = 9, height = 4.5)



