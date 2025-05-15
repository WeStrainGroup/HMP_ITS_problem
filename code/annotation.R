#!/usr/bin/env Rscript

# Set variables
project0 <- "PRJNA356769_ITS2"
project <- paste0(project0, "_")
path_in <- "/data/wangxinyu/HMP_ITS/data/PRJNA356769_ITS2/0_raw"
path_out <- "/data/wangxinyu/HMP_ITS/data/PRJNA356769_ITS2"
if(!dir.exists(path_out)) {dir.create(path_out)}
marker <- "ITS2"
pattern_F = ".1.fastq.gz"
pattern_R = ".2.fastq.gz"
maxEE <- 5
itsx_e <- 0.01
ref <- "/data/wangxinyu/ITS/Ref/Unite/sh_general_release_all_19.02.2025"
minBoot <- 50
vsearch_id = 0.985
filter_abundance <- 1/10000
filter_depth <- 5000
thread <- 100


# Load required R packages
suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(data.table)
  library(fs)
  library(dada2)
  library(Biostrings)
  library(phyloseq)
  library(microbiome)
  library(vegan)
})


# Get tool paths dynamically (via Conda environment PATH) 
seqkit <- Sys.which("seqkit")
itsx <- Sys.which("ITSx")
vsearch <- Sys.which("vsearch")

# Set ITSx environment variable
itsx_bin <- dirname(itsx)
Sys.setenv(PATH = paste(Sys.getenv("PATH"), itsx_bin, sep = ":"))

# Define function for FASTQ quality assessment and statistics
check_fq_stats <- function(F, R) {
    title <- paste0(project, "fq_stats_", sub(".*_", "", substitute(F)))
    p_in <- path_common(path_dir(c(F, R))) # Find common path
    p_out <- file.path(path_check, paste0(title, ".tsv"))
    arg_seqkit <- c("stats", "-a", "-b", "-e", "-T", "-j", thread,
                    paste0(p_in, "/*"),
                    paste(">", p_out))
    system2(seqkit, arg_seqkit)

    fa_raw_check <- read.delim(p_out, header = TRUE)     
    his1 <- ggplot(fa_raw_check, aes(x = num_seqs)) +
        geom_histogram(aes(y = after_stat(count)/sum(after_stat(count))), fill = "grey", color = "black", alpha = 0.75) +
        labs(title = title,
             x = "Num_seqs of each sample",
             y = "Proportion") +
        theme_bw()
    ggsave(file.path(path_check,  paste0(title, "1.jpg")), his1, width = 6, height = 6)

    his2 <- ggplot(fa_raw_check, aes(x = avg_len)) +
            geom_histogram(aes(y = after_stat(count)/sum(after_stat(count))), fill = "grey", color = "black", alpha = 0.75) +
            labs(title = title,
                x = "Avg_len of each sample",
                y = "Proportion") +
            theme_bw()
      ggsave(file.path(path_check,  paste0(title, "2.jpg")), his2, width = 6, height = 6)

    his3 <- ggplot(fa_raw_check, aes(x = AvgQual)) +
            geom_histogram(aes(y = after_stat(count)/sum(after_stat(count))), fill = "grey", color = "black", alpha = 0.75) +
            labs(title = title,
                x = "Avg_qual of each sample",
                y = "Proportion") +
            theme_bw()
    ggsave(file.path(path_check, paste0(title, "3.jpg")), his3, width = 6, height = 6)
}


save_quality_plot <- function(F, R) {
  title_F <- paste0(project, "quality_", substitute(F))
  title_R <- paste0(project, "quality_", substitute(R))

  p_quality_R1 <- plotQualityProfile(F[1:9]) + labs(title = title_F)
  p_quality_R2 <- plotQualityProfile(R[1:9]) + labs(title = title_R)
  
  ggsave(file.path(path_check, paste0(title_F, ".jpg")), p_quality_R1, width = 6, height = 6)
  ggsave(file.path(path_check, paste0(title_R, ".jpg")), p_quality_R2, width = 6, height = 6)
}


# Echo both command and output
options(echo = TRUE)

# --- 1 DADA2-based sequence processing pipeline ---
path1 <- file.path(path_out, "1_dada2")
if(!dir.exists(path1)) {dir.create(path1)}

# Input FASTQ files
fqF_raw <- sort(list.files(path_in, pattern = pattern_F, full.names = TRUE))
fqR_raw <- sort(list.files(path_in, pattern = pattern_R, full.names = TRUE))

length(fqF_raw)
head(fqF_raw)
length(fqR_raw)
head(fqR_raw)

path_check <- file.path(path1, "check")
if(!dir.exists(path_check)) {dir.create(path_check)}

check_fq_stats(fqF_raw, fqR_raw)
save_quality_plot(fqF_raw, fqR_raw)

# --- 1.1 Quality filtering and trimming using DADA2 ---
path_filtered <- file.path(path1, "filtered")
if(!dir.exists(path_filtered)) dir.create(path_filtered)
fqF_filtered <- file.path(path_filtered, basename(fqF_raw))
fqR_filtered <- file.path(path_filtered, basename(fqR_raw))

out <- filterAndTrim(fqF_raw, fqF_filtered, fqR_raw, fqR_filtered,
                     truncQ = 2, maxN = 0, maxEE = c(maxEE, maxEE), 
                     rm.phix = TRUE, compress = T, multithread = thread, verbose = TRUE)  

rownames(out) <- gsub(pattern_F, "", rownames(out))
head(out)

# Keep only FASTQ files with reads
fqF_filtered <- sort(list.files(path_filtered, pattern = pattern_F, full.names = TRUE))
fqR_filtered <- sort(list.files(path_filtered, pattern = pattern_R, full.names = TRUE))

# Check
check_fq_stats(fqF_filtered, fqR_filtered)
save_quality_plot(fqF_filtered, fqR_filtered)

# --- 1.2 DADA2 error rate estimation and denoising ---
# Learn the Error Rates
errF <- learnErrors(fqF_filtered, nbases = 1e9, randomize = TRUE, multithread = thread) 
errR <- learnErrors(fqR_filtered, nbases = 1e9, randomize = TRUE, multithread = thread) 

# Sample inference to the de-replicated data
dadaFs <- dada(fqF_filtered, err = errF, multithread = thread)
dadaRs <- dada(fqR_filtered, err = errR, multithread = thread)

dadaFs[[1]] 
dadaRs[[1]] 

# Plot error rates
p_error_F <- plotErrors(errF, nominalQ = TRUE) + labs(title = paste0(project, "error_F"))
p_error_R <- plotErrors(errR, nominalQ = TRUE) + labs(title = paste0(project, "error_R"))

ggsave(file.path(path_check, paste0(project, "error_F.png")), p_error_F, width = 6, height = 6) 
ggsave(file.path(path_check, paste0(project, "error_R.png")), p_error_R, width = 6, height = 6) 

# --- 1.3 Paired-end read merging and ASV table construction ---
# Merge paired reads
mergers <- mergePairs(dadaFs, fqF_filtered, dadaRs, fqR_filtered, trimOverhang = TRUE, verbose = TRUE)
head(mergers[[1]]) 

# Construct sequence (ASV) table
seqtab <- makeSequenceTable(mergers)
rownames(seqtab) <- gsub(pattern_F, "", rownames(seqtab))

dim(seqtab)
table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths

# Write output
path_result <- file.path(path1, "result")
if(!dir.exists(path_result)) dir.create(path_result)

write.csv(seqtab, file = file.path(path_result, paste0(project, "seqtab_raw.csv")))

# --- 1.4 ITS region extraction and sequence dereplication ---
path_itsx <- file.path(path1, "itsx")
if(!dir.exists(path_itsx)) dir.create(path_itsx)

seqtab2 <- as.data.frame(t(seqtab))

# Add ASV tags
seqtab2$seq <- row.names(seqtab2)
row.names(seqtab2) <- paste0(project, "ASV", seq_along(1:nrow(seqtab2)))
colnames(seqtab2)[1:6]

fasta <- DNAStringSet(seqtab2$seq)
names(fasta) <- row.names(seqtab2)
fasta
writeXStringSet(fasta, file.path(path_itsx, paste0(project, "asv_raw.fasta")), format = "fasta")

# Run ITSx
arg_ITSx <- c("-i", file.path(path_itsx, paste0(project, "asv_raw.fasta")),
              "-o", file.path(path_itsx, paste0(project, "asv_itsx")),
              "--detailed_results T",
              "-E 0.01", 
              "--cpu", thread,
              "--table T") # Output a table directly

system2(itsx, arg_ITSx)

ASV_itsx <- readDNAStringSet(file.path(path_itsx, paste0(project, "asv_itsx.", marker, ".fasta")), format = "fasta")
ASV_itsx

# Keep only sequences longer than 50 bp to prevent annotation errors
ASV_itsx_filter <- ASV_itsx[width(ASV_itsx) >= 50]
ASV_itsx_filter

table(ASV_itsx@ranges@width)
table(ASV_itsx_filter@ranges@width)

seq_itsx <- data.frame(seq_itsx = as.character(ASV_itsx_filter), ASV_itsx = names(ASV_itsx_filter))
row.names(seq_itsx) <- gsub("\\|.*", "", seq_itsx$ASV_itsx)
row.names(seq_itsx)[1:6]

# Merge based on ASV tags
seqtab3 <- merge(seqtab2, seq_itsx, by = "row.names", all.y = TRUE) # Keep only sequences with ITSx regions
dim(seqtab2)
dim(seqtab3)
head(colnames(seqtab3))
tail(colnames(seqtab3))

seqtab3 <- seqtab3 %>% select(-Row.names, -seq, -ASV_itsx)
length(seqtab3$seq_itsx)
length(unique(seqtab3$seq_itsx))

# Merge identical seq_its2 sequences
setDT(seqtab3)
head(sapply(seqtab3, class)) # Check object types
tail(sapply(seqtab3, class))

seqtab4 <- seqtab3[, lapply(.SD, sum, na.rm = TRUE), by = seq_itsx, .SDcols = where(is.numeric)]
seqtab4 <- as.data.frame(seqtab4)

# Convert sequence to row names and transpose
seqtab4 <- seqtab4 %>%
     column_to_rownames(var = "seq_itsx")

seqtab5 <- t(seqtab4)
class(seqtab5)

head(row.names(seqtab5))
tail(row.names(seqtab5))
head(colnames(seqtab5))
tail(colnames(seqtab5))

# --- 1.5 Chimera removal and final ASV table generation ---
seqtab.nochim <- removeBimeraDenovo(seqtab5, method = "consensus", multithread = thread, verbose = TRUE)

dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab5)
table(nchar(getSequences(seqtab.nochim)))

# Update sample names 
row.names(seqtab.nochim) <- gsub(pattern_F, "", row.names(seqtab.nochim))
head(row.names(seqtab.nochim))

# Write output
write.csv(seqtab.nochim, file = file.path(path_result, paste0(project, "seqtab.csv")))
rm(seqtab2, seqtab3, seqtab4, seqtab5) # Free up memory

# --- 1.6 Taxonomic assignment ---
seqs <- colnames(seqtab.nochim)  
head(seqs)

batch_size <- thread * 100
batch_size
num_batches <- ceiling(length(seqs) / batch_size) 
num_batches

taxa_list <- list()
boot_list <- list()

for (i in 1:num_batches) {
  up <- min(i * batch_size, length(seqs)) 
  down <- (i - 1) * batch_size + 1
  seqs_batch <- seqs[down:up] # Select sequences for this batch
  taxa_batch <- assignTaxonomy(seqs_batch, ref, minBoot = minBoot, outputBootstraps = TRUE, multithread = thread, tryRC = TRUE, verbose = TRUE)
  
  taxa_list[[i]] <- taxa_batch[["tax"]]
  boot_list[[i]] <- taxa_batch[["boot"]]
  
  rm(seqs_batch, taxa_batch) # Free up memory
  cat(paste0("\n", "-- batch", i, " / ", num_batches,": seq", down, "-", up ," was down! --", "\n"))
}

taxa <- do.call(rbind, taxa_list)
boot <- do.call(rbind, boot_list)

nrow(taxa)
nrow(boot)
sum(is.na(taxa[, "Species"]))
sum(is.na(taxa[, "Species"]))/nrow(taxa)
sum(is.na(taxa[, "Genus"]))
sum(is.na(taxa[, "Genus"]))/nrow(taxa)

# Write output
write.csv(taxa, file = file.path(path_result, paste0(project, "taxa.csv")))
write.csv(boot, file = file.path(path_result, paste0(project, "taxa_boot.csv")))

# --- 2 VSEARCH-based sequence clustering pipeline ---
path2 <- file.path(path_out, "2_vsearch")
if(!dir.exists(path2)) {dir.create(path2)}

# --- 2.1 ASV table and FASTA file preparation for clustering ---
# Add ASV tags
seqtab6 <- as.data.frame(t(seqtab.nochim))
seqtab6$ASV <- paste0(project, "ASV", seq_along(1:nrow(seqtab6)))

head(colnames(seqtab6))
tail(colnames(seqtab6))

# Write sequences and ASV tags to DNAStringSet
fasta <- DNAStringSet(row.names(seqtab6))
names(fasta) <- seqtab6$ASV
fasta

writeXStringSet(fasta, file.path(path2, paste0(project, "ASV_unclustered.fasta")), format = "fasta")

# --- 2.2 VSEARCH clustering execution ---
arg <- c("--cluster_fast", file.path(path2, paste0(project, "ASV_unclustered.fasta")),
         "--centroids", file.path(path2, paste0(project, "ASV_c", vsearch_id*100,".fasta")),
         "--uc", file.path(path2, paste0(project, "ASV_c", vsearch_id*100,"_index1.tsv")),
         "--otutabout", file.path(path2, paste0(project, "ASV_c", vsearch_id*100,"_index2.tsv")),
         "--id", vsearch_id,
         "--strand both",
         "--threads", thread)

system2(vsearch, arg)

# --- 2.3 ASV table reconstruction based on clustering results ---
cluster_index0 <- read.delim(file.path(path2, paste0(project, "ASV_c", vsearch_id*100, "_index1.tsv")), header = F) 
cluster_index <- cluster_index0 %>%
        filter(V1 != "C") %>% 
        mutate(V10 = ifelse(V10 == "*", V9, V10)) %>% 
        select(V9, V10) %>%
        dplyr::rename(ASV = V9, rASV = V10) 

nrow(cluster_index)
length(unique(cluster_index$rASV))

# Add index information to seqtab
seqtab7 <- merge(cluster_index, seqtab6, by = "ASV", all.y = TRUE)
seqtab7 <- seqtab7 %>%
  select(-ASV)

# Check
head(colnames(seqtab7))
head(seqtab7$rASV)

# Merge identical rASVs
setDT(seqtab7)

seqtab8 <- seqtab7[, lapply(.SD, sum, na.rm = TRUE), by = rASV, .SDcols = where(is.numeric)]
seqtab8 <- as.data.frame(seqtab8)

# Convert sequence to row names and transpose
seqtab8 <- seqtab8 %>% column_to_rownames(var = "rASV")

seqtab8[1:3, 1:3]

seqtab_cluster <- t(seqtab8)
seqtab_cluster[1:3, 1:3]
dim(seqtab_cluster)

# --- 2.4 Taxonomic table refinement for clustered ASVs ---
rASV_list <- colnames(seqtab_cluster)
length(rASV_list)

taxa2 <- merge(seqtab6, taxa, by = "row.names", all.y = TRUE)
dim(taxa2)
head(colnames(taxa2))
tail(colnames(taxa2))

taxa3 <- taxa2 %>%
  select("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

tax_cluster <- taxa3[taxa3$ASV %in% rASV_list, ]
rownames(tax_cluster) <- tax_cluster$ASV
tax_cluster <- tax_cluster %>% select(-ASV)

dim(tax_cluster)
head(tax_cluster)

# Write out
write.csv(seqtab_cluster, file.path(path2, paste0(project,"seqtab_c", vsearch_id*100,".csv")))
write.csv(tax_cluster, file.path(path2, paste0(project, "taxa_c", vsearch_id*100,".csv")))

# --- 3 Phyloseq-based community analysis pipeline ---
path3 <- file.path(path_out, "3_phyloseq")
if(!dir.exists(path3)) {dir.create(path3)}

# Construct the phyloseq subject
otu_table <- phyloseq::otu_table(as.matrix(seqtab_cluster), taxa_are_rows = FALSE)

sam_data <- data.frame(row.names = rownames(otu_table), project = rep(project0, nrow(otu_table)))
sam_data <- sample_data(sam_data)

tax_table <- phyloseq::tax_table(as.matrix(tax_cluster))

refseq <- readDNAStringSet(file.path(path2, paste0(project, "ASV_c", vsearch_id*100,".fasta")), format = "fasta")

ps <- phyloseq(otu_table,
               sam_data,
               tax_table,
               refseq)
ps
microbiome::summarize_phyloseq(ps)

# --- 3.1 Global community analysis and basic filtering ---
path_all <- file.path(path3, "all")
if(!dir.exists(path_all)) dir.create(path_all)

# --- 3.1.1 Do general filteration  --- 
# Remove ASVs with unclear phylum annotation
ps_filter <- subset_taxa(ps, Phylum != "NA" & Kingdom != "k__Eukaryota_kgd_Incertae_sedis")

ntaxa(ps) 
ntaxa(ps_filter) 

# Set ASV counts to 0 in samples where relative abundance is below filter_abundance
otu <- ps_filter@otu_table
otu[otu / rowSums(otu) < filter_abundance] <- 0 
otu_table(ps_filter) <- otu # Return otu table to ps
ps_filter <- filter_taxa(ps_filter, function(otu) sum(otu) > 1, TRUE) # Remove ASVs that are all zeros after filtering
ps_filter <- prune_samples(sample_sums(ps_filter) > 0, ps_filter) # Remove samples without any reads

ntaxa(ps)
ntaxa(ps_filter)
nsamples(ps)
nsamples(ps_filter)

# --- 3.1.2 Write output  --- 
path_all
write.csv(readcount(ps_filter), file.path(path_all, paste0(project, "readcount_all.csv")))
writeXStringSet(ps_filter@refseq, file = file.path(path_all, paste0(project, "refseq_all.fasta")), format = "fasta")
write.csv(ps_filter@tax_table, file.path(path_all, paste0(project, "taxa_all.csv")))
write.csv(ps_filter@otu_table, file.path(path_all, paste0(project,"otu_all.csv")))

# Alpha diversity
try({diversity_all <- estimate_richness(ps_filter, split = TRUE, measures = NULL)
     write.csv(diversity_all, file.path(path_all, paste0(project, "diversity_all.csv")))
     }, silent = F)


# --- 3.2 Fungal community related ---
path_fungi <- file.path(path3, "fungi")
if(!dir.exists(path_fungi)) dir.create(path_fungi)

ps_fungi <- subset_taxa(ps_filter, Kingdom == "k__Fungi")

# Check
ntaxa(ps_fungi) 
get_taxa_unique(ps_fungi, "Phylum") 
microbiome::summarize_phyloseq(ps_fungi) 

# --- 3.2.1 filter mushroom taxa from ps subject ---
# Identify mushroom taxa
mushroom_list <- read.csv(file.path(Sys.getenv("CONDA_PREFIX"), "lib/R/library/fungap/ref/mushroom_genus.csv"))
mushroom_list <- mushroom_list$Genus
mushroom_list <- paste0("g__", mushroom_list)
head(mushroom_list)

taxa_table <- as.data.frame(ps_fungi@tax_table)
taxa_table_mushroom <- taxa_table %>% filter(Genus %in% mushroom_list)
head(taxa_table_mushroom)

nrow(taxa_table) # All ASVs
nrow(taxa_table_mushroom) # ASVs of mushrooms
length(unique(taxa_table$Genus)) # All genera
length(unique(taxa_table_mushroom$Genus)) # Genera of mushrooms

# Remove identified ASVs
ASV_to_drop <- rownames(taxa_table_mushroom)
head(ASV_to_drop)
ps_fungi2 <- prune_taxa(!(taxa_names(ps_fungi) %in% ASV_to_drop), ps_fungi)

# Check
ntaxa(ps_fungi) 
ntaxa(ps_fungi2) 
microbiome::summarize_phyloseq(ps_fungi2)

# --- 3.2.2 Write output without read depth filteration  --- 
# ASV level output
path_fungi_ASV <- file.path(path_fungi, "ASV")
if(!dir.exists(path_fungi_ASV)) dir.create(path_fungi_ASV)
path_fungi_ASV

write.csv(readcount(ps_fungi2), file.path(path_fungi_ASV, paste0(project, "readcount_fungi.csv")))
writeXStringSet(ps_fungi2@refseq, file = file.path(path_fungi_ASV, paste0(project, "refseq_fungi.fasta")), format = "fasta")
write.csv(ps_fungi2@tax_table, file.path(path_fungi_ASV, paste0(project, "taxa_fungi.csv")))
write.csv(ps_fungi2@otu_table, file.path(path_fungi_ASV, paste0(project,"otu_fungi.csv")))

# Alpha diversity
try({diversity_fungi <- estimate_richness(ps_fungi2, split = TRUE, measures = NULL)
     write.csv(diversity_fungi, file.path(path_fungi_ASV, paste0(project, "diversity_fungi.csv")))
     })

# Genus level output
path_fungi_gen <- file.path(path_fungi, "genus")
if(!dir.exists(path_fungi_gen)) dir.create(path_fungi_gen)
path_fungi_gen

ps_fungi_gen_p0 <- aggregate_taxa(ps_fungi2, level = "Genus")
ps_fungi_gen_p1 <- aggregate_rare(ps_fungi2, level = "Genus", detection = 0, prevalence = 1/100)
ps_fungi_gen_p5 <- aggregate_rare(ps_fungi2, level = "Genus", detection = 0, prevalence = 5/100)
ps_fungi_gen_p10 <- aggregate_rare(ps_fungi2, level = "Genus", detection = 0, prevalence = 10/100)

otu_gen_p0_clr <- vegan::decostand(ps_fungi_gen_p0@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
otu_gen_p1_clr <- vegan::decostand(ps_fungi_gen_p1@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
otu_gen_p5_clr <- vegan::decostand(ps_fungi_gen_p5@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
otu_gen_p10_clr <- vegan::decostand(ps_fungi_gen_p10@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)

write.csv(ps_fungi_gen_p0@tax_table, file.path(path_fungi_gen, paste0(project,"taxa_fungi_gen.csv")))
write.csv(ps_fungi_gen_p0@otu_table, file.path(path_fungi_gen, paste0(project,"otu_fungi_gen_p0.csv")))
write.csv(ps_fungi_gen_p1@otu_table, file.path(path_fungi_gen, paste0(project,"otu_fungi_gen_p1.csv")))
write.csv(ps_fungi_gen_p5@otu_table, file.path(path_fungi_gen, paste0(project,"otu_fungi_gen_p5.csv")))
write.csv(ps_fungi_gen_p10@otu_table, file.path(path_fungi_gen, paste0(project,"otu_fungi_gen_p10.csv")))

write.csv(otu_gen_p0_clr, file.path(path_fungi_gen, paste0(project,"otu_fungi_gen_p0_clr.csv")))
write.csv(otu_gen_p1_clr, file.path(path_fungi_gen, paste0(project,"otu_fungi_gen_p1_clr.csv")))
write.csv(otu_gen_p5_clr, file.path(path_fungi_gen, paste0(project,"otu_fungi_gen_p5_clr.csv")))
write.csv(otu_gen_p10_clr, file.path(path_fungi_gen, paste0(project,"otu_fungi_gen_p10_clr.csv")))

# Phylum level output
path_fungi_phy <- file.path(path_fungi, "phylum")
if(!dir.exists(path_fungi_phy)) dir.create(path_fungi_phy)
path_fungi_phy

ps_fungi_phy_p0 <- aggregate_taxa(ps_fungi2, level = "Phylum")
ps_fungi_phy_p1 <- aggregate_rare(ps_fungi2, level = "Phylum", detection = 0, prevalence = 1/100)
ps_fungi_phy_p5 <- aggregate_rare(ps_fungi2, level = "Phylum", detection = 0, prevalence = 5/100)
ps_fungi_phy_p10 <- aggregate_rare(ps_fungi2, level = "Phylum", detection = 0, prevalence = 10/100)

otu_phy_p0_clr <- vegan::decostand(ps_fungi_phy_p0@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
otu_phy_p1_clr <- vegan::decostand(ps_fungi_phy_p1@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
otu_phy_p5_clr <- vegan::decostand(ps_fungi_phy_p5@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
otu_phy_p10_clr <- vegan::decostand(ps_fungi_phy_p10@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)

write.csv(ps_fungi_phy_p0@tax_table, file.path(path_fungi_phy, paste0(project,"taxa_fungi_phy.csv")))
write.csv(ps_fungi_phy_p0@otu_table, file.path(path_fungi_phy, paste0(project,"otu_fungi_phy_p0.csv")))
write.csv(ps_fungi_phy_p1@otu_table, file.path(path_fungi_phy, paste0(project,"otu_fungi_phy_p1.csv")))
write.csv(ps_fungi_phy_p5@otu_table, file.path(path_fungi_phy, paste0(project,"otu_fungi_phy_p5.csv")))
write.csv(ps_fungi_phy_p10@otu_table, file.path(path_fungi_phy, paste0(project,"otu_fungi_phy_p10.csv")))

write.csv(otu_phy_p0_clr, file.path(path_fungi_phy, paste0(project,"otu_fungi_phy_p0_clr.csv")))
write.csv(otu_phy_p1_clr, file.path(path_fungi_phy, paste0(project,"otu_fungi_phy_p1_clr.csv")))
write.csv(otu_phy_p5_clr, file.path(path_fungi_phy, paste0(project,"otu_fungi_phy_p5_clr.csv")))
write.csv(otu_phy_p10_clr, file.path(path_fungi_phy, paste0(project,"otu_fungi_phy_p10_clr.csv")))

# --- 3.2.3 Write output with read depth filteration  --- 
# Remove samples with insufficient reads
ps_fungi3 <- prune_samples(sample_sums(ps_fungi2) >= filter_depth, ps_fungi2)
ps_fungi3 <- filter_taxa(ps_fungi3, function(otu) sum(otu) > 1, TRUE) # Remove ASVs that are all zeros after filtering

# Check
ntaxa(ps_fungi2) 
ntaxa(ps_fungi3) 
nsamples(ps_fungi2) 
nsamples(ps_fungi3) 
microbiome::summarize_phyloseq(ps_fungi3)

# ASV level output
path_fungi_ASV <- file.path(path_fungi, "ASV")
if(!dir.exists(path_fungi_ASV)) dir.create(path_fungi_ASV)
path_fungi_ASV

write.csv(readcount(ps_fungi3), file.path(path_fungi_ASV, paste0(project, "readcount_fungi_filterdp.csv")))
writeXStringSet(ps_fungi3@refseq, file = file.path(path_fungi_ASV, paste0(project, "refseq_fungi_filterdp.fasta")), format = "fasta")
write.csv(ps_fungi3@tax_table, file.path(path_fungi_ASV, paste0(project, "taxa_fungi_filterdp.csv")))
write.csv(ps_fungi3@otu_table, file.path(path_fungi_ASV, paste0(project,"otu_fungi_filterdp.csv")))

# Alpha diversity
try({diversity_fungi_filterdp <- estimate_richness(ps_fungi3, split = TRUE, measures = NULL)
     write.csv(diversity_fungi_filterdp, file.path(path_fungi_ASV, paste0(project, "diversity_fungi_filterdp.csv")))
     })

# Genus level output
path_fungi_gen <- file.path(path_fungi, "genus")
if(!dir.exists(path_fungi_gen)) dir.create(path_fungi_gen)
path_fungi_gen

ps_fungi_gen_p0 <- aggregate_taxa(ps_fungi3, level = "Genus")
ps_fungi_gen_p1 <- aggregate_rare(ps_fungi3, level = "Genus", detection = 0, prevalence = 1/100)
ps_fungi_gen_p5 <- aggregate_rare(ps_fungi3, level = "Genus", detection = 0, prevalence = 5/100)
ps_fungi_gen_p10 <- aggregate_rare(ps_fungi3, level = "Genus", detection = 0, prevalence = 10/100)

otu_gen_p0_clr <- vegan::decostand(ps_fungi_gen_p0@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
otu_gen_p1_clr <- vegan::decostand(ps_fungi_gen_p1@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
otu_gen_p5_clr <- vegan::decostand(ps_fungi_gen_p5@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
otu_gen_p10_clr <- vegan::decostand(ps_fungi_gen_p10@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)

write.csv(ps_fungi_gen_p0@tax_table, file.path(path_fungi_gen, paste0(project,"taxa_fungi_gen_filterdp.csv")))
write.csv(ps_fungi_gen_p0@otu_table, file.path(path_fungi_gen, paste0(project,"otu_fungi_gen_p0_filterdp.csv")))
write.csv(ps_fungi_gen_p1@otu_table, file.path(path_fungi_gen, paste0(project,"otu_fungi_gen_p1_filterdp.csv")))
write.csv(ps_fungi_gen_p5@otu_table, file.path(path_fungi_gen, paste0(project,"otu_fungi_gen_p5_filterdp.csv")))
write.csv(ps_fungi_gen_p10@otu_table, file.path(path_fungi_gen, paste0(project,"otu_fungi_gen_p10_filterdp.csv")))

write.csv(otu_gen_p0_clr, file.path(path_fungi_gen, paste0(project,"otu_fungi_gen_p0_clr_filterdp.csv")))
write.csv(otu_gen_p1_clr, file.path(path_fungi_gen, paste0(project,"otu_fungi_gen_p1_clr_filterdp.csv")))
write.csv(otu_gen_p5_clr, file.path(path_fungi_gen, paste0(project,"otu_fungi_gen_p5_clr_filterdp.csv")))
write.csv(otu_gen_p10_clr, file.path(path_fungi_gen, paste0(project,"otu_fungi_gen_p10_clr_filterdp.csv")))

# Phylum level output
path_fungi_phy <- file.path(path_fungi, "phylum")
if(!dir.exists(path_fungi_phy)) dir.create(path_fungi_phy)
path_fungi_phy

ps_fungi_phy_p0 <- aggregate_taxa(ps_fungi3, level = "Phylum")
ps_fungi_phy_p1 <- aggregate_rare(ps_fungi3, level = "Phylum", detection = 0, prevalence = 1/100)
ps_fungi_phy_p5 <- aggregate_rare(ps_fungi3, level = "Phylum", detection = 0, prevalence = 5/100)
ps_fungi_phy_p10 <- aggregate_rare(ps_fungi3, level = "Phylum", detection = 0, prevalence = 10/100)

otu_phy_p0_clr <- vegan::decostand(ps_fungi_phy_p0@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
otu_phy_p1_clr <- vegan::decostand(ps_fungi_phy_p1@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
otu_phy_p5_clr <- vegan::decostand(ps_fungi_phy_p5@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
otu_phy_p10_clr <- vegan::decostand(ps_fungi_phy_p10@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)

write.csv(ps_fungi_phy_p0@tax_table, file.path(path_fungi_phy, paste0(project,"taxa_fungi_phy_filterdp.csv")))
write.csv(ps_fungi_phy_p0@otu_table, file.path(path_fungi_phy, paste0(project,"otu_fungi_phy_p0_filterdp.csv")))
write.csv(ps_fungi_phy_p1@otu_table, file.path(path_fungi_phy, paste0(project,"otu_fungi_phy_p1_filterdp.csv")))
write.csv(ps_fungi_phy_p5@otu_table, file.path(path_fungi_phy, paste0(project,"otu_fungi_phy_p5_filterdp.csv")))
write.csv(ps_fungi_phy_p10@otu_table, file.path(path_fungi_phy, paste0(project,"otu_fungi_phy_p10_filterdp.csv")))

write.csv(otu_phy_p0_clr, file.path(path_fungi_phy, paste0(project,"otu_fungi_phy_p0_clr_filterdp.csv")))
write.csv(otu_phy_p1_clr, file.path(path_fungi_phy, paste0(project,"otu_fungi_phy_p1_clr_filterdp.csv")))
write.csv(otu_phy_p5_clr, file.path(path_fungi_phy, paste0(project,"otu_fungi_phy_p5_clr_filterdp.csv")))
write.csv(otu_phy_p10_clr, file.path(path_fungi_phy, paste0(project,"otu_fungi_phy_p10_clr_filterdp.csv")))


# --- 3.3 Vegetable-related ---
path_vegetable <- file.path(path3, "vegetable")
if(!dir.exists(path_vegetable)) dir.create(path_vegetable)

# --- 3.3.1 Mushroom-related ---
try <- try({ps_mushroom <- prune_taxa((taxa_names(ps_fungi) %in% ASV_to_drop), ps_fungi)})

if (inherits(try, "try-error")) {
  message("Note: No mushroom reads were detected for analysis\n")
} else {
    try({
    path_mushroom <- file.path(path_vegetable, "mushroom")
    if(!dir.exists(path_mushroom)) dir.create(path_mushroom)

    # ASV level out
    path_mushroom_ASV <- file.path(path_mushroom, "ASV")
    if(!dir.exists(path_mushroom_ASV)) dir.create(path_mushroom_ASV)
    path_mushroom_ASV

    write.csv(readcount(ps_mushroom), file.path(path_mushroom_ASV, paste0(project, "readcount_mushroom.csv")))
    writeXStringSet(ps_mushroom@refseq, file = file.path(path_mushroom_ASV, paste0(project, "refseq_mushroom.fasta")), format = "fasta")
    write.csv(ps_mushroom@tax_table, file.path(path_mushroom_ASV, paste0(project, "taxa_mushroom.csv")))
    write.csv(ps_mushroom@otu_table, file.path(path_mushroom_ASV, paste0(project, "otu_mushroom.csv")))

    # Genus level out
    path_mushroom_gen <- file.path(path_mushroom, "genus")
    if(!dir.exists(path_mushroom_gen)) dir.create(path_mushroom_gen)
    path_mushroom_gen

    ps_mushroom_gen_p0 <- aggregate_taxa(ps_mushroom, level = "Genus")
    ps_mushroom_gen_p1 <- aggregate_rare(ps_mushroom, level = "Genus", detection = 0, prevalence = 1/100)
    ps_mushroom_gen_p5 <- aggregate_rare(ps_mushroom, level = "Genus", detection = 0, prevalence = 5/100)
    ps_mushroom_gen_p10 <- aggregate_rare(ps_mushroom, level = "Genus", detection = 0, prevalence = 10/100)

    otu_gen_p0_clr <- vegan::decostand(ps_mushroom_gen_p0@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
    otu_gen_p1_clr <- vegan::decostand(ps_mushroom_gen_p1@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
    otu_gen_p5_clr <- vegan::decostand(ps_mushroom_gen_p5@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
    otu_gen_p10_clr <- vegan::decostand(ps_mushroom_gen_p10@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)

    write.csv(ps_mushroom_gen_p0@tax_table, file.path(path_mushroom_gen, paste0(project, "taxa_mushroom_gen_p0.csv")))
    write.csv(ps_mushroom_gen_p1@tax_table, file.path(path_mushroom_gen, paste0(project, "taxa_mushroom_gen_p1.csv")))
    write.csv(ps_mushroom_gen_p5@tax_table, file.path(path_mushroom_gen, paste0(project, "taxa_mushroom_gen_p5.csv")))
    write.csv(ps_mushroom_gen_p10@tax_table, file.path(path_mushroom_gen, paste0(project, "taxa_mushroom_gen_p10.csv")))

    write.csv(ps_mushroom_gen_p0@otu_table, file.path(path_mushroom_gen, paste0(project, "otu_mushroom_gen_p0.csv")))
    write.csv(ps_mushroom_gen_p1@otu_table, file.path(path_mushroom_gen, paste0(project, "otu_mushroom_gen_p1.csv")))
    write.csv(ps_mushroom_gen_p5@otu_table, file.path(path_mushroom_gen, paste0(project, "otu_mushroom_gen_p5.csv")))
    write.csv(ps_mushroom_gen_p10@otu_table, file.path(path_mushroom_gen, paste0(project, "otu_mushroom_gen_p10.csv")))

    write.csv(otu_gen_p0_clr, file.path(path_mushroom_gen, paste0(project, "otu_mushroom_gen_p0_clr.csv")))
    write.csv(otu_gen_p1_clr, file.path(path_mushroom_gen, paste0(project, "otu_mushroom_gen_p1_clr.csv")))
    write.csv(otu_gen_p5_clr, file.path(path_mushroom_gen, paste0(project, "otu_mushroom_gen_p5_clr.csv")))
    write.csv(otu_gen_p10_clr, file.path(path_mushroom_gen, paste0(project, "otu_mushroom_gen_p10_clr.csv")))

    # Phylum level out
    path_mushroom_phy <- file.path(path_mushroom, "phylum")
    if(!dir.exists(path_mushroom_phy)) dir.create(path_mushroom_phy)
    path_mushroom_phy

    ps_mushroom_phy_p0 <- aggregate_taxa(ps_mushroom, level = "Phylum")
    ps_mushroom_phy_p1 <- aggregate_rare(ps_mushroom, level = "Phylum", detection = 0, prevalence = 1/100)
    ps_mushroom_phy_p5 <- aggregate_rare(ps_mushroom, level = "Phylum", detection = 0, prevalence = 5/100)
    ps_mushroom_phy_p10 <- aggregate_rare(ps_mushroom, level = "Phylum", detection = 0, prevalence = 10/100)

    otu_phy_p0_clr <- vegan::decostand(ps_mushroom_phy_p0@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
    otu_phy_p1_clr <- vegan::decostand(ps_mushroom_phy_p1@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
    otu_phy_p5_clr <- vegan::decostand(ps_mushroom_phy_p5@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
    otu_phy_p10_clr <- vegan::decostand(ps_mushroom_phy_p10@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)

    write.csv(ps_mushroom_phy_p0@tax_table, file.path(path_mushroom_phy, paste0(project, "taxa_mushroom_phy_p0.csv")))
    write.csv(ps_mushroom_phy_p1@tax_table, file.path(path_mushroom_phy, paste0(project, "taxa_mushroom_phy_p1.csv")))
    write.csv(ps_mushroom_phy_p5@tax_table, file.path(path_mushroom_phy, paste0(project, "taxa_mushroom_phy_p5.csv")))
    write.csv(ps_mushroom_phy_p10@tax_table, file.path(path_mushroom_phy, paste0(project, "taxa_mushroom_phy_p10.csv")))

    write.csv(ps_mushroom_phy_p0@otu_table, file.path(path_mushroom_phy, paste0(project, "otu_mushroom_phy_p0.csv")))
    write.csv(ps_mushroom_phy_p1@otu_table, file.path(path_mushroom_phy, paste0(project, "otu_mushroom_phy_p1.csv")))
    write.csv(ps_mushroom_phy_p5@otu_table, file.path(path_mushroom_phy, paste0(project, "otu_mushroom_phy_p5.csv")))
    write.csv(ps_mushroom_phy_p10@otu_table, file.path(path_mushroom_phy, paste0(project, "otu_mushroom_phy_p10.csv")))

    write.csv(otu_phy_p0_clr, file.path(path_mushroom_phy, paste0(project, "otu_mushroom_phy_p0_clr.csv")))
    write.csv(otu_phy_p1_clr, file.path(path_mushroom_phy, paste0(project, "otu_mushroom_phy_p1_clr.csv")))
    write.csv(otu_phy_p5_clr, file.path(path_mushroom_phy, paste0(project, "otu_mushroom_phy_p5_clr.csv")))
    write.csv(otu_phy_p10_clr, file.path(path_mushroom_phy, paste0(project, "otu_mushroom_phy_p10_clr.csv")))
  })
}


# --- 3.3.2 Plant-related ---
try <- try({ps_plant <- subset_taxa(ps_filter, Kingdom == "k__Viridiplantae")})

if (inherits(try, "try-error")) {
  message("Note: No plant reads were detected for analysis\n")
} else {
    try({
    path_plant <- file.path(path_vegetable, "plant")
    if(!dir.exists(path_plant)) dir.create(path_plant)

    # ASV level output
    path_plant_ASV <- file.path(path_plant, "ASV")
    if(!dir.exists(path_plant_ASV)) dir.create(path_plant_ASV)
    path_plant_ASV

    write.csv(readcount(ps_plant), file.path(path_plant_ASV, paste0(project, "readcount_plant.csv")))
    writeXStringSet(ps_plant@refseq, file = file.path(path_plant_ASV, paste0(project, "refseq_plant.fasta")), format = "fasta")
    write.csv(ps_plant@tax_table, file.path(path_plant_ASV, paste0(project, "taxa_plant.csv")))
    write.csv(ps_plant@otu_table, file.path(path_plant_ASV, paste0(project, "otu_plant.csv")))

    # Genus level output
    path_plant_gen <- file.path(path_plant, "genus")
    if(!dir.exists(path_plant_gen)) dir.create(path_plant_gen)
    path_plant_gen

    ps_plant_gen_p0 <- aggregate_taxa(ps_plant, level = "Genus")
    ps_plant_gen_p1 <- aggregate_rare(ps_plant, level = "Genus", detection = 0, prevalence = 1/100)
    ps_plant_gen_p5 <- aggregate_rare(ps_plant, level = "Genus", detection = 0, prevalence = 5/100)
    ps_plant_gen_p10 <- aggregate_rare(ps_plant, level = "Genus", detection = 0, prevalence = 10/100)

    otu_gen_p0_clr <- vegan::decostand(ps_plant_gen_p0@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
    otu_gen_p1_clr <- vegan::decostand(ps_plant_gen_p1@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
    otu_gen_p5_clr <- vegan::decostand(ps_plant_gen_p5@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
    otu_gen_p10_clr <- vegan::decostand(ps_plant_gen_p10@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)

    write.csv(ps_plant_gen_p0@tax_table, file.path(path_plant_gen, paste0(project, "taxa_plant_gen.csv")))
    write.csv(ps_plant_gen_p0@otu_table, file.path(path_plant_gen, paste0(project, "otu_plant_gen_p0.csv")))
    write.csv(ps_plant_gen_p1@otu_table, file.path(path_plant_gen, paste0(project, "otu_plant_gen_p1.csv")))
    write.csv(ps_plant_gen_p5@otu_table, file.path(path_plant_gen, paste0(project, "otu_plant_gen_p5.csv")))
    write.csv(ps_plant_gen_p10@otu_table, file.path(path_plant_gen, paste0(project, "otu_plant_gen_p10.csv")))

    write.csv(otu_gen_p0_clr, file.path(path_plant_gen, paste0(project, "otu_plant_gen_p0_clr.csv")))
    write.csv(otu_gen_p1_clr, file.path(path_plant_gen, paste0(project, "otu_plant_gen_p1_clr.csv")))
    write.csv(otu_gen_p5_clr, file.path(path_plant_gen, paste0(project, "otu_plant_gen_p5_clr.csv")))
    write.csv(otu_gen_p10_clr, file.path(path_plant_gen, paste0(project, "otu_plant_gen_p10_clr.csv")))

    # Phylum level output
    path_plant_phy <- file.path(path_plant, "phylum")
    if(!dir.exists(path_plant_phy)) dir.create(path_plant_phy)
    path_plant_phy

    ps_plant_phy_p0 <- aggregate_taxa(ps_plant, level = "Phylum")
    ps_plant_phy_p1 <- aggregate_rare(ps_plant, level = "Phylum", detection = 0, prevalence = 1/100)
    ps_plant_phy_p5 <- aggregate_rare(ps_plant, level = "Phylum", detection = 0, prevalence = 5/100)
    ps_plant_phy_p10 <- aggregate_rare(ps_plant, level = "Phylum", detection = 0, prevalence = 10/100)

    otu_phy_p0_clr <- vegan::decostand(ps_plant_phy_p0@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
    otu_phy_p1_clr <- vegan::decostand(ps_plant_phy_p1@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
    otu_phy_p5_clr <- vegan::decostand(ps_plant_phy_p5@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)
    otu_phy_p10_clr <- vegan::decostand(ps_plant_phy_p10@otu_table, MARGIN = 2, method = "clr", pseudocount = 1)

    write.csv(ps_plant_phy_p0@tax_table, file.path(path_plant_phy, paste0(project, "taxa_plant_phy.csv")))
    write.csv(ps_plant_phy_p0@otu_table, file.path(path_plant_phy, paste0(project, "otu_plant_phy_p0.csv")))
    write.csv(ps_plant_phy_p1@otu_table, file.path(path_plant_phy, paste0(project, "otu_plant_phy_p1.csv")))
    write.csv(ps_plant_phy_p5@otu_table, file.path(path_plant_phy, paste0(project, "otu_plant_phy_p5.csv")))
    write.csv(ps_plant_phy_p10@otu_table, file.path(path_plant_phy, paste0(project, "otu_plant_phy_p10.csv")))

    write.csv(otu_phy_p0_clr, file.path(path_plant_phy, paste0(project, "otu_plant_phy_p0_clr.csv")))
    write.csv(otu_phy_p1_clr, file.path(path_plant_phy, paste0(project, "otu_plant_phy_p1_clr.csv")))
    write.csv(otu_phy_p5_clr, file.path(path_plant_phy, paste0(project, "otu_plant_phy_p5_clr.csv")))
    write.csv(otu_phy_p10_clr, file.path(path_plant_phy, paste0(project, "otu_plant_phy_p10_clr.csv")))
  })
}

