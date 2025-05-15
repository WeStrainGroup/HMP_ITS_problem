library(Biostrings)
library(dplyr)

# --- APHunter-related ---
# GNHS
# R1
blast_R1 <- read.csv("/data/wangxinyu/HMP_ITS/aphunter/GNHS_R1/blast/per_sample_final_blast_results.csv")

# check the num of all sample
length(unique(blast_R1$query_id)) # 1633

# check the num of N/A (sample without blast hit)
sum(blast_R1$subject_id == "N/A") # 0

# check the primer
table(blast_R1$subject_id)

# check the ITS1F
blast_R1_filter <- blast_R1 %>% filter(subject_id == "ITS3_2_fwd")
table(blast_R1_filter$mismatches)

# R2 
blast_R2 <- read.csv("/data/wangxinyu/HMP_ITS/aphunter/GNHS_R2/blast/per_sample_final_blast_results.csv")

# check the num of all sample
length(unique(blast_R2$query_id)) # 1633

# check the num of N/A (sample without blast hit)
sum(blast_R2$subject_id == "N/A") # 0

# check the primer
table(blast_R2$subject_id) 

# check the ITS2
blast_R2_filter <- blast_R2 %>% filter(subject_id == "ITS4_2_rev")
table(blast_R2_filter$mismatches)


# MK-SpikeSeq
# R1
blast_R1 <- read.csv("/data/wangxinyu/HMP_ITS/aphunter/MK-SpikeSeq_R1/blast/per_sample_final_blast_results.csv")

# check the num of all sample
length(unique(blast_R1$query_id)) # 308

# check the num of N/A (sample without blast hit)
sum(blast_R1$subject_id == "N/A") # 0

# check the primer
table(blast_R1$subject_id)

# check the ITS1F
blast_R1_filter <- blast_R1 %>% filter(subject_id == "ITS1F_1_fwd")
table(blast_R1_filter$mismatches) # all 0

# R2 
blast_R2 <- read.csv("/data/wangxinyu/HMP_ITS/aphunter/MK-SpikeSeq_R2/blast/per_sample_final_blast_results.csv")

# check the num of all sample
length(unique(blast_R2$query_id)) # 308

# check the num of N/A (sample without blast hit)
sum(blast_R2$subject_id == "N/A") # 57

# check the primer
table(blast_R2$subject_id) 

# check the ITS2
blast_R2_filter <- blast_R2 %>% filter(subject_id == "ITS2_1_rev")
table(blast_R2_filter$mismatches) # all 0


# HMP
# Check the consensus sequences
consensus_R1 <- readDNAStringSet("/data/wangxinyu/HMP_ITS/aphunter/HMP_R1/blast/sample_consensus_for_blast.fasta", format = "fasta")
consensus_R2 <- readDNAStringSet("/data/wangxinyu/HMP_ITS/aphunter/HMP_R2/blast/sample_consensus_for_blast.fasta", format = "fasta")

consensus_R1 # 310
consensus_R2 # 308

# Filter sequences based on length
table(width(consensus_R1))
table(width(consensus_R2))

consensus_R1_filtered <- consensus_R1[width(consensus_R1) >= 25]
consensus_R2_filtered <- consensus_R2[width(consensus_R2) >= 25]

consensus_R1_filtered # 296
consensus_R2_filtered # 280

# Analyze BLAST results
blast_R1 <- read.csv("/data/wangxinyu/HMP_ITS/aphunter/HMP_R1/blast/per_sample_final_blast_results.csv")

# check the num of all sample
length(unique(blast_R1$query_id)) # 310

# check the num of N/A (sample without blast hit)
sum(blast_R1$subject_id == "N/A") # 72

# number of sample with significant hit
length(unique(blast_R1$query_id)) - sum(blast_R1$subject_id == "N/A") # 238

# check the primer
table(blast_R1$subject_id)

# check the ITS86F
blast_R1_filter <- blast_R1 %>% filter(subject_id == "ITS86F_2_fwd")
table(blast_R1_filter$mismatches)

# R2
blast_R2 <- read.csv("/data/wangxinyu/HMP_ITS/aphunter/HMP_R2/blast/per_sample_final_blast_results.csv")

# check the num of all sample
length(unique(blast_R2$query_id)) # 308

# check the num of N/A (sample without blast hit)
sum(blast_R2$subject_id == "N/A") # 68

# number of sample with significant hit
length(unique(blast_R2$query_id)) - sum(blast_R2$subject_id == "N/A") # 240

# check the primer
table(blast_R2$subject_id) 

# check the LR0B
blast_R2_filter <- blast_R2 %>% filter(subject_id == "LR0B_2_rev")
table(blast_R2_filter$mismatches)


# --- ASV-related ---
# Calculate the proportion of raw reads retained in the ASV table
fq_stats_raw <- read.delim("/data/wangxinyu/HMP_ITS/data/PRJNA356769_ITS2/1_dada2/check/PRJNA356769_ITS2_fq_stats_raw.tsv",
                           sep = "\t",
                           header = T)
sum_num_reads_raw <- sum(fq_stats_raw$num_seqs)
sum_num_reads_raw / 2 # 12383662

seqtab_raw <- read.csv("/data/wangxinyu/HMP_ITS/data/PRJNA356769_ITS2/1_dada2/result/PRJNA356769_ITS2_seqtab_raw.csv", row.names = 1)
sum_num_reads_asv <- sum(seqtab_raw)
sum_num_reads_asv # 10658025

round(sum_num_reads_asv / (sum_num_reads_raw / 2), 3) # 0.861

# Calculate the proportion of reads in ASVs containing ITS2 region
asv_without_its <- readDNAStringSet("/data/wangxinyu/HMP_ITS/data/PRJNA356769_ITS2/1_dada2/itsx/PRJNA356769_ITS2_asv_itsx_no_detections.fasta",
                                    format = "fasta")
asv_without_its

# Extract sequences
seq <- unname(as.character(asv_without_its))
seq

# Remove ASVs without detected ITS2 region
seqtab_filter <- seqtab_raw %>% select(-all_of(seq))
dim(seqtab_raw)
dim(seqtab_filter)

sum_num_reads_asv_its <- sum(seqtab_filter)
sum_num_reads_asv_its # 10187352

round(sum_num_reads_asv_its / (sum_num_reads_raw / 2), 3) # 0.823

# Analyze primer hits in ASVs containing ITS2 region
asv_blast <- read.delim("/data/wangxinyu/HMP_ITS/blast/Blast_HMP_ASV.tsv", sep = "\t", header = T)
head(asv_blast)
tail(asv_blast)

asv_blast_filter <- asv_blast %>% filter(!qseqid %in% names(asv_without_its))
dim(asv_blast)
dim(asv_blast_filter)
length(unique(asv_blast_filter$qseqid)) # number of ASVs with hit = 1381

table(asv_blast_filter$sseqid)

# Analyze forward primer (ITS86F) hits 
asv_blast_filter2 <- asv_blast_filter %>% filter(sseqid == "ITS86F_2_fwd")
table(asv_blast_filter2$mismatch)

length(asv_blast_filter2$qseqid) # 1198
length(asv_blast_filter2$qseqid) / length(unique(asv_blast_filter$qseqid)) # 0.867

# 算一下这些ASV所对应的seqtale_raw 中 reads数量
asv_filter2 <- asv_all[names(asv_all) %in% asv_blast_filter2$qseqid]
asv_filter2

seq <- unname(as.character(asv_filter2))

seqtab_filter2 <- seqtab_raw %>% select(all_of(seq))
dim(seqtab_raw)
dim(seqtab_filter2)
sum(seqtab_filter2) # 9444830

round(sum(seqtab_filter2) / sum(seqtab_raw), 3) # 0.886

# Analyze reverse primer (LR0B) hits
asv_blast_filter3 <- asv_blast_filter %>% filter(sseqid == "LR0B_2_rev")
table(asv_blast_filter3$mismatch)

length(unique(asv_blast_filter3$qseqid)) # 451
length(unique(asv_blast_filter3$qseqid)) / length(unique(asv_blast_filter$qseqid)) # 0.326

# 算一下这些ASV所对应的seqtale_raw 中 reads数量
asv_filter3 <- asv_all[names(asv_all) %in% asv_blast_filter3$qseqid]
asv_filter3

seq <- unname(as.character(asv_filter3))

seqtab_filter3 <- seqtab_raw %>% select(all_of(seq))
dim(seqtab_raw)
dim(seqtab_filter3)
sum(seqtab_filter3) # 6972935

round(sum(seqtab_filter3) / sum(seqtab_raw), 3) # 0.654

# Analyze reverse primer (ITS4-Seb) hits
asv_blast_filter4 <- asv_blast_filter %>% filter(sseqid == "ITS4-Seb_2_rev")
table(asv_blast_filter4$mismatch)

length(unique(asv_blast_filter4$qseqid)) # 455
length(unique(asv_blast_filter4$qseqid)) / length(unique(asv_blast_filter$qseqid)) # 0.329

# 算一下这些ASV所对应的seqtale_raw 中 reads数量
asv_filter4 <- asv_all[names(asv_all) %in% asv_blast_filter4$qseqid]
asv_filter4

seq <- unname(as.character(asv_filter4))

seqtab_filter4 <- seqtab_raw %>% select(all_of(seq))
dim(seqtab_raw)
dim(seqtab_filter4)
sum(seqtab_filter4) # 6475880

round(sum(seqtab_filter4) / sum(seqtab_raw), 3) # 0.608
