#!/usr/bin/env Rscript
options(echo = T)  

library(dada2)
cutadapt <- "/home/wangxinyu/miniconda3/envs/fq/bin/cutadapt"

# define the input folder
# HMP
# path_in <- "/data/wangxinyu/HMP_ITS/data/PRJNA356769_ITS2/0_raw_merged" 

# GNHS
# path_in <- "/data/wangxinyu/HMP_ITS/data/GNHS/raw_merged" 

# MK-SpikeSeq
path_in <- "/data/wangxinyu/HMP_ITS/data/MK-SpikeSeq/raw_merged"


# define the pattern of fastq file
pattern_F = ".1.fastq.gz" 
pattern_R = ".2.fastq.gz"

fqF_raw <- sort(list.files(path_in, pattern = pattern_F, full.names = TRUE))
fqR_raw <- sort(list.files(path_in, pattern = pattern_R, full.names = TRUE))
head(fqF_raw)
head(fqR_raw)

path_filterA <- "/data/wangxinyu/HMP_ITS/cutadapt/linker_ITS3_ITS4_MK-SpikeSeq_default"
if(!dir.exists(path_filterA)) {dir.create(path_filterA)}

# deine the output folder
fqF_filterA <- file.path(path_filterA, basename(fqF_raw))
fqR_filterA <- file.path(path_filterA, basename(fqR_raw))

# construct the cutadapt flag
# fwd <- "GCATCGATGAAGAACGCAGC" #ITS3
# rev <- "TCCTCCGCTTATTGATATGC" #ITS4

fwd <- "CTYGGTCATTTAGAGGAAGTAA" #ITS1F (modified)
rev <- "GCTGCGTTCTTCATCGATGC" #ITS2

fwd_rc <- dada2:::rc(fwd)
rev_rc <- dada2:::rc(rev)

# default
F_flag <- paste0("-a ", '"', fwd, ";min_overlap = ", round((0.9 * nchar(fwd)), 0), ";required",'...', rev_rc, ";min_overlap = ", round((0.9 * nchar(rev_rc)), 0), ";optional", '"') 
R_flag <- paste0("-A ", '"', rev, ";min_overlap = ", round((0.9 * nchar(rev)), 0), ";required",'...', fwd_rc, ";min_overlap = ", round((0.9 * nchar(fwd_rc)), 0), ";optional", '"') 

# relax 
# F_flag <- paste0("-a ", '"', fwd, ";min_overlap = ", round((0.75 * nchar(fwd)), 0), ";e=0.2;required",'...', rev_rc, ";min_overlap = ", round((0.75 * nchar(rev_rc)), 0), ";optional", '"') 
# R_flag <- paste0("-A ", '"', rev, ";min_overlap = ", round((0.75 * nchar(rev)), 0), ";e=0.2;required",'...', fwd_rc, ";min_overlap = ", round((0.75 * nchar(fwd_rc)), 0), ";optional", '"') 

F_flag 
R_flag

# run cutadapt
for(i in seq_along(fqF_raw)) {
  system2(cutadapt,
          args = c(F_flag, R_flag,
                   "-o", fqF_filterA[i], "-p", fqR_filterA[i],
                   fqF_raw[i], fqR_raw[i],
                   "--discard-untrimmed",
                   "--minimum-length 1",  # drop empty reads
                   "-j", 100))
}

