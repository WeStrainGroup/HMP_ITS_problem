
options(echo = T)  

library(dada2)
cutadapt <- "/home/wangxinyu/miniconda3/envs/fq/bin/cutadapt"

# define the input folder
path_in <- "/data/wangxinyu/HMP_ITS/data/PRJNA356769_ITS2/0_raw"

# define the pattern of fastq file
pattern_F = ".1.fastq.gz"
pattern_R = ".2.fastq.gz"

fqF_raw <- sort(list.files(path_in, pattern = pattern_F, full.names = TRUE))
fqR_raw <- sort(list.files(path_in, pattern = pattern_R, full.names = TRUE))
head(fqF_raw)
head(fqR_raw)

# output
path_filterA <- "/data/wangxinyu/HMP_ITS/cutadapt/ITS86F_LR0B_HMP_default"
if(!dir.exists(path_filterA)) {dir.create(path_filterA)}

fqF_filterA <- file.path(path_filterA, basename(fqF_raw))
fqR_filterA <- file.path(path_filterA, basename(fqR_raw))

# set the cutadatp flag
fwd <- "GTGAATCATCGAATCTTTGAA" #ITS86F
rev <- "GGTAGTCCTACCTGATTTG" # LR0B

# default
F_flag <- paste0("-g ", '"', fwd, ";min_overlap = ", round((0.9 * nchar(fwd)), 0),'"') 
R_flag <- paste0("-g ", '"', rev, ";min_overlap = ", round((0.9 * nchar(rev)), 0),'"') 

# relex
# F_flag <- paste0("-g ", '"', fwd, ";min_overlap = ", round((0.75 * nchar(fwd)), 0),";e=0.2",'"') 
# R_flag <- paste0("-g ", '"', rev, ";min_overlap = ", round((0.75 * nchar(rev)), 0),";e=0.2",'"') 

F_flag
R_flag

# run R1
for(i in seq_along(fqF_raw)) {
  system2(cutadapt,
          args = c(F_flag, 
                  "-o", fqF_filterA[i],
                   fqF_raw[i],
                   "--discard-untrimmed",
                   "--minimum-length 1",  # drop empty reads
                   "-j", 100))
}


# run R2
for(i in seq_along(fqR_raw)) {
  system2(cutadapt,
          args = c(R_flag, 
                  "-o", fqR_filterA[i],
                   fqR_raw[i],
                   "--discard-untrimmed",
                   "--max-n 0",
                   "--minimum-length 1",  # drop empty reads
                   "-j", 100))
}
