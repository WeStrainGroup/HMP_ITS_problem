# Randomly select 250 paired-end files from GNHS dataset and merge them into a single file
cd /data/wangxinyu/HMP_ITS/data/GNHS/raw_merged
find /data/wangxinyu/ITS/Project/GNHS/dada/raw -name "*1.fq.gz" | shuf | head -n 250 | tee selected_R1.txt | sed 's/.1.fq.gz/.2.fq.gz/' > selected_R2.txt && \
xargs -a selected_R1.txt gzip -dc | gzip > GNHS_250_merged.1.fastq.gz && \
xargs -a selected_R2.txt gzip -dc | gzip > GNHS_250_merged.2.fastq.gz

# Merge all PE300 sequecning sample from PRJEB36435
cd /data/wangxinyu/HMP_ITS/data/MK-SpikeSeq/raw_merged
ls /data/huangkailang/project/fungi_ITS/rawdata/PE/PRJEB36435b_ITS1/*.1.fastq.gz | sort | tee selected_R1.txt | xargs zcat | gzip > MK-SpikeSeq_merged.1.fastq.gz
ls /data/huangkailang/project/fungi_ITS/rawdata/PE/PRJEB36435b_ITS1/*.2.fastq.gz | sort | tee selected_R2.txt | xargs zcat | gzip > MK-SpikeSeq_merged.2.fastq.gz


# Run cutadapt with different parameters for various datasets
# Process HMP dataset with default parameters using linked ITS3-ITS4 primers
nohup Rscript /data/wangxinyu/HMP_ITS/code/cutadapt_linked.R > /data/wangxinyu/HMP_ITS/cutadapt/linker_ITS3_ITS4_HMP_default.log 2>&1 &
# Process GNHS dataset with default parameters using linked ITS3-ITS4 primers
nohup Rscript /data/wangxinyu/HMP_ITS/code/cutadapt_linked.R > /data/wangxinyu/HMP_ITS/cutadapt/linker_ITS3_ITS4_GNHS_default.log 2>&1 &
# Process MK-SpikeSeq dataset with default parameters using linked ITS1F-ITS2 primers
nohup Rscript /data/wangxinyu/HMP_ITS/code/cutadapt_linked.R > /data/wangxinyu/HMP_ITS/cutadapt/linker_ITS1F_ITS2_MK-SpikeSeq_default.log 2>&1 &
nohup Rscript /data/wangxinyu/HMP_ITS/code/cutadapt_linked.R > /data/wangxinyu/HMP_ITS/cutadapt/linker_ITS1F_ITS2_MK-SpikeSeq_spilted_default.log 2>&1 &


# Process HMP dataset with relaxed parameters using ITS3-ITS4 primers
nohup Rscript /data/wangxinyu/HMP_ITS/code/cutadapt_linked.R > /data/wangxinyu/HMP_ITS/cutadapt/linker_ITS3_ITS4_HMP_relax.log 2>&1 &


# Process HMP dataset with default parameters using ITS86F and LR0B primers
nohup Rscript /data/wangxinyu/HMP_ITS/code/cutadapt_single.R > /data/wangxinyu/HMP_ITS/cutadapt/ITS86F_LR0B_HMP_default.log 2>&1 &
# Process HMP dataset with relex parameters using ITS86F and LR0B primers
nohup Rscript /data/wangxinyu/HMP_ITS/code/cutadapt_single.R > /data/wangxinyu/HMP_ITS/cutadapt/ITS86F_LR0B_HMP_relex.log 2>&1 &

# Generate sequence statistics using SeqKit
seqkit stats /data/wangxinyu/HMP_ITS/cutadapt/ITS86F_LR0B_HMP_default/*.1.fastq.gz -a -b -e -T -j 100 > /data/wangxinyu/HMP_ITS/cutadapt/fq_stats_F_cutITS86F_default.tsv
seqkit stats /data/wangxinyu/HMP_ITS/cutadapt/ITS86F_LR0B_HMP_default/*.2.fastq.gz -a -b -e -T -j 100 > /data/wangxinyu/HMP_ITS/cutadapt/fq_stats_R_cutLR0B_default.tsv

seqkit stats /data/wangxinyu/HMP_ITS/cutadapt/ITS86F_LR0B_HMP_relax/*.1.fastq.gz -a -b -e -T -j 100 > /data/wangxinyu/HMP_ITS/cutadapt/fq_stats_F_cutITS86F_relax.tsv
seqkit stats /data/wangxinyu/HMP_ITS/cutadapt/ITS86F_LR0B_HMP_relax/*.2.fastq.gz -a -b -e -T -j 100 > /data/wangxinyu/HMP_ITS/cutadapt/fq_stats_R_cutLR0B_relax.tsv


# Run APHunter for primer analysis
# GNHS dataset
aphunter -i /data/wangxinyu/ITS/Project/GNHS/dada/raw -o /data/wangxinyu/HMP_ITS/aphunter/GNHS_R1 -s .1.fq.gz -t 100
aphunter -i /data/wangxinyu/ITS/Project/GNHS/dada/raw -o /data/wangxinyu/HMP_ITS/aphunter/GNHS_R2 -s .2.fq.gz -t 100

# MK-SpikeSeq dataset
aphunter -i /data/huangkailang/project/fungi_ITS/rawdata/PE/PRJEB36435b_ITS1 -o /data/wangxinyu/HMP_ITS/aphunter/MK-SpikeSeq_R1 -s .1.fastq.gz -t 100
aphunter -i /data/huangkailang/project/fungi_ITS/rawdata/PE/PRJEB36435b_ITS1 -o /data/wangxinyu/HMP_ITS/aphunter/MK-SpikeSeq_R2 -s .2.fastq.gz -t 100

# HMP dataset
aphunter -i /data/wangxinyu/HMP_ITS/data/PRJNA356769_ITS2/0_raw -o /data/wangxinyu/HMP_ITS/aphunter/HMP_R1 -s .1.fastq.gz -t 100
aphunter -i /data/wangxinyu/HMP_ITS/data/PRJNA356769_ITS2/0_raw -o /data/wangxinyu/HMP_ITS/aphunter/HMP_R2 -s .2.fastq.gz -t 100

# Blast
makeblastdb -in /data/wangxinyu/HMP_ITS/blast/its_primer_db/its_primer_db.fasta -dbtype nucl -out /data/wangxinyu/HMP_ITS/blast/its_primer_db/its_primer_db

bash /data/wangxinyu/HMP_ITS/code/blast.sh