#!/bin/bash

# 设置变量
QUERY="/data/wangxinyu/HMP_ITS/data/PRJNA356769_ITS2/1_dada2/itsx/PRJNA356769_ITS2_asv_raw.fasta"
OUTPUT="/data/wangxinyu/HMP_ITS/blast/Blast_HMP_ASV.tsv"
DB="/data/wangxinyu/HMP_ITS/blast/its_primer_db/its_primer_db"

# 创建带有表头的输出文件
echo -e "qseqid\tsseqid\tpident\tlength\tqcovs\tmismatch\tgapopen\tqlen\tqstart\tqend\tslen\tsstart\tsend\tsstrand\tevalue\tbitscore" > "$OUTPUT"

# 运行blastn
blastn -query "$QUERY" \
       -db "$DB" \
       -out "$OUTPUT.tmp" \
       -outfmt '6 qseqid sseqid pident length qcovs mismatch gapopen qlen qstart qend slen sstart send sstrand evalue bitscore' \
       -task blastn-short \
       -evalue 0.01 \
       -num_threads 25

# 将blast结果追加到输出文件
cat "$OUTPUT.tmp" >> "$OUTPUT"

# 处理未匹配的序列
comm -23 \
  <(grep '^>' "$QUERY" | sed 's/^>//' | sort) \
  <(cut -f1 "$OUTPUT.tmp" | sort | uniq) \
| awk -v OFS='\t' '{print $1, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"}' \
>> "$OUTPUT"

# 删除临时文件
rm "$OUTPUT.tmp"