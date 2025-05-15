# README

# Issues in ITS2 Sequencing Data from the Human Microbiome Project

## Introduction

While the Human Microbiome Project (HMP) significantly advanced our understanding of human microbiome, our recent revisiting of the fecal ITS2 sequencing data from HMP reveals two notable issues. First, technical variation in sequencing length or reads trimming introduces batch effects in downstream analyses. Second, using newly developed pipelines, we discovered a mismatch between the reported and actual amplification primers, which has significant implications for methodology standardization and cross-study comparisons. 

## Content

The HMP raw data was downloaded from PRJ356769. The annotation and analysis results are stored in the data folder. We divided the samples by read length into two groups: PRJ356769a (PE250) and PRJ356769b (PE300). The meta.csv file contains sample read types with more than 5,000 non-mushroom fungal reads after filtering. Results from Cutadapt, APHunter, and BLAST analyses are stored in their respective folders. All annotation and analysis code can be found in the code folder. 

## Contact

If you have any questions, please leave a comment in the Issues section or email the authors: Xinyu Wang ([wangxinyu30@westlake.edu.cn](mailto:wangxinyu30@westlake.edu.cn)) or Wenhao Zhou ([zhouwenhao@westlake.edu.cn](mailto:zhouwenhao@westlake.edu.cn)).