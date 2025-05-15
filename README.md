# README

# Issues in ITS2 Sequencing Data from the Human Microbiome Project

## Introduction

While the Human Microbiome Project (HMP) significantly advanced our understanding of the human microbiome, our revisiting of its fecal ITS2 sequencing data reveals two significant issues. First, technical variations in sequencing or data processing introduce batch effects in downstream analyses. Second, we discovered a mismatch between the reported and actual amplification primers, which requires careful consideration when interpreting, reusing, or referencing HMP ITS2 data and methodology.

## Content

The HMP raw data was downloaded from PRJ356769. The annotation and analysis results are stored in the data folder. We divided the samples by read length into two groups: PRJ356769a (PE250) and PRJ356769b (PE300). The meta.csv file contains sample read types with more than 5,000 non-mushroom fungal reads after filtering. Results from Cutadapt, APHunter, and BLAST analyses are stored in their respective folders. All annotation and analysis code can be found in the code folder. 

## Contact

If you have any questions, please leave a comment in the Issues section or email the authors: Xinyu Wang ([wangxinyu30@westlake.edu.cn](mailto:wangxinyu30@westlake.edu.cn)) or Wenhao Zhou ([zhouwenhao@westlake.edu.cn](mailto:zhouwenhao@westlake.edu.cn)).