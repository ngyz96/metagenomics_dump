#!/bin/bash
module load anaconda2020/python3 
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_microbes

qiime itsxpress trim-pair-output-unmerged \
--i-per-sample-sequences /Soil-Microbiomes/Tembusu-Microbiome/qiime2_its/its_paired_end.qza \
--p-region ITS2 \
--p-taxa F \
--p-cluster-id 1.0 \
--p-threads 12 \
--o-trimmed /Soil-Microbiomes/Tembusu-Microbiome/qiime2_its/its_paired_end_trimmed_unmerged.qza
