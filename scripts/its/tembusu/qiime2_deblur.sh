#!/bin/bash

module load anaconda2020/python3
eval "$(conda shell.bash hook)"
conda activate /usr/local_sbs/conda_microbes
PWD="/home/bmb/Soil-Microbiomes/Tembusu-Microbiome/qiime2_its"
qimme quality filter q-score \
--i-demux $PWD/its_paired_end.qza \
--o-filtered-sequences its_paired_end_filtered.qza \
--o-filter-stats its_paired_end_filter_stats.qza


