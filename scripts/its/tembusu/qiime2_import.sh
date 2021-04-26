#!/bin/bash

module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_microbes

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path /home/bmb/Soil-Microbiomes/Tembusu-Microbiome/qiime2_its/its_manifest.tsv \
--output-path /home/bmb/Soil-Microbiomes/Tembusu-Microbiome/qiime2_its/its_paired_end.qza \
--input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
--i-data /home/bmb/Soil-Microbiomes/Tembusu-Microbiome/qiime2_its/its_paired_end.qza \
--o-visualization /home/bmb/Soil-Microbiomes/Tembusu-Microbiome/qiime2_its/its_paired_end.qzv
