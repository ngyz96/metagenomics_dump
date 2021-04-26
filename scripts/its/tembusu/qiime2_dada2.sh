#!/bin/bash

module load anaconda2020/python3
eval "$(conda shell.bash hook)"
conda activate /usr/local_sbs/conda_microbes
WD="/home/bmb/Soil-Microbiomes/Tembusu-Microbiome/qiime2_its"
OUTDIR="/home/bmb/Soil-Microbiomes/Tembusu-Microbiome/qiime2_its/dada2"
qiime dada2 denoise-paired \
--i-demultiplexed-seqs $WD/its_paired_end_trimmed_unmerged.qza \
--p-trunc-len-f 0 \
--p-trunc-len-r 0 \
--p-n-threads 12 \
--o-table $OUTDIR/dada2_table_trimmed.qza \
--o-representative-sequences $OUTDIR/dada2_rep_seq_trimmed.qza \
--o-denoising-stats $OUTDIR/dada2_denoising_stats_trimmed.qza

qiime metadata tabulate \
--m-input-file $OUTDIR/dada2_denoising_stats_trimmed.qza \
--o-visualization $OUTDIR/dada2_denoising_stats_trimmed.qzv

qiime feature-table summarize \
--i-table $OUTDIR/dada2_table_trimmed.qza \
--o-visualization $OUTDIR/dada2_table_trimmed.qzv \
--m-sample-metadata-file $WD/tembusu_its_metadata.tsv \

qiime feature-table tabulate-seqs \
--i-data $OUTDIR/dada2_rep_seq_trimmed.qza \
--o-visualization $OUTDIR/dada2_rep_seq_trimmed.qzv
