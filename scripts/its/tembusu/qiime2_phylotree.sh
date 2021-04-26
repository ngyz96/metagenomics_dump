#!/bin/bash
module load anaconda2020/python3 
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_microbes
WD="/home/bmb/Soil-Microbiomes/Tembusu-Microbiome/qiime2_its/dada2"

qiime phylogeny align-to-tree-mafft-iqtree \
--i-sequences $WD/dada2_rep_seq_trimmed.qza \
--o-alignment $WD/dada2_alignments_trimmed.qza \
--o-masked-alignment $WD/dada2_masked_alignments_trimmed.qza \
--o-tree $WD/dada2_tree_trimmed.qza \
--o-rooted-tree $WD/dada2_rooted_tree_trimmed.qza \
--p-n-threads 4
