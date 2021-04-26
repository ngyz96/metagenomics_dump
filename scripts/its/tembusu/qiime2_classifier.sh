#!/bin/bash
module load anaconda2020/python3 
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_microbes
unitepath="/home/bmb/Soil-Microbiomes/databases/sh_qiime_release_04.02.2020"
WD="/home/bmb/Soil-Microbiomes/Tembusu-Microbiome/qiime2_its"

qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path $unitepath/sh_refs_qiime_ver8_dynamic_04.02.2020.fasta \
--output-path $WD/unite_trimmed.qza

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path $unitepath/sh_taxonomy_qiime_ver8_dynamic_04.02.2020.txt \
--output-path $WD/unite_trimmed_taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads $WD/unite_trimmed.qza \
--i-reference-taxonomy $WD/unite_trimmed_taxonomy.qza \
--o-classifier $WD/unite_naive_bayes_trimmed_classifier.qza

qiime feature-classifier classify-sklearn \
--i-reads $WD/dada2/dada2_rep_seq_trimmed.qza \
--i-classifier $WD/unite_naive_bayes_trimmed_classifier.qza \
--o-classification $WD/dada2/dada2_unite_taxonomy_nb_trimmed.qza

qiime metadata tabulate \
--m-input-file $WD/dada2/dada2_unite_taxonomy_nb_trimmed.qza \
--o-visualization $WD/dada2/dada2_unite_taxonomy_nb_trimmed.qzv
