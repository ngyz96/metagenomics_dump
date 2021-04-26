module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_microbes
WD="/home/bmb/Soil-Microbiomes/Tembusu-Microbiome/qiime2_its"

qiime tools export \
--input-path $WD/dada2/dada2_table_trimmed.qza \
--output-path $WD/dada2/export_data

qiime tools export \
--input-path $WD/dada2/dada2_rep_seq_trimmed.qza \
--output-path $WD/dada2/export_data

qiime tools export \
--input-path $WD/dada2/dada2_unite_taxonomy_nb_trimmed.qza \
--output-path $WD/dada2/export_data
