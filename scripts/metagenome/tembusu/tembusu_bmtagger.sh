module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_metawrap

bmtool -d /home/bmb/Soil-Microbiomes/bmtagger_index/F_fragrans_modified_CanuHiC.faa -o /home/bmb/Soil-Microbiomes/bmtagger_index/F_fragrans_modified_CanuHiC.bitmask
srprism mkindex -i /home/bmb/Soil-Microbiomes/bmtagger_index/F_fragrans_modified_CanuHiC.faa -o /home/bmb/Soil-Microbiomes/bmtagger_index/F_fragrans_modified_CanuHiC.faa.srprism -M 12000
