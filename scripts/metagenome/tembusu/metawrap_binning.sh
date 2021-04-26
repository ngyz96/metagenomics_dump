#!/bin/bash

module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_metawrap

SCRATCH="/scratch/yng072/metaspades"
ASSEMBLY="/scratch/yng072/metaspades/B_CD_megahit/final.contigs.fa"
cd $SCRATCH

metawrap binning -o INITIAL_BINNING -t 60 -m 100 -a ${ASSEMBLY} --metabat2 --maxbin2 --concoct B_CD_1.fastq B_CD_2.fastq 
