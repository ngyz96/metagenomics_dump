#!/bin/bash
        
module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate tensorflow_gpuenv
conda activate /usr/local_sbs/conda_assemblers
cd /scratch/yng072/dip/megahit/contigs/
quast.py -m 500 -t 40 -o /scratch/yng072/dip/megahit/quast/ `cat /scratch/yng072/dip/megahit/contigs/contig_list`
