#!/bin/bash

module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_metawrap

SCRATCH="/scratch/yng072/metaspades"

cd $SCRATCH

metawrap bin_refinement -o BIN_REFINEMENT -t 60 -m 100 \
-A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins/ -C INITIAL_BINNING/concoct_bins/ -c 70 -x 10 
