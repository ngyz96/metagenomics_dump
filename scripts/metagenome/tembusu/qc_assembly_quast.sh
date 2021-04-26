module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"

conda activate tensorflow_gpuenv
conda activate /usr/local_sbs/conda_assemblers

SCRATCHDIR="/scratch/yng072/metaspades/B_CD_megahit"
cd $SCRATCHDIR 

quast.py -m 500 -t 10 -o ./ final.contigs.fa 
