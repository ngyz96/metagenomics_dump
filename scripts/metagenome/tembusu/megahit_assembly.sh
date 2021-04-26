module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_microbes

SCRATCHDIR="/scratch/yng072/metaspades"
cd $SCRATCHDIR

megahit -1 B_CD_1.fastq -2 B_CD_2.fastq --presets meta-large -o ./B_CD_megahit \
-m 400000000000 -t 24 --min-contig-len 1000
