module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_metawrap

megahit="/scratch/yng072/metaspades/B_CD_megahit/final.contigs.fa"
metaspades="/scratch/yng072/metaspades/B_CD_assembly/final_assembly.fasta"
WDIR="/scratch/yng072/metaspades" 
cd $WDIR

metawrap kraken2 -o kraken2_out -t 60 $megahit $metaspades B_CD_1.fastq B_CD_2.fastq
