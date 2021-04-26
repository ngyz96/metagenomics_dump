module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_metawrap

SCRATCHDIR="/scratch/yng072/metaspades"
cd $SCRATCHDIR 

for f in `ls -1 *_1.fastq | sed 's/_1.fastq//'`
do
	metawrap assembly -1 ${f}_1.fastq -2 ${f}_2.fastq -m 400 -t 48 --metaspades -o ${f}_assembly 	
done
