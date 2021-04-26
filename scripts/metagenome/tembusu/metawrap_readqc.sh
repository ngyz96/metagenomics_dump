module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_metawrap

INDIR="/home/bmb/Soil-Microbiomes/Tembusu-Microbiome/metagenome/X401SC20122148-Z01-F001/fastp_data"
OUTDIR="/home/bmb/Soil-Microbiomes/Tembusu-Microbiome/metagenome/X401SC20122148-Z01-F001/post_qc"
for f in `ls -1 $INDIR/*_1.fq | sed 's/_1.fq//'`
do
	filename=`echo $f | cut -d "/" -f9`
	samplename=`echo $filename | cut -d "_" -f1,2,3,4`
	metawrap read_qc -1 ${f}_1.fq -2 ${f}_2.fq -o ${OUTDIR}/${samplename} -t 48 -x F_fragrans_modified_CanuHiC \
	--skip-trimming --skip-pre-qc-report 	
done
