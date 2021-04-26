#1/bin/bash

if [ $# == 0 ]; then
        echo "Usage: ./bwamem.sh -a ASSEMBLY -r READS [-o OUTDIR -t THREADS]"
	echo "Runs bwa mem on all paired end fastq files on their corresponding assembly"
        echo "		-a: Full file path for directory with assembly"
	echo "		-r: Full file path for directory with fastq files"
        echo "		-o: Full file path for output directory (Default=ASSEMBLY)"
        echo "		-t: Number of threads to use (Default=24)"
        exit 1
fi

while getopts "a:r:o:t:" opt
do
        case ${opt} in
        a) ASSEMBLY=${OPTARG} ;;
        r) READS=${OPTARG} ;;
	o) OUTDIR=${OPTARG} ;;
        t) THREADS=${OPTARG} ;;
        esac
done

: ${ASSEMBLY:?Missing -a} ${READS:?Missing -r} 
: ${THREADS:=24} ${OUTDIR:=$ASSEMBLY}

module load samtools/1.10
module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
source activate tensorflow_gpuenv
conda activate /usr/local_sbs/conda_aligners

if [ ! -d "$OUTDIR" ]; then
	mkdir -p $OUTDIR
fi

echo "*************************************************************************************************************"
echo "*****************************************    Starting bwa mem...    *****************************************"
echo "*************************************************************************************************************"
echo " " 
cd $READS

for FILE in `ls -1 *_1.fq.gz | sed 's/_1.fq.gz//'`
do
	SAMPLE=`echo $FILE | cut -d "_" -f1,2`
	if [ ! -f "${ASSEMBLY}/${SAMPLE}.contigs.fa" ]; then
		echo "Assembly missing for $SAMPLE!"
		continue
	fi
	if [ -f "${FILE}_2.fq.gz" ]; then
		echo "Running bwamem for $SAMPLE in paired end mode"
		bwa mem ${ASSEMBLY}/${SAMPLE}.contigs.fa ${FILE}_1.fq.gz ${FILE}_2.fq.gz -t ${THREADS} | samtools sort -@${THREADS} -o ${OUTDIR}/${SAMPLE}.bam
	else
		echo "Read 2 missing! Running bwamem in $SAMPLE in single end mode"
		bwa mem ${ASSEMBLY}/${SAMPLE}.contigs.fa ${FILE}_1.fq.gz -t ${THREADS} | samtools sort -@${THREADS} -o ${OUTDIR}/${SAMPLE}.bam
        fi
	echo " "
	echo "assembly done for $SAMPLE"
done
