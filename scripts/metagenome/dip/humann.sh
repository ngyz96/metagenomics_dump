#1/bin/bash

if [ $# == 0 ]; then
        echo "Usage: ./humann.sh -i INDIR [-o OUTDIR -t THREADS]"
	echo "Runs HUMAnN3 on all QC'ed fastq files in a directory (cat paired end fastq files prior to use)"
	echo "This script runs with the default ChocoPhlan pangenome database and Uniref50 database for use with soil samples. Please change the script and humann-config as required."
        echo "		-i: Full file path for directory with fastq files"
        echo "		-o: Full file path for output directory (Default=INDIR)"
        echo "		-t: Number of threads to use (Default=24)"
        exit 1
fi

while getopts "i:o:t:" opt
do
        case ${opt} in
        i) INDIR=${OPTARG} ;;
        o) OUTDIR=${OPTARG} ;;
        t) THREADS=${OPTARG} ;;
        esac
done

: ${INDIR:?Missing -i}
: ${THREADS:=24} ${OUTDIR:=$INDIR}


module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/biobakery

if [ ! -d "$OUTDIR" ]; then
	mkdir -p $OUTDIR
fi

echo "*************************************************************************************************************"
echo "*****************************************    Starting HUMAnN3...    *****************************************"
echo "*************************************************************************************************************"
echo " " 
cd $INDIR

for FILE in `ls -1 *.fq.gz | sed 's/.fq.gz//'`
do
        SAMPLE=`echo $FILE | cut -d "_" -f1,2`
        echo "Running HUMAnN3 for $SAMPLE"
	humann --input ${FILE}.fq.gz --output ${OUTDIR}/${SAMPLE}_humann_output --threads $THREADS
	echo " "
	echo "HUMAnN3 done for $SAMPLE"
done
