#!/bin/bash

if [ $# == 0 ]; then
        echo "Usage: ./subsample.sh -i INDIR [-n NUM_READS]"
	echo "Runs seqtk on fastq files in a directory"
        echo "		-i: Full file path for directory with fastq files"
        echo "		-n: number of reads after subsampling(Default=1000000)"
        exit 1
fi

while getopts "i:n:" opt
do
        case ${opt} in
        i) INDIR=${OPTARG} ;;
        n) NUM=${OPTARG} ;;
        esac
done

: ${INDIR:?Missing -i}
: ${NUM:=1000000}


module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_genepredictors

echo "**************************************************************************************************************"
echo "*****************************************    Starting sampling...    *****************************************"
echo "**************************************************************************************************************"
echo " " 
cd $INDIR

for FILE in `ls -1 *_1.fq.gz | sed 's/_1.fq.gz//'`
do
        SAMPLE=`echo $FILE | cut -d "_" -f1,2`
        echo "Subsampling to $NUM reads for $SAMPLE"
	seqtk sample -s42 ${FILE}_1.fq.gz $NUM > ${SAMPLE}_subsample_1.fq
	seqtk sample -s42 ${FILE}_2.fq.gz $NUM > ${SAMPLE}_subsample_2.fq
	echo " "
	echo "Sampling done for $SAMPLE"
	echo " "
done

echo "Compressing all files..."
gzip *.fq
echo "Compressing done! Terminating..."
