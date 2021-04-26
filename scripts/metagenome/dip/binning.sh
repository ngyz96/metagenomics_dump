#1/bin/bash

if [ $# == 0 ]; then
        echo "Usage: ./bwamem.sh -a ASSEMBLY -r READS [-o OUTDIR -t THREADS]"
	echo "Runs metaWRAP binning module on all paired end fastq files on their corresponding assembly"
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

module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_metawrap

if [ ! -d "$OUTDIR" ]; then
	mkdir -p $OUTDIR
fi

echo "**************************************************************************************************************"
echo "*****************************************    Starting metaWRAP...    *****************************************"
echo "**************************************************************************************************************"
echo " " 
cd $READS

echo "Unzipping files and writing new files for metawrap module..."
echo " "
for OLD in `ls -1 *.fq.gz | sed 's/.fq.gz//'`
do
	gunzip -c ${OLD}.fq.gz > ${OLD}.fastq
done 
echo "done with preprocessing... moving on to binning..."
echo " "

for FILE in `ls -1 *_1.fastq | sed 's/_1.fastq//'`
do
	SAMPLE=`echo $FILE | cut -d "_" -f1,2`
	if [ ! -f "${ASSEMBLY}/${SAMPLE}.contigs.fa" ]; then
		echo "Assembly missing for $SAMPLE!"
		continue
	fi
	if [ -f "${FILE}_2.fastq" ]; then
		echo "Running metaWRAP binning for $SAMPLE"
		echo "Creating dir for $SAMPLE"
		if [ ! -d "${OUTDIR}/${SAMPLE}" ]; then
			mkdir -p ${OUTDIR}/${SAMPLE}
		fi
		metawrap binning -o ${OUTDIR}/${SAMPLE} -t $THREADS -m 240 -a ${ASSEMBLY}/${SAMPLE}.contigs.fa \
		--metabat2 --maxbin2 --concoct ${FILE}_1.fastq ${FILE}_2.fastq
		echo " "
		echo "binning done for $SAMPLE"
	fi
done

echo "binning done for all samples, removing unzipped files..."
rm *.fastq
echo "clean up done! terminating..."



