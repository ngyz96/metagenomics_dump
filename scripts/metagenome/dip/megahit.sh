#1/bin/bash

if [ $# == 0 ]; then
        echo "Usage: ./megahit.sh -i INDIR [-o OUTDIR -t THREADS]"
	echo "Runs megahit on all paired end fastq files in a directory (files ending with _{1,2}.fq.gz)"
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
conda activate /usr/local_sbs/conda_microbes

if [ ! -d "$OUTDIR" ]; then
	mkdir -p $OUTDIR
fi

echo "*************************************************************************************************************"
echo "*****************************************    Starting MEGAHIT...    *****************************************"
echo "*************************************************************************************************************"
echo " " 
cd $INDIR

for FILE in `ls -1 *_1.fq.gz | sed 's/_1.fq.gz//'`
do
        SAMPLE=`echo $FILE | cut -d "_" -f1,2`
        echo "Running MEGAHIT for $SAMPLE"
	megahit -1 ${FILE}_1.fq.gz -2 ${FILE}_2.fq.gz \
	--presets meta-large -o ./${SAMPLE}_megahit --out-prefix ${SAMPLE} \
	-t $THREADS --min-contig-len 500 -m 400000000000
	echo " "
	echo "assembly done for $SAMPLE"
done
