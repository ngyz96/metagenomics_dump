#1/bin/bash

if [ $# == 0 ]; then
        echo "Usage: ./kraken2.sh -i INDIR -o OUTDIR [-d DATABASE -t THREADS -l LEVEL -b]"
	echo "Runs Kraken2 on all paired end fastq files in a directory (files ending with _{1,2}.fq.gz)"
        echo "		-i: Full file path for directory with fastq files"
        echo "		-o: Full file path for output directory (Output folders will be created for each sample)"
        echo "		-d: Full file path for Kraken2 database located at /home/bmb/Soil-Microbiomes/databases (default=/home/bmb/Soil-Microbiomes/databases/kraken2)"
        echo "		-t: Number of threads to use (Default=24)"
	echo "		-b: Run Bracken on Kraken2 output (Default to 150 bp Bracken database, change in script if required)"
        echo "		-l: level to run Bracken on (Default = P, can any taxonomical level)"
	exit 1
fi

while getopts "i:o:d:t:l:b" opt
do
        case ${opt} in
        i) INDIR=${OPTARG} ;;
        o) OUTDIR=${OPTARG} ;;
        d) DATABASE=${OPTARG} ;;
        t) THREADS=${OPTARG} ;;
        b) BRACKEN=true ;;
	l) LEVEL=${OPTARG} ;;
        esac
done

: ${INDIR:?Missing -i} ${OUTDIR:?Missing -o}
: ${THREADS:=24} ${BRACKEN:=false} ${DATABASE:="/home/bmb/Soil-Microbiomes/databases/kraken2"} ${LEVEL:="P"}

READ_LEN=150
THRESHOLD=10

module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_microbes

if [ ! -d "$OUTDIR" ]; then
	mkdir -p $OUTDIR
fi

echo "*************************"
echo "** Starting Kraken2... **"
echo "*************************"
cd $INDIR

for FILE in `ls -1 *_1.fq.gz | sed 's/_1.fq.gz//'`
do
        SAMPLE=`echo $FILE | cut -d "_" -f1,2`
        echo "Running Kraken2 for $SAMPLE"
	if [ ! -d ${OUTDIR}/${SAMPLE}_output ]; then
		echo "Creating $SAMPLE output directory"
		mkdir -p ${OUTDIR}/${SAMPLE}_output;
	fi
        kraken2 --db ${DATABASE} --threads ${THREADS} --report ${OUTDIR}/${SAMPLE}_output/${SAMPLE}.kreport \
	--paired ${FILE}_1.fq.gz ${FILE}_2.fq.gz --output ${OUTDIR}/${SAMPLE}_output/${SAMPLE}.kraken 
        echo "Kraken2 done for $SAMPLE"
        if $BRACKEN
        then
                echo "Running Bracken on $SAMPLE Kraken2 output"
                bracken -d ${DATABASE} -i ${OUTDIR}/${SAMPLE}_output/${SAMPLE}.kreport \
		-o ${OUTDIR}/${SAMPLE}_output/${SAMPLE}.bracken -r ${READ_LEN} -l ${LEVEL} -t ${THRESHOLD}
        fi
done
