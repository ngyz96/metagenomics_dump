#1/bin/bash

if [ $# == 0 ]; then
        echo "Usage: ./kraken2.sh -i INDIR -o OUTDIR -d DATABASE [-tb]"
	echo "Runs Kraken2 on all fasta files in a directory"
        echo "		-i: Full file path for directory with fasta files"
        echo "		-o: Full file path for output directory (Output folders will be created for each sample)"
        echo "		-d: Full file path for Kraken2 database located at /home/bmb/Soil-Microbiomes/databases (default=/home/bmb/Soil-Microbiomes/databases/kraken2)"
        echo "		-t: Number of threads to use (Default=24)"
	echo "		-b: Run Bracken on Kraken2 output (Default to 150 bp Bracken database and Genus level, change in script if required)"
        exit 1
fi

while getopts "i:o:d:t:b" opt
do
        case ${opt} in
        i) INDIR=${OPTARG} ;;
        o) OUTDIR=${OPTARG} ;;
        d) DATABASE=${OPTARG} ;;
        t) THREADS=${OPTARG} ;;
        b) BRACKEN=true ;;
        esac
done

: ${INDIR:?Missing -i} ${OUTDIR:?Missing -o}
: ${THREADS:=24} ${BRACKEN:=false} ${DATABASE:="/home/bmb/Soil-Microbiomes/databases/kraken2"}

READ_LEN=150
LEVEL="G"
THRESHOLD=10

module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_microbes

if [ ! -d "$OUTDIR" ]; then
	mkdir -p $OUTDIR
fi

echo "***************************************************************************"
echo "************************** Running Kraken2... *****************************"
echo "***************************************************************************"

for FILE in `ls -1 *.fa | sed 's/.fa//'`
do
        SAMPLE=`echo $FILE | cut -d "." -f1`
        echo "Running Kraken2 for $SAMPLE"
	if [ ! -d ${OUTDIR}/${SAMPLE}_output ]; then
		echo "Creating $SAMPLE output directory"
		mkdir -p ${OUTDIR}/${SAMPLE}_output;
	fi
        kraken2 --db ${DATABASE} --threads ${THREADS} --report ${OUTDIR}/${SAMPLE}_output/${SAMPLE}.kreport \
	--output ${OUTDIR}/${SAMPLE}_output/${SAMPLE}.kraken ${FILE}.fa
        echo "Kraken2 done for $SAMPLE"
        if $BRACKEN
        then
                echo "Running Bracken on $SAMPLE Kraken2 output"
                bracken -d ${DATABASE} -i ${OUTDIR}/${SAMPLE}_output/${SAMPLE}.kreport \
		-o ${OUTDIR}/${SAMPLE}_output/${SAMPLE}.bracken -r ${READ_LEN} -l ${LEVEL} -t ${THRESHOLD}
        fi
done
