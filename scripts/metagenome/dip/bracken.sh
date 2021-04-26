#1/bin/bash

if [ $# == 0 ]; then
        echo "Usage: ./bracken.sh -i INDIR -o OUTDIR -d DATABASE [-tb]"
	echo "Runs Bracken on kreport files generated from Kraken/Kraken2"
        echo "		-i: Full file path for directory with kreport files"
        echo "		-o: Full file path for output directory (Default = INDIR/bracken_out)"
        echo "		-d: Full file path for Bracken database located at /home/bmb/Soil-Microbiomes/databases (Default=/home/bmb/Soil-Microbiomes/databases/kraken2)"
        echo "		-t: Number of threads to use (Default=24)"
	echo "		-l: taxonomic level to run Bracken for (Default=S, can be K, P, C, O, F, G, S)"
        echo "		-r: length of reads that Kraken/Kraken2 was generated from (Default=150)"
	exit 1
fi

while getopts "i:o:d:t:l:r:" opt
do
        case ${opt} in
        i) INDIR=${OPTARG} ;;
        o) OUTDIR=${OPTARG} ;;
        d) DATABASE=${OPTARG} ;;
        t) THREADS=${OPTARG} ;;
        l) LEVEL=${OPTARG} ;;
	r) READLEN=${OPTARG} ;;
        esac
done

: ${INDIR:?Missing -i}
: ${THREADS:=24} ${BRACKEN:=false} ${DATABASE:="/home/bmb/Soil-Microbiomes/databases/kraken2"}
: ${OUTDIR:="${INDIR}/bracken_out"} ${LEVEL:="S"} ${READLEN:=150}

module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_microbes

if [ ! -d "$OUTDIR" ]; then
	mkdir -p $OUTDIR
fi
THRESHOLD=10
echo "************************************************************************************************************************"
echo "************************************************ Starting Bracken... ***************************************************"
echo "************************************************************************************************************************"
cd $INDIR

for FILE in `ls -1 *.kreport | sed 's/.kreport//'`
do
        echo "Running Bracken on $FILE Kraken2 output"
        bracken -d ${DATABASE} -i ${INDIR}/${FILE}.kreport \
		-o ${OUTDIR}/${FILE}.bracken -r ${READLEN} -l ${LEVEL} -t ${THRESHOLD}
done
