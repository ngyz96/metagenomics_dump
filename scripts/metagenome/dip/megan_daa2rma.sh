#1/bin/bash

if [ $# == 0 ]; then
        echo "Usage: ./diamond.sh -i INDIR -o OUTDIR -d MEGANMAP [-t]"
	echo "Runs daa2rma on paired daa files in a directory and output rma files for MEGAN"
        echo "		-i: Full file path for directory with paired daa"
        echo "		-o: Full file path for output directory (Default = INDIR/rma_out)"
        echo "		-d: Full file path for DIAMOND database located at /home/bmb/Soil-Microbiomes/databases (default=/home/bmb/Soil-Microbiomes/databases/megan_map/megan-map-Jan2021.db)"
        echo "		-t: Number of threads to use (Default=24)"
	echo "		-a: Algorithm for the LCA. (Default=naive, only accepts naive, weighted or longReads)"
        exit 1
fi

while getopts "i:o:d:t:a:" opt
do
        case ${opt} in
        i) INDIR=${OPTARG} ;;
        o) OUTDIR=${OPTARG} ;;
        d) DATABASE=${OPTARG} ;;
        t) THREADS=${OPTARG} ;;
	a) ALGO=${OPTARG}
        esac
done

: ${INDIR:?Missing -i}
: ${THREADS:=24} ${DATABASE:="/home/bmb/Soil-Microbiomes/databases/megan_map/megan-map-Jan2021.db"} ${OUTDIR:="${INDIR}/rma_out"} ${ALGO:="naive"}

if [ ! -d "$OUTDIR" ]; then
	echo "output directory missing, creating directory..."
	mkdir -p $OUTDIR
fi

echo "******************************************************************************************************"
echo "************************************* Starting daa2rma... ********************************************"
echo "******************************************************************************************************"
cd $INDIR

for FILE in `ls -1 *_1.daa | sed 's/_1.daa//'`
do
	echo "Running daa2rma for $SAMPLE"
        /usr/local_sbs/MEGAN6/tools/daa2rma -i ${FILE}_1.daa ${FILE}_2.daa -o ${OUTDIR}/${FILE}.rma -p true \
	-alg ${ALGO} -mdb ${DATABASE} -t ${THREADS} -v 
        echo "DIAMOND done for $SAMPLE"
done
