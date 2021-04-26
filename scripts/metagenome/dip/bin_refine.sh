#1/bin/bash

if [ $# == 0 ]; then
        echo "Usage: ./bin_refine.sh -i INDIR [-o OUTDIR -1 ALGO1 -2 ALGO2 -3 ALGO3 -t THREADS -c COMPLETION -x CONTAMINATION]"
	echo "Runs metaWRAP bin refinement module on all binning folders produced by the metaWRAP binning module"
        echo "		-i: file path for the input directory will all the binning folders from metawrap binning module"
        echo "		-o: Full file path for output directory (INDIR/BIN_REFINEMENT)"
	echo "		-1: Prefix of the first binning algorithm (DEFAULT=maxbin2)"
	echo "		-2: Prefix of the second binning algorithm (DEFAULT=metabat2)"
	echo "		-3: Prefix of the third binning algorithm (DEFAULT=concoct)"
        echo "		-t: Number of threads to be used (DEFAULT=24)"
	echo "		-c: Completion percentage (DEFAULT=90)"
	echo "		-x: Contamination percentage (DEFAULT=5)"
        exit 1
fi

while getopts "i:o:1:2:3:t:c:x:" opt
do
        case ${opt} in
        i) INDIR=${OPTARG} ;;
	o) OUTDIR=${OPTARG} ;;
        t) THREADS=${OPTARG} ;;
	c) COMPLETION=${OPTARG} ;;
	x) CONTAMINATION=${OPTARG} ;;
	1) ALGO1=${OPTARG} ;;
	2) ALGO2=${OPTARG} ;;
	3) ALGO3=${OPTARG} ;;
        esac
done

: ${INDIR:?Missing -i}
: ${THREADS:=24} ${OUTDIR:="${INDIR}/BIN_REFINEMENT"} ${COMPLETION:=90} ${CONTAMINATION:=5}
: ${ALGO1:="maxbin2"} ${ALGO2:="metabat2"} ${ALGO3:="concoct"}

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

cd $INDIR
for DIR in `ls -d1 *_*`
do
	if [ ! -d "${OUTDIR}/${DIR}" ]; then
		mkdir -p ${OUTDIR}/${DIR}
	fi
	metawrap bin_refinement -o ${OUTDIR}/${DIR} -t $THREADS \
	-A ${INDIR}/${DIR}/${ALGO1}_bins -B ${INDIR}/${DIR}/${ALGO2}_bins -C ${INDIR}/${DIR}/${ALGO3}_bins \
	-c $COMPLETION -x $CONTAMINATION
done

echo "all binning refinement done... exiting..."
