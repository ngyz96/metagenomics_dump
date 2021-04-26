#1/bin/bash

if [ $# == 0 ]; then
        echo "Usage: ./diamond.sh -i INDIR -o OUTDIR -d DATABASE [-t]"
	echo "Runs DIAMOND on all fastq files in a directory and output daa files (Now in paired end mode!)"
        echo "		-i: Full file path for directory with fastq files"
        echo "		-o: Full file path for output directory (Default = INDIR/diamond_out)"
        echo "		-d: Full file path for DIAMOND database located at /home/bmb/Soil-Microbiomes/databases (default=/home/bmb/Soil-Microbiomes/databases/diamond2021/nr.dmnd)"
        echo "		-t: Number of threads to use (Default=24)"
        exit 1
fi

while getopts "i:o:d:t:" opt
do
        case ${opt} in
        i) INDIR=${OPTARG} ;;
        o) OUTDIR=${OPTARG} ;;
        d) DATABASE=${OPTARG} ;;
        t) THREADS=${OPTARG} ;;
        esac
done

: ${INDIR:?Missing -i}
: ${THREADS:=24} ${DATABASE:="/home/bmb/Soil-Microbiomes/databases/diamond2021/nr.dmnd"} ${OUTDIR:="${INDIR}/diamond_out"}

module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_microbes

if [ ! -d "$OUTDIR" ]; then
	echo "output directory missing, creating directory..."
	mkdir -p $OUTDIR
fi

echo "*************************"
echo "** Starting DIAMOND... **"
echo "*************************"
cd $INDIR

for FILE in `ls -1 *_1.fq.gz | sed 's/_1.fq.gz//'`
do
        SAMPLE=`echo $FILE | cut -d "_" -f1,2`
	echo "Running DIAMOND for $SAMPLE"
        diamond blastx --threads $THREADS --db $DATABASE --query ${FILE}_1.fq.gz ${FILE}_2.fq.gz --log --out ${OUTDIR}/${SAMPLE}.daa --outfmt 100 -k 10 -b15.0 -c1
        echo "DIAMOND done for $SAMPLE"
done
