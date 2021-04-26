#!/bin/bash

if [ $# == 0 ]; then
	echo "Usage: ./bowtie2_host_removal.sh -i INDIR -o OUTDIR -s SPECIES [-tkh]"
	echo "		-i: Full file path for directory with fastq files"
	echo "		-o: Full file path for output directory"
	echo "		-s: species name for Bowtie2 index located at /home/bmb/Soil-Microbiomes/bowtie2_index/ (exact match)"
	echo "		-t: Number of threads to use (Default=24)"
	echo "		-k: Whether to keep intermediate sam/bam files"
	echo "		-h: Whether to remove human contaminants"
	exit 1
fi

while getopts "i:o:s:t:hk" opt
do
	case ${opt} in 
	i) INDIR=${OPTARG} ;;
	o) OUTDIR=${OPTARG} ;;
	s) SPECIES=${OPTARG} ;;
	t) THREADS=${OPTARG} ;;
	k) KEEP=true ;;
	h) HUMAN=true ;;
	esac
done
echo $INDIR $OUTDIR $SPECIES $THREADS
: ${SPECIES:?Missing -s} ${INDIR:?Missing -i} ${OUTDIR:?Missing -o}
: ${THREADS:=24} ${KEEP:=false} ${HUMAN:=false}

module load anaconda2020/python3
eval "$(/usr/local/anaconda3-2020/bin/conda shell.bash hook)"
conda activate /usr/local_sbs/conda_microbes

BWT2DIR="/home/bmb/Soil-Microbiomes/bowtie2_index/${SPECIES}"
BWT2HUMAN="/home/bmb/Soil-Microbiomes/bowtie2_index/GRCh38_noalt_as"
cd $INDIR
echo "Removing $SPECIES contaminants"
for file in `ls -1 *_1.fq.gz | sed 's/_1.fq.gz//'`
do
	sample=`echo $file | cut -d "_" -f1,2`
	echo "removing $SPECIES contaminants for $sample"
	bowtie2 -1 ${file}_1.fq.gz -2 ${file}_2.fq.gz -x $BWT2DIR \
	-p $THREADS -S ${sample}_all.sam
	samtools view -bS ${sample}_all.sam > ${sample}_all.bam
	samtools view -b -f 12 -F 256 ${sample}_all.bam > ${sample}_unmapped.bam
	samtools sort -n -m 60G -@ $THREADS ${sample}_unmapped.bam -o ${sample}_unmapped_sorted.bam
	samtools fastq -@ $THREADS ${sample}_unmapped_sorted.bam \
	-1 $OUTDIR/${sample}_host_removed_1.fq.gz \
 	-2 $OUTDIR/${sample}_host_removed_2.fq.gz \
 	-0 /dev/null -s /dev/null -n
	echo "$SPECIES removed for $sample"
	if ! $KEEP
	then 
		echo "Removing intermediate sam/bam files"
		rm ${sample}_all.sam ${sample}_all.bam ${sample}_unmapped.bam ${sample}_unmapped_sorted.bam
	fi
done

if $HUMAN
then
	echo "Removing Human contaminants"
	cd $OUTDIR
	for file in `ls -1 *_host_removed_1.fq.gz | sed 's/_host_removed_1.fq.gz//'`
	do 
		echo "removing human contaminants for $file"
		bowtie2 -1 ${file}_host_removed_1.fq.gz -2 ${file}_host_removed_2.fq.gz -x $BWT2HUMAN \
        	-p $THREADS -S ${file}_host_removed_all.sam
        	samtools view -bS ${file}_host_removed_all.sam > ${file}_host_removed_all.bam
        	samtools view -b -f 12 -F 256 ${file}_host_removed_all.bam > ${file}_host_removed_unmapped.bam
        	samtools sort -n -m 60G -@ $THREADS ${file}_host_removed_unmapped.bam -o ${file}_host_removed_unmapped_sorted.bam
        	samtools fastq -@ $THREADS ${file}_host_removed_unmapped_sorted.bam \
        	-1 ${file}_human_host_removed_1.fq.gz \
        	-2 ${file}_human_host_removed_2.fq.gz \
        	-0 /dev/null -s /dev/null -n
		echo "Human contaminants done for $file!"
        	if ! $KEEP
        	then
                	echo "Removing intermediate sam/bam files"
                	rm ${file}_host_removed_all.sam ${file}_host_removed_all.bam \
			${file}_host_removed_unmapped.bam ${file}_host_removed_unmapped_sorted.bam
        	fi
	done
fi
echo "QC done!"
