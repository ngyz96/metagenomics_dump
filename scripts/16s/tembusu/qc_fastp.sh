#!/bin/bash

module load sbs/fastp/0.20.0
INDIR="Soil-Microbiomes/Tembusu-Microbiome/16S/result/00.RawData"
OUTDIR="Soil-Microbiomes/Tembusu-Microbiome/16S/fastp_data/raw"

for file in `ls -1 $INDIR/**/*_1.fq.gz | grep 'raw' | sed 's/_1.fq.gz//'`
do
	filename=`echo $file | cut -d '/' -f10`
	fastp --in1 ${file}_1.fq.gz --in2 ${file}_2.fq.gz \
	--out1 ${OUTDIR}/${filename}_cleaned_1.fq.gz --out2 ${OUTDIR}/${filename}_cleaned_2.fq.gz \
	-w 8 -j ${OUTDIR}/${filename}_fastp.json -h ${OUTDIR}/${filename}_fastp.html
done
