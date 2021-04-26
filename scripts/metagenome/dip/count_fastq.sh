#!/bin/bash

for FILE in `ls -1 *_1.fq.gz | sed 's/_1.fq.gz//'`
do
	READF=`echo $(zcat ${FILE}_1.fq.gz| wc -l)/4|bc`
	READR=`echo $(zcat ${FILE}_2.fq.gz| wc -l)/4|bc`
	echo "$FILE	$READF	$READR" >> log.tsv
done
