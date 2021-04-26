#/bin/bash

for DIR in `ls -1`
do
	CONCOCT=`ls -1 ${DIR}/concoct_bins/ | wc -l`
	MAXBIN2=`ls -1 ${DIR}/maxbin2_bins/ | wc -l`
	METABAT2=`ls -1 ${DIR}/metabat2_bins/ | wc -l`
	
	echo "$DIR	$(($CONCOCT-1))	$MAXBIN2	$(($METABAT2-1))" >> bin_counts.tsv
done
