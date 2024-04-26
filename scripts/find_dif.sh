#!/bin/bash
#20% hsp len min plus minimum 20 bp alignment

infile=$1
threads=$2

touch dif_hits.stranded

makeblastdb -in $infile -dbtype nucl 1> /dev/null
blastn -task blastn-short -outfmt '6 std qlen slen qseq sseq' -num_threads $threads -dust no -query ${BEAV_DIR}/databases/dif_db.fna -db $infile -qcov_hsp_perc 20 | awk '$4 >= 20' | awk '$11 <= 0.01' | while read line; do
	contig=`echo -e "$line" | cut -f 2`
	queryname=`echo -e "$line" | cut -f 1`
	startpos=`echo -e "$line" | cut -f 9`
	endpos=`echo -e "$line" | cut -f 10`
	evalue=`echo -e "$line" | cut -f 11`
	if [[ "$endpos" -gt "$startpos" ]]; then
		strand="+1"
	else
		strand="-1"
		temppos="$startpos"
		startpos="$endpos"
		endpos="$temppos"
	fi
	echo -e "$contig	$startpos	$endpos	$strand	$queryname	$evalue" >> dif_hits.stranded 
done 

touch dif.table

cat dif_hits.stranded | sort -k1,1 -k2,2g | bedtools merge -c 4,5,6 -o collapse,collapse,collapse | while read curborderline; do
	contig=`echo -e "$curborderline" | cut -f 1`
	startpos=`echo -e "$curborderline" | cut -f 2`
	endpos=`echo -e "$curborderline" | cut -f 3`
	
	difevalues=`echo -e "$curborderline" | cut -f 6`
	difnames=`echo -e "$curborderline" | cut -f 5`
	difstrand=`echo -e "$curborderline" | cut -f 4`
	
	evallist=`echo -e "$difevalues" | sed 's/,/\n/g'`
	namelist=`echo -e "$difnames" | sed 's/,/\n/g'`
	strandlist=`echo -e "$difstrand" |  tr ',' '\n'`
	
	newname=`echo "$namelist" | paste - <(echo "$evallist") | paste - <(echo "$strandlist") | sort -k2,2g | head -n1 `
	newstrand=`echo -e "$newname" | cut -f 3`
	newannot=`echo -e "$newname" | cut -f 1-2 | tr '\t' ';'`
	echo -e "$contig	$startpos	$endpos	100	dif	$newstrand">> dif.table

done

rm dif_hits.stranded

rm ${infile}.ndb  ${infile}.nhr  ${infile}.nin  ${infile}.njs  ${infile}.not  ${infile}.nsq  ${infile}.ntf  ${infile}.nto
