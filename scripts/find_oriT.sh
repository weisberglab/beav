#!/bin/bash
#this strategy: https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkad964/7335748?rss=1&login=false
#origin of oriT seqs: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10123127/#sup1
#removed several duplicate and short oriTs
#20% hsp len min plus minimum 20 bp alignment

infile=$1
threads=$2

touch oriT_hits.stranded

makeblastdb -in $infile -dbtype nucl

blastn -task blastn-short -outfmt '6 std qlen slen qseq sseq' -num_threads $threads -dust no -query ${BEAV_DIR}/databases/oriT_db.fna -db $infile -qcov_hsp_perc 20 | awk '$4 >= 20' | while read line; do
	contig=`echo -e "$line" | cut -f 2`
	queryname=`echo -e "$line" | cut -f 1`
	startpos=`echo -e "$line" | cut -f 9`
	endpos=`echo -e "$line" | cut -f 10`
	evalue=`echo -e "$line" | cut -f 11`
	if [[ "$endpos" -gt "$startpos" ]]; then
		strand="+"
	else
		strand="-"
		temppos="$startpos"
		startpos="$endpos"
		endpos="$temppos"
	fi
	echo -e "$contig	$startpos	$endpos	$strand	$queryname	$evalue" >> oriT_hits.stranded 
done 

touch oriT.table

cat oriT_hits.stranded | sort -k1,1 -k2,2g | bedtools merge -c 4,5,6 -o collapse,collapse,collapse | while read curborderline; do
	contig=`echo -e "$curborderline" | cut -f 1`
	startpos=`echo -e "$curborderline" | cut -f 2`
	endpos=`echo -e "$curborderline" | cut -f 3`
	
	oriTevalues=`echo -e "$curborderline" | cut -f 6`
	oriTnames=`echo -e "$curborderline" | cut -f 5`
	oriTstrand=`echo -e "$curborderline" | cut -f 4`
	
	evallist=`echo -e "$oriTevalues" | sed 's/,/\n/g'`
	namelist=`echo -e "$oriTnames" | sed 's/,/\n/g'`
	strandlist=`echo -e "$oriTstrand" |  tr ',' '\n'`
	
	newname=`echo "$namelist" | paste - <(echo "$evallist") | paste - <(echo "$strandlist") | sort -k2,2g | head -n1 `
	newstrand=`echo -e "$newname" | cut -f 3`
	newannot=`echo -e "$newname" | cut -f 1-2 | tr '\t' ';'`
	echo -e "$contig	oriT	$startpos	$endpos	$newstrand	Reference:$newannot" >> oriT.table
done

rm oriT_hits.stranded

rm ${infile}.ndb  ${infile}.nhr  ${infile}.nin  ${infile}.njs  ${infile}.not  ${infile}.nsq  ${infile}.ntf  ${infile}.nto
