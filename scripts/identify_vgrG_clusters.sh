#!/bin/bash

strain=$1

hmmsearch --tblout vgrG.TIGR01646.1.table.out ${BEAV_DIR}/models/TIGR01646.1.hmm ./bakta/${strain}.faa > /dev/null
hmmsearch --tblout vgrG.TIGR03361.1.table.out ${BEAV_DIR}/models/TIGR03361.1.hmm ./bakta/${strain}.faa > /dev/null

cat vgrG.TIGR*.out | grep -v '#' | sed 's/\s.*//g' | sort | uniq | while read vgrGlocus; do
	operon=`grep "^$vgrGlocus	" operon-mapper_results/list_of_operons.table | cut -f 2`
	othergenes=`grep "	$operon$" operon-mapper_results/list_of_operons.table | cut -f 1`
	
	curout=""
	while read curgene; do
		found=`grep "ID=$curgene" ./bakta/${strain}.gff3 | grep 'CDS' | cut -f 1,4,5,7,9 | sed 's/ID=.*locus_tag=//g;s/;.*product=/\t/g;s/;.*gene=/\t/g;s/;.*//g'`
		firstpart=`echo -e "$found" | cut -f 1-5 -d '	'`
		productpart=`echo -e "$found" | cut -f 6 -d '	'`
		genepart=`echo -e "$found" | cut -f 7 -d '	'`
		found="$firstpart	$genepart	$productpart"
		curout="${curout}${curgene}	$found\n"
	done < <(echo "$othergenes")
	
	
	echo -e "vgrG-like $vgrGlocus	${operon}:"
	
	strand=`echo -e "$curout" | head -n 1 | cut -f 5`
	if [[ "$strand" == "+" ]]; then
		echo -e "$curout" | head -n-1
	else
		echo -e "$curout" | head -n-1 | tac
	fi
	
	echo -e ""
done

rm vgrG.TIGR01646.1.table.out 
rm vgrG.TIGR03361.1.table.out
