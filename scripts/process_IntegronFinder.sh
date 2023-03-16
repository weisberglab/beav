#!/bin/bash

strain="$1"
#strain="AcBa1656-2"
touch ${strain}.integronloci
touch ${strain}.integronborders

cat Results_Integron_Finder_${strain}/${strain}.integrons | grep -v '#\|ID_integron' | cut -f 1,2,4,5,6,9,11 | while read line; do
    integron=`echo -e "$line" | cut -f 1`
    contig=`echo -e "$line" | cut -f 2`
    startpos=`echo -e "$line" | cut -f 3`
    endpos=`echo -e "$line" | cut -f 4`
    strand=`echo -e "$line" | cut -f 5`
    locustype=`echo -e "$line" | cut -f 6`
    completeness=`echo -e "$line" | cut -f 7`
    locustag=`echo -e "$contig	$startpos	$endpos	$integron	$completeness	$strand	$locustype" | bedtools intersect -F 0.5 -a - -b ${strain}.gff3 -wb | grep 'CDS' | cut -f 16 | sed 's/^.*;locus_tag=//g;s/;.*//g'`
	if [[ ! -z $locustag  ]]; then
		echo -e "$locustag	IntegronFinder:within $integron;integron_type:$completeness;function:$locustype" | sed 's/;function:protein$//g' >> ${strain}.integronloci
	
	elif [[ "$locustype" == "attC" || "$locustype" == Pint* ]]; then
    	echo -e "$locustype	$contig	$startpos	$endpos	$strand	IntegronFinder:within $integron;integron_type:$completeness" >> ${strain}.integronborders
			
	fi
done
