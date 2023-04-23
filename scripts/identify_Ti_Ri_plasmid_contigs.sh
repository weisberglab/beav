#!/usr/bin/env bash

infile=$1
strain=`echo -e "$infile" | sed 's/.fna//g'`

inref=`ls -1 $BEAV_DIR/databases/oncogenic_plasmids/Weisberg2022PhilTransB/*.sketch $BEAV_DIR/databases/oncogenic_plasmids/burr/*.sketch | tr '\n' ','`

echo "Identifying oncogenic (Ti/Ri) plasmid contigs"

comparesketch.sh in=$infile ref=$inref minwkid=0.5 out=${strain}.oncogenic_plasmid_type.comparesketch.out
comparesketch.sh in=$infile ref=$inref mode=sequence minwkid=0.5 out=${strain}.oncogenic_plasmid_contigs.comparesketch.out 
echo ""
echo ""
echo "oncogenic plasmid: (best hit)" | tee ${strain}.oncogenicplasmid_final.out
grep 'WKID' -A 1 ${strain}.oncogenic_plasmid_type.comparesketch.out | tail -n 1 | sed 's/^.*\s//g' | sed 's/__/\t/g' | tee -a ${strain}.oncogenic_plasmid_final.out
echo "" | tee -a ${strain}.oncogenic_plasmid_final.out
echo "oncogenic plasmid contigs: (best hit)" | tee -a ${strain}.oncogenic_plasmid_final.out
grep 'Query' ${strain}.oncogenic_plasmid_contigs.comparesketch.out -A2 | grep -v 'WKID' | sed 's/ \[.*/\t/g' | sed 's/Query: //g' | grep -v '^--$' | grep -v '^$' | sed 's/^[0-9]\+.[0-9]\+.*\s//g' | paste -d " " - - | grep -v 'No hits.' | tee -a ${strain}.oncogenic_plasmid_final.out | tee ${strain}.oncogenic_plasmid_final.out.contiglist
