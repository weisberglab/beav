#!/usr/bin/env bash

infile=$1
input=`basename $infile`
strain=`echo -e "$input" | sed 's/.fna//g'`
inref=`ls -1 $BEAV_DIR/databases/oncogenic_plasmids/Weisberg2022PhilTransB/*.sketch $BEAV_DIR/databases/oncogenic_plasmids/burr/*.sketch | tr '\n' ','`

echo "Identifying oncogenic (Ti/Ri) plasmid contigs"

comparesketch.sh in=$input ref=$inref minwkid=0.5 out=${strain}.oncogenic_plasmid_type.sketch.out 2>&1 | tee oncogenicplasmid.log
comparesketch.sh in=$input ref=$inref mode=sequence minwkid=0.5 out=${strain}.oncogenic_plasmid_contigs.sketch.out 2>&1 | tee --append oncogenicplasmid.log 
echo ""
echo ""
echo "oncogenic plasmid type: (best hit reference strain)" | tee ${strain}.oncogenic_plasmid_final.out | tee --append oncogenicplasmid.log
grep 'WKID' -A 1 ${strain}.oncogenic_plasmid_type.sketch.out | tail -n 1 | sed 's/^.*\s//g' | sed 's/__/\t/g' | tee -a ${strain}.oncogenic_plasmid_final.out | tee --append oncogenicplasmid.log
echo "" | tee -a ${strain}.oncogenic_plasmid_final.out | tee --append oncogenicplasmid.log
echo "oncogenic plasmid-like contigs: (best hit reference strain)" | tee -a ${strain}.oncogenic_plasmid_final.out | tee --append oncogenicplasmid.log
echo "contig	length	weighted_k-mer_id	k-mer_id	plasmid_type	best_hit_ref_strain" | tee -a ${strain}.oncogenic_plasmid_final.out | tee --append oncogenicplasmid.log

grep 'WKID' ${strain}.oncogenic_plasmid_contigs.sketch.out -A1 -B1 | grep -v '^--$' | grep -v 'WKID' | sed 's/^Query: //g' | sed 's/\s.*Bases: /___BASES___/g' | sed 's/\s.*File:.*$//g' | paste -d " " - - | sed 's/\s\+/\t/g' | cut -f 1,2,3,14 | sed 's/___BASES___/\t/g' | sed 's/__\(\S\+\)$/\t\1/g' | tee -a ${strain}.oncogenic_plasmid_final.out | tee ${strain}.oncogenic_plasmid_final.out.contiglist | tee --append oncogenicplasmid.log
