#!/bin/bash
strain=$1
agrobacteria=$2

echo -e "************************************************************************************************"
echo -e "* SUMMARY:                                                                                     *"
echo -e "************************************************************************************************"

echo -e "SECRETION SYTEMS (MacSyFinder):"
cat ${strain}/${strain}_macsyfinder_TXSS/best_solution.tsv | cut -f 1,6 | sed 's/_[0-9]\+$//g' | sed 's/\t\S\+_/\t/g' | uniq | grep -v '#\|^$' | sed "s/${strain}.gbff_//g"
echo -e ""

echo -e "MOBILE GENETIC ELEMENTS:"
echo -e ""

echo -e "INTEGRATIVE CONJUGATIVE/MOBILIZABLE ELEMENTS (ICEs/IMEs; TIGER2):"
cat ${strain}/${strain}_TIGER2_final.table.out
echo -e ""
echo -e "(Note: only includes monopartite ICEs/IMEs. Polypartite elements are not shown)" 
echo -e ""

echo -e "PHAGE (DBSCAN-SWA):" 
cat ${strain}/prophage.table
echo -e ""

echo -e "INTEGRONS (IntegronFinder):"
cat ${strain}/integron.table
echo -e ""

echo -e "MGE DEFENSE SYSTEMS (DefenseFinder):"
echo -e "locus_tag  system"
cat ${strain}/${strain}_defensefinder.tsv.table
echo -e ""

echo -e "VIRULENCE/SYMBIOSIS LOCI:"
echo -e ""

echo -e "PREDICTED SECRETED EFFECTORS:"
echo -e "contig	start	end	locus_tag  product	gene"
grep 'Type III effector' ${strain}/${strain}.gff3 | grep 'CDS' | cut -f 1,4,5,9 | sed 's/ID=/\t/g;s/;Name=/\t/g;s/;locus_tag=.*//g'

echo -e ""
echo -e "(Note: not a de novo analysis, only homology to known, named effectors)"
echo -e ""

echo -e "NON-AGROBACTERIUM PHYTOPATHOGEN VIRULENCE LOCI:"
echo -e "contig	start	end	locus_tag  product	gene"
grep 'att locus\|fas locus\|FasDF\|gene=tomA\|gene=nec1\|gene=txt' ${strain}/${strain}.gff3 | grep 'CDS' | cut -f 1,4,5,9 | sed 's/ID=/\t/g;s/;Name=/\t/g;s/;locus_tag=.*;gene=/\t/g' 
echo -e ""

echo -e "NITROGEN FIXATION LOCI:"
echo -e "contig	start	end	locus_tag  product	gene"
grep 'Nodulation protein Nod\|Nitrogen fixation protein ' ${strain}/${strain}.gff3 | grep 'CDS' | cut -f 1,4,5,9 | sed 's/ID=/\t/g;s/;Name=/\t/g;s/;locus_tag=.*//g'
echo -e ""


if [[ $agrobacteria == "Agrobacterium" ]]; then
	echo -e ""
	echo -e ""
	echo -e "************************************************************************************************"
	echo -e "* AGROBACTERIUM SUMMARY:                                                                       *"
	echo -e "************************************************************************************************"
	echo -e ""

	echo -e "AGROBACTERIUM LINEAGE:"
	cat ${strain}/${strain}.agrobacteria_taxonomy.out
	echo -e ""
	
	echo -e "PREDICTED ONCOGENIC PLASMID:"
	cat ${strain}/${strain}.oncogenic_plasmid_final.out
	echo -e ""
	
	echo -e "T-DNA BORDERS/OVERDRIVE:"
	echo -e "contig	border	start	end"
	grep 'T-DNA\|overdrive' ${strain}/borders.table
	echo -e ""

	echo -e "TRA BOX/VIR BOX:"
	echo -e "contig	border	start	end"
	grep 'box' ${strain}/borders.table
	echo -e ""

	echo -e "VIR GENES:"
	echo -e "contig	start	end	locus_tag  product	gene"
	grep 'vir locus\|GALLS\|virF\|virH\|virK' ${strain}/${strain}.gff3 | grep 'CDS' | cut -f 1,4,5,9 | sed 's/ID=/\t/g;s/;Name=/\t/g;s/;locus_tag=.*//g'
	echo -e ""
	
	echo -e "OPINE SYNTHESIS:"
	echo -e "contig	start	end	locus_tag  product	gene"
	grep 'opine synthase\|nopaline synthase' ${strain}/${strain}.gff3 | grep 'CDS' | cut -f 1,4,5,9 | sed 's/ID=/\t/g;s/;Name=/\t/g;s/;locus_tag=.*//g'
	echo -e ""
	
	echo -e "AGROCIN84 SYNTHESIS:"
	echo -e "contig	start	end	locus_tag  product	gene"
	grep 'gene=agn' ${strain}/${strain}.gff3 | grep 'CDS' | cut -f 1,4,5,9 | sed 's/ID=/\t/g;s/;Name=/\t/g;s/;locus_tag=.*//g' 
	echo -e ""
fi

echo -e "************************************************************************************************"

