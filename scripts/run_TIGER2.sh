#!/usr/bin/env bash
#TIGER2 ICE finder

strain=$1

blastdb=$2
#/nfs7/BPP/Weisberg_Lab/projects/annotation_pipeline/databases/blast/allagro.fna

cpus=$3
#8

echo -e "TIGER2: preparing input"

#for bakta version, unset pythonpath and activate bakta conda

mkdir TIGER2
cd TIGER2
cp ../bakta/${strain}.fna ./genome.fa
echo -e "357	Bacteria sp. $strain	Bacteria;;;;;;; 11	$strain" > genome.tax
mkdir protein
cp ../bakta/${strain}.faa protein/protein.faa
cp ../bakta/${strain}.gff3 protein/protein.gff
oldlocus=`head -n1 protein/protein.faa | sed 's/>//g;s/_[0-9]\+.*//g'`
sed -i "s/$oldlocus/$strain/g" protein/protein*

echo -e "TIGER2: running Islander"
perl $BEAV_DIR/software/TIGER/bin/islander.pl -tax genome.tax -cpu $cpus genome.fa &> islander.log
#refseq genomic is very slow. make a lineage specific blast db (blastdb made using -parse_seqids)
echo -e "TIGER2: running Tiger"
perl $BEAV_DIR/software/TIGER/bin/tiger.pl -verbose -cpu $cpus -db $blastdb -fasta genome.fa &> tiger.log
echo -e "TIGER2: running Typing"
perl $BEAV_DIR/software/TIGER/bin/typing.pl genome.island.nonoverlap.gff &> typing.log
echo -e "TIGER2: running Resolve"
perl $BEAV_DIR/software/TIGER/bin/resolve.pl mixed lenient > resolved.log 
#does not work well#perl $BEAV_DIR/software/TIGER/bin/typing.pl resolve3.gff &> typing_final.log 

echo -e "TIGER2: parsing output"

declare -A aatable
aatable=( ["C"]="Cys" 
		  ["D"]="Asp"
		  ["S"]="Ser"
		  ["Q"]="Gln"
		  ["K"]="Lys"
		  ["I"]="Ile"
		  ["P"]="Pro"
		  ["T"]="Thr"
		  ["F"]="Phe"
		  ["N"]="Asn"
		  ["G"]="Gly"
		  ["H"]="His"
		  ["L"]="Leu"
		  ["R"]="Arg"
		  ["W"]="Trp"
		  ["A"]="Ala"
		  ["V"]="Val"
		  ["E"]="Glu"
		  ["Y"]="Tyr"
		  ["M"]="Met" )
cd ..
touch ${strain}_TIGER2_final.table.out 
while read line; do
	replicon=`echo -e "$line" |  cut -f 1`
	islandtype=`echo -e "$line" |  cut -f 2`
	startpos=`echo -e "$line" |  cut -f 4`
	endpos=`echo -e "$line" |  cut -f 5`

	target=`echo -e "$line" |  cut -f 9 | sed 's/;coord=.*//g;s/^.*target=//g'`
	targetlen=`expr length "$target"`
	if [[ "$targetlen" == 1 ]]; then
		target="tRNA_${aatable[$target]}"
	fi
	leftborderseq=`echo -e "$line" |  cut -f 9 | sed 's/^.*isleLseq=//g;s/;.*//g' | sed 's/[A-Z]//g'`
	rightborderseq=`echo -e "$line" |  cut -f 9 | sed 's/^.*isleRseq=//g;s/;.*//g' | sed 's/[A-Z]//g'`
	targetseq=`echo -e "$line" |  cut -f 9 | sed 's/^.*unintSeq=//g;s/;.*//g' | sed 's/[A-Z]//g'`
	icelength=`echo -e "$line" |  cut -f 9 | sed 's/^.*;len=//g;s/;.*//g'`
	tandem=`echo -e "$line" |  cut -f 9 | sed 's/^.*;tandem=//g;s/;.*//g' | sed 's/1,1/no/g;s/\([0-9]\),/count=\1,/g;s/,/,pos=/g'`
	
	#get loci for all integrases
	integraseposlist=`echo -e "$line" |  cut -f 9 | sed 's/^.*;ints=//g;s/;.*//g'`
	curloci=""
	while read intlocus; do
		inttype=`echo -e "$intlocus" | sed 's/\..*//g'` 
		intposstart=`echo -e "$intlocus" | sed 's/^.*://g;s/-.*//g'` 
		intposend=`echo -e "$intlocus" | sed 's/.*-//g'` 
		intlocus=`echo "$replicon	$intposstart	$intposend" | bedtools intersect -nonamecheck -a stdin -b ./${strain}.nofasta.gff3 -wb -f 0.90  | grep 'CDS' | sed 's/^.*ID=//g;s/;.*//g' | tr '\n' ',' | sed 's/,$//g'`
		curloci="$curloci,$inttype:$intlocus"
	done < <( echo -e "$integraseposlist" | tr ',' '\n' | sort -k2,2g -t : )  
	curloci=`echo -e "$curloci" | sed 's/^,//g'`
	echo -e "$replicon	$startpos	$endpos	TIGER2:integrative element (ICE/IME);target=$target;attL=$leftborderseq;attR=$rightborderseq;attB=$targetseq;length=$icelength;integrase=$curloci;tandem=$tandem;pred_model=$islandtype" >> ${strain}_TIGER2_final.table.out
done < ./TIGER2/resolve3.gff

echo -e "TIGER2: Done"
