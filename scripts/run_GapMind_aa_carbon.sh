#!/usr/bin/env bash

strain=$1
numCPUs=$2

GapMindPath="$BEAV_DIR/software/PaperBLAST/"

echo -e "GapMind: Running AA biosynthesis pathway analysis"

mkdir GapMind
$GapMindPath/bin/buildorgs.pl -out GapMind/orgs -orgs "file:${strain}.faa:$strain"

#AA biosynthesis
# Use bin/gapsearch.pl to compare the queries (one file per pathway) to
# the proteins.
#usearch
#$GapMindPath/bin/gapsearch.pl -orgs GapMind/orgs -set aa -out GapMind/aa.hits -nCPU $numCPUs
#diamond
diamond makedb --in GapMind/orgs.faa -d GapMind/orgs.aa.dmnd
$GapMindPath/bin/gapsearch.pl -diamond -orgs GapMind/orgs -set aa -out GapMind/aa.hits -nCPU $numCPUs

# Use bin/gaprevsearch.pl to see if the candidates are similar to
# characterized proteins that have other functions.
#usearch
#$GapMindPath/bin/gaprevsearch.pl -orgs GapMind/orgs -hits GapMind/aa.hits -curated $GapMindPath/tmp/path.aa/curated.faa.udb -out GapMind/aa.revhits -nCPU $numCPUs
#diamond
$GapMindPath/bin/gaprevsearch.pl -diamond -orgs GapMind/orgs -hits GapMind/aa.hits -curated $GapMindPath/tmp/path.aa/curated.faa.dmnd -out GapMind/aa.revhits -nCPU $numCPUs

# Use bin/gapsummary.pl to score all the candidates and steps and to
# find the best-scoring path for each pathway. This produces three
# tables: one for each pathway and rule (including the rule "all" for
# the entire pathway); one for each step; and one for each candidate for
# each step.
$GapMindPath/bin/gapsummary.pl -orgs GapMind/orgs -set aa -hits GapMind/aa.hits -rev GapMind/aa.revhits -out GapMind/aa.sum

## Optionally, use bin/checkGapRequirements.pl to check dependencies
## between pathways. The output table will list any warnings.
$GapMindPath/bin/checkGapRequirements.pl -org GapMind -stepsDb $GapMindPath/tmp/path.aa/steps.db -results . -set aa -out GapMind/aa.sum.warn

## Optionally, use orgsVsMarkers.pl to compare the genome to organisms
## with known gaps (for amino acid biosynthesis only).
#$GapMindPath/bin/orgsVsMarkers.pl -orgs GapMind/orgs -vs $GapMindPath/gaps/aa/aa.known.gaps.markers.faa -out GapMind/aa.sum.knownsim

## Optionally, use buildGapsDb.pl to combine the results into a sqlite3
## database.
##bin/buildGapsDb.pl -gaps GapMind/aa.sum -requirements GapMind/aa.sum.warn -steps GapMind/path.aa/steps.db -out GapMind/aa.sum.db -markersim GapMind/aa.sum.knownsim

#process and combine GapMind AA biosynthesis results
awk '$6 > 0' GapMind/aa.sum.steps | cut -f 4- | cut -f 1,2,4,5 | grep -v '      ' | sort -k3 -r | awk -F"\t" '!seen[$1, $4]++' | sed 's/\(\S\+\)\t\(\S\+\)\t\(\S\+\)$/\3:\2:\1/g' | awk 'BEGIN{FS="\t"; OFS=FS}; { arr[$2] = arr[$2] == ""? $1 : arr[$2] "," $1 } END {for (i in arr) print i, arr[i] }' | sed 's/:/\t/g' | sort -k1 | grep -v '^locusId' | while read line; do
    locus=`echo -e "$line" | cut -f 1`
    confidence=`echo -e "$line" | cut -f 2 | sed 's/0/low/g;s/1/medium/g;s/2/high/g'`
    geneid=`echo -e "$line" | cut -f 3`
    aa=`echo -e "$line" | cut -f 4`
    echo -e "$locus	GapMind:$aa biosynthesis; gene=$geneid; $confidence confidence" >> GapMind/combined_GapMind_results.tab
done

echo -e "GapMind: Running small molecule carbon catabolism analysis"
#carbon metabolism
$GapMindPath/bin/gapsearch.pl -diamond -orgs GapMind/orgs -set carbon -out GapMind/carbon.hits -nCPU $numCPUs
#$GapMindPath/bin/gapsearch.pl -orgs GapMind/orgs -set carbon -out GapMind/carbon.hits -nCPU $numCPUs
$GapMindPath/bin/gaprevsearch.pl -diamond -orgs GapMind/orgs -hits GapMind/carbon.hits -curated $GapMindPath/tmp/path.carbon/curated.faa.dmnd -out GapMind/carbon.revhits -nCPU $numCPUs
#$GapMindPath/bin/gaprevsearch.pl -orgs GapMind/orgs -hits GapMind/carbon.hits -curated $GapMindPath/tmp/path.carbon/curated.faa.udb -out GapMind/carbon.revhits -nCPU $numCPUs
$GapMindPath/bin/gapsummary.pl -orgs GapMind/orgs -set carbon -hits GapMind/carbon.hits -rev GapMind/carbon.revhits -out GapMind/carbon.sum
$GapMindPath/bin/checkGapRequirements.pl -org GapMind -stepsDb $GapMindPath/tmp/path.carbon/steps.db -results . -set carbon -out GapMind/carbon.sum.warn

#process and combine GapMind carbon catabolism results
awk '$6 > 0' GapMind/carbon.sum.steps | cut -f 4- | cut -f 1,2,4,5 | grep -v '		' | sort -k3 -r | awk -F"\t" '!seen[$1, $4]++' | sed 's/\(\S\+\)\t\(\S\+\)\t\(\S\+\)$/\3:\2:\1/g' | awk 'BEGIN{FS="\t"; OFS=FS}; { arr[$2] = arr[$2] == ""? $1 : arr[$2] "," $1 } END {for (i in arr) print i, arr[i] }' | sed 's/:/\t/g' | sort -k1 | grep -v '^locusId' | while read line; do
	locus=`echo -e "$line" | cut -f 1`
	confidence=`echo -e "$line" | cut -f 2 | sed 's/0/low/g;s/1/medium/g;s/2/high/g'`
	geneid=`echo -e "$line" | cut -f 3`
	carbon=`echo -e "$line" | cut -f 4`
	echo -e "$locus	GapMind:$carbon catabolism; gene=$geneid; $confidence confidence" >> GapMind/combined_GapMind_results.tab
done

echo -e "GapMind: Done"
