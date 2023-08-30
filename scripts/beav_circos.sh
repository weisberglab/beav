#!/bin/bash

beav_dir=$1

# Address to main beav gbk file
GBK=`ls $1/*gbk`

# Contig list
CONTIG=`cat ${beav_dir}/*oncogenic_plasmid_final.out.contiglist | cut -f 1 | head -n 1`

# Conditionally run python script
if [ -z $CONTIG ]
then
    python3 beav_circos.py $GBK
else
    python3 beav_circos.py $GBK $CONTIG
fi
