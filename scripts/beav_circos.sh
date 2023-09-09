#!/bin/bash

beav_dir=$1

# Address to main beav gbk file
GBK=`ls $1/*gbk`

# Contig list
CONTIG=`cat ${beav_dir}/*oncogenic_plasmid_final.out.contiglist | cut -f 1 | tr '\n' ' '`

# Conditionally run python script
if [ -z "$CONTIG" ]
then
    echo "python3 beav_circos.py --input $GBK"
    python3 beav_circos.py --input $GBK
else
    echo "python3 beav_circos.py --input $GBK --contig $CONTIG"
    python3 beav_circos.py --input $GBK --contig $CONTIG
fi
