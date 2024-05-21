#!/bin/bash

beav_dir=$1

# Address to main beav gbk file
GBK=`ls $1/*_final.gbk`

# Contig list
CONTIG=`cat ${beav_dir}/*oncogenic_plasmid_final.out.contiglist | cut -f 1 | tr '\n' ' '`

# Conditionally run python script
if [ -z "$CONTIG" ]
then
    echo "python3 beav_circos.py --input $GBK"
    python3 $BEAV_DIR/scripts/beav_circos.py --input $GBK
else
    # Get plasmid type
    pTi=`cat $1/*oncogenic_plasmid_final.out | head -2 | tail -1 | cut -f1`
    echo "python3 beav_circos.py --input $GBK --pTi $CONTIG --plasmid $pTi"
    python3 $BEAV_DIR/scripts/beav_circos.py --input $GBK --pTi $CONTIG --plasmid $pTi
fi
