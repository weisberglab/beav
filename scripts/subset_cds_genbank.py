#!/bin/python
import sys
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
strain = sys.argv[1]

cds_path = f"./{strain}_cds_subset.gbk"
cds_record = []
for record in SeqIO.parse(f"./{strain}.gbff","gb"):
    has_genes = False
    for feature in record.features:
        if feature.type == "CDS":
            has_genes = True
            break 
    if has_genes:
        cds_record.append(record)
cds_output = open(cds_path, 'w')
SeqIO.write(cds_record,cds_output,"genbank")
cds_output.close()            
