#!/usr/bin/env python
from Bio import SeqIO
import glob, os

faa_filename = "allgembase.faa"
topo_filename = "allgembase.topo"

output_handle = open(faa_filename, "w")
output_handletopo = open(topo_filename, "w")


for genbank_file in glob.glob("*.gbff"):
    input_handle  = open(genbank_file, "r")
    for seq_record in SeqIO.parse(input_handle, "genbank") :
        output_handletopo.write("%s_%s : circular\n" % (
            genbank_file.replace(".gbk",""),
            seq_record.name
        ))
        #print "Dealing with GenBank record %s" % seq_record.id
        for seq_feature in seq_record.features :
            if seq_feature.type=="CDS" and 'translation' in seq_feature.qualifiers:
                shortlocus = seq_feature.qualifiers['locus_tag'][0].replace("_","=")
                output_handle.write(">%s_%s_%s %s %s\n%s\n" % (
                   genbank_file.replace(".gbk",""),
                   seq_record.name,
                   shortlocus,
                   seq_feature.qualifiers['locus_tag'][0],
                   seq_feature.qualifiers['product'][0],
                   seq_feature.qualifiers['translation'][0]))
    input_handle.close()

output_handle.close()
