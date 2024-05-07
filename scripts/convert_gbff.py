#!/bin/python

import sys
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from BCBio import GFF
strain = sys.argv[1]
gbk_filename = f"{strain}/bakta/{strain}.gbff"

fna_filename = f"{strain}.fna"
faa_filename = f"{strain}.faa"
gff3_filename = f"{strain}.gff3"
ffn_filename = f"{strain}.ffn"

gbk_input_handle = open(gbk_filename, "r")
#gbk to fna
fna_output_handle = open(fna_filename, "w")
for seq_record in SeqIO.parse(gbk_input_handle, "genbank"):
#       print ("Dealing with GenBank record %s" % seq_record.id)
        fna_output_handle.write(">%s %s\n%s\n" % (
                seq_record.id,
                seq_record.description,
                seq_record.seq))
fna_output_handle.close()
gbk_input_handle.close()

#gbk to faa
gbk_input_handle = open(gbk_filename, "r")
faa_output_handle = open(faa_filename, "w")
for seq_record in SeqIO.parse(gbk_input_handle, "genbank") :
#       print ("Dealing with GenBank record %s" % seq_record.id)
        for seq_feature in seq_record.features:
                if seq_feature.type=="CDS" and 'translation' in seq_feature.qualifiers and 'pseudogene' not in seq_feature.qualifiers:
                        assert len(seq_feature.qualifiers['translation'])==1
                        faa_output_handle.write(">%s %s\n%s\n" % (
                                seq_feature.qualifiers['locus_tag'][0],
                                seq_feature.qualifiers['product'][0],
                                seq_feature.qualifiers['translation'][0]))

faa_output_handle.close()
gbk_input_handle.close()

#gbk to gff3
gbk_input_handle = open(gbk_filename, "r")
gff3_output_handle = open(gff3_filename, "w")
GFF.write(SeqIO.parse(gbk_input_handle, "genbank"), gff3_output_handle)
gff3_output_handle.close()
gbk_input_handle.close()

#gbk to ffn
gbk_input_handle = open(gbk_filename, "r")
ffn_input_handle = open(ffn_filename, 'w')
for seq_record in SeqIO.parse(gbk_input_handle, "genbank") :
#        print ("Dealing with GenBank record %s" % seq_record.id)
    for seq_feature in seq_record.features:
        if seq_feature.type=="CDS":
            product = seq_feature.qualifiers.get("product")[0]
            gene = seq_feature.qualifiers.get("gene")
            if gene is not None:
                gene = seq_feature.qualifiers.get("gene")[0]
            else:
                gene = ""
            locus = seq_feature.qualifiers.get("locus_tag")
            if locus is not None:
                 locus = seq_feature.qualifiers.get("locus_tag")[0]
            else:
                 locus = ""
            nucleotide = seq_feature.extract(seq_record).seq
            locus = seq_feature.qualifiers.get("locus_tag")[0]
            ffn_input_handle.write(">%s %s %s\n%s\n" % (
                                locus,
                                product,
                                gene,
                                nucleotide))
ffn_input_handle.close()
gbk_input_handle.close()
