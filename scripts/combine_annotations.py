#!/bin/python

import sys
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
#from BCBio import GFF
import os
import os.path
strain = sys.argv[1]

#Check if input exists
macsyfinder_path = f"./{strain}/macsyfinder.tsv.table"
defensefinder_path = f"./{strain}/{strain}_defensefinder.tsv.table"
hmmdb_path = f"./{strain}/{strain}_uniq_borders.table"
dbscan_path = f"./{strain}/prophage.table"
antismash_path = f"./{strain}/{strain}_antismash.table.beav.subset"
tiger_path = f"./{strain}/{strain}_TIGER2_final.table.out"
integron_path = f"./{strain}/integron.table"
integron_gene_path = f"./{strain}/integron_gene.table"
gapmind_path = f"./{strain}/GapMind/combined_GapMind_results.tab"
operon_path = f"./{strain}/operon-mapper_results/list_of_operons.table"

outfile_path = f"./{strain}/{strain}_final.gbk"
faa_path = f"./{strain}/{strain}_final.faa"
ffn_path = f"./{strain}/{strain}_final.ffn"

#read in annotation tables as dictionary
locus_dict = {}
if os.path.isfile(macsyfinder_path) == True:
    with open (f"./{strain}/macsyfinder.tsv.table", "r") as f:
        for line in f:
            (key,values) = line.split()
            locus_dict[key] = values 


defense_dict = {}
if os.path.isfile(defensefinder_path) == True:
    with open (f"./{strain}/{strain}_defensefinder.tsv.table", 'r') as k:
        for line in k:
            (key,values) = line.split()
            defense_dict[key] = values

gapmind_dict = {}
if os.path.isfile(gapmind_path) == True:
    with open (gapmind_path, 'r') as gp:
        for line in gp:
            (key,values) = line.rstrip().split("\t")
            gapmind_dict[key] = values

operon_dict = {}
if os.path.isfile(operon_path) == True:
    with open (operon_path, 'r') as opp:
        for line in opp:
            (key,values) = line.rstrip().split("\t")
            operon_dict[key] = values

hmmdb = {}
if os.path.isfile(hmmdb_path) == True:
    with open(f"./{strain}/{strain}_uniq_borders.table", 'r') as hmmtable_file:
        for l in hmmtable_file:
            replicon,start,end,annot,strand = l.strip('\n').split('\t')
            if replicon in hmmdb:
                hmmdb[replicon].append((replicon,start,end,annot,strand))
            else:
                hmmdb[replicon]=[(replicon,start,end,annot,strand)]

dbscan = {}
if os.path.isfile(dbscan_path) == True:     
    with open (f"./{strain}/prophage.table", 'r') as prophage_file:
        for l in prophage_file:
            replicon,start,end,annot,category = l.strip().split('\t')
            if replicon in dbscan:
                dbscan[replicon].append((start,end,annot,category))
            else:
                dbscan[replicon] = [(start,end,annot,category)] 

antismash_dict = {}
if os.path.isfile(antismash_path) == True:
    with open (f"./{strain}/{strain}_antismash.table.beav.subset", 'r') as antismash_file:
        for line in antismash_file:
            locus,cluster,function,nrps,domain = line.rstrip('\n').split('\t')
            if locus in antismash_dict:
                antismash_dict[locus].append((cluster,function,nrps,domain))
            else:
                antismash_dict[locus] = [(cluster,function,nrps,domain)]


tiger_dict = {}
if os.path.isfile(tiger_path) == True:          
    with open (f"./{strain}/{strain}_TIGER2_final.table.out", 'r') as tiger_file:
        for line in tiger_file:
            replicon,start,end,annot = line.strip().split('\t')
            if replicon in tiger_dict:
                tiger_dict[replicon].append((start,end,annot))
            else:
                tiger_dict[replicon] = [(start,end,annot)]
integron = {}
if os.path.isfile(integron_path) == True:                
    with open (f"./{strain}/integron.table", 'r') as integron_table:
        for l in integron_table:
            replicon,element,start,end,strand,complete = l.strip().split('\t')
            if replicon in integron:
                integron[replicon].append((element,start,end,strand,complete))
            else:
                integron[replicon] = [(element,start,end,strand,complete)]

integron_gene = {}
if os.path.isfile(integron_gene_path) == True:
    protein = "IntI"
    with open (f"./{strain}/integron_gene.table", 'r') as integron_gene_file:
        for locus in integron_gene_file:
            key = locus.strip()
            integron_gene[key] = protein       
#Parse genbank and match locus tags to annotation                        
new_records = []
for record in SeqIO.parse(f"./{strain}/{strain}.gbff","gb"):
    accession = record.annotations.get("accessions")
    for feature in record.features:
        locus_tags = feature.qualifiers.get("locus_tag")
        if locus_tags is not None and feature.type == "CDS":
            if locus_tags[0] in locus_dict:
                if "note" in feature.qualifiers:
                    feature.qualifiers["note"].append("MacSyFinder: " + locus_dict[locus_tags[0]])
                else:
                    feature.qualifiers["note"] = ["MacSyFinder: " + locus_dict[locus_tags[0]]]

            if locus_tags[0] in defense_dict:
                if "note" in feature.qualifiers:
                    feature.qualifiers["note"].append("DefenseFinder: " + defense_dict[locus_tags[0]])
                else: 
                    feature.qualifiers["note"] = ["DefenseFinder: " + defense_dict[locus_tags[0]]]
                    
            if locus_tags[0] in antismash_dict:
                cluster = antismash_dict[locus_tags[0]][0][0]
                gene_product = antismash_dict[locus_tags[0]][0][1]
                NRPS_PKS = antismash_dict[locus_tags[0]][0][2]
                domains = antismash_dict[locus_tags[0]][0][3]
                if cluster != "":
                    if  "note" in feature.qualifiers:
                        feature.qualifiers["note"].append("antiSMASH cluster: " +  cluster)
                    else:
                        feature.qualifiers["note"] = ["antiSMASH cluster: " + cluster]
                if gene_product != "":
                    if "note" in feature.qualifiers:
                        feature.qualifiers["note"].append("antiSMASH gene product: " +  gene_product)
                    else:
                        feature.qualifiers["note"] = ["antiSMASH gene product: " + gene_product]

                if NRPS_PKS != "":
                    if "note" in feature.qualifiers:
                        feature.qualifiers["note"].append("antiSMASH NRPS/PKS: " +  NRPS_PKS)
                    else:
                        feature.qualifiers["note"] = ["antiSMASH NRPS/PKS: " + NRPS_PKS]

                if domains != "":
                    if "note" in feature.qualifiers:
                        feature.qualifiers["note"].append("antiSMASH domain: " +  domains)
                    else:
                        feature.qualifiers["note"] = ["antiSMASH domain: " + domains]

            if locus_tags[0] in integron_gene:
                if "note" in feature.qualifiers:
                    feature.qualifiers["note"].append("IntegronFinder: " + integron_gene[locus_tags[0]])
                else:
                    feature.qualifiers["note"] = ["IntegronFinder: " + integron_gene[locus_tags[0]]]

            if locus_tags[0] in gapmind_dict:
                if "note" in feature.qualifiers:
                    feature.qualifiers["note"].append(gapmind_dict[locus_tags[0]])
                else:
                    feature.qualifiers["note"] = [gapmind_dict[locus_tags[0]]]
            
            if locus_tags[0] in operon_dict:
                feature.qualifiers["operon"] = operon_dict[locus_tags[0]]

#adding new features
    if record.id in hmmdb:
        for current_border in hmmdb[record.id]:
            new_feat = SeqFeature((FeatureLocation(int(current_border[1]), int(current_border[2]), strand = int(current_border[4]))),type="misc_feature", qualifiers= {"note": [current_border[3]], "inference" : "BEAV"})
            record.features.append(new_feat)
        
    if record.id in integron:
        for integrons in integron[record.id]:
            integron_newfeat = SeqFeature(FeatureLocation(int(integrons[1]), int(integrons[2])), type="misc_feature", qualifiers= {"note": [integrons[0]], "inference" : "MacSyFinder TXSS models"}, strand = int(integrons[3]))
            integron_newfeat.qualifiers["note"].append("Integron Status: " +integrons[4])
            record.features.append(integron_newfeat)

    if record.id in tiger_dict:
        for ice in tiger_dict[record.id]:
            tiger_newfeat = SeqFeature(FeatureLocation(int(ice[0]),int(ice[1])), type="mobile_element", qualifiers={"mobile_element_type": "integrative element", "note": ice[2] , "inference" : "TIGER2"})
            record.features.append(tiger_newfeat)
    for number in accession:
        if number in dbscan:
            for phage in dbscan[number]:
                newfeat = SeqFeature(FeatureLocation(int(phage[0]),int(phage[1])), type="mobile_element", qualifiers={"mobile_element_type": "phage", "note": phage[2] + " " + phage[3] , "inference" : "DBSCAN-SWA"})
                record.features.append(newfeat)
    record.features.sort(key=lambda x: x.location.start,reverse=False)  
    new_records.append(record)
output_gbk_handle = open(outfile_path, 'w')
SeqIO.write(new_records,output_gbk_handle, "genbank")
output_gbk_handle.close()

#writing output into different formats
output_faa_handle = open(faa_path, 'w')
for record in new_records:
        for feature in record.features:
                if feature.type=="CDS":
                        product = feature.qualifiers.get("product")[0]
                        gene = feature.qualifiers.get("gene")
                        if gene is not None:
                                gene = feature.qualifiers.get("gene")[0]
                        else:
                                gene = ""

                        translation = feature.qualifiers.get("translation")
                        if translation is not None:
                                translation = feature.qualifiers.get("translation")[0]
                        else:
                                translation = ""
                        output_faa_handle.write(">%s %s %s\n%s\n" % (
                                feature.qualifiers['locus_tag'][0],
                                product,
                                gene,
                                translation))
output_faa_handle.close()


output_ffn_handle = open(ffn_path, 'w')

for record in new_records:
        for feature in record.features:
                if feature.type=="CDS":
                        product = feature.qualifiers.get("product")[0]
                        gene = feature.qualifiers.get("gene")
                        if gene is not None:
                                gene = feature.qualifiers.get("gene")[0]
                        else:
                                gene = ""

                        locus = feature.qualifiers.get("locus_tag")
                        if locus is not None:
                                locus = feature.qualifiers.get("locus_tag")[0]
                        else:
                                locus = ""
                        nucleotide = feature.extract(record).seq
                        locus = feature.qualifiers.get("locus_tag")[0]
                        output_ffn_handle.write(">%s %s %s\n%s\n" % (
                                locus,
                                product,
                                gene,
                                nucleotide))
output_ffn_handle.close()


