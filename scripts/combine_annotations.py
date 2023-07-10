#!/bin/python

import sys
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
import os
import os.path
strain = sys.argv[1]

macsyfinder_path = f"./{strain}/macsyfinder.tsv.table"
defensefinder_path = f"./{strain}/{strain}_defensefinder.tsv.table"
hmmdb_path = f"./{strain}/{strain}_borders.table"
dbscan_path = f"./{strain}/prophage.table"
antismash_path = f"./{strain}/antismash_locustags.table"
tiger_path = f"./{strain}/{strain}_TIGER2_final.table.out"
integron_path = f"./{strain}/integron.table"
integron_gene_path = f"./{strain}/integron_gene.table"
gapmind_path = f"./{strain}/GapMind/combined_GapMind_results.tab"

outfile_path = f"./{strain}/{strain}_final.gbk"

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

hmmdb = {}
if os.path.isfile(hmmdb_path) == True:
    with open(f"./{strain}/{strain}_borders.table", 'r') as hmmtable_file:
        for l in hmmtable_file:
            replicon,start,end,value,annot,strand = l.strip().split('\t')
            if replicon in hmmdb:
                hmmdb[replicon].append((replicon,start,end,value,annot,strand))
            else:
                hmmdb[replicon]=[(replicon,start,end,value,annot,strand)]

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
    with open (f"./{strain}/antismash_locustags.table", 'r') as antismash_file:
        for line in antismash_file:
            (locus,annot) = line.split()
            antismash_dict[locus] = annot

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
                        

new_records = []
for record in SeqIO.parse(f"./{strain}/{strain}.gbff","gb"):
    for feature in record.features:
        locus_tags = feature.qualifiers.get("locus_tag")
        if locus_tags is not None and feature.type == "CDS":
            if locus_tags[0] in locus_dict:
                if "note" in locus_dict:
                    feature.qualifiers["note"].append("MacSyFinder: " + locus_dict[locus_tags[0]])
                else:
                    feature.qualifiers["note"] = "MacSyFinder: " + locus_dict[locus_tags[0]]

            if locus_tags[0] in defense_dict:
                if "note" in defense_dict:
                    feature.qualifiers["note"].append("DefenseFinder: " + defense_dict[locus_tags[0]])
                else: 
                    feature.qualifiers["note"] = "DefenseFinder: " + defense_dict[locus_tags[0]]

            if locus_tags[0] in antismash_dict:
                if "note" in antismash_dict:
                    feature.qualifiers["note"].append("antiSMASH: " + antismash_dict[locus_tags[0]])
                else:
                    feature.qualifiers["note"] = "antiSMASH: " + antismash_dict[locus_tags[0]]

            if locus_tags[0] in integron_gene:
                if "note" in integron_gene:
                    feature.qualifiers["note"].append("IntegronFinder: " + integron_gene[locus_tags[0]])
                else:
                    feature.qualifiers["note"] = "IntegronFinder: " + integron_gene[locus_tags[0]]

            if locus_tags[0] in gapmind_dict:
                if "note" in integron_gene:
                    feature.qualifiers["note"].append(gapmind_dict[locus_tags[0]])
                else:
                    feature.qualifiers["note"] = gapmind_dict[locus_tags[0]]


    if record.id in hmmdb:
        for current_border in hmmdb[record.id]:
            new_feat = SeqFeature(FeatureLocation(int(current_border[1]), int(current_border[2])),type="misc_feature", qualifiers= {"note": [current_border[4]], "inference" : "fuzznuc pattern or hmm"}, strand = int(current_border[5]))
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
    if record.id in dbscan:
        for phage in dbscan[record.id]:
            newfeat = SeqFeature(FeatureLocation(int(phage[0]),int(phage[1])), type="mobile_element",qualifiers={"mobile_element_type": "phage", "note": phage[2] + " " + phage[3] , "inference" : "DBSCAN-SWA"})
            record.features.append(newfeat)
    record.features.sort(key=lambda x: x.location.start,reverse=False)  
    new_records.append(record)
out_handle = open(outfile_path, 'w')
SeqIO.write(new_records,out_handle, "genbank")
