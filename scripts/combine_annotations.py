#!/bin/python

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation


locus_dict = {}
with open ('./AS1D4/macsyfinder.tsv.table') as f:
	for line in f:
		(key,values) = line.split()
		locus_dict[key] = values 

defense_dict = {}
with open ('./AS1D4/AS1D4_defensefinder.tsv.table') as k:
	for line in k:
		(key,values) = line.split()
		defense_dict[key] = values
#table lines look like: pTi  T-DNA_left_border 1967645 1967615
hmmdb = {}
with open('./AS1D4/borders.table', 'r') as hmmtable_file:
	for l in hmmtable_file:
		replicon,annot,start,end = l.strip().split('\t')
		if replicon in hmmdb:
			hmmdb[replicon].append((annot,start,end))
		else:
			hmmdb[replicon]=[(annot,start,end)]
		
dbscan = {}
with open ('./AS1D4/prophage.table', 'r') as prophage_file:
	for l in prophage_file:
		replicon,start,end,annot,category = l.strip().split('\t')
		if replicon in dbscan:
			dbscan[replicon].append((start,end,annot,category))
		else:
			dbscan[replicon] = [(start,end,annot,category)]	

antismash_dict = {}
with open ('./AS1D4/antismash_locustags.table') as antismash_file:
	for line in antismash_file:
		(locus,annot) = line.split()
		antismash_dict[locus] = annot
			
tiger_dict = {}
with open ('./AS1D4/AS1D4_TIGER2_final.table.out') as tiger_file:
	for line in tiger_file:
		replicon,start,end,annot = line.strip().split('\t')
		if replicon in tiger_dict:
			tiger_dict[replicon].append((start,end,annot))
		else:
			tiger_dict[replicon] = [(start,end,annot)]
#integron = {}
#with open ('./AS1D4/integron.table', 'r') as integron_table:
#	for l in integron_table:
#		replicon,element,start,end,strand,complete = l.strip().split('\t')
#		if replicon in integron:
#			integron[replicon].append((element,start,end,strand,complete))
#		else:
#			integron[replicon] = [(element,start,end,strand,complete)]
#integron_gene = {}
#protein = "IntI"
#with open ('./AS1D4/integron_gene.table') as integron_gene_file:
#	for locus in integron_gene_file:
#		key = locus.strip()
#		integron_gene[key] = protein
		
						

new_records = []
for record in SeqIO.parse("./AS1D4/AS1D4.gbff","gb"):
	for feature in record.features:
		locus_tags = feature.qualifiers.get("locus_tag")
		if locus_tags is not None and feature.type == "CDS":
				#print(locus_tags)
			if locus_tags[0] in locus_dict:
				if "note" in locus_dict:
					feature.qualifiers["note"].append("Macsyfinder: " + locus_dict[locus_tags[0]])
				else:
					feature.qualifiers["note"] = "Macsyfinder: " + locus_dict[locus_tags[0]]

			if locus_tags[0] in defense_dict:
				if "note" in defense_dict:
					feature.qualifiers["note"].append("Defensefinder: " + defense_dict[locus_tags[0]])
				else: 
					feature.qualifiers["note"] = "Defensefinder: " + defense_dict[locus_tags[0]]

			if locus_tags[0] in antismash_dict:
				if "note" in antismash_dict:
					feature.qualifiers["note"].append("Antismash: " + antismash_dict[locus_tags[0]])
				else:
					feature.qualifiers["note"] = "Antismash: " + antismash_dict[locus_tags[0]]

		#	if locus_tags[0] in integron_gene:
				if "note" in feature.qualifiers:
		#			feature.qualifiers["note"].append("Integronfinder: " + integron_gene[locus_tags[0]])
				else: 
					feature.qualifiers["note"] = "Integronfinder: " + integron_gene[locus_tags[0]]	
	if record.id in hmmdb:
		for current_border in hmmdb[record.id]:
			start = int(current_border[1])
			end = int(current_border[2])
			feature_location = FeatureLocation(start,end)
			feat_type = "misc_feature"
			#feat_id = current_border[0]
			new_feat = SeqFeature(feature_location,type=feat_type,)
			new_feat.qualifiers["note"] = current_border[0]
			record.features.append(new_feat)
	
#	if record.id in integron:
#		for integrons in integron[record.id]:
#			integron_newfeat = SeqFeature(FeatureLocation(int(integrons[1]), int(integrons[2])), type="misc_feature", qualifiers= {"note": ["Integronfinder: " +integrons[0]]}, strand = int(integrons[3]))
#			integron_newfeat.qualifiers["note"].append("Integron Status: " +integrons[4])

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
out_handle = open("test_output.gbk", 'w')
SeqIO.write(new_records,out_handle, "genbank")








