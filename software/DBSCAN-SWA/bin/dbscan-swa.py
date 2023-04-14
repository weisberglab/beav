import os,re,argparse
from os.path import dirname, abspath
import codecs
import numpy as np
from sklearn.cluster import DBSCAN
from Bio import Entrez, SeqIO
from collections import Counter
import threading
import time
import sys
import multiprocessing
import json
def inputinstructions():
	return """
Usage: DBSCAN-SWA [options]
--input <file name>        : Query phage file path: FASTA or Multi-Fasta or GenBank file
--output <folder name>     : Output folder in which results will be stored
--prefix <prefix>          : default: bac:
--evalue <x>               : maximal E-value of searching for homology virus proteins from viral UniProt TrEML database. default:1e-7
--min_protein_num <x>      : optional,the minimal number of proteins forming a phage cluster in DBSCAN, default:6
--protein_number <x>       : optional,the number of expanding proteins when finding prophage att sites, default:10
--add_annotation <options> : optional,1.PGPD: a phage genome and protein database,2.phage_path:specified phage genome to detect whether the phage infects the query bacteria
							3.none:no phage annotation. default:PGPD
--per <x>                  : Minimal % percentage of hit proteins on hit prophage region(default:30)
--idn <x>                  : Minimal % identity of hit region on hit prophage region by making blastn(default:70)
--cov <x>                  : Minimal % coverage of hit region on hit prophage region by making blastn(default:30)
--thread_num <x>           : the number of threads(default:10)
"""

def pro_distance(pro1, pro2):
	pro_list = list(pro1)+list(pro2)
	pro_list = sorted(list(map(int,pro_list)))
	two_pro_list = [[np.min(pro1),np.max(pro1)],[np.min(pro2),np.max(pro2)]]
	two_pro_list.sort(key=lambda x:x[0])
	if two_pro_list[0][1]>two_pro_list[1][0]:
		distance = 1
	else:
		#distance = abs(np.mean(pro1)-np.mean(pro2))
		distance = pro_list[2]-pro_list[1]
	#print(distance)
	return distance

def dbscan(phage_like_protein_list,outdir,prefix):
	cluster_list_file = os.path.join(outdir, prefix+"_dbscan_cluster.txt")
	cluster_protein_file = os.path.join(outdir, prefix+"_dbscan_cluster_phage_protein.txt")
	cluster_fw = codecs.open(cluster_list_file, "w", encoding="utf-8")
	phage_fw = codecs.open(cluster_protein_file, "w", encoding="utf-8")
	phage_fw.write('dbscan_region_id\tregion_location\tregion_key_proteins\tprotein_id\thit_uniprot_id\tidentity\tevalue\n')
	cluster_fw.write('dbscan_region_id\tregion_location\n')
	protein_name_list = [item[0] for item in phage_like_protein_list]
	protein_start_list1 = [int(item[1]) for item in phage_like_protein_list]
	protein_end_list1 = [int(item[2]) for item in phage_like_protein_list]
	#print(protein_start_list)
	protein_start_list = np.array(protein_start_list1)[:, np.newaxis]
	protein_end_list = np.array(protein_end_list1)[:, np.newaxis]
	# scaler = StandardScaler()
	# protein_start_list = scaler.fit_transform(protein_start_list)
	# protein_end_list = scaler.fit_transform(protein_end_list)
	protein_location_list = np.hstack([protein_start_list,protein_end_list]) 
	model = DBSCAN(eps=3000,min_samples=int(min_protein_num),metric=lambda a, b: pro_distance(a, b))
	prediction = model.fit_predict(protein_location_list)
	protein_index_list = np.argsort(protein_start_list1)
	#print(protein_location_list)
	regions_id = []
	region_details = {}
	save_dbscan_region = [list(map(int,phage_like_protein_list[protein_index_list[0]][1:3]))]
	#print(save_dbscan_region)
	counter = 0
	for index in protein_index_list:
		# print(index)
		if prediction[index] != -1:
			if prediction[index] not in regions_id:
				regions_id.append(prediction[index])
			protein_info = phage_like_protein_list[index]
			orf_name = protein_info[0]
			start = min(int(protein_info[1]),int(protein_info[2]))
			end = max(int(protein_info[1]),int(protein_info[2]))
			if start-save_dbscan_region[-1][1]<=3000:
				save_dbscan_region[-1][1] = end
				if counter not in region_details.keys():
					region_details.update({counter:[[],[]]})
			else:
				counter = counter+1
				save_dbscan_region.append([start,end])
				region_details.update({counter:[[],[]]})
			
			if len(orf_name.split('|'))==4:
				key_protein = orf_name.split('|')[-1].split(',')
			else:
				key_protein = ['NA']
			# if prediction[index] not in region_details.keys():
			# 	region_details.update({prediction[index]:[[],[]]})
			region_details[counter][0] = region_details[counter][0]+key_protein
			region_details[counter][0] = list(set(region_details[counter][0]))
			if 'NA' in region_details[counter][0]:
				region_details[counter][0].remove('NA')
			region_details[counter][1].append([start,end,orf_name]+protein_info[3:])

	for region_id in region_details.keys():
		items = region_details[region_id][1]
		key_proteins = region_details[region_id][0]
		sorted(items,key=(lambda x:int(x[0])))
		region_start = items[0][0]
		region_end = items[-1][1]
		region_id = str(region_id)
		cluster_fw.write(region_id+'\t'+str(region_start)+':'+str(region_end)+'\n')
		for item in items:
			start = item[0]
			end = item[1]
			orf_name = item[2]
			if len(key_proteins)>0:
				phage_fw.write(str(region_id) +'\t'+ str(region_start)+":"+str(region_end)+'\t'+','.join(key_proteins)+'\t'+'\t'.join(item[2:])+"\n")
			else:
				phage_fw.write(str(region_id) +'\t'+ str(region_start)+":"+str(region_end)+'\t'+'NA'+'\t'+'\t'.join(item[2:])+"\n")

def mkdir(dir):
	command = "mkdir -p "+dir
	os.system(command)

def prophage_window(blast_file):#n: a window of 60 proteins
	with open(blast_file) as f:
		contents = f.readlines()
	proteins_key = []
	proteins_local = []
	for protein in contents:
		protein = protein.strip().split('\t')
		pro_title = protein[0]
		pro_titles = pro_title.split('|')
		if len(pro_titles)==5:#flag=0 no keyword
			flag = 1
		else:
			flag = 0
		phage_title = protein[1]
		identity = protein[2]
		evalue = protein[-2]
		if pro_title not in proteins_local:
			proteins_local.append(pro_title)
			item = {'name':pro_title,'key':flag,'phage':[]}
			item['phage'].append([phage_title,identity,evalue])
			proteins_key.append(item)
		else:
			item['phage'].append([phage_title,identity,evalue])
	return proteins_key

def get_cluster_window(proteins_key,outdir,k):	
	num_region = 0
	num_proteins = len(proteins_key)	
	cluster_win = os.path.join(outdir,'cluster2_gb.txt')
	cluster_pro = os.path.join(outdir,'cluster2_pro.txt')
	f1 = open(cluster_win,'w')
	f2 = open(cluster_pro,'w')
	s = 0
	flag = 0
	while s<num_proteins-k:
		temp_pros = proteins_key[s:s+k]
		temp_protein = []
		keys = []
		counter =0 
		for pro in temp_pros:
			counter = counter + pro['key']
			if pro['key']==1:
				#print(pro['name'].split('|')[-1].split(','))
				temp_protein.append({'name':pro['name'],'phage':pro['phage']})
				keyword = pro['name'].split('|')[-1].split(',')
				for key in keyword:
					if key not in keys:
						keys.append(key)
		if counter>=6:
			start = proteins_key[s]
			end = proteins_key[s+k-1]
			start_loc = start['name'].split('|')[2]
			end_loc = end['name'].split('|')[2]
			start_location = re.findall("\d+\.?\d*",str(start_loc))[0]
			end_location = re.findall("\d+\.?\d*",str(end_loc))[1]
			f1.write(str(num_region)+'\t'+str(start_location)+'\t'+str(end_location)+'\t'+','.join(keys)+'\n')
			#detail.txt
			for key,temp in temp_protein.items():
				pros = temp['name'].split('|')
				pro_id = pros[1]
				for phages in temp['phage']:
					phage = phages[0].rstrip('|').split('|')
					phage_inf = phage[0]
					phage_pro_id = phage[-1]
					identity = phages[1]
					evalue = phages[2]
					f2.write(str(num_region)+'\t'+str(start_location)+'\t'+str(end_location)+'\t'+pro_id+'\t'+phage_pro_id+'\t'+phage_inf+'\t'+identity+'\t'+evalue+'\n')			
			num_region = num_region+1
			flag = 1
		else:
			flag = 0
		if flag==0:
			s = s+1
		else:
			s =s+k
	f1.close()
	f2.close()

def get_faa_protein(protein_file):#choose these proteins that are not predicted in cluster
	with open(protein_file) as f:
		contents = f.read()
	proteins = contents.strip().split('>')[1:]
	proteins_key = {}
	for protein in proteins:
		
		pro_title = protein.strip().split('\n')[0]
		sequence = protein.strip().split('\n')[1].strip()

		pro_titles = pro_title.split('|')
		proid = pro_titles[1]
		item = {proid:sequence}
		proteins_key.update(item)
	return proteins_key

def predict_prophage_swa(protein_file,blastp_file,outdir,k,prefix):#choose these proteins that are not predicted in cluster
	cluster_win = os.path.join(outdir,prefix+'_swa_cluster.txt')
	cluster_pro = os.path.join(outdir,prefix+'_swa_cluster_phage_protein.txt')
	f1 = open(cluster_win,'w')
	f2 = open(cluster_pro,'w')
	f1.write('swa_region_id\tregion_location\n')
	f2.write('swa_region_id\tregion_location\tregion_key_proteins\tprotein_id\thit_uniprot_id\tidentity\tcoverage\n')
	bac_protein_homo_dict = {}
	with open(blastp_file) as f:
		contents = f.readlines()
	for line in contents:
		line = line.strip().split('\t')
		bac_protein_id = line[0]
		hit_uniprot_id = line[1]
		bac_protein_homo_dict.update({bac_protein_id:[hit_uniprot_id,line[2],line[-2]]})
	
	with open(protein_file) as f:
		contents = f.read()
	proteins = contents.strip().split('>')[1:]
	proteins_key = []
	for protein in proteins:
		pro_title = protein.strip().split('\n')[0]
		pro_titles = pro_title.split('|')
		if len(pro_titles)==4:#flag=0 no keyword
			flag = 1
		else:
			flag = 0
		item = {'name':pro_title,'key':flag}
		proteins_key.append(item)
	num_region = 0
	num_proteins = len(proteins_key)	
	s = 0
	flag = 0
	while s<num_proteins-k:
		temp_pros = proteins_key[s:s+k]
		temp_protein = []
		keys = []
		counter =0
		bounder  = []
		for pro in temp_pros:
			if pro['key']==1:
				bounder.append(temp_pros.index(pro))
				# counter = counter + pro['key']
				#print(pro['name'].split('|')[-1].split(','))
				# temp_protein.append({'name':pro['name'],'sequence':pro['sequence']})
				keyword = pro['name'].split('|')[-1].split(',')
				if not ((len(keyword)==1) and (keyword[0] == 'tRNA')):
					counter = counter + pro['key']
				for key in keyword:
					if key not in keys:
						keys.append(key)
		if counter>=6:
			# start = proteins_key[s]
			# end = proteins_key[s+k-1]
			start = temp_pros[bounder[0]]
			end = temp_pros[bounder[-1]]
			start_loc = start['name'].split('|')[2]
			end_loc = end['name'].split('|')[2]
			start_location = start_loc.split('_')[0]
			end_location = end_loc.split('_')[1]
			#detail.txt
			if bounder[-1]==len(temp_pros)-1:
				temp_pros1 = temp_pros[bounder[0]:]
			else:
				temp_pros1 = temp_pros[bounder[0]:bounder[-1]+1]
			for temp in temp_pros1:
				pro_id = temp['name']
				if pro_id in bac_protein_homo_dict.keys():
					f2.write(str(num_region)+'\t'+str(start_location)+':'+str(end_location)+'\t'+','.join(keys)+'\t'+pro_id+'\t'+'\t'.join(bac_protein_homo_dict[pro_id])+'\n')
				else:
					f2.write(str(num_region)+'\t'+str(start_location)+':'+str(end_location)+'\t'+','.join(keys)+'\t'+pro_id+'\t'+'\t'.join(['NA']*3)+'\n')
			f1.write(str(num_region)+'\t'+str(start_location)+':'+str(end_location)+'\n')		
			num_region = num_region+1
			flag = 1
		else:
			flag = 0
		if flag==0:
			s = s+1
		else:
			s =s+k
	f1.close()
	f2.close()

def GetFaaSequenc(fileName,saveFaaPath,prefix,add_genome_id='no'):   #parse special protein from genbank files in phaster
	special_pros = ['capsid','head','plate','tail','coat','portal','holin','integrase','transposase','terminase','protease','lysis','bacteriocin','tRNA']
	records = SeqIO.parse(fileName, "gb")
	fileID = fileName.split('/')[-1]
	#fileID = fileName.strip('.gb')
	counter = 0
	# outFileName = os.path.join(saveFaaPath,fileID + '.faa')
	outFileName = os.path.join(saveFaaPath,prefix+'_protein.faa')
	outfiledef = os.path.join(saveFaaPath,prefix+'_protein_def')
	savefile = open(outFileName, 'w')
	savefile_protein = open(outfiledef,'w')
	for record in records:
		for feature in record.features:
			if (feature.type == 'CDS') or (feature.type == 'tRNA') or (feature.type == 'tmRNA'):                
				location = feature.location
				if str(location).find('+') != -1:
					direction = '+'
				elif str(location).find('-') != -1:
					direction = '-'
				locations = re.findall("\d+\.?\d*",str(location))
				min_start = locations[0]
				max_end = locations[1]
				for loc in location:
					if int(loc)<int(min_start):
						min_start = loc
					if int(loc)>int(max_end):
						max_end = loc
				location = str(min_start)+"_"+str(max_end)+'_'+direction
				counter = counter+1
				if 'product' in feature.qualifiers:
					product = feature.qualifiers['product'][0]
					key_pros = []
					for special in special_pros:
						if special in product:
							key_pros.append(special)      
					if 'protein_id' in feature.qualifiers:
						proteinId = feature.qualifiers['protein_id'][0]
					else:
						if 'locus_tag' in feature.qualifiers:
							proteinId = feature.qualifiers['locus_tag'][0]
						else:
							proteinId = 'unknown'
					if 'translation' in feature.qualifiers:
						translation = feature.qualifiers['translation'][0]
						if add_genome_id=='no':
							if len(key_pros)>0:
								savefile.write('>ref' + '|' + str(proteinId) + '|' + str(location) + '|' +','.join(key_pros)+'\n')
								savefile_protein.write(proteinId+'\t'+str(location)+'\t'+product+'\t'+','.join(key_pros)+'\n')
							else:
								savefile.write('>ref' + '|' + str(proteinId) + '|' + str(location) +'\n')
								savefile_protein.write(proteinId+'\t'+str(location)+'\t'+product+'\n')
						else:
							savefile.write('>'+add_genome_id + '|' + str(proteinId) + '|' + str(location) +'\n')
					# savefile.write(">"+fileID+ '_'+str(counter)+'\n')                   
						if translation[-1] == '\n':
							savefile.write(translation)
						else:
							savefile.write(translation + '\n')                  
	savefile.close()
	savefile_protein.close()

def diamond_blastp(file,outfile,database,format,evalue,diamond_thread_num=20):
	#num_threads = 20
	script = "diamond blastp -d "+database+" -q "+file+" -f "+str(format)+" -e "+str(evalue)+" -o "+outfile+" -p "+str(diamond_thread_num)+" --max-target-seqs 1"
	print(script)
	os.system(script)

def blastp(file,outfile,database,format,evalue):
	num_threads = 20
	script = "blastp -db "+database+" -query "+file+" -outfmt "+str(format)+" -evalue "+str(evalue)+" -out "+outfile+" -num_threads "+str(num_threads)+" -max_target_seqs 1"

	os.system(script)

def blastn(file,outfile,database,format,evalue):
	num_threads = 20
	format = 0
	script = "blastn -db "+database+" -query "+file+" -outfmt "+str(format)+" -evalue "+str(evalue)+" -out "+outfile+" -num_threads "+str(num_threads)+" -num_alignments 1 -word_size 11"
	os.system(script)

def blastn_file(file,outfile,subject_file,format,evalue):
	num_threads = 20
	format = 0
	script = "blastn -subject "+subject_file+" -query "+file+" -outfmt "+str(format)+" -evalue "+str(evalue)+" -out "+outfile+" -word_size 11"
	print(script)
	os.system(script)

def blastp_file(file,outfile,subject_file,format,evalue):
	num_threads = 20
	script = "blastp -subject "+subject_file+" -query "+file+" -outfmt "+str(format)+" -evalue "+str(evalue)+" -out "+outfile+" -max_target_seqs 1"
	os.system(script)

def get_length(file):#get the genome length
	if is_faa(file):
		with open(file) as f:
			contents = f.readlines()[1:]
		length = len('\n'.join(contents).strip())
	else:
		with open(file) as f:
			contents = f.readlines()[0]
		length = contents.split()[2]
	return str(length)

def get_protein_def(file):#save as dict
	with open(file) as f:
		contents = f.readlines()	
	proteins_dict = {}
	for line in contents:
		line = line.strip().split('\t')		
		if len(line)==4:
			proteins_dict.update({line[0]:{'prodef':line[2],'key':line[-1]}})
		else:
			proteins_dict.update({line[0]:{'prodef':line[2],'key':'NA'}})
	return proteins_dict

	
	# for line in contents:
	# 	line = line.strip().split('\t')		
	# 	if len(line)==3:
	# 		proteins_dict.update({line[0]:{'prodef':line[1],'key':line[-1]}})
	# 	else:
	# 		proteins_dict.update({line[0]:{'prodef':line[1],'key':'NA'}})

def get_uniprot_def1(file,protein_id):
	file = os.path.join(root_path,'db','profiles','uniprot_species.txt')
	with open(file) as f:
		contents = f.read()
	start = contents.find(protein_id)
	end = contents.find(protein_id,start)
	protein_species = contents[start:end].strip().split('\t')[-1]
	return protein_species

def get_uniprot_def(file):
	file = os.path.join(root_path,'db','profiles','uniprot_species.txt')
	with open(file) as f:
		contents = f.readlines()
	proteins_species = {}
	for line in contents:
		line = line.strip().split('\t')
		if len(line) > 0:
			proid = line[0]
			species = line[-1]
			if 'virus' in species:
				species = species[0:species.find('virus')]+'virus'
			if 'phage' in species:
				species = species[0:species.find('phage')]+'phage'
			proteins_species.update({proid:species})
	return proteins_species

def dict_to_file(dict,outfile):
	f = open(outfile,'w')
	json_str = json.dumps(dict,indent=4)
	f.write(json_str)
	f.close()

def GetFnaSequence(fileName,outFileName):
    # fileID = fileName.split('/')[-1].strip('.gb')
    # outFileName = os.path.join(os.path.dirname(fileName),fileID+'.fasta')
    handle = open(fileName)
    SeqIO.convert(handle, 'genbank',outFileName, 'fasta')

def get_faa_protein_fasta(pro_file):
	with open(pro_file) as f:
		contents = f.read().strip().split('>')
	pro_sequence_dict = {}
	for cds in contents[1:]:
		head = cds.strip().split('\n')[0]
		sequence = cds.strip().split('\n')[-1]
		cds_start = head.strip().split('_')[-3]
		cds_end = head.strip().split('_')[-2]
		flag = head.strip().split('_')[-1]
		if flag =='+':
			proid = cds_start+'..'+cds_end
		else:
			proid = "complement("+cds_start+'..'+cds_end+")"
		pro_sequence_dict.update({proid:sequence})
	return pro_sequence_dict

def parse_blastn_0(file):
	with open(file) as f:
		contents = f.read()
	hit_genome = contents.split('>')
	if len(hit_genome)>=2:
		hit_genome1 = hit_genome[1].split('\n')[0].strip()
		hit_id = hit_genome1.split()[0]
		hit_def = ' '.join(hit_genome1.split()[1:])
	else:
		hit_id = 'NA'
		hit_def = 'NA'
	blastn_homo = []
	if len(hit_genome)>=2:
		blastn_details = hit_genome[1].split('Score')[1:]
		for each_align in blastn_details:
			parameters = each_align.split('Strand')[0]
			evalue = get_num(parameters.split('\n')[0])[2]
			identity = get_num(parameters.split('\n')[1])[2]+'%'
			alignment_length = int(get_num(parameters.split('\n')[1])[0])
			blastn_homo.append({'length':alignment_length,'blastn_identity':identity,'blastn_evalue':evalue})
		blastn_homo.sort(key=lambda pair:pair['length'])
	return hit_id,hit_def,blastn_homo

def get_str_index(string,substring):#get all substring index,return list
	lis = [m.start() for m in re.finditer(substring,string)]
	return lis

def get_region_cluster(file,file_detail):
	with open(file) as f:
		contents = f.readlines()
	with open(file_detail) as f:
		contents_detail = f.read()

def get_num(string):
	strlist = re.findall("\d+\.?\d*",str(string))
	return strlist

def predict_prophage_region_orf(resufile,phaster_dir):
	if os.path.exists(resufile):			
		with open(resufile) as f:
			phaster_contents = f.readlines()
	else:
		phaster_contents = []
	if '' in phaster_contents:
		phaster_contents.remove('')
	if len(phaster_contents)>0:
		for line_phaster in phaster_contents:
			line_phaster = line_phaster.strip().split('\t')
			start = line_phaster[1]
			end = line_phaster[2]
			region_dir = os.path.join(phaster_dir,'prophage_'+start+'_'+end)
			fna_file = os.path.join(region_dir,'prophage_region.fna')
			faa_file = os.path.join(region_dir,'prophage_region.faa')
			if not os.path.exists(region_dir+'.faa'):
				predict_orf(fna_file,region_dir)

def get_inf1(file):
	global strain_id,strain_def,type
	with open(file) as f:
		contents = f.read().strip()
	if contents[0]=='>':
		type = 'fasta'
		record = SeqIO.read(file,"fasta")
		strain_id = record.id
		strain_def = record.description
	else:
		type = "genbank"
		record = SeqIO.read(file,"genbank")
		strain_id = record.id
		strain_def = record.description
	return strain_id,strain_def,type

def get_inf(file,outdir,prefix='bac'):
	with open(file) as f:
		contents = f.read().strip()
	if contents[0]=='>':
		type ='fasta'
		try:
			record = SeqIO.read(file,"fasta")
			genome_id = record.id
			strain_def = record.description
			strain_info_dict = {genome_id.split()[0].split('.')[0]:strain_def}
		except:
			strain_info_dict,strain_sequence_dict  = get_strain_info(file,outdir,prefix)
	else:
		type="genbank"
		strain_info_dict = split_genbank(file,outdir,prefix)		
		# record = SeqIO.read(file,"genbank")
		# genome_id = record.id
		# strain_def = record.description
		# strain_info_dict = {genome_id.split()[0].split('.')[0]:strain_def}
	
	strain_inf_dict_file = os.path.join(outdir,'strain_inf_dict')
	# np.save(strain_inf_dict_file,strain_info_dict)
	save_dict(strain_info_dict,strain_inf_dict_file)
	return strain_info_dict,type

def get_protein_to_position_genbank(pro_file,start,end,outfile,method,strain_id):
	f_result = open(outfile,'w')
	with open(pro_file) as f:
		contents = f.read().strip().split('>ref')
	all_pro_sequence_list = []
	counter = 0
	for line in contents[1:]:
		pro_start = line.split('\n')[0].split('|')[2].split('_')[-3]
		pro_end = line.split('\n')[0].split('|')[2].split('_')[-2]
		pro_id = line.split('\n')[0].split('|')[1]
		direction = line.split('\n')[0].split('|')[2].split('_')[-1]
		sequence = ''.join(line.split('\n')[1:]).strip()
		#print(sequence)
		if (min(int(start),int(end)) <= min(int(pro_start),int(pro_end))) and (max(int(start),int(end)) >= max(int(pro_start),int(pro_end))):
			counter = counter+1
			f_result.write('>'+strain_id+'|'+str(start)+':'+str(end)+'|'+line.split('\n')[0].strip().strip('|')+'|'+method+'\n')
			f_result.write(sequence.strip()+'\n')
			all_pro_sequence_list.append('>'+strain_id+'|'+str(start)+':'+str(end)+'|'+line.split('\n')[0].strip().strip('|')+'|'+method)
			all_pro_sequence_list.append(sequence)
	f_result.close()
	return all_pro_sequence_list,counter

def get_protein_to_position_fasta(pro_file,start,end,outfile,method):
	f_result = open(outfile,'w')
	with open(pro_file) as f:
		contents = f.read().strip().split('>')
	for line in contents[1:]:
		pro_start = line.split('\n')[0].split('_')[-3]
		pro_end = line.split('\n')[0].split('_')[-2]
		direction = line.split('\n')[0].split('_')[-1]
		sequence = '\n'.join(line.split('\n')[1:])
		#print(pro_start+'\t'+pro_end+'\t'+start+'\t'+end+'\n')
		if (min(int(start),int(end)) <= min(int(pro_start),int(pro_end))) and (max(int(start),int(end)) >= max(int(pro_start),int(pro_end))):
			f_result.write('>'+strain_id+'|'+str(start)+':'+str(end)+'|'+pro_start+'_'+pro_end+'_'+direction+'|'+method+'\n')
			f_result.write(sequence.strip()+'\n')
	f_result.close()

def get_att_6(out_blastn_file,save_file):
	f_result = open(save_file,'w')
	f_result.write('prophage_region\tattL\tattR\tatt_length\tbit_score\n')
	position = []
	bit_score = 0
	info = out_blastn_file.split('/')[-1]
	up_start = info.split('_')[0].split(':')[0]
	down_start = info.split('_')[1].split(':')[0]
	prophage_region = str(int(info.split('_')[0].split(':')[1])+1)+'-'+str(int(down_start)-1)
	if os.path.exists(out_blastn_file):
		if os.path.getsize(out_blastn_file)>0:
			with open(out_blastn_file) as f:
				contents = f.readlines()
			for line in contents:
				line = line.strip().split('\t')
				hit_length = line[3]
				query_start = line[6]
				query_end = line[7]
				subject_start = line[8]
				subject_end = line[9]
				mismatch = line[4]
				query_start = min(int(query_start),int(query_end))
				query_end = max(int(query_start),int(query_end))
				subject_start = min(int(subject_start),int(subject_end))
				subject_end = max(int(subject_start),int(subject_end))
				query_start = int(up_start)+query_start
				query_end = int(up_start)+query_end
				subject_start = int(down_start)+subject_start
				subject_end = int(down_start)+subject_end
				if int(hit_length)>=12 and int(mismatch)==0:
					position.append([prophage_region,str(query_start)+'-'+str(query_end),str(subject_start)+'-'+str(subject_end),hit_length,line[-1]])
	
	if len(position)>0:
		sorted(position,key=(lambda x:float(x[-1])),reverse=True)
		for each in position:
			f_result.write('\t'.join(each)+'\n')
	f_result.close()

def get_att_0(out_blastn_file,save_file):
	f_result = open(save_file,'w')
	f_result.write('prophage_region\tattL\tattR\tatt_length\tbit_score\n')
	position = []
	bit_score = 0
	info = out_blastn_file.split('/')[-1]
	up_start = info.split('_')[0].split(':')[0]
	down_start = info.split('_')[1].split(':')[0]
	#prophage_region = str(int(info.split('_')[0].split(':')[1])+1)+':'+str(int(down_start)-1)
	prophage_region = ':'.join(out_blastn_file.split('/')[-2].split('_'))
	if os.path.exists(out_blastn_file):
		if os.path.getsize(out_blastn_file)>0:
			with open(out_blastn_file) as f:
				contents = f.read().strip()
			# if len(contents)<2:
			# 	return 0
			#details = contents[1]
			up_pro_name = contents.split('Query=')[1].split('\n')[0].strip()
			down_pro_name = contents.split('Subject=')[1].split('\n')[0].strip()			
			hit_details = contents.split('Score =')
			for line in hit_details[1:]:
				bit_score = get_num(line.split('\n')[0])[0]
				evalue = get_num(line.split('\n')[0])[2]
				hit_length_t = get_num(line.split('\n')[1])[0]
				hit_length_f = get_num(line.split('\n')[1])[1]
				mismatch = int(hit_length_f)-int(hit_length_t)
				hit_sequence = '\n'.join(line.split('\n')[4:7]).strip()
				query_start = get_num(line.split('\n')[4])[0]
				query_end = get_num(line.split('\n')[4])[1]
				subject_start = get_num(line.split('\n')[6])[0]
				subject_end = get_num(line.split('\n')[6])[1]
				query_start = min(int(query_start),int(query_end))
				query_end = max(int(query_start),int(query_end))
				subject_start = min(int(subject_start),int(subject_end))
				subject_end = max(int(subject_start),int(subject_end))
				query_start = int(up_start)+query_start-1
				query_end = int(up_start)+query_end-1
				subject_start = int(down_start)+subject_start-1
				subject_end = int(down_start)+subject_end-1
				subject_end = subject_start+int(hit_length_t)-1
				#print(get_num(line.split('\n')[6]))
				# print(subject_end)
				if int(hit_length_t)>=12 and int(mismatch)==0:
					#print(subject_start,subject_end)
					position.append([prophage_region,str(query_start)+':'+str(query_end),str(subject_start)+':'+str(subject_end),hit_length_t,bit_score,hit_sequence])	
	if len(position)>0:
		sorted(position,key=(lambda x:float(x[-2])),reverse=True)
		for each in position:
			f_result.write('\t'.join(each[0:-1])+'\n'+each[-1]+'\n\n')
	f_result.close()

def getUpStreamProt_fasta(faa_file,region_start,num_prot,genome_id):
    up_stream_prot_temp = []
    hitFlag = 0
    titles = []
    index = 0
    with open(faa_file, 'r')as fin:
        content = fin.read()
        elems = content.split('>')
        if '' in elems:
            elems.remove('')
        prot_count = 0
        for elem in elems[::-1]:
            cur_genome_id = '_'.join(elem.split('\n')[0].strip().split('_')[0:-3])
            #if (genome_id in cur_genome_id) or (cur_genome_id in genome_id):
            if len(elem.split('\n')[0].split('_'))>3:
                cur_prot_start = elem.split('\n')[0].split('_')[-3]
                cur_prot_end = elem.split('\n')[0].split('_')[-2]
                if int(cur_prot_end) < int(region_start):
                    if prot_count == num_prot:
                        break
                    else:
                        prot_count = prot_count + 1
                        prot_size = str(int((int(cur_prot_end) - int(cur_prot_start) + 1) / 3))
                        prot_info = '%s-%s|%saa' % (cur_prot_start, cur_prot_end, prot_size)
                        titles.append(prot_size+'aa|up_'+str(index)+'|'+elem.strip().split('\n')[0])
                        index = index+1
                        up_stream_prot_temp.append('>'+elem.strip())
        titles = titles[::-1]
        if len(up_stream_prot_temp)>0:
            num_prior_prot = len(up_stream_prot_temp)
            if len(up_stream_prot_temp)==num_prot:
                up_stream_prot_list = up_stream_prot_temp
            else:
                titles = ['NA']*(num_prot-len(up_stream_prot_temp))+titles
                up_stream_prot_list = up_stream_prot_temp
        else:
            titles = ['NA']*num_prot
            up_stream_prot_list = ['NA']*num_prot
    return [up_stream_prot_list,titles]

def getDownStreamProt_fasta(faa_file,region_end,num_prot,genome_id):
    down_stream_prot_list = []
    titles = []
    index = 0
    with open(faa_file,'r')as fin:
        content = fin.read()
        elems = content.split('>')
        if '' in elems:
            elems.remove('')
        prot_count = 0
        for elem in elems:
            # cur_genome_id = elem.split('\n')[0].split('_')[-0].split('.')[0]
            cur_genome_id = '_'.join(elem.split('\n')[0].strip().split('_')[0:-3])
            # if (genome_id in cur_genome_id) or (cur_genome_id in genome_id):
            if len(elem.split('\n')[0].split('_'))>3:
                cur_prot_start = elem.split('\n')[0].split('_')[-3]
                cur_prot_end = elem.split('\n')[0].split('_')[-2]
                if int(cur_prot_start) > int(region_end):
                    if prot_count == num_prot:
                        break
                    else:
                        prot_count = prot_count + 1
                        prot_size = str(int((int(cur_prot_end) - int(cur_prot_start) + 1) / 3))
                        titles.append(prot_size+'aa|down_'+str(index)+'|'+elem.strip().split('\n')[0])
                        index = index+1
                        prot_info = '%s-%s|%saa'%(cur_prot_start,cur_prot_end,prot_size)
                        down_stream_prot_list.append('>'+elem.strip())
        if len(down_stream_prot_list)>0:
            num_prior_prot = len(down_stream_prot_list)
            if len(down_stream_prot_list)<num_prot:
                titles = titles+['NA']*(num_prot-len(down_stream_prot_list))
                down_stream_prot_list = down_stream_prot_list
        else:
            down_stream_prot_list = ['NA']*num_prot
            titles = ['NA']*num_prot
    return [down_stream_prot_list,titles]

def getUpStreamProt_gb(faa_file,region_start,num_prot,genome_id):
    up_stream_prot_temp = []
    hitFlag = 0
    titles = []
    with open(faa_file, 'r')as fin:
        content = fin.read()
        elems = content.split('>ref')
        if '' in elems:
            elems.remove('')
        prot_count = 0
        for elem in elems[::-1]:
            defLine = elem.split('\n')[0]
            cur_genome_id = defLine.split('|')[1].split('.')[0]
            #if (genome_id in cur_genome_id) or (cur_genome_id in genome_id):
            if len(defLine.split('|'))>=3:
                prot_location = defLine.split('|')[2]
                replace_list = ['>', '<']
                for replace_item in replace_list:
                    prot_location = prot_location.replace(replace_item, '')
                if 'join' in prot_location:
                    prot_location = prot_location.replace('join{', '')
                    prot_location = prot_location.replace('}', '')
                    prot_location_list = prot_location.split(',')
                    cur_prot_start = min(int(prot_location_list[0].split(':')[0].split('[')[1]),
                                         int(prot_location_list[0].split(':')[1].split(']')[0]),
                                         int(prot_location_list[-1].split(':')[0].split('[')[1]),
                                         int(prot_location_list[-1].split(':')[1].split(']')[0])
                                         )
                    cur_prot_end = max(int(prot_location_list[0].split(':')[0].split('[')[1]),
                                         int(prot_location_list[0].split(':')[1].split(']')[0]),
                                         int(prot_location_list[-1].split(':')[0].split('[')[1]),
                                         int(prot_location_list[-1].split(':')[1].split(']')[0]))
                else:
                    cur_prot_start = min(int(get_num(prot_location)[0]),
                    int(get_num(prot_location)[-1]))
                    cur_prot_end = max(int(get_num(prot_location)[0]),
                    int(get_num(prot_location)[-1]))
                if int(cur_prot_end) < int(region_start):
                    if prot_count == num_prot:
                        break
                    else:
                        prot_count = prot_count + 1
                        prot_size = str(int((int(cur_prot_end) - int(cur_prot_start) + 1) / 3))
                        prot_info = '%s-%s|%saa' % (cur_prot_start, cur_prot_end, prot_size)
                        titles.append('ref'+elem.strip().split('\n')[0])
                        up_stream_prot_temp.append('>ref'+elem.strip())
        titles = titles[::-1]
        if len(up_stream_prot_temp)>0:
            num_prior_prot = len(up_stream_prot_temp)
            if len(up_stream_prot_temp)==num_prot:
                up_stream_prot_list = up_stream_prot_temp
            else:
                titles = ['NA']*(num_prot-len(up_stream_prot_temp))+titles
                up_stream_prot_list = up_stream_prot_temp
        else:
            titles = ['NA']*num_prot
            up_stream_prot_list = ['NA']*num_prot
    return [up_stream_prot_list,titles]

def getDownStreamProt_gb(faa_file,region_end,num_prot,genome_id):
    down_stream_prot_list = []
    titles = []
    with open(faa_file,'r')as fin:
        content = fin.read()
        elems = content.split('>ref')
        if '' in elems:
            elems.remove('')
        prot_count = 0
        for elem in elems:
            # cur_genome_id = elem.split('\n')[0].split('_')[-0].split('.')[0]
            defLine = elem.split('\n')[0]
            #cur_genome_id = defLine.split('|')[1].split('.')[0]
            # if (genome_id in cur_genome_id) or (cur_genome_id in genome_id):
            if len(defLine.split('|'))>=3:
                prot_location = defLine.split('|')[2]
                replace_list = ['>', '<']
                for replace_item in replace_list:
                    prot_location = prot_location.replace(replace_item, '')
                if 'join' in prot_location:
                    prot_location = prot_location.replace('join{', '')
                    prot_location = prot_location.replace('}', '')
                    prot_location_list = prot_location.split(',')
                    cur_prot_start = min(int(prot_location_list[0].split(':')[0].split('[')[1]),
                                         int(prot_location_list[0].split(':')[1].split(']')[0]),
                                         int(prot_location_list[-1].split(':')[0].split('[')[1]),
                                         int(prot_location_list[-1].split(':')[1].split(']')[0])
                                         )
                    cur_prot_end = max(int(prot_location_list[0].split(':')[0].split('[')[1]),
                                         int(prot_location_list[0].split(':')[1].split(']')[0]),
                                         int(prot_location_list[-1].split(':')[0].split('[')[1]),
                                         int(prot_location_list[-1].split(':')[1].split(']')[0]))
                else:
                    cur_prot_start = min(int(get_num(prot_location)[0]),
                    int(get_num(prot_location)[-1]))
                    cur_prot_end = max(int(get_num(prot_location)[0]),
                    int(get_num(prot_location)[-1]))
                if int(cur_prot_start) > int(region_end):
                    if prot_count == num_prot:
                        break
                    else:
                        prot_count = prot_count + 1
                        prot_size = str(int((int(cur_prot_end) - int(cur_prot_start) + 1) / 3))
                        titles.append('ref'+elem.strip().split('\n')[0])
                        prot_info = '%s-%s|%saa'%(cur_prot_start,cur_prot_end,prot_size)
                        down_stream_prot_list.append('>ref'+elem.strip())
        if len(down_stream_prot_list)>0:
            num_prior_prot = len(down_stream_prot_list)
            if len(down_stream_prot_list)<num_prot:
                titles = titles+['NA']*(num_prot-len(down_stream_prot_list))
                down_stream_prot_list = down_stream_prot_list
        else:
            down_stream_prot_list = ['NA']*num_prot
            titles = ['NA']*num_prot
    return [down_stream_prot_list,titles]

def getProt_region_gb(faa_file,region_start,region_end):
	prot_list = []
	titles = []
	with open(faa_file,'r')as fin:
		content = fin.read()
		elems = content.split('>ref')
		if '' in elems:
			elems.remove('')
		prot_count = 0
		for elem in elems:           
			cur_prot_start = elem.split('\n')[0].split('|')[2].split('_')[0]
			cur_prot_end = elem.split('\n')[0].split('|')[2].split('_')[1]
			if (int(cur_prot_start) >= int(region_start)) and (int(cur_prot_end) <= int(region_end)):
				prot_count = prot_count + 1
				prot_size = str(int((int(cur_prot_end) - int(cur_prot_start) + 1) / 3))
				titles.append('ref'+elem.strip().split('\n')[0])
				index = index+1
				prot_info = '%s-%s|%saa'%(cur_prot_start,cur_prot_end,prot_size)
				prot_list.append('>ref'+elem.strip())
	return [prot_list,titles]

def get_acc(filename):
    with open(filename) as f:
        acc = f.readlines()[0].split()
    if is_faa(filename):
        acc = acc[0].strip('>').split('.')[0]
    else:
        acc = acc[1].strip()
    return acc

def short_blastn(file1,file2,outfile):
	command = "blastn -query "+file1+" -subject "+file2+" -out "+outfile+" -task blastn-short -evalue 1000"
	os.system(command)

def identify_att_before(usrfile,resufile,outdir,type,faa_file):
	genome_id = strain_id
	if os.path.exists(resufile):
		with open(resufile) as f:
			contents = f.readlines()
		if len(contents)==1:
			return 0
		with open(usrfile) as f:
			sequence = f.readlines()[1:]
		bac_sequence = ""
		for line in sequence:
			bac_sequence = bac_sequence+line.strip()
		for region in contents[1:]:
			head_summary = region.strip().split('\t')
			region_start = head_summary[3]
			region_end = head_summary[4]
			key_protein = head_summary[-1]
			if 'integrase' not in key_protein:
				continue
			up_start = 0
			up_down = 0
			down_start = 0
			down_end = 0
			#just get the region up and down protein
			if type=='fasta':
				up_pro_sequence,up_pro_name = getUpStreamProt_fasta(faa_file,region_start,1,genome_id)
				down_pro_sequence,down_pro_name = getDownStreamProt_fasta(faa_file,region_end,1,genome_id)
				if up_pro_name[-1]!='NA':
					up_start = int(up_pro_name[-1].split('_')[-3])
				else:
					up_start = 1
				if down_pro_name[0]!='NA':
					down_end = int(down_pro_name[0].split('_')[-2])
				else:
					down_end = len(bac_sequence)
				down_start = int(region_end)+1
				up_end = int(region_start)-1
			else:
				up_pro_sequence,up_pro_name = getUpStreamProt_gb(faa_file,region_start,1,genome_id)
				down_pro_sequence,down_pro_name = getDownStreamProt_gb(faa_file,region_end,1,genome_id)
				if up_pro_name[-1]!='NA':
					up_start = int(get_num(up_pro_name[-1].split('|')[2])[0])
				else:
					up_start = 1
				if down_pro_name[0]!='NA':
					down_end = int(get_num(down_pro_name[0].split('|')[2])[-1])
				else:
					down_end = len(bac_sequence)
				down_start = int(region_end)+1
				up_end = int(region_start)-1
			upstream_sequence = bac_sequence[int(up_start)-1:int(up_end)]
			downstream_sequence = bac_sequence[int(down_start)-1:int(down_end)]
			if (upstream_sequence == '' or downstream_sequence == ''):
				continue
			region_dir = os.path.join(outdir,region_start+'_'+region_end)
			mkdir(region_dir)
			down_file = os.path.join(region_dir,'down_sequence')
			up_file = os.path.join(region_dir,'up_sequence')
			with open(up_file,'w') as f:
				f.write('>'+region_start+':'+region_end+'|'+str(up_start)+'_'+str(up_end)+' '+up_pro_name[-1]+'\n'+upstream_sequence)			
			with open(down_file,'w') as f:
				f.write('>'+region_start+':'+region_end+'|'+str(down_start)+'_'+str(down_end)+' '+down_pro_name[0]+'\n'+downstream_sequence)			
			out_blastn_file = os.path.join(region_dir,str(up_start)+':'+str(up_end)+'_'+str(down_start)+':'+str(down_end)+'_blastn_result')
			short_blastn(up_file,down_file,out_blastn_file)
			out_att_file = os.path.join(region_dir,'att_info.txt')
			get_att_0(out_blastn_file,out_att_file)

def identify_att(strain_id,bac_fna_file,bac_faa_file,prophage_region_detail_file,outdir,prefix):
	if os.path.exists(prophage_region_detail_file):
		with open(prophage_region_detail_file) as f:
			contents = f.read().split('>prophage')
		if len(contents)<2:
			return 0
		with open(bac_fna_file) as f:
			sequence = f.readlines()[1:]
		bac_sequence = ""
		for line in sequence:
			bac_sequence = bac_sequence+line.strip()
		for region in contents[1:]:
			head_summary = region.strip().split('\n')[1].split('\t')
			
			region_start = head_summary[3]
			region_end = head_summary[4]
			key_protein = head_summary[-1]
			region_proteins = region.strip().split('\n')[2:]
			if 'integrase' in key_protein:
				#continue
				for region_protein in region_proteins[::-1]:
					c_pro_keyword = region_protein.strip().split('\t')[2]
					if 'integrase' in c_pro_keyword:
						key_integrase_pro_id = region_protein.strip().split('\t')[0]
						break

				key_integrase_pro_start = key_integrase_pro_id.split('|')[1].split('_')[0]
				key_integrase_pro_end = key_integrase_pro_id.split('|')[1].split('_')[1]
				up_start = 0
				up_down = 0
				down_start = 0
				down_end = 0
				#just get the region up and down 10 proteins to calculate att sites and use integrase as anchor			
				up_pro_sequence,up_pro_name = getUpStreamProt_gb(bac_faa_file,region_start,int(att_pro_num),strain_id)
				down_pro_sequence,down_pro_name = getDownStreamProt_gb(bac_faa_file,region_end,int(att_pro_num),strain_id)
				while 'NA' in up_pro_name:
					up_pro_name.remove('NA')
				up_start = 1
				if len(up_pro_name)>0:
					if up_pro_name[0]!='NA':
						up_start = int(up_pro_name[0].split('|')[2].split('_')[0])			
				while 'NA' in down_pro_name:
					down_pro_name.remove('NA')
				down_end = len(bac_sequence)
				if len(down_pro_name)>0:
					if down_pro_name[-1]!='NA':
						down_end = int(down_pro_name[-1].split('|')[2].split('_')[0])
				if int(key_integrase_pro_end)<len(bac_sequence):
					down_start = int(key_integrase_pro_end)+1
				else:
					down_start = len(bac_sequence)
				if int(key_integrase_pro_start)>0:
					up_end = int(key_integrase_pro_start)-1
				else:
					up_end = 1
				
				upstream_sequence = bac_sequence[int(up_start)-1:int(up_end)]
				downstream_sequence = bac_sequence[int(down_start)-1:int(down_end)]
				if (upstream_sequence == '' or downstream_sequence == ''):
					continue
				region_dir = os.path.join(outdir,region_start+'_'+region_end)
				mkdir(region_dir)
				down_file = os.path.join(region_dir,'down_sequence')
				up_file = os.path.join(region_dir,'up_sequence')
				with open(up_file,'w') as f:
					f.write('>'+region_start+':'+region_end+'|'+str(up_start)+'_'+str(up_end)+'\n'+upstream_sequence)			
				with open(down_file,'w') as f:
					f.write('>'+region_start+':'+region_end+'|'+str(down_start)+'_'+str(down_end)+'\n'+downstream_sequence)			
				out_blastn_file = os.path.join(region_dir,str(up_start)+':'+str(up_end)+'_'+str(down_start)+':'+str(down_end)+'_blastn_result')
				short_blastn(up_file,down_file,out_blastn_file)
				out_att_file = os.path.join(region_dir,'att_info.txt')
				get_att_0(out_blastn_file,out_att_file)

def decide_boundary(strain_id,bac_fna_file,bac_faa_file,prophage_region_detail_file,att_dir,save_prophage_file,save_prophage_summary_file,save_prophage_protein_file,save_prophage_nucl_file,bac_protein_blastp_phage_db_file,bac_protein_def_file,prefix):
	bac_protein_homo_dict = {}
	with open(bac_protein_blastp_phage_db_file) as f:
		contents = f.readlines()
	for line in contents:
		line = line.strip().split('\t')
		bac_protein_id = line[0]
		hit_uniprot_id = line[1]
		bac_protein_homo_dict.update({bac_protein_id:[hit_uniprot_id,line[2],line[-2]]})

	uniprot_proteins_species = get_uniprot_def("")
	bac_protein_def_dict = get_protein_def(bac_protein_def_file)
	f_save = open(save_prophage_file,'w')
	f_save_summary = open(save_prophage_summary_file,'w')
	f_save_summary.write('bacteria_id\tbacteria_def\tgenome_size\tprophage_start\tprophage_end\tkey_proteins\tbest_hit_species\tCDS_number\tattl_region\tattr_region\n')
	f_save_prophage_protein = open(save_prophage_protein_file,'w')
	f_save_prophage_nucl = open(save_prophage_nucl_file,'w')
	#boundary locating
	if os.path.exists(prophage_region_detail_file):
		with open(prophage_region_detail_file) as f:
			contents = f.read().split('>prophage')
		f_save.write(contents[0].strip()+'\n')
		with open(bac_fna_file) as f:
			sequence = f.readlines()[1:]
		bac_sequence = ""
		for line in sequence:
			bac_sequence = bac_sequence+line.strip()

		new_prophage_regions = []
		bac_info = contents[1].split('\n')[1].split('\t')[0:3]
		for region in contents[1:]:
			head_summary = region.strip().split('\n')[1].split('\t')
			
			region_start = head_summary[3]
			region_end = head_summary[4]
			#print(region_start+':'+region_end)
			key_protein = head_summary[-1]
			region_proteins = region.strip().split('\n')[2:]			
			region_dir = os.path.join(att_dir,region_start+'_'+region_end)
			out_att_file = os.path.join(region_dir,'att_info.txt')	
			write_flag = 0
			if os.path.exists(out_att_file):	
				with open(out_att_file) as f:
					att_infos = f.readlines()
				if len(att_infos)>1:
					write_flag = 1
					best_att = att_infos[1]
					prophage_start = best_att.split('\t')[0].split(':')[0]
					prophage_end = best_att.split('\t')[0].split(':')[1]
					attl_start = best_att.split('\t')[1].split(':')[0]
					attl_end = best_att.split('\t')[1].split(':')[1]
					attl_sequence = bac_sequence[int(attl_start)-1:int(attl_end)]
					attr_start = best_att.split('\t')[2].split(':')[0]
					attr_end = best_att.split('\t')[2].split(':')[1]
					attr_sequence = bac_sequence[int(attr_start)-1:int(attr_end)]
					
					if int(attl_start)<int(prophage_start):
						prophage_start = attl_start				
					if int(attr_end)>int(prophage_end):
						prophage_end = attr_end
					new_prophage_regions.append([int(region_start),int(region_end),[[attl_start+':'+attl_end,attr_start+':'+attr_end]]])
			#if write_flag==0:
			new_prophage_regions.append([int(region_start),int(region_end),[]])
		
		new_prophage_regions.sort(key=lambda x:x[0])
		#print(new_prophage_regions)
		new_merge_region = []
		for region in new_prophage_regions:
			region_start = region[0]
			region_end = region[1]
			if (len(new_merge_region)==0) or (region_start-new_merge_region[-1][1])>3000:
				new_merge_region.append(region)
			else:
				new_merge_region[-1][1] = max(region_end,new_merge_region[-1][1])
				if len(region[2])>0:
					new_merge_region[-1][2].append(region[2][0])
		#print(new_merge_region)

		for region_index,region in enumerate(new_merge_region):
			prophage_start = str(region[0])
			prophage_end = str(region[1])
			att_infos = region[2]
			if len(att_infos)>0:				 
				wide_attl = att_infos[0][0]
				wide_attr = att_infos[-1][1]
			else:
				wide_attl = 'NA'
				wide_attr = 'NA'
			c_save_prophage_protein_file = os.path.join(att_dir,prophage_start+'_'+prophage_end+'_prophage.faa')
			region_pro_sequence_list,prophage_pro_num = get_protein_to_position_genbank(bac_faa_file,prophage_start,prophage_end,c_save_prophage_protein_file,'DBSCAN-SWA',bac_info[0])
			
			key_words = [item.split('|')[-2] for item in region_pro_sequence_list[0::2] if len(item.split('|'))==6]
			key_words = ','.join(list(set(','.join(key_words).split(','))))
					
			all_hit_uniprot_species = [uniprot_proteins_species[bac_protein_homo_dict['ref|'+'|'.join(item.split('|')[2:-1])][0].split('|')[1]] for item in region_pro_sequence_list[0::2] if 'ref|'+'|'.join(item.split('|')[2:-1]) in bac_protein_homo_dict.keys()]
			while 'NA' in all_hit_uniprot_species:
				all_hit_uniprot_species.remove('NA')
			try:
				per = list(Counter(all_hit_uniprot_species).most_common(1)[0])[1]
				per = round(float(per)/float(len(all_hit_uniprot_species))*100,2)
				hit_mf = max(Counter(all_hit_uniprot_species),key=Counter(all_hit_uniprot_species).get)+"("+str(per)+"%)"
				if per<10:
					continue
			except:
				continue
			for index,region_pro in enumerate(region_pro_sequence_list[0::2]):
				f_save_prophage_protein.write(region_pro+'|'+str(prophage_pro_num)+'\n'+region_pro_sequence_list[2*index+1].strip()+'\n')
				f_save_prophage_protein.flush()
					
			f_save_prophage_nucl.write('>'+bac_info[0]+'|'+prophage_start+':'+prophage_end+'|DBSCAN-SWA\n'+bac_sequence[int(prophage_start)-1:int(prophage_end)]+'\n')
			f_save_prophage_nucl.flush()
			
			f_save.write('>prophage '+str(region_index+1)+'\n')
			
			f_save.write('\t'.join(bac_info).strip()+'\t'+'\t'.join([prophage_start,prophage_end,key_words,hit_mf,str(prophage_pro_num),wide_attl,wide_attr])+'\n')
			f_save.flush()
			f_save_summary.write('\t'.join(bac_info).strip()+'\t'+'\t'.join([prophage_start,prophage_end,key_words,hit_mf,str(prophage_pro_num),wide_attl,wide_attr])+'\n')	
			f_save_summary.flush()
			attl_flag = len(att_infos)
			attr_flag = len(att_infos)
			for index,region_protein in enumerate(region_pro_sequence_list[0::2]):
				flag = 0
				c_protein_start = region_protein.split('|')[3].split('_')[0]
				c_protein_end = region_protein.split('|')[3].split('_')[1]
				try:
					next_protein_end = region_pro_sequence_list[index+1].split('|')[3].split('_')[0]
				except:
					next_protein_end = c_protein_end
				if len(region_protein.split('|'))==6:
					key_protein = region_protein.split('|')[4]
				else:
					key_protein = 'NA'
				
				protein_id = region_protein.split('|')[2]
				protein_def = bac_protein_def_dict[protein_id]['prodef']
				if 'ref|'+'|'.join(region_protein.split('|')[2:-1]) in bac_protein_homo_dict.keys():
					hit_uniprot = bac_protein_homo_dict['ref|'+'|'.join(region_protein.split('|')[2:-1])]
					hit_uniprot_id = hit_uniprot[0].split('|')[1]
					hit_uniprot_def = uniprot_proteins_species[hit_uniprot_id]
					save_record = ['|'.join(region_protein.strip('>').split('|')[2:-1]),protein_def,key_protein,hit_uniprot_id,hit_uniprot_def]+hit_uniprot[1:]
				else:
					save_record = ['|'.join(region_protein.strip('>').split('|')[2:-1]),protein_def,key_protein]+['NA']*4
				#print(save_record)
				if wide_attl != 'NA':
					wide_attr_end = wide_attr.split(':')[1]
					for att in att_infos:
						attl_start = att[0].split(':')[0]
						attl_end = att[0].split(':')[1]
						attr_start = att[1].split(':')[0]
						attr_end = att[1].split(':')[1]
						attl_sequence = bac_sequence[int(attl_start)-1:int(attl_end)]
						attr_sequence = bac_sequence[int(attr_start)-1:int(attr_end)]
						if attl_flag>0:
							if int(attl_start)<=int(c_protein_start):
								f_save.write(attl_start+':'+attl_end+'\tattL\t'+attl_sequence+'\t'+'\t'.join(['NA']*4)+'\n')
								attl_flag = attl_flag-1
								f_save.write('\t'.join(save_record)+'\n')
								flag = 1
						if attr_flag > 0:	
							if (int(attr_end)<=int(c_protein_end)):
								f_save.write('\t'.join(save_record)+'\n')
								flag = 1
								f_save.write(attr_start+':'+attr_end+'\tattR\t'+attr_sequence+'\t'+'\t'.join(['NA']*4)+'\n')
								attr_flag = attr_flag-1
						
						if index==len(region_pro_sequence_list[0::2])-1:	
					
							if int(attr_end) > int(region_pro_sequence_list[0::2][-1].strip().split('|')[3].split('_')[1]):
								f_save.write('\t'.join(save_record)+'\n')
								flag = 1
								f_save.write(attr_start+':'+attr_end+'\tattR\t'+attr_sequence+'\t'+'\t'.join(['NA']*4)+'\n')
								attr_flag = 0
								break
				if flag ==0:
					f_save.write('\t'.join(save_record)+'\n')
			
	f_save_prophage_protein.close()
	f_save_prophage_nucl.close()
	f_save.close()

def get_strain_info(file,outdir,prefix):
	multi_dir = os.path.join(outdir,'results')
	mkdir(multi_dir)
	
	#fasta or multi-fasta
	strain_file = os.path.join(outdir,'strain_inf.txt')
	f_result = open(strain_file,'w')
	strain_inf_dict = {}
	strain_sequence_dict = {}
	with open(file) as f:
		contents = f.read().strip()
	counter  =0
	if '\n>' in contents:		
		for strain in contents.split('\n>'):
			strain_title = strain.split('\n')[0].strip()
			genome_id = strain_title.split()[0].strip('>')
			strain_def = strain_title.strip('>')
			strain_inf_dict.update({genome_id.split('.')[0]:strain_def})
			f_result.write(genome_id+'\t'+strain_def+'\n')
			f_result.flush()
			sequence = ''.join(strain.split('\n')[1:]).strip()
			strain_sequence_dict.update({genome_id.split('.')[0]:sequence})
			c_outdir = os.path.join(multi_dir,prefix+'_'+str(counter))
			c_prefix = prefix+'_'+str(counter)
			mkdir(c_outdir)
			c_inputfile_bac = os.path.join(c_outdir,prefix+'_'+str(counter)+'.fna')
			with open(c_inputfile_bac,'w') as f:
				f.write('>'+strain_title.strip()+'\n'+sequence.strip())
			counter = counter+1
	else:
		strain_title = contents.split('\n')[0].strip()
		genome_id = strain_title.split()[0].strip('>')
		strain_def = strain_title.strip('>')
		strain_inf_dict.update({genome_id.split('.')[0]:strain_def})
		f_result.write(genome_id+'\t'+strain_def+'\n')
		f_result.flush()
		sequence = ''.join(contents.split('\n')[1:]).strip()
		strain_sequence_dict.update({genome_id.split('.')[0]:sequence})
	f_result.close()
	return strain_inf_dict,strain_sequence_dict

def annotate_bacteria_fasta(fasta_file,outdir,prefix,kingdom):
	command = "prokka %s --outdir %s --prefix %s --kingdom %s --force"%(fasta_file,outdir,prefix,kingdom)
	print(command)
	os.system(command)

def get_bac_protein(inputfile_bac,protein_dir,strain_type,save_prefix,add_genome_id='no',kingdom='Bacteria'):
	if strain_type == 'fasta':
		#protein_dir = os.path.join(outdir,'protein_annotation')
		annotate_bacteria_fasta(inputfile_bac,protein_dir,save_prefix,kingdom) #prokka
		
		bac_gb_file = os.path.join(protein_dir,save_prefix+'.gbk')
		if not os.path.exists(bac_gb_file):
			bac_gb_file = os.path.join(protein_dir,save_prefix+'.gbf')
		if add_genome_id == 'no':
			GetFaaSequenc(bac_gb_file,protein_dir,save_prefix)
		else:
			GetFaaSequenc(bac_gb_file,protein_dir,save_prefix,add_genome_id)
	else:
		if add_genome_id == 'no':
			GetFaaSequenc(inputfile_bac,protein_dir,save_prefix) #get protein from gb file
		else:
			GetFaaSequenc(inputfile_bac,protein_dir,save_prefix,add_genome_id)
		bac_fna_file = os.path.join(protein_dir,save_prefix+'.fna')
		GetFnaSequence(inputfile_bac,bac_fna_file)

def get_phage_like_gene(bac_protein_file,phage_gene_annotation_dir,prefix,diamond_thread_num):
	blastp_file = os.path.join(phage_gene_annotation_dir,prefix+'_blastp_uniprot.txt')
	database = os.path.join(root_path,'db','database','uniprot.dmnd')
	diamond_blastp(bac_protein_file,blastp_file,database,6,str(blastp_evalue),diamond_thread_num)
	if os.path.exists(blastp_file):
		if os.path.getsize(blastp_file)>0:
			phage_like_protein_list = []
			with open(blastp_file) as f:
				contents = f.readlines()
			phage_like_protein_list = []
			for line in contents:
				line = line.strip().split('\t')
				bac_pro_id = line[0]
				pro_start = bac_pro_id.split('|')[2].split('_')[0]
				pro_end = bac_pro_id.split('|')[2].split('_')[1]
				hit_pro_id = line[1]
				evalue = line[-2]
				identity = line[2]
				phage_like_protein_list.append([bac_pro_id,pro_start,pro_end,hit_pro_id,identity,evalue])
			return phage_like_protein_list
	print('No phage or phage-like gene was identified in the query genome!')
	return []

def predict_prophage_dbscan_swa(strain_id,bac_fna_file,bac_faa_file,dbscan_region_file,swa_region_file,bac_protein_def_file,outdir,prefix):
	strain_def = strain_inf_dict[strain_id]
	with open(dbscan_region_file) as f:
		dbscan_contents = f.readlines()
	with open(swa_region_file) as f:
		swa_contents = f.readlines()
	
	intervals = []
	regions_inf = {}
	contents = dbscan_contents+swa_contents[1:]
	if len(contents)==1:
		print('0 prophage region was detected in the query bacterial genome!')
		sys.exit(1)
	bac_protein_def_dict = get_protein_def(bac_protein_def_file)
	uniprot_proteins_species = get_uniprot_def("")
	with open(bac_fna_file) as f:
		bac_fna_sequence = ''.join(f.read().strip().split('\n')[1:]).strip()
	bac_length = len(bac_fna_sequence)
	for line in contents[1:]:
		line = line.strip().split('\t')
		region_start = line[1].split(':')[0]
		region_end = line[1].split(':')[1]
		proid = line[3]
		key_proteins = line[2]
		pro = bac_protein_def_dict[proid.split('|')[1]]
		prodef = pro['prodef']
		key_protein = pro['key']
		if line[4]!='NA':
			hit_uniprot_id = line[4].split('|')[1]
			hit_uniprot_species = uniprot_proteins_species[hit_uniprot_id]
		else:
			hit_uniprot_id = 'NA'
			hit_uniprot_species = 'NA'
		identity = line[5]
		evalue = line[6]
		if [int(region_start),int(region_end)] not in intervals:
			intervals.append([int(region_start),int(region_end)])
			regions_inf.update({region_start+":"+region_end:[]})
			homo = {'proid':proid.strip('ref|'),'prodef':prodef,'key':key_protein,"blastp_homo_proid":hit_uniprot_id,"blastp_homo_def":hit_uniprot_species,'identity':identity,'evalue':evalue}
			regions_inf[region_start+":"+region_end].append(homo)
		else:
			homo = {'proid':proid.strip('ref|'),'prodef':prodef,'key':key_protein,"blastp_homo_proid":hit_uniprot_id,"blastp_homo_def":hit_uniprot_species,'identity':identity,'evalue':evalue}
			regions_inf[region_start+":"+region_end].append(homo)
		
	merge_inf = []
	intervals.sort(key=lambda x: x[0])		
	merged = []
	blastp_homo_species = []
	for interval in intervals:
		if not merged or int(merged[-1][1]) < int(interval[0]):
			merged.append(interval)
			start = merged[-1][0]
			end = merged[-1][1]
			hit_inf = regions_inf[str(interval[0])+":"+str(interval[1])]					
			merge_inf.append({"place":str(start),'end':str(end),'label':'DBSCAN-SWA','blastp_homo':hit_inf})
			species = []
			for homo in hit_inf:
				species.append(homo['blastp_homo_def'])
			blastp_homo_species.append({"place":str(start),'end':str(end),"blastp_phage":species})
		else:
			hit_inf_last = merge_inf[-1]['blastp_homo']
			merged[-1][1] = max(merged[-1][1],interval[1])
			start = merged[-1][0]
			end = merged[-1][1]
			hit_inf = regions_inf[str(interval[0])+":"+str(interval[1])]				
			for hit in hit_inf:
				if hit not in hit_inf_last:
					hit_inf_last.append(hit)
			merge_inf[-1]['end'] = str(end)
			merge_inf[-1]['blastp_homo'] = hit_inf_last
			blastp_homo_species[-1]['end'] = str(end)				
			for homo in hit_inf:
				if isinstance(homo,dict):
					blastp_homo_species[-1]['blastp_phage'].append(homo['blastp_homo_def'])
					if 'NA' in blastp_homo_species[-1]['blastp_phage']:
						blastp_homo_species[-1]['blastp_phage'].remove('NA') 
	#get the best hit species according to frequency
	# blastp_besthit = os.path.join(outdir,prefix+'_besthit_species_by_uniprot_db')
	# f = open(blastp_besthit,'w')
	region_blastp_besthit = {}
	for item in blastp_homo_species:
		start = item['place']
		end = item['end']
		per = list(Counter(item['blastp_phage']).most_common(1)[0])[1]
		per = round(float(per)/float(len(item['blastp_phage']))*100,2)
		hit_mf = max(Counter(item['blastp_phage']),key=Counter(item['blastp_phage']).get)+"("+str(per)+"%)"
		# f.write(str(start)+'\t'+str(end)+'\t'+str(length)+'\t'+hit_mf+'\n')
		region_blastp_besthit.update({str(start)+":"+str(end):hit_mf})
	# f.close()

	#record the predicted prophage region
	save_prophage_region_file = os.path.join(outdir,prefix+'_prophage_region.txt')
	save_prophage_region_detail_file = os.path.join(outdir,prefix+'_prophage_region_details.txt')
	f_save_region = open(save_prophage_region_file,'w')
	
	f_save_region.write('bacteria_id\tgenome_size\tprophage_start\tprophage_end\tpredict_method\tbest_hit_species\tkey_proteins\n')
	f_save_region_detail = open(save_prophage_region_detail_file,'w')
	f_save_region_detail.write('''The following contents displays predicted prophage regions
first line of each prophage describes the prophage information and the following lines describe the proteins and homology proteins in uniprot database
prophage_protein_ID\tprophage_protein_product\tkey_proteins\thit_protein_id\thit_species\tidentity\tevalue\n
''')
	save_prophage_protein_file = os.path.join(outdir,prefix+'_prophage.faa')
	save_prophage_nucl_file = os.path.join(outdir,prefix+'_prophage.fna')
	f_save_protein = open(save_prophage_protein_file,'w')
	f_save_nucl = open(save_prophage_nucl_file,'w')
	merge_inf1 = []
	prophage_counter = 1
	for inf in merge_inf:
		keys = []
		inf.update({"blastp_taxonomy":region_blastp_besthit[inf['place']+":"+inf['end']]})
		merge_inf1.append(inf)
		keys = []
		for homo in inf['blastp_homo']:
			if isinstance(homo,dict) and (homo['key'] not in keys):
				if homo['key']!='NA':
					keys.append(homo['key'])
		region_dir = os.path.join(outdir,'prophage_'+inf['place']+"_"+inf['end'])
		mkdir(region_dir)
		fna_sequence_region = bac_fna_sequence[int(inf['place'])-1:int(inf['end'])]
		fna_sequence_file = os.path.join(region_dir,'prophage_region.fna')
		fna_w = open(fna_sequence_file,'w')
		fna_w.write('>'+strain_id+'|'+inf['place']+":"+inf['end']+'|DBSCAN-SWA'+'\n'+fna_sequence_region+'\n')
		fna_w.close()
		f_save_nucl.write('>'+strain_id+'|'+inf['place']+":"+inf['end']+'|DBSCAN-SWA'+'\n'+fna_sequence_region+'\n')
		f_save_nucl.flush()
		region_faa_file = os.path.join(region_dir,'prophage_region.faa')
		region_pro_sequence_list,prophage_pro_num = get_protein_to_position_genbank(bac_faa_file,inf['place'],inf['end'],region_faa_file,'DBSCAN-SWA',strain_id)
		for pro_index,this_pro in enumerate(region_pro_sequence_list[0::2]):
			# print(this_pro)
			#print(region_pro_sequence_list[pro_index*2+1])
			f_save_protein.write(this_pro+'|'+str(prophage_pro_num)+'\n'+region_pro_sequence_list[2*pro_index+1]+'\n')
			f_save_protein.flush()
		if len(keys) ==0:
			key_protein = 'NA'
		else:
			key_protein = ','.join(keys)
		f_save_region.write(strain_id+'\t'+strain_def+'\t'+str(bac_length)+'\t'+str(inf['place'])+"\t"+str(inf['end'])+'\t'+inf['label']+'\t'+inf['blastp_taxonomy']+'\t'+key_protein+'\n')
		f_save_region_detail.write('>prophage '+str(prophage_counter)+'\n')
		prophage_counter = prophage_counter+1
		f_save_region_detail.write(strain_id+'\t'+strain_def+'\t'+str(bac_length)+'\t'+str(inf['place'])+"\t"+str(inf['end'])+'\t'+inf['label']+'\t'+inf['blastp_taxonomy']+'\t'+key_protein+'\n')
			
		for homo in inf['blastp_homo']:
			if isinstance(homo,dict):
				f_save_region_detail.write(homo['proid']+'\t'+homo['prodef']+'\t'+homo['key']+'\t'+homo['blastp_homo_proid']+'\t'+homo['blastp_homo_def']+'\t'+homo['identity']+'\t'+homo['evalue']+'\n')
	f_save_protein.close()
	f_save_region_detail.close()
	f_save_nucl.close()
	f_save_region.close()

def bac_blastn_phagedb(prophage_file,outfile):
	phage_database = os.path.join(root_path,'db','database','phage_nucl','phage_nucl_db')
	format = 6
	evalue = 0.01
	script = "blastn -db "+phage_database+" -query "+prophage_file+" -outfmt "+str(format)+" -evalue "+str(evalue)+" -out "+outfile+" -num_threads 20 -word_size 11 -perc_identity 70 -reward 1 -penalty -2 -gapopen 0 -gapextend 0 -max_target_seqs 1000000"
	os.system(script)

def bac_blastn_phage(prophage_file,phage_file,outfile):
	format = 6
	evalue = 0.01
	command = "blastn -query "+prophage_file+" -subject "+phage_file+" -outfmt "+str(format)+" -evalue "+str(evalue)+" -out "+outfile+" -word_size 11 -perc_identity 70 -reward 1 -penalty -2 -gapopen 0 -gapextend 0"
	print(command)
	os.system(command)

def diamond_blastp_nomax(file,outfile,database):
	num_threads = 20
	format = 6
	evalue = 1
	identity = 0.4
	coverage = 0.7
	script = "diamond blastp -d "+database+" -q "+file+" -f "+str(format)+" -e "+str(evalue)+" -o "+outfile+" -p "+str(num_threads)+" --id "+str(identity)+" --query-cover "+str(coverage)+" -k 1000000"
	os.system(script)

def bac_blastp_phagedb(prophage_file,outfile):
	num_threads = 20
	format = 6
	database = os.path.join(root_path,'db','database','phage_protein_db')
	diamond_blastp_nomax(prophage_file,outfile,database)

def bac_blastp_phage(prophage_file,phage_file,outfile):
	num_threads = 20
	evalue = 0.01
	format = 6
	command = "blastp -query "+prophage_file+" -subject "+phage_file+" -outfmt "+str(format)+" -evalue "+str(evalue)+" -out "+outfile
	print(command)
	os.system(command)

def combine(lst):
	lst.sort(key=lambda obj:int(obj[0]))
	templist = lst[0]
	combinedList = []
	if len(lst)==1:
		return lst
	for newlist in lst[1:]:
		newListmin = min(int(newlist[0]),int(newlist[1]))
		newListmax = max(int(newlist[0]), int(newlist[1]))
		tempListmax = max(int(templist[0]), int(templist[1]))
		if int(newListmin) <= int(tempListmax):
			templist[1] = str(max(int(tempListmax),int(newListmax)))
		else:
			combinedList.append(templist)
			templist = newlist
	combinedList.append(templist)
	return combinedList

def save_dict(my_dict,save_file):
	js = json.dumps(my_dict) 
	file = open(save_file, 'w') 
	file.write(js) 
	file.close()  

def load_dict(save_file):
	file = open(save_file, 'r') 
	protein_js = file.read()
	file.close()
	my_dict = json.loads(protein_js)
	return my_dict

def parse_prophage_blastp_phagedb_6(outfile,outdir):
	result_dict = {}
	outfile_dict = os.path.join(outdir,'prophage_blastp_phage_result_dict')
	resufile = os.path.join(outdir,'prophage_blastp_phage_result.txt')
	f_result = open(resufile,'w')
	f_result.write('bac_id\tbac_def\thit_prophage\thit_phage_id\thit_phage_def\tprophage_pro_num\thit_pro_percent\thit_pro_num\n')
	if os.path.exists(outfile):
		with open(outfile) as f:
			contents = f.readlines()
		for line in contents:
			line = line.strip().split('\t')
			prophage_pro_id = line[0]
			phage_pro = line[1]
			# if phage_type=='fasta':
			# 	phage_id = '_'.join(phage_pro.split('_')[0:-3]).split('.')[0]
			# else:
			phage_id = phage_pro.split('|')[0].split('.')[0]
			bac_id = prophage_pro_id.split('|')[0].split('.')[0]
			method = prophage_pro_id.split('|')[-2]
			prophage_region_pro_num = float(prophage_pro_id.split('|')[-1])
			prophage_id =  '|'.join(prophage_pro_id.split('|')[0:2])+'|'+'|'.join(prophage_pro_id.split('|')[-2:])
			start = prophage_id.split('|')[1].split(':')[0]
			end = prophage_id.split('|')[1].split(':')[1]
			if bac_id not in result_dict.keys():
				result_dict.update({bac_id:{}})
			if phage_id not in result_dict[bac_id].keys():
				result_dict[bac_id].update({phage_id:{}})
			if method not in result_dict[bac_id][phage_id].keys():
				result_dict[bac_id][phage_id].update({method:{}})
			if start+':'+end not in result_dict[bac_id][phage_id][method].keys():
				result_dict[bac_id][phage_id][method].update({start+':'+end:[]})
			result_dict[bac_id][phage_id][method][start+':'+end].append([prophage_pro_id,phage_pro,prophage_region_pro_num])
	result_dict1 = {}
	for bac_id in result_dict.keys():
		try:
			bac_def = strain_inf_dict[bac_id.split('.')[0]]
		except:
			bac_def = bac_id
		for phage_id in result_dict[bac_id].keys():
			phage_def = phage_inf_dict[phage_id]
			for method,regions_infos in result_dict[bac_id][phage_id].items():
				for region,hit_infos in regions_infos.items():
					region_pro_num = hit_infos[0][-1]
					prophage_hit_pro_num = len(list(set([item[0] for item in hit_infos])))
					percent = prophage_hit_pro_num/float(region_pro_num)
					if percent>=float(per)/100:
						if bac_id not in result_dict1.keys():
							result_dict1.update({bac_id:{}})
						if phage_id not in result_dict1[bac_id].keys():
							result_dict1[bac_id].update({phage_id:{}})
						if method not in result_dict1[bac_id][phage_id].keys():
							result_dict1[bac_id][phage_id].update({method:{}})
						result_dict1[bac_id][phage_id][method].update({region:[]})
						result_dict1[bac_id][phage_id][method][region] = [percent,len(hit_infos),region_pro_num,hit_infos]
					prophage_id = '|'.join(hit_infos[0][0].split('|')[0:2]+hit_infos[0][0].split('|')[-2:])
					f_result.write('\t'.join([bac_id,bac_def,prophage_id,phage_id,phage_def,str(region_pro_num),str(percent),str(len(hit_infos))])+'\n')
					f_result.flush()
	# np.save(outfile_dict,result_dict1)
	save_dict(result_dict1,outfile_dict)
	f_result.close()

def parse_prophage_blastn_phagedb_6(outfile,outdir):
	with open(outfile) as f:
		contents = f.readlines()
	result_dict = {}
	outfile_dict = os.path.join(outdir,'prophage_blastn_phage_result_dict')
	resufile = os.path.join(outdir,'prophage_blastn_phage_result.txt')
	f_result = open(resufile,'w')
	f_result.write('bac_id\tbac_def\tprophage_id\thit_phage_id\thit_phage_def\tprophage_length\thit_length\tidentity\tcoverage\n')
	for line in contents:
		line = line.strip().split('\t')
		prophage_id = line[0]
		phage_id = line[1].split('.')[0]
		bac_id = prophage_id.split('|')[0].split('.')[0]
		method = prophage_id.split('|')[-1]
		prophage_length = abs(int(prophage_id.split('|')[1].split(':')[1])-int(prophage_id.split('|')[1].split(':')[0])+1)
		identity = line[2]
		alignment_length = int(line[3])
		mismatch = int(line[4])
		bit_score = float(line[-1])
		start = prophage_id.split('|')[1].split(':')[0]
		end = prophage_id.split('|')[1].split(':')[1]
		query_start = line[-6]
		query_end = line[-5]
		hit_start = line[-4]
		hit_end = line[-3]
		if float(identity) >= float(idn):#mismatch<=2
			if bac_id not in result_dict.keys():
				result_dict.update({bac_id:{}})
			#coverage = alignment_length/float(prophage_length)
			if phage_id not in result_dict[bac_id].keys():
				result_dict[bac_id].update({phage_id:{}})
			if method not in result_dict[bac_id][phage_id].keys():
				result_dict[bac_id][phage_id].update({method:{}})
			if start+':'+end not in result_dict[bac_id][phage_id][method].keys():
				result_dict[bac_id][phage_id][method].update({start+':'+end:[]})
			result_dict[bac_id][phage_id][method][start+':'+end].append([prophage_id,identity,[query_start,query_end],bit_score])
	
	result_dict1 = {}
	for bac_id in result_dict.keys():
		try:
			bac_def = strain_inf_dict[bac_id.split('.')[0]]
		except:
			bac_def = bac_id
		for phage_id in result_dict[bac_id].keys():
			phage_def = phage_inf_dict[phage_id]
			for method in result_dict[bac_id][phage_id].keys():
				identitys = []
				for region in result_dict[bac_id][phage_id][method].keys():
					#print(region)
					hit_infos = result_dict[bac_id][phage_id][method][region]
					hit_prophage = hit_infos[0][0]
					prophage_length = abs(int(region.split(':')[1])-int(region.split(':')[0]))+1
					hit_region = [x[2] for x in hit_infos]
					identitys = [float(x[1]) for x in hit_infos]
					identity = np.mean(identitys)
					merge_region = combine(hit_region)
					hit_length = 0
					for region1 in merge_region:
						region_length = abs(int(region1[0])-int(region1[1]))+1
						hit_length = hit_length+region_length
					coverage = float(hit_length)/prophage_length
					if coverage>=float(cov)/100:
						if bac_id not in result_dict1.keys():
							result_dict1.update({bac_id:{}})
						if phage_id not in result_dict1[bac_id].keys():
							result_dict1[bac_id].update({phage_id:{}})
						if method not in result_dict1[bac_id][phage_id].keys():
							result_dict1[bac_id][phage_id].update({method:{}})
						result_dict1[bac_id][phage_id][method].update({region:[]})
						result_dict1[bac_id][phage_id][method][region] = [prophage_id,identity,coverage,hit_infos]
					f_result.write('\t'.join([bac_id,bac_def,hit_prophage,phage_id,phage_def,str(prophage_length),str(hit_length),str(identity),str(coverage)])+'\n')
					f_result.flush()
	# np.save(outfile_dict,result_dict1)
	save_dict(result_dict1,outfile_dict)
	f_result.close()

def filter_identity_coverage(blastp_file,filter_file,target):
	f_result = open(filter_file,'w')
	with open(blastp_file) as f:
		contents = f.readlines()
	for line in contents:
		line = line.strip().split('\t')
		query_pro = line[0]
		hit_pro = line[1]
		identity = line[2]
		if target == 'query':
			query_start = line[-6]
			query_end = line[-5]
			hit_length = abs(int(query_start)-int(query_end))+1
			pro_length = int((abs(int(query_pro.split('|')[3].split('_')[-3])-int(query_pro.split('|')[3].split('_')[-2]))+1)/3)
		else:
			hit_start = line[-4]
			hit_end = line[-3]
			hit_length = abs(int(hit_start)-int(hit_end))+1
			pro_length = int((abs(int(hit_pro.split('|')[3].split('_')[-3])-int(hit_pro.split('|')[3].split('_')[-2]))+1)/3)
		coverage = float(hit_length)/pro_length
		if (coverage>=0.7) and ((float(identity)/100)>=0.4):
			f_result.write('\t'.join(line)+'\t'+str(coverage)+'\n')
			f_result.flush()
	f_result.close()

def prophage_annotate_phagedb(prophage_region_file,prophage_protein_file,prophage_nucl_file,outdir,prefix):
	outfile_blastp = os.path.join(outdir,'prophage_blastp_phage')
	outfile_blastn = os.path.join(outdir,'prophage_blastn_phage')
	if os.path.exists(prophage_protein_file) or os.path.exists(prophage_nucl_file):
		m1 = threading.Thread(target=bac_blastn_phagedb,args=(prophage_nucl_file,outfile_blastn))
		m2 = threading.Thread(target=bac_blastp_phagedb,args=(prophage_protein_file,outfile_blastp))
		m1.start()
		m2.start()
		m1.join()
		m2.join()
	else:
		print('prophage protein file or nucleotide file was not created, no prophage was detected!!')
		return(0)
	time.sleep(1)
	filter_file = os.path.join(outdir,'prophage_blastp_phage_identity_0.4_coverage_0.7')
	filter_identity_coverage(outfile_blastp,filter_file,'query')
	time.sleep(1)
	if os.path.exists(outfile_blastn):
		m1 = threading.Thread(target=parse_prophage_blastn_phagedb_6,args=(outfile_blastn,outdir))
	else:
		print('prophage nucleotide file was not created, no prophage was detected!!')
	if os.path.exists(outfile_blastp):
		m2 = threading.Thread(target=parse_prophage_blastp_phagedb_6,args=(filter_file,outdir))
	else:
		print('prophage protein file was not created, no prophage was detected!!')	
	if os.path.exists(outfile_blastn):	
		m1.start()
	if os.path.exists(outfile_blastp):
		m2.start()
	if os.path.exists(outfile_blastn):
		m1.join()
	if os.path.exists(outfile_blastp):
		m2.join()

	outfile_blastn_dict_file = os.path.join(outdir,'prophage_blastn_phage_result_dict.npy')
	outfile_blastp_dict_file = os.path.join(outdir,'prophage_blastp_phage_result_dict.npy')
	
	try:
		get_bac_phage_feature(prophage_region_file,outfile_blastp_dict_file,outfile_blastn_dict_file,outdir,prefix)
	except:
		print('no infecting phage was detected!')

def get_bac_phage_feature(prophage_region_file,prophage_blastp_dict_file,prophage_blastn_dict_file,outdir,prefix):
	all_prophage_info_dict = {}
	with open(prophage_region_file) as f:
		contents = f.read().strip().split('>prophage')
	for line in contents[1:]:
		line = line.strip().split('\n')[1].split('\t')
		bac_def = line[1]
		prophage_region = line[3]+':'+line[4]
		bac_id = line[0]
		all_prophage_info_dict.update({bac_id+'_'+prophage_region:line})
	#print(all_prophage_info_dict)
	prophage_blastp_result_file = os.path.join(outdir,"prophage_blastp_phage_result.txt")
	prophage_blastn_result_file = os.path.join(outdir,"prophage_blastn_phage_result.txt")
	prophage_blastp_original_file = os.path.join(outdir,'prophage_blastp_phage_identity_0.4_coverage_0.7')
	prophage_blastn_original_file = os.path.join(outdir,'prophage_blastn_phage')
	
	prophage_file_blastp_dict = {}
	if os.path.exists(prophage_blastp_result_file):
		with open(prophage_blastp_result_file) as f:
			contents = f.readlines()
		for line in contents[1:]:
			line = line.strip().split('\t')
			bac_id = line[0].split('.')[0]
			bac_def = line[1]
			phage_id = line[3].split('.')[0]
			hit_prophage  = line[2]
			method = hit_prophage.split('|')[-2]#
			region = hit_prophage.split('|')[1]
			percent = line[-2]
			if phage_id not in prophage_file_blastp_dict.keys():
				prophage_file_blastp_dict.update({phage_id:{}})
			if bac_id not in prophage_file_blastp_dict[phage_id].keys():
				prophage_file_blastp_dict[phage_id].update({bac_id:{}})
			if method not in prophage_file_blastp_dict[phage_id][bac_id].keys():
				prophage_file_blastp_dict[phage_id][bac_id].update({method:{}})
			if region not in prophage_file_blastp_dict[phage_id][bac_id][method].keys():
				prophage_file_blastp_dict[phage_id][bac_id][method].update({region:float(percent)})

	prophage_file_blastn_dict = {}
	if os.path.exists(prophage_blastn_result_file):
		with open(prophage_blastn_result_file) as f:
			contents = f.readlines()
		for line in contents[1:]:
			line = line.strip().split('\t')
			bac_id = line[0].split('.')[0]
			phage_id = line[3].split('.')[0]
			hit_prophage  = line[2]
			method = hit_prophage.split('|')[-1]
			region = hit_prophage.split('|')[1]
			identity = line[-2]
			coverage = line[-1]
			coverage = round(float(coverage),3)
			if phage_id not in prophage_file_blastn_dict.keys():
				prophage_file_blastn_dict.update({phage_id:{}})
			if bac_id not in prophage_file_blastn_dict[phage_id].keys():
				prophage_file_blastn_dict[phage_id].update({bac_id:{}})
			if method not in prophage_file_blastn_dict[phage_id][bac_id].keys():
				prophage_file_blastn_dict[phage_id][bac_id].update({method:{}})
			if region not in prophage_file_blastn_dict[phage_id][bac_id][method].keys():
				prophage_file_blastn_dict[phage_id][bac_id][method].update({region:[]})
			prophage_file_blastn_dict[phage_id][bac_id][method][region].append([float(identity),float(coverage)])
	
	prophage_blastp_hit_dict = {}
	with open(prophage_blastp_original_file) as f:
	    contents = f.readlines()
	for line in contents:
		line = line.strip().split('\t')
		phage_pro_id = line[1]
		prophage_pro_id = line[0]
		identity = line[2]
		coverage = line[-1]
		coverage = round(float(coverage),3)
		evalue = line[-3]
		bac_id = prophage_pro_id.split('|')[0].split('.')[0]
		method = prophage_pro_id.split('|')[-2]
		#supposing id no |
		if '|' in phage_pro_id:
		    if len(phage_pro_id.split('|'))>=3:
		        phage_id = '|'.join(phage_pro_id.split('|')[0:-2]).split('.')[0]
		    else:
		        phage_id = '_'.join(phage_pro_id.split('_')[0:-3]).split('.')[0]
		else:
		    phage_id = '_'.join(phage_pro_id.split('_')[0:-3]).split('.')[0]
		region = prophage_pro_id.split('|')[1]
		if phage_id not in prophage_blastp_hit_dict.keys():
		    prophage_blastp_hit_dict.update({phage_id:{}})
		if bac_id not in prophage_blastp_hit_dict[phage_id].keys():
		    prophage_blastp_hit_dict[phage_id].update({bac_id:{}})
		if method not in prophage_blastp_hit_dict[phage_id][bac_id].keys():
		    prophage_blastp_hit_dict[phage_id][bac_id].update({method:{}})
		if region not in prophage_blastp_hit_dict[phage_id][bac_id][method].keys():
		   prophage_blastp_hit_dict[phage_id][bac_id][method].update({region:[]})
		prophage_blastp_hit_dict[phage_id][bac_id][method][region].append({'phage_id':phage_id,
			'bac_id':bac_id,
			'phage_pro_id':phage_pro_id,
			'prophage_pro_id':prophage_pro_id,
			'identity':identity,
			'coverage':coverage,
			'evalue':evalue})
	
	#print(prophage_blastp_hit_dict)
	prophage_blastn_hit_dict = {}
	with open(prophage_blastn_original_file) as f:
		contents = f.readlines()
	for line in contents:
		line = line.strip().split('\t')
		phage_id = line[1].split('.')[0]
		prophage_id = line[0]
		identity = line[2]
		evalue = line[-2]
		hit_prophage_start = line[-6]
		hit_prophage_end = line[-5]
		phage_hit_start = line[-4]
		phage_hit_end = line[-3]
		phage_hit_region = phage_hit_start+":"+phage_hit_end
		prophage_hit_region = hit_prophage_start+":"+hit_prophage_end
		alignment_length = abs(int(hit_prophage_start)-int(hit_prophage_end))+1
		bac_id = prophage_id.split('|')[0].split('.')[0]
		method = prophage_id.split('|')[2]
		region = prophage_id.split('|')[1]
		if float(identity)>=70:#web shows the original contents
			if phage_id not in prophage_blastn_hit_dict.keys():
				prophage_blastn_hit_dict.update({phage_id:{}})
			if bac_id not in prophage_blastn_hit_dict[phage_id].keys():
				prophage_blastn_hit_dict[phage_id].update({bac_id:{}})
			if method not in prophage_blastn_hit_dict[phage_id][bac_id].keys():
				prophage_blastn_hit_dict[phage_id][bac_id].update({method:{}})
			if region not in prophage_blastn_hit_dict[phage_id][bac_id][method].keys():
				prophage_blastn_hit_dict[phage_id][bac_id][method].update({region:[]})
			prophage_blastn_hit_dict[phage_id][bac_id][method][region].append({'phage_id':phage_id,
			'bac_id':bac_id,
			'prophage_id':prophage_id,
			'identity':identity,
			'evalue':evalue,
			'prophage_hit_region':prophage_hit_region,
			'alignment_length':alignment_length,
			'phage_hit_region':phage_hit_region
			})

	prophage_blastp_dict = {}
	prophage_blastn_dict = {}
	try:
		prophage_blastp_dict = load_dict(prophage_blastp_dict_file)
	except:
		pass
	try:
		prophage_blastn_dict = load_dict(prophage_blastn_dict_file)
	except:
		pass
	# print(prophage_blastp_dict)
	# print(prophage_file_blastp_dict)
	save_prophage_hit_phage_dict = {}	
	result_file = os.path.join(outdir,prefix+'_prophage_annotate_phage.txt')
	f_prophage_result = open(result_file,'w')
	f_prophage_result.write('bac_id\tbac_def\tgenome_size\tprophage_start\tprophage_end\tmethod\tbest_hit_phage_species(turnout)\tprophage_key_proteins\thit_phage_id\thit_phage_def\tprophage_homolog_percent\tprophage_alignment_identity\tprophage_alignment_coverage\n')
	result_detail_file = os.path.join(outdir,prefix+'_prophage_annotate_phage_details.txt')
	f_prophage_result_detail = open(result_detail_file,'w')
	#f_prophage_result_detail.write('bac_id\tbac_def\tgenome_size\tprophage_start\tprophage_end\tmethod\tbest_hit_phage_species(turnout)\tprophage_key_proteins\thit_phage_id\thit_phage_def\tprophage_homolog_percent\tprophage_alignment_identity\tprophage_alignment_coverage\n')

	all_query_list = list(set(list(prophage_blastp_dict.keys())+list(prophage_blastn_dict.keys())))
	for contig_id in all_query_list:
		contig_id = contig_id.split('.')[0]	
		if contig_id in strain_inf_dict.keys():
			contig_def = strain_inf_dict[contig_id]
		else:
			conti_def = query_id
		if contig_id in prophage_blastp_dict.keys():
			blastp_hit_list = list(prophage_blastp_dict[contig_id].keys())
		else:
			blastp_hit_list = []
		if contig_id in prophage_blastn_dict.keys():
			blastn_hit_list = list(prophage_blastn_dict[contig_id].keys())
		else:
			blastn_hit_list = []
		hit_phage_list = list(set(blastp_hit_list+blastn_hit_list))
		print(hit_phage_list)
		for phage_id in hit_phage_list:
			phage_id = phage_id.split('.')[0]
			try:
				phage_def = phage_inf_dict[phage_id]
			except:
				phage_def = phage_id
			if contig_id in prophage_blastp_dict.keys():
				if phage_id in prophage_blastp_dict[contig_id].keys():
					blastp_prophage_info = prophage_blastp_dict[contig_id][phage_id]
					print(blastp_prophage_info)
					if contig_id not in save_prophage_hit_phage_dict.keys():
						save_prophage_hit_phage_dict.update({contig_id:{}})
					for method in blastp_prophage_info.keys():
						for prophage_region,hit_info in blastp_prophage_info[method].items():
							if prophage_region not in save_prophage_hit_phage_dict[contig_id].keys():
								save_prophage_hit_phage_dict[contig_id].update({prophage_region:[]})
							prophage_homology_percent = prophage_file_blastp_dict[phage_id][contig_id][method][prophage_region]
							if (prophage_region=='2652209:2659269') and (phage_id=='KU052037'):
								print(hit_info)
								print(prophage_homology_percent)
							try:
								#print(contig_id,phage_id,method,prophage_region)
								prophage_identity =prophage_file_blastn_dict[phage_id][contig_id][method][prophage_region][0][0]
								prophage_coverage = prophage_file_blastn_dict[phage_id][contig_id][method][prophage_region][0][1]
							except:
								prophage_identity = 0
								prophage_coverage = 0
							save_prophage_hit_phage_dict[contig_id][prophage_region].append([phage_id,phage_def,prophage_homology_percent,prophage_identity,prophage_coverage]) 

			if contig_id in prophage_blastn_dict.keys():
				if phage_id in prophage_blastn_dict[contig_id].keys():
					blastn_prophage_info = prophage_blastn_dict[contig_id][phage_id]
					contig_id = contig_id.split('.')[0]
					if contig_id not in save_prophage_hit_phage_dict.keys():
						save_prophage_hit_phage_dict.update({contig_id:{}})
					for method in blastn_prophage_info.keys():
						for region,hit_infos in blastn_prophage_info[method].items():
							if region not in save_prophage_hit_phage_dict[contig_id].keys():
								save_prophage_hit_phage_dict[contig_id].update({region:[]})
								region_identity = hit_infos[1]
								region_coverage = hit_infos[2]
								try:
									prophage_homology_percent = prophage_file_blastp_dict[phage_id][contig_id][method][region]
								except:
									#print(phage_id,contig_id,method,prophage_region)
									prophage_homology_percent = 0
								# if (region=='1288022:1331541') and (phage_id=='GQ421471'):
								# 	print(hit_infos)
								# 	print(prophage_homology_percent)
								save_prophage_hit_phage_dict[contig_id][region].append([phage_id,phage_def,prophage_homology_percent,region_identity,region_coverage])
	counter = 1
	#print(save_prophage_hit_phage_dict)
	for contig_id in save_prophage_hit_phage_dict.keys():		
		for prophage_region in save_prophage_hit_phage_dict[contig_id].keys():
			for hit_phage_info in save_prophage_hit_phage_dict[contig_id][prophage_region]:
				if contig_id+'_'+prophage_region not in all_prophage_info_dict.keys():
					continue
				prophage_info = all_prophage_info_dict[contig_id+'_'+prophage_region]
				phage_id = hit_phage_info[0]
				#print(phage_id,contig_id,prophage_region)

				try:
					blastp_infos = prophage_blastp_hit_dict[phage_id][contig_id]['DBSCAN-SWA'][prophage_region] 
				except:
					blastp_infos = []
				try:
					blastn_infos = prophage_blastn_hit_dict[phage_id][contig_id]['DBSCAN-SWA'][prophage_region]
				except:
					blastn_infos = []
				f_prophage_result.write('\t'.join(prophage_info)+'\t'+'\t'.join(list(map(str,hit_phage_info)))+'\n')
				f_prophage_result_detail.write('>PHIS '+str(counter)+'\n')
				counter = counter+1
				f_prophage_result_detail.write('\t'.join(prophage_info)+'\t'+'\t'.join(list(map(str,hit_phage_info)))+'\n')
				f_prophage_result_detail.write(''.join(['#blastp']*10)+'\n')
				for blastp_info in blastp_infos:
					f_prophage_result_detail.write('\t'.join([
						blastp_info['bac_id'],
						blastp_info['prophage_pro_id'],
						blastp_info['phage_id'],
						blastp_info['phage_pro_id'],
						str(blastp_info['identity']),
						str(blastp_info['coverage']),
						str(blastp_info['evalue'])])+'\n')
				f_prophage_result_detail.write(''.join(['#blastn']*10)+'\n')
				for blastn_info in blastn_infos:
					f_prophage_result_detail.write('\t'.join([
						blastn_info['bac_id'],
						blastn_info['prophage_id'],
						blastn_info['prophage_hit_region'],
						blastn_info['phage_id'],
						blastn_info['phage_hit_region'],
						str(blastn_info['alignment_length']),
	
						str(blastn_info['identity']),
						str(blastn_info['evalue'])]
						)+'\n')
	f_prophage_result.close()
	f_prophage_result_detail.close()

def prophage_annotate_phage(prophage_region_file,prophage_protein_file,prophage_nucl_file,phage_pro_file,fasta_phage_file,outdir,prefix):
	#step2:prophage blast phage database
	outfile_blastp = os.path.join(outdir,'prophage_blastp_phage')
	outfile_blastn = os.path.join(outdir,'prophage_blastn_phage')
	m1 = threading.Thread(target=bac_blastp_phage,args=(prophage_protein_file,phage_pro_file,outfile_blastp))
	m2 = threading.Thread(target=bac_blastn_phage,args=(prophage_nucl_file,fasta_phage_file,outfile_blastn))
	m1.start()
	m2.start()
	m1.join()
	m2.join()

	filter_file = os.path.join(outdir,'prophage_blastp_phage_identity_0.4_coverage_0.7')
	filter_identity_coverage(outfile_blastp,filter_file,'query')
	
	#step3:parse the blast result
	try:
		parse_prophage_blastp_phagedb_6(filter_file,outdir)
	except:
		print('parse blastp file:%s Error!'%outfile_blastp)
	try:
		parse_prophage_blastn_phagedb_6(outfile_blastn,outdir)
	except:
		print('parse blastn file:%s Error!'%outfile_blastn)
	
	#step4:get the feature
	prophage_blastp_dict = os.path.join(outdir,'prophage_blastp_phage_result_dict.npy')
	prophage_blastn_dict = os.path.join(outdir,'prophage_blastn_phage_result_dict.npy')
	try:
		get_bac_phage_feature(prophage_region_file,prophage_blastp_dict,prophage_blastn_dict,outdir,prefix)
	except:
		print('no phage detected!')

def	annotate_prophage_region(strain_id,prophage_region_file,prophage_protein_file,prophage_nucl_file,annotate_outdir,add_annotation,prefix):
	if add_annotation == 'PGPD':
		prophage_annotate_phagedb(prophage_region_file,prophage_protein_file,prophage_nucl_file,annotate_outdir,prefix)	
	else:
		#user input phage
		if os.path.exists(add_annotation):
			with open(add_annotation) as f:
				contents = f.read().strip()
			if contents[0]=='>':
				phage_nucl_file = add_annotation
				phage_id = contents.split('\n')[0].split()[0].split('.')[0].strip('>')
				try:
					get_bac_protein(add_annotation,annotate_outdir,'fasta','phage',phage_id,'Viruses')
				except:
					print('parse %s error!'%add_annotation)
					sys.exit(1)
			else:
				record = SeqIO.read(add_annotation,"genbank")
				phage_id = record.id
				phage_id = phage_id.split('.')[0]
				phage_nucl_file = os.path.join(annotate_outdir,'phage.fna')
				try:
					get_bac_protein(add_annotation,annotate_outdir,'genbank','phage',phage_id,'Viruses')
				except:
					print('parse %s error!'%add_annotation)
					sys.exit(1)
			phage_pro_file = os.path.join(annotate_outdir,'phage_protein.faa')
			prophage_annotate_phage(prophage_region_file,prophage_protein_file,prophage_nucl_file,phage_pro_file,phage_nucl_file,annotate_outdir,prefix)
		else:
			print('%s does not exist'%add_annotation)
			sys.exit(1)

def visualize(save_prophage_file,save_prophage_nucl_file,save_prophage_protein_file,prophage_region_annotate_dir,add_annotation,templeate_file,save_html_file,prefix):
	prophage_nucl_dict = {}
	with open(save_prophage_nucl_file) as f:
		contents = f.read().strip().split('>')	
	for line in contents[1:]:
		prophage_id = line.split('\n')[0].strip().split('|')[1]
		prophage_nucl_dict.update({prophage_id:'>'+line.strip()})
	
	prophage_protein_dict = {}
	with open(save_prophage_protein_file) as f:
		contents = f.read().strip().split('>')	
	for line in contents[1:]:
		prophage_id = line.split('\n')[0].strip().split('|')[1]
		pro_id = '|'.join(line.split('\n')[0].strip().split('|')[2:-2])
		if prophage_id not in prophage_protein_dict.keys():
			prophage_protein_dict.update({prophage_id:{}})
		prophage_protein_dict[prophage_id].update({pro_id:'>'+''.join(line).strip()})
	
	with open(templeate_file) as f:
		templeates = f.read().strip()
	f_save = open(save_html_file,'w')
	f_save.write(templeates+'\n')
	with open(save_prophage_file) as f:
		contents = f.read().strip().split('>prophage')
	
	genome_size = contents[1].split('\n')[1].split('\t')[2]
	f_save.write('''<plasmid sequencelength=%s plasmidheight=600 plasmidwidth=800 style='border:0px;'>
			'''%genome_size)
	f_save.write('''<plasmidtrack trackstyle="fill:#ccc" width="5" radius="160"></plasmidtrack>
			<plasmidtrack trackstyle="fill:rgba(225,225,225,0.5)" radius=150>
			''')
	counter = 0
	for line in contents[1:]:
		prophage_info = line.split('\n')[1].split('\t')
		strain_id = prophage_info[0]
		strain_def = prophage_info[1]
		genome_size = prophage_info[2]
		prophage_start = prophage_info[3]
		prophage_end = prophage_info[4]
		attl_region = prophage_info[-2]
		attr_region = prophage_info[-1]
		key_word = prophage_info[-5]
		if attl_region!='NA':
			attl_start = attl_region.split(':')[0]
			attl_end = attl_region.split(':')[1]
			attr_start = attr_region.split(':')[0]
			attr_end = attr_region.split(':')[1]
			attl_sequence = [item.split('\t')[2] for item in line.strip().split('\n')[2:] if item.split('\t')[1]=='attL']
			attr_sequence = [item.split('\t')[2] for item in line.strip().split('\n')[2:] if item.split('\t')[1]=='attR']
		if counter==0:
			f_save.write('''<tracklabel text='%s' style="font-size:15px;font-weight:bold"></tracklabel>
  <tracklabel text='%s' vadjust=-250 style="font-size:17px;font-weight:bold"></tracklabel>
				'''%(strain_id,strain_def))
			f_save.write('''<trackscale interval='10000' style='stroke:#ccc'></trackscale>
  <trackscale interval='50000' direction='in' style='stroke:#ccc'></trackscale>
  <trackscale interval='400000' direction='in' showlabels=1 labelstyle='fill:#999;stroke:none;text-anchor:middle;alignment-baseline:middle;font-size:10px'></trackscale>
			''')
		counter = counter+1
		if attl_region!='NA':
			f_save.write('''<trackmarker start='%s' end='%s'  markerstyle='stroke:#CD5C5C;stroke-dasharray:2,2;stroke-width:2px' wadjust="30">
        <markerlabel text='%s' vadjust="30" style='fill:#CD5C5C'></markerlabel>
      </trackmarker>
			'''%(attl_start,attl_end,'attL:'+attl_sequence[0]))	
			f_save.write('''<trackmarker start='%s' end='%s'  markerstyle='stroke:#CD5C5C;stroke-dasharray:2,2;stroke-width:2px' wadjust="50">
        <markerlabel text='%s' vadjust="50" style='fill:#CD5C5C'></markerlabel>
      </trackmarker>
			'''%(attr_start,attr_end,'attR:'+attr_sequence[0]))
		#print(key_word)
		if key_word.strip()=='':
			f_save.write('''<trackmarker start='%s' end='%s' markerstyle='fill:#7CEDDB'>
      <markerlabel text='%s' vadjust=20 style='fill:#7CEDDB'></markerlabel>
    </trackmarker>
			'''%(prophage_start,prophage_end,str(counter)))
		else:
			f_save.write('''<trackmarker start='%s' end='%s' markerstyle='fill:#7CEDDB'>
      <markerlabel text='%s' vadjust=20 style='fill:#AD1A0D'></markerlabel>
    </trackmarker>
			'''%(prophage_start,prophage_end,str(counter)))
	f_save.write('</plasmidtrack>\n</plasmid>\n</div>\n</div>\n')
	f_save.write('''<div id="CpG" style="margin:50px;font-size: 12px;">
  <table class="display" class="table table-striped table-hover table-bordered" style="border:0px;"> 
<thead>
      <tr style="text-align:center;">
      <th>Prophage_ID</th>
      <th>Prophage_region</th>
      <th>prophage_CDS_num</th>
      <th>Prophage_key_proteins</th>
      <th>Prophage_best_hit_uniprot</th>
      <th>CDS_details</th>
      </tr>
    </thead>
			''')	
	counter = 0
	prophage_id_dict = {}
	for line in contents[1:]:
		counter = counter+1
		prophage_info = line.split('\n')[1].split('\t')
		strain_id = prophage_info[0]
		strain_def = prophage_info[1]
		genome_size = prophage_info[2]
		prophage_start = prophage_info[3]
		prophage_end = prophage_info[4]
		attl_region = prophage_info[-2]
		attr_region = prophage_info[-1]
		key_word = prophage_info[-5]
		prophage_id_dict.update({strain_id+'|'+prophage_start+':'+prophage_end:str(counter)})
		if key_word.strip()=='':
			key_word = 'NA'
		best_hit_uniprot = prophage_info[-4]
		prophage_pro_num = prophage_info[-3]
		if attl_region!='NA':
			attl_start = attl_region.split(':')[0]
			attl_end = attl_region.split(':')[1]
			attr_start = attr_region.split(':')[0]
			attr_end = attr_region.split(':')[1]
		prophage_nucl_sequence = prophage_nucl_dict[prophage_start+':'+prophage_end]
		prophage_protein_sequence = '\n'.join(list(prophage_protein_dict[prophage_start+':'+prophage_end].values()))
		f_save.write('''<tr style="text-align:center;">     
      <td>
      <a data-toggle="modal" data-target="#%s" style="cursor: pointer;">%s <span class="glyphicon glyphicon-info-sign"></span></a>
      <div class="modal fade" id="%s" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true">
    <div class="modal-dialog">
    <div class="modal-content" style="width:700px;text-align: left;">
      <div class="modal-header">
        <button type="button" class="close" data-dismiss="modal" aria-hidden="true">
        </button>
        <!-- <h4 class="modal-title" id="myModalLabel">Prophage region Sequence</h4> -->
      </div>
      <div class="modal-body">
        <h4>%s region DNA</h4>
       <p style="word-wrap:break-word;word-break:break-all;">
       <pre>%s</pre>
       </p>
        <h4>%s region protein</h4>
       <p style="word-wrap:break-word;word-break:break-all;">
       <pre>%s</pre>
       </p>
      </div>
       <div class="modal-footer">
        <button type="button" class="btn btn-default" data-dismiss="modal">close
        </button>
      </div>
    </div>
    </div>
  </div>
  </td>
			'''%(str(counter),'DBSCAN-SWA_'+str(counter),str(counter),'DBSCAN-SWA_'+str(counter),prophage_nucl_sequence,'DBSCAN-SWA_'+str(counter),prophage_protein_sequence))
		
		f_save.write('''
		<td>%s</td>
      <td>%s</td>
      <td>%s</td>
      <td>%s</td>
			'''%(prophage_start+':'+prophage_end,str(prophage_pro_num),key_word,best_hit_uniprot))
		
		f_save.write('''
		<td>
      <button data-toggle="modal" data-target="#detail_%s" style="background-color: #AD2417;color:white;width:50px;height:25px;border-radius: 5px;">view</button>
		<div class="modal fade" id="detail_%s" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true">
    <div class="modal-dialog">
    <div class="modal-content" style="width:900px;">
      <div class="modal-header">
        <button type="button" class="close" data-dismiss="modal" aria-hidden="true"></button>
        <h4 class="modal-title" id="myModalLabel">DBSCAN-SWA_%s CDS details</h4>
      </div>
      <div class="modal-body">
      <p style="color:#FFA07A;text-align: left;">The bacterium proteins that are colored denote the protein is present at specific phage-related keywords (such as 'capsid', 'head', 'integrase', 'plate', 'tail', 'fiber', 'coat', 'transposase', 'portal', 'terminase', 'protease' or 'lysin' and 'tRNA')</p>
        <table class="display" class="table table-striped table-hover table-bordered" style="text-align: center;">
          <thead>
          <tr>
            <th>bacterium_Protein</th>
            <th>bacterium_info</th>
            <th>key_word</th>
            <th>hit_uniprot_id</th>
            <th>hit_protein_info</th>
            <th>Identity(%%)</th>
            <th style="word-break: keep-all;white-space:nowrap;">E-value</th>
          </tr>
          </thead>
          <tbody>
          '''%(str(counter),str(counter),str(counter)))
		
		for prophage_protein in line.strip().split('\n')[2:]:
			prophage_region = prophage_protein.strip().split('\t')
			if prophage_region[2]!='NA':
				style_color = '#DD5044'
			else:
				style_color = 'black'
			f_save.write('''   
          <tr style="color:%s">
            <td style="word-break: break-all; word-wrap:break-word;">%s</td>
            <td style="word-break: break-all; word-wrap:break-word;">%s</td>
            <td style="word-break: break-all; word-wrap:break-word;">%s</td>
            <td>%s</td>
            <td>%s</td>
            <td>%s</td>
            <td style="word-break: keep-all;white-space:nowrap;">%s</td>
          </tr>
			'''%(style_color,prophage_region[0],prophage_region[1],prophage_region[2],prophage_region[3],prophage_region[4],prophage_region[5],prophage_region[6]))
		f_save.write('''
        	</tbody>
        </table>
        </div>
        </div>
        </div>
        </td>
        </tr>''')
	
	f_save.write('''
		</tbody>
        </table>
        </div>''')
	
	if add_annotation!='none':
		f_save.write('''<p style='font-size:20px;text-align:center;'><strong>Phage Annotation</strong></p>
			<hr style="color:black;">
			<div id="CpG" style="margin:50px;font-size: 12px;">
  <table class="display" class="table table-striped table-hover table-bordered" style="border:0px;"> 
<thead>
      <tr style="text-align:center;">
      <th>Prophage_ID</th>
      <th>Prophage_region</th>
      <th>Hit_Phage</th>
      <th>Prophage_homology_percent</th>
      <th>Prophage_alignment_identity(%)</th>
      <th>Prophage_alignment_coverage</th>
      <th>Hit_info</th>
      </tr>
    </thead>
    <tbody>
			''')
		prophage_hit_phage_file = os.path.join(prophage_region_annotate_dir,prefix+'_prophage_annotate_phage_details.txt')
		if os.path.exists(prophage_hit_phage_file):
			with open(prophage_hit_phage_file) as f:
				contents = f.read().strip().split('>PHIS')		
			counter = 0
			if len(contents)<2:
				return 0
			for phis in contents[1:]:
				phis = phis.strip().split('\n')
				phis_line = phis[1].strip().split('\t')
				bac_id = phis_line[0]
				hit_phage_id = phis_line[-5]
				hit_phage_def = phis_line[-4]
				prophage_id = bac_id+'|'+':'.join(phis_line[3:5])
				prophage_counter = prophage_id_dict[prophage_id]
				counter = counter+1
				f_save.write('''
				<tr>
					<td>DBSCAN-SWA_%s</td>
					<td>%s</td>
	      			<td><a href="https://www.ncbi.nlm.nih.gov/nuccore/%s">%s</a></td>
	      			<td>%s</td>
	      			<td>%s</td>
	      			<td>%s</td>
	      			'''%(prophage_counter,prophage_id,hit_phage_id,hit_phage_def,phis_line[-3],phis_line[-2],phis_line[-1]))
				if len(phis)==2:
					f_save.write('''
					<td style="color:#F05555;font-size: 15px;">
					<span class="glyphicon glyphicon-remove-sign" aria-hidden="true"></span>
					</td>
					</tr>
					''')
				else:
					f_save.write('''
	<td style="font-size: 13px;">
	<button data-toggle="modal" data-target="#DBSCAN-SWA_%s" style="background-color: #AD2417;color:white;width:50px;height:25px;border-radius: 5px;">detail</button>
	<div class="modal fade" id="DBSCAN-SWA_%s" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true">
	    <div class="modal-dialog">
	    <div class="modal-content" style="width:800px;text-align:center;">
	      <div class="modal-header">
	        <button type="button" class="close" data-dismiss="modal" aria-hidden="true"></button>
	        <h4 class="modal-title" id="myModalLabel">Blastp and Blastn of %s VS %s</h4>
	      </div>
	      <div class="modal-body">
	        <div class="btn-group" data-toggle="buttons">
	        <table>
	        <tr>
	        <td id="blastp--DBSCAN-SWA_%s" style="font-size: 20px;border:2px solid #D8F2AF;text-align: center;" onclick="show(this.id)">
	         <button type="button" class="btn btn-default" style="vertical-align: center;font-size: 15px;padding-top: 10px;padding-bottom: 10px;">
	           blastp results
	        </button>
	        </td>       
	          <td id="blastn--DBSCAN-SWA_%s" style="font-size: 15px;border:2px solid #D8F2AF;text-align: center;" onclick="show(this.id)">
	          <button type="button" class="btn btn-default" style="vertical-align: center;font-size: 15px;padding-top: 10px;padding-bottom: 10px;">
	         blastn results
	        </button>
	        </td>
	       
	        </tr>
	        
	      </table>
	    </div>'''%(str(counter),str(counter),bac_id,hit_phage_id,str(counter),str(counter)))
					f_save.write('''
			<div id="blastp--DBSCAN-SWA_%s_result">
	        <table class="display" class="table table-striped table-hover table-bordered" style="text-align: center;">
	          <thead>
	          <tr>
	            <th>Prophage_Protein</th>
	            <th>Hit_Phage_Protein</th>
	            <th>Identity(%%)</th>
	            <th>Coverage</th>
	            <th style="word-break: keep-all;white-space:nowrap;">E-value</th>
	          </tr>
	          </thead>
	          <tbody>
						'''%(str(counter)))
					for blastp_info in phis[3:]:
						if '#blastn' in blastp_info:
							break
						blastp_info = blastp_info.strip().split('\t')
						f_save.write('''
			<tr>
	            <td style="word-break: break-all; word-wrap:break-word;">%s</td>
	            <td style="word-break: break-all; word-wrap:break-word;">%s</td>
	            <td>%s</td>
	            <td>%s</td>
	            <td style="word-break: keep-all;white-space:nowrap;">%s</td>
	        </tr>
							'''%(blastp_info[1],blastp_info[3],blastp_info[4],blastp_info[5],blastp_info[6]))
					f_save.write('''
							</tbody>
	       				 </table>
	      				</div>
							''')
					f_save.write('''
			<div id="blastn--DBSCAN-SWA_%s_result">      
	        <table class="display" class="table table-striped table-hover table-bordered">
	          <thead>
	          <tr>
	            <th>Hit_Prophage_Region</th>
	            <th>Hit_Phage_Region</th>
	            <th>Alignment Length</th>
	            <th>Identity</th>
	            <th>E-value</th>
	          </tr>
	          </thead>
	          <tbody>
						'''%(str(counter)))
					for blastn_info in '\n'.join(phis).split('#blastn')[-1].split('\n'):
						if blastn_info.strip()!='':
							blastn_info = blastn_info.strip().split('\t')
							f_save.write('''
			<tr>
	            <td style="word-break: break-all; word-wrap:break-word;">%s|%s</td>
	            <td>%s|%s</td>
	            <td>%s</td>            
	            <td>%s</td>
	            <td style="word-break: keep-all;white-space:nowrap;">%s</td>
	        </tr>
				'''%(blastn_info[1],blastn_info[2],hit_phage_id,blastn_info[4],blastn_info[5],blastn_info[6],blastn_info[7]))
					f_save.write('''
							</tbody>
	       				 </table>
	      				</div>
							''')
					f_save.write('''
	      </div>      
	      <div class="modal-footer">
	        <button type="button" class="btn btn-default" data-dismiss="modal">close
	        </button>        
	      </div>
	    </div><!-- /.modal-content -->
	  </div><!-- /.modal -->
	</div>
	</td>
	</tr>
						''')
		f_save.write('''
			</tbody>
			</table>
			</div>
			''')
	f_save.write('''
		</body>
		</html>''')

def predict_prophage(strain_id,inputfile_bac,outdir,add_annotation,prefix,diamond_thread_num):
	start = time.time()
	
	#step1:get the proteins in bacterial genome,annotate genome using prokka if in fasta format, else extract proteins
	print('step1:extract or predict the proteins in bacterial genome')
	protein_dir = os.path.join(outdir,'protein_annotation')
	mkdir(protein_dir)
	print(protein_dir)
	get_bac_protein(inputfile_bac,protein_dir,type,prefix)
	
	bac_protein_file = os.path.join(protein_dir,prefix+'_protein.faa')
	bac_protein_annotation_file = os.path.join(protein_dir,prefix+'_protein_def')
	
	#step2:identify phage or phage-like genes
	print('step2:identify phage or phage-like genes')
	phage_gene_annotation_dir = os.path.join(outdir,'phage_gene_cluster')
	mkdir(phage_gene_annotation_dir)
	phage_like_protein_list = get_phage_like_gene(bac_protein_file,phage_gene_annotation_dir,prefix,diamond_thread_num)
	bac_protein_blastp_phage_db_file = os.path.join(phage_gene_annotation_dir,prefix+'_blastp_uniprot.txt')
	
	#step3: execute dbscan,SWA
	print('step3:predict prophage regions by DBSCAN and SWA')
	m1 = threading.Thread(target=dbscan,args=(phage_like_protein_list,phage_gene_annotation_dir,prefix))
	m2 = threading.Thread(target=predict_prophage_swa,args=(bac_protein_file,bac_protein_blastp_phage_db_file,phage_gene_annotation_dir,60,prefix))
	m1.start()
	m2.start()
	m1.join()
	m2.join()
	dbscan_region_file = os.path.join(phage_gene_annotation_dir,prefix+'_dbscan_cluster_phage_protein.txt')
	swa_region_file = os.path.join(phage_gene_annotation_dir,prefix+'_swa_cluster_phage_protein.txt')
	
	#step4:merge prophage region
	print('step4:predict prophage regions based on step3')
	prophage_outdir = os.path.join(outdir,'prophage_region')
	mkdir(prophage_outdir)
	bac_fna_file = os.path.join(protein_dir,prefix+'.fna')
	predict_prophage_dbscan_swa(strain_id,bac_fna_file,bac_protein_file,dbscan_region_file,swa_region_file,bac_protein_annotation_file,prophage_outdir,prefix)
	prophage_region_file = os.path.join(prophage_outdir,prefix+'_prophage_region.txt')
	prophage_region_detail_file = os.path.join(prophage_outdir,prefix+'_prophage_region_details.txt')
	
	save_prophage_protein_file = os.path.join(prophage_outdir,prefix+'_prophage.faa')
	save_prophage_nucl_file = os.path.join(prophage_outdir,prefix+'_prophage.fna')

	#step5:predict attachment sites and boundary locating
	print('step5:predict attachment sites and determine the boundaries of the prophage region')
	att_outdir = os.path.join(prophage_outdir,'att_detection')
	mkdir(att_outdir)
	save_prophage_file = os.path.join(outdir,prefix+'_DBSCAN-SWA_prophage.txt')
	save_prophage_summary_file = os.path.join(outdir,prefix+'_DBSCAN-SWA_prophage_summary.txt')
	save_prophage_nucl_file = os.path.join(outdir,prefix+'_DBSCAN-SWA_prophage.fna')
	save_prophage_protein_file = os.path.join(outdir,prefix+'_DBSCAN-SWA_prophage.faa')
	identify_att(strain_id,bac_fna_file,bac_protein_file,prophage_region_detail_file,att_outdir,prefix)
	decide_boundary(strain_id,bac_fna_file,bac_protein_file,prophage_region_detail_file,att_outdir,save_prophage_file,save_prophage_summary_file,save_prophage_protein_file,save_prophage_nucl_file,bac_protein_blastp_phage_db_file,bac_protein_annotation_file,prefix)
	
	#step7:annotate prophage region
	prophage_region_annotate_dir = os.path.join(outdir,'prophage_annotation')
	mkdir(prophage_region_annotate_dir)
	if add_annotation!='none':
		mkdir(prophage_region_annotate_dir)
		global phage_inf_dict,phage_type
		if add_annotation=='PGPD':
			phage_type = 'fasta'
			phage_inf_dict = load_dict(os.path.join(root_path,'db','profiles','phage_inf_dict.txt'))
		else:
			phage_inf_dict,phage_type = get_inf(add_annotation,prophage_region_annotate_dir,prefix)
		print('step7:start to annotate the predicted prophage regions')
		annotate_prophage_region(strain_id,save_prophage_file,save_prophage_protein_file,save_prophage_nucl_file,prophage_region_annotate_dir,add_annotation,prefix)
	else:
		print('No more annoation!')
	
	#step8:visualizations
	print('start to visualize the results!')
	templeate_file = os.path.join(root_path,'static',"templeate.html")
	static_dir = os.path.join(root_path,'static')
	command = "cp -r "+static_dir+" "+outdir
	os.system(command)
	save_html_file = os.path.join(outdir,prefix+'_prophage_visualization.html')
	visualize(save_prophage_file,save_prophage_nucl_file,save_prophage_protein_file,prophage_region_annotate_dir,add_annotation,templeate_file,save_html_file,prefix)
	
	end = time.time()
	run_time = str(float((end-start))/60)
	print('Running time: %s minutes! Finished prediction in %s'%(run_time,outdir))

class MyThread(threading.Thread):
	def __init__(self,strain_id,inputfile_bac,c_outdir,add_annotation,c_prefix,diamond_thread_num):
		threading.Thread.__init__(self)
		self.strain_id = strain_id
		self.inputfile_bac = inputfile_bac
		self.c_outdir = c_outdir
		self.add_annotation = add_annotation
		self.c_prefix = c_prefix
		self.diamond_thread_num = diamond_thread_num
	def run(self):
		predict_prophage(self.strain_id,self.inputfile_bac,self.c_outdir,self.add_annotation,self.c_prefix,self.diamond_thread_num)

def getFaaFromGB(gb_file_path,faa_file_path):
	records = SeqIO.parse(gb_file_path, "gb")
	if os.path.exists(faa_file_path):
		os.remove(faa_file_path)
	savefile = open(faa_file_path, 'w')
	for record in records:
		fileID = record.id
		prot_num = 0
		for feature in record.features:
			if feature.type == 'CDS':
				prot_num = prot_num + 1
				location = feature.location
				if str(location).find('+') != -1:
					direction = '+'
				elif str(location).find('-') != -1:
					direction = '-1'
				if feature.type == 'CDS':
					if 'product' in feature.qualifiers:
						product = feature.qualifiers['product'][0]
						if ' ' in product:
							product = product.replace(' ','_')
						else:
							product = 'unkown'
					else:
						product = 'unknown'
				else:
					product = 'unkown'
				if 'protein_id' in feature.qualifiers:
					proteinId = feature.qualifiers['protein_id'][0]
				else:
					if 'inference' in feature.qualifiers:
						strInference = str(feature.qualifiers['inference'])
						if 'RefSeq' in strInference:
							proteinId = strInference.split('RefSeq:')[1].rstrip(']').rstrip('\'')
						elif 'SwissProt' in strInference:
							proteinId = strInference.split('SwissProt:')[1].rstrip(']').rstrip('\'')
						else:
							proteinId = 'unknown'
					else:
					    proteinId = 'unknown'
				if 'translation' in feature.qualifiers:
					translation = feature.qualifiers['translation'][0]
				else:
					translation = 'unkown'
				savefile.write('>ref' + '|'+fileID+'|'+proteinId+'|'+str(location)+'|'+str(product)+ '|'+str(prot_num)+'\n')
				if translation[-1] == '\n':
					savefile.write(translation)
				else:
					savefile.write(translation + '\n')
	savefile.close()

def getFnaFromGB(gb_file_path,fa_file_path):
	handle = open(gb_file_path)
	if os.path.exists(fa_file_path):
		os.remove(fa_file_path)
	SeqIO.convert(handle, 'genbank', fa_file_path, 'fasta')

def split_genbank(gb_file,outdir,prefix):
	strain_info_dict = {}
	
	multi_dir = os.path.join(outdir,'results')
	records = SeqIO.parse(gb_file, "gb")
	strain_proteins = ''
	counter = 0
	for record in records:
		strain_id = record.id
		strain_def = record.description
		strain_info_dict.update({strain_id.split()[0].split('.')[0]:strain_def})
		c_outdir = os.path.join(multi_dir,prefix+'_'+str(counter))
		c_prefix = prefix+'_'+str(counter)
		mkdir(c_outdir)
		contig_file = os.path.join(c_outdir,prefix+'_'+str(counter)+'.gb')
		SeqIO.write(record,contig_file, "gb")
		contig_faa_file = os.path.join(c_outdir,prefix+'_'+str(counter)+'.faa')
		contig_fa_file = os.path.join(c_outdir,prefix+'_'+str(counter)+'.fa')
		if not os.path.exists(contig_faa_file):
			getFaaFromGB(contig_file, contig_faa_file)
		if not os.path.exists(contig_fa_file):
			(contig_file, contig_fa_file)
		counter = counter+1
	return strain_info_dict


if __name__=='__main__':
	multiprocessing.freeze_support()
	parser = argparse.ArgumentParser()
	parser.add_argument('--input',help='Path of the input file.\n')
	parser.add_argument('--output', help='Path of the output folder.\n')
	parser.add_argument('--prefix', help='prefix of the prediction files.\n')	
	parser.add_argument('--add_annotation', help='optional,1.PGPD,2.phage_path:specified phage genome to detect whether the phage infects the query bacteria,3.none.\n')
	parser.add_argument('--evalue', help='optional,evalue cutoff when anotating .\n')
	parser.add_argument('--per', help='optional,homology percent cutoff when anotating .\n')
	parser.add_argument('--idn', help='optional,blastn identity cutoff when anotating .\n')
	parser.add_argument('--cov', help='optional,blastn coverage cutoff when anotating .\n')
	parser.add_argument('--protein_number', help='optional,the protein number of expanding when finding prophage boundaries .\n')
	parser.add_argument('--min_protein_num', help='optional,the protein number of expanding when finding prophage boundaries .\n')
	parser.add_argument('--thread_num', help='optional,the number of threads for multiple contigs.\n')
	parser.add_argument('--diamond_thread_num', help='optional,the number of threads for diamond blastp.\n')
	parser.add_argument('--h',help='''print help information''')
	args = parser.parse_args()
	global root_path
	if args.h:
		print(inputinstructions())
		sys.exit(-1)
	if args.input:
		inputfile_bac = args.input
	else:
		print('please input bacterial genome in fasta or GenBank format!')
		print(inputinstructions())
		sys.exit(1)
	
	if args.output:
		outdir = args.output
	else:
		print('please input the output directory path')
		print(inputinstructions())
		sys.exit(1)
	if not os.path.exists(outdir):
		mkdir(outdir)
	global prefix
	if args.prefix:
		prefix = args.prefix
	else:
		prefix = 'bac'
		print('Warning:The program will specify the prefix of the ouput files bac!If you want to use another prefix,please use the parameter:--prefix <prefix>!')
	
	if args.add_annotation:
		add_annotation = args.add_annotation
	else:
		print('Warning:The parameter add_annotation default:PGPD,If you want to run fast, please set add_annotation to none')
		add_annotation = 'PGPD'
	
	global blastp_evalue
	if args.evalue:
		blastp_evalue = args.evalue
		if float(blastp_evalue)>1e-7:
			print('Warning:The evalue %s can not be larger than 1e-4!'%blastp_evalue)
			blastp_evalue = 1e-7
			print('The evalue of homology search for virus uniprot database has been set to 1e-4')
	else:
		blastp_evalue = 1e-7
	global per,idn,cov
	if args.per:
		per = args.per
	else:
		per = 30
	if args.idn:
		idn = args.idn
	else:
		idn = 70
	
	if args.cov:
		cov = args.cov
	else:
		cov = 30
	thread_num = 10
	if args.thread_num:
		thread_num = args.thread_num
	else:
		thread_num = 10
    global diamond_thread_num
	if args.diamond_thread_num:
		diamond_thread_num = args.diamond_thread_num
	else:
		diamond_thread_num = 20
	
	global att_pro_num
	if args.protein_number:
		att_pro_num = args.protein_number
		if int(att_pro_num)<1:
			print('Warning:the protein number of expanding can not be smaller than 1')
			att_pro_num = 10
			print('the protein number of expanding has been set to 10')
	else:
		att_pro_num = 10
	global min_protein_num
	if args.min_protein_num:
		min_protein_num = args.min_protein_num
		if int(min_protein_num)<6:
			print('Warning:the protein number of shaping a cluster can not be smaller than 6')
			min_protein_num = 6
			print('the protein number of expanding has been set to 6')
	else:
		min_protein_num = 6
	root_path = dirname(dirname(abspath(__file__)))
	database = os.path.join(root_path,'db')
	flag = 0
	if os.path.exists(database):
		if ('profiles' in os.listdir(database)) and ('database' in os.listdir(database)):
			flag = 1
	if flag ==0:
		print('You are running DBSCAN-SWA first time, the database will be downloaded in %s'%database)
		command = "wget -P %s http://www.microbiome-bigdata.com/PHISDetector/static/download/DBSCAN-SWA/db.tar.gz"%database
		os.system(command)
		gz_db_file = os.path.join(database,'db.tar.gz')
		if os.path.exists(gz_db_file):
			command = "tar -zxvf %s -C %s"%(gz_db_file,root_path)
			os.system(command)
			command = "rm -f "+gz_db_file
			os.system(command)
		else:
			print('download error!please download from http://www.microbiome-bigdata.com/PHISDetector/static/download/DBSCAN-SWA/db.tar.gz')
			sys.exit(-1)
	
	global strain_inf_dict,type
	strain_inf_dict,type = get_inf(inputfile_bac,outdir,prefix)
	if len(list(strain_inf_dict.keys()))==1:
		strain_id = list(strain_inf_dict.keys())[0]
		predict_prophage(strain_id,inputfile_bac,outdir,add_annotation,prefix,diamond_thread_num)
	else:
		counter = 0
		tsk = []
		multi_dir = os.path.join(outdir,'results')
		mkdir(multi_dir)
		print(strain_inf_dict)
		for strain_id in strain_inf_dict.keys():
			c_outdir = os.path.join(multi_dir,prefix+'_'+str(counter))
			c_prefix = prefix+'_'+str(counter)
			mkdir(c_outdir)
			print(counter)
			c_inputfile_bac = os.path.join(c_outdir,prefix+'_'+str(counter)+'.gb')
			c_inputfile_bac_protein = os.path.join(c_outdir,prefix+'_'+str(counter)+'.faa')
			if not os.path.exists(c_inputfile_bac):
				c_inputfile_bac = os.path.join(c_outdir,prefix+'_'+str(counter)+'.fna')
			if not os.path.exists(c_inputfile_bac_protein):
				c_inputfile_bac = os.path.join(c_outdir,prefix+'_'+str(counter)+'.fna')
			else:
				if os.path.getsize(c_inputfile_bac_protein)==0:
					c_inputfile_bac = os.path.join(c_outdir,prefix+'_'+str(counter)+'.fna')
			# with open(c_inputfile_bac,'w') as f:
			# 	f.write('>'+strain_id+' '+strain_inf_dict[strain_id]+'\n'+strain_sequence_dict[strain_id]+'\n')
			pro_file = os.path.join(c_outdir,'protein_annotation',c_prefix+'_protein.faa')
			if os.path.getsize(c_inputfile_bac)>0:
				tsk.append(MyThread(strain_id,c_inputfile_bac,c_outdir,add_annotation,c_prefix,diamond_thread_num))			
			counter = counter+1
		for t in tsk:
			t.start()
			#time.sleep(1)
			while True:
				if (len(threading.enumerate()) <= int(thread_num)):
					break
		for t in tsk:
			t.join()
		save_file = os.path.join(outdir,prefix+'_DBSCAN-SWA_prophage_summary.txt')
		f_save = open(save_file,'w')
		f_save.write('prefix\tbacteria_id\tbacteria_def\tgenome_size\tprophage_start\tprophage_end\tkey_proteins\tbest_hit_species\tCDS_number\tattl_region\tattr_region\n')
		f_save.flush()
		save_file_detail = os.path.join(outdir,prefix+'_DBSCAN-SWA_prophage.txt')
		f_save_detail = open(save_file_detail,'w')
		f_save_detail.write('''The following contents displays predicted prophage regions
first line of each prophage describes the prophage information and the following lines describe the proteins and homology proteins in uniprot database
contig_id\tcontig_def\tgenome_size\tprophage_start\tprophage_end\tkey_proteins\tbest_hit_species\tCDS_number\tattl_region\tattr_region\n
''')
		save_file_pro = os.path.join(outdir,prefix+'_DBSCAN-SWA_prophage.faa')
		f_save_pro = open(save_file_pro,'w')
		save_file_nucl = os.path.join(outdir,prefix+'_DBSCAN-SWA_prophage.fna')
		f_save_nucl = open(save_file_nucl,'w')
		
		for counter in range(0,len(list(strain_inf_dict.keys()))):
			c_outdir = os.path.join(multi_dir,prefix+'_'+str(counter))
			c_prefix = prefix+'_'+str(counter)
			c_prophage_file = os.path.join(c_outdir,c_prefix+'_DBSCAN-SWA_prophage_summary.txt')
			if os.path.exists(c_prophage_file):
				with open(c_prophage_file) as f:
					c_contents = f.readlines()
				if len(c_contents)>1:
					for line in c_contents[1:]:
						f_save.write(prefix+'\t'+line.strip()+'\n')
						f_save.flush()
					c_detail_file = os.path.join(c_outdir,c_prefix+'_DBSCAN-SWA_prophage.txt')
					with open(c_detail_file) as f:
						c_contents_detail = f.read().strip().split('>prophage')
					for prophage_detail in c_contents_detail[1:]:
						f_save_detail.write('>prophage '+str(counter)+'_'+prophage_detail.split('\n')[0].strip()+'\n')
						f_save_detail.write('\n'.join(prophage_detail.split('\n')[1:]).strip()+'\n')
						f_save_detail.flush()
					c_faa_file = os.path.join(c_outdir,c_prefix+'_DBSCAN-SWA_prophage.faa')
					c_fna_file = os.path.join(c_outdir,c_prefix+'_DBSCAN-SWA_prophage.fna')
					with open(c_faa_file) as f:
						c_contents_faa = f.read().strip()
					with open(c_fna_file) as f:
						c_contents_fna = f.read().strip()
					f_save_pro.write(c_contents_faa+'\n')
					f_save_pro.flush()
					f_save_nucl.write(c_contents_fna+'\n')
					f_save_nucl.flush()	
		f_save.close()
		f_save_detail.close()
		f_save_pro.close()
		f_save_nucl.close()
