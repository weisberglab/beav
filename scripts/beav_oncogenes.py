#!/usr/bin/python3

# This script enlists Ti plasmid oncogenes to be used in beav_circos.py

# Dictionaries to define gene of interest in Ti Plasmid

# agrocinopine transport/catabolism
agrocinopine_dict = {
    'accA': 'agrocinopine ABC transporter, substrate binding protein accA',
    'accB': 'agrocinopine ABC transporter, nucleotide binding/ATPase accB',
    'accC': 'agrocinopine ABC transporter, nucleotide binding/ATPase accC',
    'accD': 'agrocinopine ABC transporter, membrane spanning protein accD',
    'accE': 'agrocinopine ABC transporter, membrane spanning protein accE',
    'accF_2': 'agrocinopine phosphodiesterase accF_2',
    'accF': 'agrocinopine phosphodiesterase accF',
    'accG_2': 'arabinose phosphate phosphatase accG_2',
    'accG': 'arabinose phosphate phosphatase accG',
    'accR': 'Repressor of the acc and arc operons accR'
}

# opine transport/catabolism
opine_tracat_dict = {
    'ocd': 'Ornithine cyclodeaminase',
    'nocQ': 'nopaline ABC transporter permease NocQ',
    'nocT': 'nopaline ABC transporter substrate-binding protein NocT',
    'nocP': 'nopaline ABC transporter ATP-binding protein NocP',
    'nocR': 'Regulatory protein NocR',
    'ooxA': 'Opine oxidase subunit A',
    'ooxB': 'Opine oxidase subunit B',
    'occT': 'Octopine-binding periplasmic protein',
    'occP': 'Octopine permease ATP-binding protein P',
    'occM': 'Octopine transport system permease protein OccM',
    'occQ': 'Octopine transport system permease protein OccQ',
    'occR': 'Octopine catabolism/uptake operon regulatory protein OccR'
}

# opine synthase genes
opine_synth_dict = {
    'acs': 'agrocinopine synthase',
    'ags': 'agropine synthase',
    'chsA': 'chrysopine synthase',
    'cus': 'cucumopine synthase',
    'mas1': 'mannopine synthase 1',
    'mas2': 'mannopine synthase 2',
    'nos': 'nopaline synthase',
    'ocs': 'octopine synthase',
    'sus': 'succinamopine synthase',
    'vis': 'vitopine synthase',
    'unk': 'unknown opine synthase'
}

# virulence genes
vir_dict = {
    'virA': 'vir locus protein VirA',
    'virB10': 'vir locus Type IV secretion system VirB10',
    'virB11': 'vir locus Type IV secretion system VirB11',
    'virB1': 'vir locus Type IV secretion system VirB1',
    'virB2': 'vir locus Type IV secretion system VirB2',
    'virB3': 'vir locus Type IV secretion system VirB3',
    'virB4': 'vir locus Type IV secretion system VirB4',
    'virB5': 'vir locus Type IV secretion system VirB5',
    'virB6': 'vir locus Type IV secretion system VirB6',
    'virB7': 'vir locus Type IV secretion system VirB7',
    'virB8': 'vir locus Type IV secretion system VirB8',
    'virB9': 'vir locus Type IV secretion system VirB9',
    'virC1': 'vir locus protein VirC1',
    'virC2': 'vir locus protein VirC2',
    'virD1': 'vir locus protein VirD1',
    'virD2': 'vir locus protein VirD2',
    'virD3': 'vir locus protein VirD3',
    'virD3': 'Protein VirD3',
    'virD4': 'vir locus protein VirD4',
    'virD5': 'vir locus protein VirD5',
    'virE1': 'vir locus protein VirE1',
    'virE2': 'vir locus protein VirE2',
    'virE3': 'vir locus protein VirE3',
    'virF': 'exported virulence protein virF',
    'virG': 'vir locus protein VirG',
    'virH1': 'P-450 monooxygenase virH1',
    'virH2': 'P-450 monoxygenase virH2',
    'virK': 'virA/G regulated gene virK',
    'GALLS': 'vir protein GALLS'
}

# plasmid conjugation tra genes
tra_dict = {
    'traA': 'conjugal transfer protein TraA',
    'traB': 'conjugal transfer protein TraB',
    'traC': 'conjugal transfer protein TraC',
    'traD': 'conjugal transfer protein TraD',
    'traF': 'conjugal transfer protein TraF',
    'traG': 'conjugal transfer protein TraG',
    'traH': 'conjugal transfer protein TraH',
    'traI': 'conjugal transfer protein TraI',
    'traM': 'conjugal transfer protein TraM',
    'traR': 'conjugal transfer protein TraR'
}

# plasmid conjugation trb genes
trb_dict = {
    'trbB': 'conjugal transfer protein TrbB',
    'trbC': 'conjugal transfer protein TrbC',
    'trbD': 'conjugal transfer protein TrbD',
    'trbE': 'conjugal transfer protein TrbE',
    'trbF': 'conjugal transfer protein TrbF',
    'trbG': 'conjugal transfer protein TrbG',
    'trbH': 'conjugal transfer protein TrbH',
    'trbI': 'conjugal transfer protein TrbI',
    'trbJ': 'conjugal transfer protein TrbJ',
    'trbK': 'conjugal transfer protein TrbK',
    'trbL': 'conjugal transfer protein TrbL'
}

# replication gene
rep_dict = {
    'repA': 'plasmid partitioning protein RepA',
    'repB': 'plasmid partitioning protein RepB',
    'repC': 'plasmid replication protein RepC'
}

# T-DNA/oncogene genes
oncogene_dict = {
    'iaaH': 'indoleacetamide hydrolase',
    'iaaM': 'Tryptophan 2-monooxygenase',
    'ipt': 'Adenylate dimethylallyltransferase',
    'tzs': 'Adenylate dimethylallyltransferase'
}

oncogene_list = [
    'Tryptophan 2-monooxygenase',
    'Adenylate dimethylallyltransferase',
    '6a protein',
    'Uncharacterized protein 6',
    '6b protein',
    "Gene 4' protein",
    'C protein',
    'D protein',
    'E protein',
    '5 protein',
    'RolB-RolC domain-containing protein',
    'T-DNA oncoprotein'
]
