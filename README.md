# beav - a bacterial genome and mobile element annotation pipeline
beav: Bacteria/Element Annotation reVamped

**beav** is a command line resource that steam-lines bacteria genome annotation by parsing the output of multiple existing prediction software and organizing these regions of interest into a final easy-to-read output. This tool takes an integrative approach to automate and provide comprehensive functional genome annotations. **Annotated features include secretion systems, integrons, anti-phage defense systems, integrative conjugative elements, amino acid biosynthesis, and carbon metabolism pathways, prophage, biosynthetic gene clusters. The agrobacterium option annotates biovar, Ti/Ri plasmids, T-DNA borders, virboxes, traboxes, pipboxes, nodboxes, hrpboxes, and ttsboxes**. 

# **Installation**
It is highly encouraged and recommended to use conda to install BEAV.

**Prerequisites:
Usearch must be installed and present in the environment: https://www.drive5.com/usearch/**

# From conda (Recommended) 
```
conda create -n beav
conda activate beav
conda install beav
```
# Alternative: From source

Prerequisites: 

Bakta

IntegronFinder

MacSyFinder

DefenseFinder

TIGER2

GapMind

DBSCAN-SWA

antiSMASH

EMBOSS (Fuzznuc)

HMMER (Nhmmer)
```
git clone https://github.com/weisberglab/beav.git
```

**If installed from source, DBSCAN-SWA, TIGER2, and PaperBLAST needs to be installed in the software folder. Then an environment variable needs to be set to BEAV_DIR and point to the software folder.**

# Install databases (Required for both)


```
beav_db
```

# Usage
```
usage: beav [--input INPUT] [--output OUPUT_DIRECTORY] [--strain STRAIN] [--bakta_arguments BAKTA_ARGUMENTS] [--tiger_arguments TIGER_ARGUMENTS][--agrobacterium AGROBACTERIUM] [--skip_macsyfinder] [--skip_integronfinder][--skip_defensefinder] [--skip_tiger] [--skip_gapmind][--skip_dcscan-swa] [--skip_antismash] [--help] [--threads THREADS]
    BEAV- Bacterial Element Annotation reVamped
    Input/Output: 
        --input, -i
                Input file in fasta nucleotide format (Required)
        --output
                Output directory (default: current working directory)
        --strain
                Strain name (default: input file prefix)
        --bakta_arguments
                Additional arguments and database options specific to Bakta 
        --antismash_arguments
                Additional arguments and database options specific to antiSMASH (Default: \"$antismash_args\") 
        --tiger_blast_database
                Path to a reference genome blast database for TIGER2 ICE analysis (Required unless --skip_tiger is used)
    Options:
        --agrobacterium
                Agrobacterium specific tools that identify biovar/species group, Ti/Ri plasmid, T-DNA borders, virboxes and traboxes
        --skip_macsyfinder
                Skip detection and annotation of secretion systems
        --skip_integronfinder
                Skip detection and annotation of integrons 
        --skip_defensefinder
                Skip detection and annotation of anti-phage defense systems 
        --skip_tiger
                Skip detection and annotation of integrative conjugative elements (ICEs)
        --skip_gapmind
                Skip detection of amino acid biosynthesis and carbon metabolism pathways
        --skip_dbscan-swa
                Skip detection and annotation of prophage
        --skip_antismash
                Skip detection and annotation of biosynthetic gene clusters
    General:
        --help, -h
                Show BEAV help message
        --threads, -t
                Number of CPU threads
```
# Options

**--antismash_arguments**

Additional antiSMASH arguments can be input into antiSMASH using the --antismash_arguments option. This allows for full usage of antiSMASH and additional databases.

**--tiger_blast_database**

In order to run TIGER2, users must provide a path to a blast database of reference genomes using the --tiger_blast_database option. 

**--bakta_arguments**

Additional Bakta arguments can be input into Bakta using the --bakta_arguments option.  This allows for full usage of Bakta and additional databases.

**--agrobacterium**

The --agrobacterium option uses fuzznuc pattern matching, HMMER models, and several unique scripts to provide agrobacterium-specific annotation. 

**--skip-PROGRAM**

The skip options allow for specified programs to be skipped if the annotation is not needed or required programs are not installed. 

# Examples
**Minimum run**

```
beav --input /path/to/file/test.fna --threads 8 --skip_tiger
```

**Complex run**

```
beav --input /path/to/file/test.fna --threads 8 --bakta_arguments --db /path/to/alternative-data-bases/bakta-1.7/ --tiger_blast_database /path/to/databases/blast/allagro.fna --agrobacterium --skip_integronfinder
```



