# beav - a bacterial genome and mobile element annotation pipeline
beav: Bacteria/Element Annotation reVamped

**beav** is a command line tool that streamlines bacterial genome and mobile genetic element annotation. It combines multiple annotation tools, automating the process of running, parsing, and combining the results into a single easy-to-read output. Annotated features include secretion systems, anti-phage defense systems, integrative & conjugative/mobilizable elements, integrons, prophage regions, amino acid biosynthesis pathways, small carbon metabolite catabolism pathways, and biosynthetic gene clusters.

The **beav** pipeline also includes several tools and databases that enhance the annotation of plant associated microbes, including phytopathogens and symbionts. Custom bakta databases provide correct gene names and annotations for phytopathogen virulence genes, effectors, and genes important for mutualist symbiosis. Other tools annotate promoter elements such as . 

An optional <i>Agrobacterium</i>-specific pipeline identifies the presence of Ti and Ri plasmids and classifies them under the Weisberg et al. 2020 scheme. It also annotates Ti/Ri plasmid elements including T-DNA borders, overdrive, virbox, trabox, and other binding sites, and determines the biovar and genomospecies of the input strain. Virulence and T-DNA genes, including opine synthase and transport/catabolism loci, are also correctly named and annotated.

# **Installation**
The **beav** pipeline requires a number of programs and databases be installed. Therefore, it is highly encouraged and recommended to use conda to install **beav** and all of its dependencies.

**Prerequisites:
Prior to installing beav, the usearch program must be installed and present in the environment and PATH variable. It can be downloaded from https://www.drive5.com/usearch/**

# From conda (Recommended) 
It is recommended to use either conda with libmamba or mamba to install beav as this will greatly speed up the time solving the environment.

instructions for conda:
```
conda create -n beav
conda env update -n beav beav
conda install beav
```
alternative instructions using mamba:
```
conda create -n beav
mamba env update -n beav beav
```

The conda environment can then be activated using:
```
conda activate beav
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

# Install all databases 


```
conda activate beav 
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

Required if running TIGER. Users must provide a path to a blast database of reference genomes using the --tiger_blast_database option. 

**--bakta_arguments**

Additional arguments can be passed to bakta using the --bakta_arguments option.

**--agrobacterium**

The --agrobacterium option activates an additional pipeline to provide agrobacterium-specific annotation. 

**--skip-PROGRAM**

The skip options allow for specified programs to be skipped if the annotation is not needed or required programs are not installed. 

# Examples
**Minimal run**

```
beav --input /path/to/file/test.fna --threads 8 --skip_tiger
```

**Standard run**
```
beav --input /path/to/file/test.fna --threads 8 --tiger_blast_database /path/to/databases/blast/refseq_genomic.fna
```

**Complex run**

```
beav --input /path/to/file/test.fna --threads 8 --bakta_arguments --db /path/to/alternative-data-bases/bakta-1.7/ --tiger_blast_database /path/to/databases/blast/allagro.fna --agrobacterium --skip_integronfinder
```



