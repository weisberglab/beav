#!/bin/bash
echo -e "downloading and setting up databases for dependencies"

echo -e "BAKTA"

baktadbpath="$1"
baktadbtype="$2"

#check if bakta db variable exists and db is downloaded
if [[ -z "$BAKTA_DB" ]] && [[ -f $BAKTA_DB/version.json ]]; then
	echo -e "BAKTA_DB already set and installed"
else
	if [[ ! -z $baktadbpath ]]; then
		echo -e "Please ensure that the first argument to this script is the path to store the bakta db"
		exit 1
	elif [[ ! ($baktadbtype == "full" || $baktadbtype == "light") ]]; then
		echo -e "Please ensure that the second argument to this script is either full or light"
		exit 1
	fi
	bakta_db download --output $baktadbpath --type $baktadbtype
fi

echo -e ""

#macsyfinder
echo -e "MACSYFINDER"
macsydata install --user --upgrade TXSScan
echo -e ""

#defense_finder
echo -e "DEFENSE FINDER"

defense-finder update


echo -e ""
echo -e "ANTISMASH"
echo -e "not for now"
#download-antismash-databases

echo -e ""
#TIGER2
echo -e "TIGER2"
if [[ ! -f $BEAV_DIR/software/TIGER/db/Pfam-A.hmm ]]; then
	echo -e "requirement Pfam-A.hmm file does not exist"
	echo -e "Downloading and extracting Pfam-A hmms"
	curl -L -O https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz --output $BEAV_DIR/software/TIGER/db/Pfam-A.hmm.gz
	gunzip $BEAV_DIR/software/TIGER/db/Pfam-A.hmm.gz
else
	echo -e "Pfam-A.hmm already exists"	
fi

echo -e ""
echo -e "DBSCAN-SWA"
echo -e "DBSCAN-SWA will download a database (~600 Mb) on first run"

echo -e ""
#GapMind
echo -e "GAPMIND"
#make user install PaperBLAST in BEAV_DIR/software/
#then run script to download and format dbs
#gapmind aa carbon and format

path_to_usearch=$(which usearch)
#if [[ -f $BEAV_DIR/software/PaperBLAST/bin/usearch ]]; then
#	echo -e "usearch is already added to the correct directory"
if [[ ! -x $path_to_usearch ]]; then 
	echo -e "Please download the usearch executable from https://www.drive5.com/usearch/download.html"
	echo -e "Extract this file, rename it to usearch, make it executable (chmod +x usearch) and place it in a location in your PATH before running beav"
	echo -e "Then re-run this setup script to format GapMind databases"
	exit 1
fi

mkdir $BEAV_DIR/software/PaperBLAST/fbrowse_data
mkdir $BEAV_DIR/software/PaperBLAST/private
mkdir $BEAV_DIR/software/PaperBLAST/tmp
mkdir $BEAV_DIR/software/PaperBLAST/tmp/path.aa
mkdir $BEAV_DIR/software/PaperBLAST/tmp/path.carbon

curl -L -O https://papers.genomics.lbl.gov/tmp/path.aa/curated.faa --output $BEAV_DIR/software/PaperBLAST/tmp/path.aa/curated.faa
curl -L -O https://papers.genomics.lbl.gov/tmp/path.aa/curated.db --output $BEAV_DIR/software/PaperBLAST/tmp/path.aa/curated.db
curl -L -O https://papers.genomics.lbl.gov/tmp/path.aa/steps.db --output $BEAV_DIR/software/PaperBLAST/tmp/path.aa/steps.db
curl -L -O https://papers.genomics.lbl.gov/tmp/path.carbon/curated.faa --output $BEAV_DIR/software/PaperBLAST/tmp/path.carbon/curated.faa
curl -L -O https://papers.genomics.lbl.gov/tmp/path.carbon/curated.db --output $BEAV_DIR/software/PaperBLAST/tmp/path.carbon/curated.db
curl -L -O https://papers.genomics.lbl.gov/tmp/path.carbon/steps.db --output $BEAV_DIR/software/PaperBLAST/tmp/path.carbon/steps.db

$BEAV_DIR/software/PaperBLAST/bin/extractHmms.pl $BEAV_DIR/software/PaperBLAST/tmp/path.aa/steps.db $BEAV_DIR/software/PaperBLAST/tmp/path.aa
$BEAV_DIR/software/PaperBLAST/bin/extractHmms.pl $BEAV_DIR/software/PaperBLAST/tmp/path.carbon/steps.db $BEAV_DIR/software/PaperBLAST/tmp/path.carbon

usearch -makeudb_ublast $BEAV_DIR/software/PaperBLAST/tmp/path.aa/curated.faa -output $BEAV_DIR/software/PaperBLAST/tmp/path.aa/curated.faa.udb
formatdb -p T -o T -i $BEAV_DIR/software/PaperBLAST/tmp/path.aa/curated.faa

usearch -makeudb_ublast $BEAV_DIR/software/PaperBLAST/tmp/path.aa/curated.faa -output $BEAV_DIR/software/PaperBLAST/tmp/path.aa/curated.faa.udb
formatdb -p T -o T -i $BEAV_DIR/software/PaperBLAST/tmp/path.aa/curated.faa


echo -e ""
echo -e "DONE"
echo -e ""

echo -e "Dont forget to set the BAKTA_DB environment variable to point to $baktadbpath"
echo -e "alternatively, provide --db $baktadbpath as an argument to bakta"
