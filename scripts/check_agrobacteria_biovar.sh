#!/usr/bin/env bash

strain=$1
cpus=$2

ls -1 $BEAV_DIR/databases/agrobacterium_taxa/*fna > genomelist

fastANI -t $cpus -q ${strain}.fna --rl genomelist -o ${strain}.agro_ani.out  
bestani=`cat ${strain}.agro_ani.out | awk '$3 >= 95' | cut -f 1-3 | sed 's/^.*\///g;s/\.fna//g' | sed 's/__/\t/g'`

touch ${strain}.agrobacteria_taxonomy.out

echo -e "strain	biovar	species_group	reference	%ANI" > ${strain}.agrobacteria_taxonomy.out
echo -e "$bestani" >> ${strain}.agrobacteria_taxonomy.out

if [[ -z "$bestani" ]]; then
	echo -e "No matching species in database, closest ANI hit:" >> ${strain}.agrobacteria_taxonomy.out
	head -n 1 ${strain}.agro_ani.out | cut -f 1-3 | sed 's/^.*\///g;s/\.fna//g' | sed 's/__/\t/g' >> ${strain}.agrobacteria_taxonomy.out
fi
