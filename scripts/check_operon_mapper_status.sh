#!/usr/bin/env bash
URL=`tail -n1 operon-mapper_results_url`
submitstatus=`head -n1 operon-mapper_results_url`
jobnum=`echo -e "$URL" | sed 's/^.*out_//g;s/\.html$//g'`

#maximum wait time for operon-mapper results
max_wait_time=120

check_job () {
	response=$(curl -L -s -w "%{http_code}" $URL)
	#http_code=200
	http_code=$(tail -n1 <<< "$response")  # get the last line
	content=$(sed '$ d' <<< "$response")   # get all but the last line which contains the status code
	#content=`cat out_590595.html.done.html` 

	if [[ "$http_code" -ne 200 ]] ; then
		echo -e "Error contacting the operon-mapper server with error code ${http_code}."
		echo -e "skipping operon-mapper parsing."
		exit 0
	else
		#check if job is finished
		if echo "$content" | grep -qF 'Download' ; then
			#the job completed successfully
			echo -e "Operon-mapper webserver Job completed. Downloading results files."
			genepairsurl=`echo -e "$content" | grep 'operonic_gene_pairs_' | sed 's/^.*href="//g;s/" TARGET.*//g' | sed 's#^\.\.#http://biocomputo.ibt.unam.mx/operon_mapper#g'`
			operonsurl=`echo -e "$content" | grep 'list_of_operons_' | sed 's/^.*href="//g;s/" TARGET.*//g' | sed 's#^\.\.#http://biocomputo.ibt.unam.mx/operon_mapper#g'`
			cogsurl=`echo -e "$content" | grep 'predicted_COGs_' | sed 's/^.*href="//g;s/" TARGET.*//g' | sed 's#^\.\.#http://biocomputo.ibt.unam.mx/operon_mapper#g'`
			funcdescurl=`echo -e "$content" | grep 'functional_descriptions_' | sed 's/^.*href="//g;s/" TARGET.*//g' | sed 's#^\.\.#http://biocomputo.ibt.unam.mx/operon_mapper#g'`
			#echo -e "$genepairsurl"
			#echo -e "$operonsurl"
			#echo -e "$cogsurl"
			#echo -e "$funcdescurl"
				
			mkdir -p operon-mapper_results
			curl -s -L -o operon-mapper_results/operonic_gene_pairs $genepairsurl
			curl -s -L -o operon-mapper_results/list_of_operons $operonsurl
			curl -s -L -o operon-mapper_results/predicted_COGs $cogsurl
			curl -s -L -o operon-mapper_results/functional_descriptions $funcdescurl

			#parse the operon-mapper output and reformat table file
			if [[ -f operon-mapper_results/list_of_operons ]]; then
				tail -n+2 operon-mapper_results/list_of_operons | cut -f 1,2 | sed 's/\.[pr]01$//g' | awk -F'\t' 'BEGIN{OFS = FS; i=0;}{if ( NF == 2 ){print $2,i;}else{i="operon_"$1;}}' > operon-mapper_results/list_of_operons.table
			fi
			mv operon-mapper_results_url ./operon-mapper_results/
			mv testresponse.html ./operon-mapper_results/
			echo "Done operon-mapper parsing."
			exit 0
		elif echo "$content" | grep -iqF 'error' ; then
			echo -e "Error in operon-mapper run, skipping output parsing." 
			exit 0
		else
			echo "Job is still running. waiting 1 minute."
			sleep 1m
		fi
	fi
}


if [[ "$submitstatus" -ne 200 ]]; then
	echo -e "Initial operon-mapper job submission was not successful. Skipping."
	exit 0
fi

wait_count=0
while true; do
	check_job
	wait_count=$((wait_count+1))
	if [[ $wait_count -ge $max_wait_time ]]; then
		echo -e "Exceeded maximum job wait time ($max_wait_time minutes), skipping operon-mapper analysis."
		exit 0
	fi
done

