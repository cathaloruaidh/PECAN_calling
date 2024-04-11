#! /bin/bash


##################################################
# as previous script, collapse RD vs PR/SR
##################################################

SAMPLE=${RESULTS_DIR}/${SAMPLE_NAME}

rm -f ${LOG_DIR}/${1}.log




ADD=""
if [[ ${GENO} ]]
then
	ADD=".gq${GQ}"
fi


ALL_RAW=${SAMPLE}.ALL${ADD}.raw_details.txt
if
[[ -f ${ALL_RAW} ]]
then
	rm -f ${ALL_RAW}
	touch ${ALL_RAW}
fi

### BETWEEN METHODS MERGING

for TYPE in BND DEL DUP INS INV
do

	RD_ALL=${SAMPLE}.RD${ADD}.${TYPE}.bed
	RD_CUT=${SAMPLE}.RD${ADD}.${TYPE}.cut.bed
	PRSR_ALL=${SAMPLE}.PRSR${ADD}.${TYPE}.bed
	PRSR_CUT=${SAMPLE}.PRSR${ADD}.${TYPE}.cut.bed
	
	ALL_FILE=${SAMPLE}.ALL${ADD}.${TYPE}.bed
	ALL_CUT=${SAMPLE}.ALL${ADD}.${TYPE}.cut.bed


	cut -f 1-4 ${RD_ALL} > ${RD_CUT}
	cut -f 1-4 ${PRSR_ALL} > ${PRSR_CUT}
	

	if [[ ( ! -s ${RD_CUT} ) && ( ! -s ${PRSR_CUT} ) ]]
	then
		echo "Both RD and PRSR files are empty for ${TYPE}!" >> ${LOG_DIR}/${1}.log 
		touch ${ALL_FILE}

	# check if RD file is empty
	elif [[ ! -s ${RD_CUT} ]]
	then
		echo "The RD file is empty for ${TYPE}!" >> ${LOG_DIR}/${1}.log
		awk -v OFS=$'\t' '{$5=$4 ; gsub(/PRSR/, "ALL", $4) ; print $1,$2,$3,$4,$5}' ${PRSR_CUT} > ${ALL_FILE}
		
	# check if PRSR file is empty
	elif [[ ! -s ${PRSR_CUT} ]]
	then
		echo "The PRSR file is empty for ${TYPE}!" >> ${LOG_DIR}/${1}.log
		awk -v OFS=$'\t' '{$5=$4 ; gsub(/RD/, "ALL", $4) ; print $1,$2,$3,$4,$5}' ${RD_CUT} > ${ALL_FILE}

	else
		# Find all sites in A and B that overlap, print both coordinates
		echo "Temp Union for ${TYPE} - RD" >> ${LOG_DIR}/${1}.log
		bedtools intersect -a ${RD_CUT} -b ${PRSR_CUT} -f ${OVLP_FRAC} -r -wa -wb | \
			tr '\t' ' ' > ${ALL_FILE}.tmp


		# Get the intersection of PRSR and RD calls
		# For each line, merge the two regions - one from each caller
		echo "Union for ${TYPE} - RD" >> ${LOG_DIR}/${1}.log
		while IFS="" read -r LINE || [ -n "$LINE" ]
		do 
			REGION=$(echo $LINE | sed 's/ /\n/4' | tr ' ' '\t' | bedtools sort -i - | bedtools merge -i - -c 4 -o distinct)
			CALLS=$(echo ${REGION} | cut -f4)

			RD_TAG=$(echo ${CALLS} | tr ',' '\n' | grep 'RD' | cut -f3 -d_ | cut -f1-2 -d:)
			PRSR_TAG=$(echo ${CALLS} | tr ',' '\n' | grep 'PRSR' | cut -f3 -d_ | cut -f3-4 -d:)

			TAG="${RD_TAG}:${PRSR_TAG}"

			echo ${REGION} | awk -v OFS=$'\t' -v TAG=${TAG} -v TYPE=${TYPE} '{print $1,$2,$3,"ALL_" TYPE "_" TAG "_" $1 ":" $2 "-" $3,$4}'
		done < ${ALL_FILE}.tmp > ${ALL_FILE}
		
		
		# Original calls for each caller from union
		echo "Common Each ${TYPE} - RD and PRSR" >> ${LOG_DIR}/${1}.log
		bedtools intersect -a ${PRSR_CUT} -b ${RD_CUT} -f ${OVLP_FRAC} -r -wa > ${ALL_FILE}.PRSR.tmp 
		bedtools intersect -a ${RD_CUT} -b ${PRSR_CUT} -f ${OVLP_FRAC} -r -wa > ${ALL_FILE}.RD.tmp


		echo "Singletons ${TYPE} - RD" >> ${LOG_DIR}/${1}.log
		grep -v -f ${ALL_FILE}.RD.tmp ${RD_CUT} > ${RD_CUT}.no_ovlp.tmp
		while IFS="" read -r LINE || [ -n "$LINE" ]
		do 
			CALLS=$(echo ${LINE} | cut -f4)

			RD_TAG=$(echo ${CALLS} | tr ',' '\n' | grep 'RD' | cut -f3 -d_ | cut -f1-2 -d:)

			TAG="${RD_TAG}:LA:MA"

			echo ${LINE} | awk -v OFS=$'\t' -v TAG=${TAG} -v TYPE=${TYPE} '{print $1,$2,$3,"RD_" TYPE "_" TAG "_" $1 ":" $2 "-" $3,$4}'
		done < ${RD_CUT}.no_ovlp.tmp >> ${ALL_FILE}


		echo "Singletons ${TYPE} - PRSR" >> ${LOG_DIR}/${1}.log
		grep -v -f ${ALL_FILE}.PRSR.tmp ${PRSR_CUT} > ${PRSR_CUT}.no_ovlp.tmp
		while IFS="" read -r LINE || [ -n "$LINE" ]
		do 
			CALLS=$(echo ${LINE} | cut -f4)

			PRSR_TAG=$(echo ${CALLS} | tr ',' '\n' | grep 'PRSR' | cut -f3 -d_ | cut -f3-4 -d:)

			TAG="CA:EA:${PRSR_TAG}"

			echo ${LINE} | awk -v OFS=$'\t' -v TAG=${TAG} -v TYPE=${TYPE} '{print $1,$2,$3,"PRSR_" TYPE "_" TAG "_" $1 ":" $2 "-" $3,$4}'
		done < ${PRSR_CUT}.no_ovlp.tmp  >> ${ALL_FILE}


		
		bedtools sort -i ${ALL_FILE} > ${ALL_FILE}.sort
		mv ${ALL_FILE}.sort ${ALL_FILE}

		# cleanup
		rm ${SAMPLE}*.tmp  
	fi

	# Ensure that all regions are contained with the Reference Genome

	echo "Final sort, and done ${TYPE}" >> ${LOG_DIR}/${1}.log
	bedtools intersect -a <( cut -f 1-4 ${ALL_FILE} ) -b ${INCLUDE_FILE} | bedtools sort -i - > ${ALL_CUT}



	# get raw calls and genotyping details for each CNV

	if [[ -s ${ALL_FILE} ]]
	then 
		while IFS='' read -r LINE || [[ -n "${LINE}" ]] 
		do
			ID=$( echo ${LINE} | cut -f4 -d ' ' )
			CALLS=$( echo ${LINE} | cut -f5 -d ' ' )
			echo ${CALLS} | tr ',' '\n' | \
				grep -w -f - <( cat ${PRSR_ALL} ${RD_ALL} ) | cut -f5 | tr ',' '\n' | \
				grep -w -f - <( cat ${SAMPLE}.cnvpytor${ADD}.${TYPE}.collapse.bed ${SAMPLE}.erds${ADD}.${TYPE}.collapse.bed ${SAMPLE}.lumpy${ADD}.${TYPE}.collapse.bed ${SAMPLE}.manta${ADD}.${TYPE}.collapse.bed ) | cut -f5 | tr ',' '\n' | \
				grep -w -f - <( cat ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}*gt${ADD}.txt ) | \
				awk -v ID="${ID}" -v OFS="\t" '{ print ID, $0}' \
			>> ${ALL_RAW}
		done < ${ALL_FILE} >> ${ALL_RAW}
	fi


	echo "****************************************" >> ${LOG_DIR}/${1}.log
	
done


