#! /bin/bash

##################################################
# collapse calls within caller types, so 
# (ERDS vs CNVpytor) and then (LUMPY vs Manta)
# (this is a little messy ... )
##################################################



SAMPLE=${RESULTS_DIR}/${SAMPLE_NAME}

rm -f ${LOG_DIR}/${1}.log



ADD=""
if [[ ${GENO} ]]
then
	ADD=".gq${GQ}"
fi



for TYPE in BND DEL DUP INS INV
do

	CNVPYTOR_COLLAPSE=${SAMPLE}.cnvpytor${ADD}.${TYPE}.collapse.bed
	ERDS_COLLAPSE=${SAMPLE}.erds${ADD}.${TYPE}.collapse.bed
	LUMPY_COLLAPSE=${SAMPLE}.lumpy${ADD}.${TYPE}.collapse.bed
	MANTA_COLLAPSE=${SAMPLE}.manta${ADD}.${TYPE}.collapse.bed

	RD_FILE=${SAMPLE}.RD${ADD}.${TYPE}.bed
	PRSR_FILE=${SAMPLE}.PRSR${ADD}.${TYPE}.bed


	### WITHIN METHOD MERGING

	## CNVPYTOR_COLLAPSE AND ERDS_COLLAPSE

	# check if both files are empty
	if [[ ( ! -s ${CNVPYTOR_COLLAPSE} ) && ( ! -s ${ERDS_COLLAPSE} ) ]]
	then
		echo "Both cnvpytor and erds files are empty for ${TYPE}!" >> ${LOG_DIR}/${1}.log
		touch ${RD_FILE}

	# check if cnvpytor file is empty
	elif [[ ! -s ${CNVPYTOR_COLLAPSE} ]]
	then
		echo "The cnvpytor file is empty for ${TYPE}!" >> ${LOG_DIR}/${1}.log
		while IFS="" read -r LINE || [ -n "$LINE" ]
		do 
			CALLS=$(echo ${LINE} | cut -f4)

			E_STAT=$(echo ${CALLS} | cut -f3 -d '_' | head -n1)


			if [[ -z "${E_STAT}" ]]
			then
				E_TAG="EA"
			else
				E_TAG="E${E_STAT}"
			fi

			TAG="CA:${E_TAG}:LA:MA"

			echo ${LINE} | awk -v OFS=$'\t' -v TAG=${TAG} -v TYPE=${TYPE} '{print $1,$2,$3,"RD_" TYPE "_" TAG "_" $1 ":" $2 "-" $3,$4}'
		done < ${ERDS_COLLAPSE} > ${RD_FILE}
		
	# check if erds file is empty
	elif [[ ! -s ${ERDS_COLLAPSE} ]]
	then
		echo "The erds file is empty for ${TYPE}!." >> ${LOG_DIR}/${1}.log
		while IFS="" read -r LINE || [ -n "$LINE" ]
		do 
			CALLS=$(echo ${LINE} | cut -f4)

			C_STAT=$(echo ${CALLS} | cut -f3 -d '_' | head -n1)


			if [[ -z "${C_STAT}" ]]
			then
				C_TAG="CA"
			else
				C_TAG="C${C_STAT}"
			fi

			TAG="${C_TAG}:EA:LA:MA"

			echo ${LINE} | awk -v OFS=$'\t' -v TAG=${TAG} -v TYPE=${TYPE} '{print $1,$2,$3,"RD_" TYPE "_" TAG "_" $1 ":" $2 "-" $3,$4}'
		done < ${CNVPYTOR_COLLAPSE} > ${RD_FILE}

	else
		# Find all sites in A and B that overlap, print both coordinates
		echo "Temp Union for ${TYPE} - RD" >> ${LOG_DIR}/${1}.log
		bedtools intersect -a ${CNVPYTOR_COLLAPSE} -b ${ERDS_COLLAPSE} -f ${OVLP_FRAC} -r -wa -wb | \
			tr '\t' ' ' > ${RD_FILE}.tmp

		# For each line, merge the two regions - one from each caller
		echo "Union for ${TYPE} - RD" >> ${LOG_DIR}/${1}.log
		while IFS="" read -r LINE || [ -n "$LINE" ]
		do 
			REGION=$(echo $LINE | sed 's/ /\n/5' | tr ' ' '\t' | bedtools sort -i - | bedtools merge -i - -c 4 -o distinct)
			CALLS=$(echo ${REGION} | cut -f4)

			C_STAT=$(echo ${CALLS} | tr ',' '\n' | grep 'cnvpytor' | cut -f3 -d '_' | head -n1)
			E_STAT=$(echo ${CALLS} | tr ',' '\n' | grep 'erds' | cut -f3 -d '_' | head -n1)


			if [[ -z "${C_STAT}" ]]
			then
				C_TAG="CA"
			else
				C_TAG="C${C_STAT}"
			fi


			if [[ -z "${E_STAT}" ]]
			then
				E_TAG="EA"
			else
				E_TAG="E${E_STAT}"
			fi


			TAG="${C_TAG}:${E_TAG}:LA:MA"

			echo ${REGION} | awk -v OFS=$'\t' -v TAG=${TAG} -v TYPE=${TYPE} '{print $1,$2,$3,"RD_" TYPE "_" TAG "_" $1 ":" $2 "-" $3,$4}'
		done < ${RD_FILE}.tmp > ${RD_FILE}
		

		# Original calls for each caller from union
		echo "Common Each ${TYPE} - RD" >> ${LOG_DIR}/${1}.log
		bedtools intersect -a ${CNVPYTOR_COLLAPSE} -b ${ERDS_COLLAPSE} -f ${OVLP_FRAC} -r -wa > ${RD_FILE}.cnvpytor.tmp 
		bedtools intersect -a ${ERDS_COLLAPSE} -b ${CNVPYTOR_COLLAPSE} -f ${OVLP_FRAC} -r -wa > ${RD_FILE}.erds.tmp

		echo "Singletons ${TYPE} - RD" >> ${LOG_DIR}/${1}.log
		while IFS="" read -r LINE || [ -n "$LINE" ]
		do 
			CALLS=$(echo ${LINE} | cut -f4)

			C_STAT=$(echo ${CALLS} | tr ',' '\n' | grep 'cnvpytor' | cut -f3 -d '_' | head -n1)


			if [[ -z "${C_STAT}" ]]
			then
				C_TAG="CA"
			else
				C_TAG="C${C_STAT}"
			fi

			TAG="${C_TAG}:EA:LA:MA"

			echo ${LINE} | awk -v OFS=$'\t' -v TAG=${TAG} -v TYPE=${TYPE} '{print $1,$2,$3,"RD_" TYPE "_" TAG "_" $1 ":" $2 "-" $3,$4}'
		done < <( grep -v -f ${RD_FILE}.cnvpytor.tmp ${CNVPYTOR_COLLAPSE} ) >> ${RD_FILE}


		while IFS="" read -r LINE || [ -n "$LINE" ]
		do 
			CALLS=$(echo ${LINE} | cut -f4)

			E_STAT=$(echo ${CALLS} | tr ',' '\n' | grep 'erds' | cut -f3 -d '_' | head -n1)


			if [[ -z "${E_STAT}" ]]
			then
				E_TAG="EA"
			else
				E_TAG="E${E_STAT}"
			fi

			TAG="CA:${E_TAG}:LA:MA"

			echo ${LINE} | awk -v OFS=$'\t' -v TAG=${TAG} -v TYPE=${TYPE} '{print $1,$2,$3,"RD_" TYPE "_" TAG "_" $1 ":" $2 "-" $3,$4}'
		done < <( grep -v -f ${RD_FILE}.erds.tmp ${ERDS_COLLAPSE})  >> ${RD_FILE}
		
		bedtools sort -i ${RD_FILE} > ${RD_FILE}.sort
		mv ${RD_FILE}.sort ${RD_FILE}

		# cleanup
		rm ${RESULTS_DIR}/${SAMPLE_NAME}*.tmp
	fi



	## LUMPY_COLLAPSE AND MANTA_COLLAPSE

	# Check if both files are empty
	if [[ ( ! -s ${LUMPY_COLLAPSE} ) && ( ! -s ${MANTA_COLLAPSE} ) ]]
	then
		echo "Both lumpy and manta files are empty for ${TYPE}!" >> ${LOG_DIR}/${1}.log
		touch ${PRSR_FILE}

	# check if lumpy file is empty
	elif [[ ! -s ${LUMPY_COLLAPSE} ]]
	then
		echo "The lumpy file is empty for ${TYPE}!." >> ${LOG_DIR}/${1}.log
		while IFS="" read -r LINE || [ -n "$LINE" ]
		do 
			CALLS=$(echo ${LINE} | cut -f4)

			M_STAT=$(echo ${CALLS} | cut -f3 -d '_' | head -n1)

			if [[ -z "${M_STAT}" ]]
			then
				M_TAG="MA"
			else
				M_TAG="M${M_STAT}"
			fi

			TAG="CA:EA:LA:${M_TAG}"

			echo ${LINE} | awk -v OFS=$'\t' -v TAG=${TAG} -v TYPE=${TYPE} '{print $1,$2,$3,"PRSR_" TYPE "_" TAG "_" $1 ":" $2 "-" $3,$4}'
		done < ${MANTA_COLLAPSE} > ${PRSR_FILE}

	# check if manta file is empty
	elif [[ ! -s ${MANTA_COLLAPSE}  ]]
	then
		echo "The manta file is empty for ${TYPE}!." >> ${LOG_DIR}/${1}.log
		while IFS="" read -r LINE || [ -n "$LINE" ]
		do 
			CALLS=$(echo ${LINE} | cut -f4)

			L_STAT=$(echo ${CALLS} | cut -f3 -d '_' | head -n1)

			if [[ -z "${L_STAT}" ]]
			then
				L_TAG="LA"
			else
				L_TAG="L${L_STAT}"
			fi

			TAG="CA:EA:${L_TAG}:MA"

			echo ${LINE} | awk -v OFS=$'\t' -v TAG=${TAG} -v TYPE=${TYPE} '{print $1,$2,$3,"PRSR_" TYPE "_" TAG "_" $1 ":" $2 "-" $3,$4}'
		done < ${LUMPY_COLLAPSE} > ${PRSR_FILE}

	else
		echo "Temp Union for ${TYPE} - PRSR" >> ${LOG_DIR}/${1}.log
		bedtools intersect -a ${LUMPY_COLLAPSE} -b ${MANTA_COLLAPSE} -f ${OVLP_FRAC} -r -wa -wb | \
			tr '\t' ' ' > \
			${PRSR_FILE}.tmp

		echo "Union for ${TYPE} - PRSR" >> ${LOG_DIR}/${1}.log


		while IFS="" read -r LINE || [ -n "$LINE" ]
		do 
			REGION=$(echo $LINE | sed 's/ /\n/5' | tr ' ' '\t' | bedtools sort -i - | bedtools merge -i - -c 4 -o distinct)
			CALLS=$(echo ${REGION} | cut -f4)

			L_STAT=$(echo ${CALLS} | tr ',' '\n' | grep 'lumpy' | cut -f3 -d '_' | head -n1)
			M_STAT=$(echo ${CALLS} | tr ',' '\n' | grep 'manta' | cut -f3 -d '_' | head -n1)

			if [[ -z "${L_STAT}" ]]
			then
				L_TAG="LA"
			else
				L_TAG="L${L_STAT}"
			fi

			if [[ -z "${M_STAT}" ]]
			then
				M_TAG="MA"
			else
				M_TAG="M${M_STAT}"
			fi

			TAG="CA:EA:${L_TAG}:${M_TAG}"

			echo ${REGION} | awk -v OFS=$'\t' -v TAG=${TAG} -v TYPE=${TYPE} '{print $1, $2, $3,"PRSR_" TYPE "_" TAG "_" $1 ":" $2 "-" $3, $4}'
		done < ${PRSR_FILE}.tmp > ${PRSR_FILE}



		bedtools intersect -a ${LUMPY_COLLAPSE} -b ${MANTA_COLLAPSE} -f ${OVLP_FRAC} -r -wa > ${PRSR_FILE}.lumpy.tmp
		bedtools intersect -a ${MANTA_COLLAPSE} -b ${LUMPY_COLLAPSE} -f ${OVLP_FRAC} -r -wa > ${PRSR_FILE}.manta.tmp


		while IFS="" read -r LINE || [ -n "$LINE" ]
		do 
			CALLS=$(echo ${LINE} | cut -f4)

			L_STAT=$(echo ${CALLS} | tr ',' '\n' | grep 'lumpy' | cut -f3 -d '_' | head -n1)


			if [[ -z "${L_STAT}" ]]
			then
				L_TAG="LA"
			else
				L_TAG="L${L_STAT}"
			fi

			TAG="CA:EA:${L_TAG}:MA"

			echo ${LINE} | awk -v OFS=$'\t' -v TAG=${TAG} -v TYPE=${TYPE} '{print $1,$2,$3,"PRSR_" TYPE "_" TAG "_" $1 ":" $2 "-" $3,$4}'
		done < <( grep -v -f ${PRSR_FILE}.lumpy.tmp ${LUMPY_COLLAPSE} ) >> ${PRSR_FILE}


		while IFS="" read -r LINE || [ -n "$LINE" ]
		do 
			CALLS=$(echo ${LINE} | cut -f4)

			M_STAT=$(echo ${CALLS} | tr ',' '\n' | grep 'manta' | cut -f3 -d '_' | head -n1)


			if [[ -z "${M_STAT}" ]]
			then
				M_TAG="MA"
			else
				M_TAG="M${M_STAT}"
			fi

			TAG="CA:EA:LA:${M_TAG}"

			echo ${LINE} | awk -v OFS=$'\t' -v TAG=${TAG} -v TYPE=${TYPE} '{print $1,$2,$3,"PRSR_" TYPE "_" TAG "_" $1 ":" $2 "-" $3,$4}'
		done < <( grep -v -f ${PRSR_FILE}.manta.tmp ${MANTA_COLLAPSE}) >> ${PRSR_FILE}

		bedtools sort -i ${PRSR_FILE} > ${PRSR_FILE}.sort
		mv ${PRSR_FILE}.sort ${PRSR_FILE}
		
		# cleanup
		rm ${RESULTS_DIR}/${SAMPLE_NAME}*.tmp
	fi

	echo " " >> ${LOG_DIR}/${1}.log
done



