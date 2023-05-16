#! /bin/bash

##################################################
# reomve variants that overlap repeat or low
# complexity regions (RLCR) 
##################################################


SAMPLE=${RESULTS_DIR}/${SAMPLE_NAME}


# get RLCR reference file
RLCR_BED=${CURR_SCRIPT_DIR}/RLCRs.${REF_BUILD}.noRepeatMasker.centromeres.bed


rm -f ${LOG_DIR}/${1}.log



ADD=""
if [[ ${GENO} ]]
then
	ADD=".gq${GQ}"
fi



for TYPE in BND DEL DUP INS INV
do

	FILE_SAMPLE_NAME=${SAMPLE}.ALL${ADD}.${TYPE}


	# first check if file is empty
	if [[ ! -s ${FILE_SAMPLE_NAME}.bed ]]
	then
		echo -e "Filtered ${SAMPLE_NAME}, type ${TYPE} - input file empty\n\n\n" >> ${LOG_DIR}/${1}.log 
		touch ${FILE_SAMPLE_NAME}.RLCR_removed.bed
	else
		echo "Filtered ${SAMPLE_NAME}, type ${TYPE} - file found" >> ${LOG_DIR}/${1}.log
		echo -e "\tGetting raw overlaps" >> ${LOG_DIR}/${1}.log


		# No proportional overlap between RLCR segment and CNV - any overlap will do. 
		bedtools intersect -a ${FILE_SAMPLE_NAME}.bed -b ${RLCR_BED} -e -wao | \
			awk -v SAMPLE_NAME=${FILE_SAMPLE_NAME} '{if($5 != ".") print > SAMPLE_NAME "-" $4 ".ovlp"}' 


		# now, sum the total overlap with all RLCR features
		echo -e "\tSumming overlaps and getting fraction of total" >> ${LOG_DIR}/${1}.log
		find ${RESULTS_DIR} -name "${SAMPLE_NAME}.ALL${ADD}.${TYPE}*.ovlp" > ${FILE_SAMPLE_NAME}.list.txt
		if [[ $( wc -l < ${FILE_SAMPLE_NAME}.list.txt ) -ne 0 ]]
		then

			while IFS="" read -r FILE || [[ -n "${FILE}" ]] 
			do
				OVERLAP=$(awk '{print $NF}' ${FILE} | paste -sd+ | bc)
				CNV=$(head -n1 ${FILE} | awk '{print $4}' )
				SIZE=$(head -n1 ${FILE} | awk '{print $3-$2}' )
				FRACTION=$(echo "${OVERLAP}/${SIZE}" | bc -l)
				echo -e "${CNV}\t${SIZE}\t${OVERLAP}\t${FRACTION}"
			done < ${FILE_SAMPLE_NAME}.list.txt  > ${FILE_SAMPLE_NAME}.RLCR_fraction.txt


			find ${RESULTS_DIR} -name "${SAMPLE_NAME}.ALL${ADD}.${TYPE}*.ovlp" -exec rm {} \;



			# remove calls where overlap fraction is above selected threshold
			echo -e "\tExtract calls that are sufficiently distinct from RLCR" >> ${LOG_DIR}/${1}.log
			rm -f ${FILE_SAMPLE_NAME}.RLCR_bad.txt
			awk -v FRAC=${RLCR_FRAC} '{if($4 > FRAC) print $1}' ${FILE_SAMPLE_NAME}.RLCR_fraction.txt > ${FILE_SAMPLE_NAME}.RLCR_bad.txt
			
			if [[ -s ${FILE_SAMPLE_NAME}.RLCR_bad.txt ]]
			then
				echo -e "\tRemove the bad CNV calls" >> ${LOG_DIR}/${1}.log
				grep -f ${FILE_SAMPLE_NAME}.RLCR_bad.txt ${FILE_SAMPLE_NAME}.bed > ${FILE_SAMPLE_NAME}.RLCR_bad.bed
				bedtools subtract -a ${FILE_SAMPLE_NAME}.bed -b ${FILE_SAMPLE_NAME}.RLCR_bad.bed -f 0.99999 -r > ${FILE_SAMPLE_NAME}.RLCR_removed.bed
			
				#rm ${FILE_SAMPLE_NAME}.RLCR_bad.txt ${FILE_SAMPLE_NAME}.RLCR_bad.bed
				echo -e "Filtered ${SAMPLE_NAME}, type ${TYPE}\n\n\n" >> ${LOG_DIR}/${1}.log
			else
				cp ${FILE_SAMPLE_NAME}.bed ${FILE_SAMPLE_NAME}.RLCR_removed.bed
				echo -e "Filtered ${SAMPLE_NAME}, type ${TYPE} - no overlap with RLCR\n\n\n" >> ${LOG_DIR}/${1}.log
			fi
			


		else
			cp ${FILE_SAMPLE_NAME}.bed ${FILE_SAMPLE_NAME}.RLCR_removed.bed
			echo -e "Filtered ${SAMPLE_NAME}, type ${TYPE} - no overlap with RLCR\n\n\n" >> ${LOG_DIR}/${1}.log
		fi

	fi

	echo "**************************************************" >> ${LOG_DIR}/${1}.log

done



