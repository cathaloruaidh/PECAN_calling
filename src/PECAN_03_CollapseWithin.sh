#! /bin/bash

##################################################
# collapse the raw output of each caller, after
# genotyping (if selected)
##################################################


rm -f ${LOG_DIR}/${1}.log



# if genotyping, add the GQ filter to filename
ADD=""
if [[ ${GENO} ]]
then
	ADD=".gq${GQ}"
fi



for CALLER in cnvpytor erds lumpy manta ; 
do 
	for TYPE in BND DEL DUP INS INV
	do 
		PREFIX=${RESULTS_DIR}/${SAMPLE_NAME}.${CALLER}${ADD}.${TYPE}


		# no data for CNV type, continue
		if [[  ( ! -s ${PREFIX}.bed ) || ( ! -f ${PREFIX}.bed )  ]]
		then
			echo "No input file for ${PREFIX}, or file is empty. " >> ${LOG_DIR}/${1}.log
			cp ${PREFIX}.bed ${PREFIX}.collapse.bed
			continue
		fi


		# split by chromosome
		rm -f ${PREFIX}.chr*.bed
		awk -v PREFIX=${PREFIX} '{print > PREFIX "." $1 ".bed"}' ${PREFIX}.bed


		# collapse per chromosome 
		# (could parallelise this to speed up)
		for CHROM in $(cut -f1 ${PREFIX}.bed | sort | uniq)
		do
			${SCRIPT_DIR}/PECAN_MakeSets.sh ${PREFIX}.${CHROM} ${CALLER}_${TYPE} ${COLLAPSE_FRAC} 
		done


		# if collapsed files exist per chromosome, combine into one file
		# and delete chromosome files
		if [[ $(cat ${PREFIX}.chr*.collapse.bed 2> /dev/null | wc -l) -ne 0 ]]
		then
			cat ${PREFIX}.chr*.collapse.bed  | bedtools sort -i - > ${PREFIX}.collapse.bed
			rm ${PREFIX}.chr?.*bed ${PREFIX}.chr??.*bed

		fi


		echo "Collapsed ${PREFIX}" >> ${LOG_DIR}/${1}.log
	done
done



