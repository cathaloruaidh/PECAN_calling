#! /bin/bash



PREFIX=${1}

NAME=${2}

if [[ ! -z "${3}"  ]]
then
	FRACTION=${3}
else
	echo "Error, no fraction of overlap passed. Exiting ..."
	exit 1
fi



# find CNVs that overlap with others in the input file (exccept themselves).
# Output the details to an overlap (ovlp) file
bedtools intersect -a ${PREFIX}.bed -b ${PREFIX}.bed -f ${FRACTION} -r -wao | \
	awk -v PREFIX=${PREFIX} '{if( $4 != $8) print > PREFIX "-" $4 ".ovlp"}'



# if there's no overlap anywhere, done!
if [[ $(ls ${PREFIX}-*.ovlp 2> /dev/null | wc -l) -eq 0 ]]
then
	awk -v OFS="\t" '{ print $1,$2,$3,$4,$4}' ${PREFIX}.bed > ${PREFIX}.collapse.bed
	exit 0
fi



# get the IDs of all CNVs that have overlap with something else
cat ${PREFIX}-*.ovlp | cut -f 4,8 | tr '\t' '\n' | sort | uniq > ${PREFIX}.names.txt



N_SETS=0
rm -f ${PREFIX}.SET*txt



# loop over al CNV IDs, making the sets
for CNV in $( cat ${PREFIX}.names.txt ) 
do
	CHROM=$(echo ${CNV} | cut -f4 -d_)

	unset SET_NUM
	unset RES


	# check if the CNV is already in one of the overlap files 
	# and if so, get the corresponding SET number
	for i in $(seq 1 ${N_SETS})
	do
		grep -w -f <( cat <( echo "${CNV}") <( cut -f8 ${PREFIX}-${CNV}.ovlp 2> /dev/null ) ) ${PREFIX}.SET_${i}.txt 2>&1 1> /dev/null
		RES=$?

		if [[ ${RES} -eq 0 ]]
		then
			SET_NUM=${i}
		fi
	done


	# If the CNV is in the set, add it:
	if [[ -n "${SET_NUM}" && ${SET_NUM} -le ${N_SETS} ]]
	then
		SET_NAME=${PREFIX}.SET_${SET_NUM}.txt

		if [[ -f ${PREFIX}-${CNV}.ovlp ]] 
		then  
			cat ${SET_NAME} <( echo ${CNV} ) <(cut -f 8 ${PREFIX}-${CNV}.ovlp) | sort | uniq > ${SET_NAME}.tmp
		else 
			cat ${SET_NAME} <( echo ${CNV} ) | sort | uniq > ${SET_NAME}.tmp
		fi 

		mv ${SET_NAME}.tmp ${SET_NAME}
		

	# If the CNV is NOT in any set, create new one:
	else

		N_SETS=$(( N_SETS + 1 ))
		SET_NAME=${PREFIX}.SET_${N_SETS}.txt

		if [[ -f ${PREFIX}-${CNV}.ovlp ]]
		then 
			# Include the CNV name!!!
			cat <( echo "${CNV}" ) <( cut -f 8 ${PREFIX}-${CNV}.ovlp ) | sort | uniq > ${SET_NAME}
		else
			echo ${CNV} > ${SET_NAME}
		fi
	fi
done



# Merge all the calls in the sets, to create the collapse file
for FILE in ${PREFIX}.SET_*.txt
do
	if [[ -s ${FILE} ]] 
	then
		grep -w -f ${FILE} ${PREFIX}.bed | bedtools sort -i - | bedtools merge -i - -c 4 -o collapse \
			| awk -F $'\t' -v OFS=$'\t' -v PREFIX=${NAME}_M '{print $1, $2, $3, PREFIX "_" $1 ":" $2 "-" $3, $4}'
	fi
done | bedtools sort -i - | uniq > ${PREFIX}.set_merged.bed


CNV_PREFIX=$( head -n1 ${PREFIX}.bed | cut -f4 | cut -f1 -d: )


cat ${PREFIX}.set_merged.bed <( grep -w -v -f ${PREFIX}.names.txt ${PREFIX}.bed ) | awk -v OFS="\t" '{if($5 == "") $5 = $4 ; print $0}' | bedtools sort -i - > ${PREFIX}.collapse.bed



rm ${PREFIX}.SET_*.txt ${PREFIX}-*ovlp ${PREFIX}.names.txt


echo "Collapsed file ${PREFIX}"




