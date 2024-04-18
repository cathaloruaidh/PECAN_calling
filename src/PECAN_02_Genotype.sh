#! /bin/bash


rm -f ${LOG_DIR}/${1}.log


# create SV2 temp directory
GT_TEMP_DIR=${RESULTS_DIR}/GENOTYPE/tmp

mkdir ${GT_TEMP_DIR}


# extract calls from VCF, >1kbp for RD callers only
bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\n' ${RESULTS_DIR}/CNVPYTOR/${SAMPLE_NAME}.cnvpytor.vcf.gz 2> /dev/null | awk '$3-$2 > 1000' > ${RESULTS_DIR}/${SAMPLE_NAME}.cnvpytor.txt
bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\n' ${RESULTS_DIR}/ERDS/${SAMPLE_NAME}.erds.vcf 2> /dev/null | awk '$3-$2 > 1000' > ${RESULTS_DIR}/${SAMPLE_NAME}.erds.txt
bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\n' ${RESULTS_DIR}/LUMPY/${SAMPLE_NAME}.lumpy.vcf.gz 2> /dev/null > ${RESULTS_DIR}/${SAMPLE_NAME}.lumpy.txt
bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\n' ${RESULTS_DIR}/MANTA/results/variants/diploidSV.vcf.gz 2> /dev/null > ${RESULTS_DIR}/${SAMPLE_NAME}.manta.txt



# convert to BED and add tag
for CALLER in cnvpytor erds lumpy manta
do
	for TYPE in DEL DUP
	do

		grep ${TYPE} ${RESULTS_DIR}/${SAMPLE_NAME}.${CALLER}.txt | \
			awk -v OFS=$'\t' -v PREFIX=${CALLER}_${TYPE}_S '{print $1, $2-1, $3-1, PREFIX "_" $1 ":" $2-1 "-" $3-1 }' > \
			${RESULTS_DIR}/${SAMPLE_NAME}.${CALLER}.${TYPE}.bed

	done



	if [[ "${FAST}" == true ]]
	then
		log "Fast mode selected: ignoring BND, INS, and INV types." 3
		touch ${RESULTS_DIR}/${SAMPLE_NAME}.${CALLER}.BND.bed
		touch ${RESULTS_DIR}/${SAMPLE_NAME}.${CALLER}.INS.bed
		touch ${RESULTS_DIR}/${SAMPLE_NAME}.${CALLER}.INV.bed

	else
		for TYPE in INV
		do

			grep ${TYPE} ${RESULTS_DIR}/${SAMPLE_NAME}.${CALLER}.txt | \
				awk -v OFS=$'\t' -v PREFIX=${CALLER}_${TYPE}_S '{print $1, $2-1, $3-1, PREFIX "_" $1 ":" $2-1 "-" $3-1 }' > \
				${RESULTS_DIR}/${SAMPLE_NAME}.${CALLER}.${TYPE}.bed

		done

		# BND and INS END coordinates are given as START+1
		for TYPE in BND INS
		do

			grep ${TYPE} ${RESULTS_DIR}/${SAMPLE_NAME}.${CALLER}.txt | \
				awk -v OFS=$'\t' -v PREFIX=${CALLER}_${TYPE}_S '{print $1, $2-1, $2, PREFIX "_" $1 ":" $2-1 "-" $3-1 }' > \
				${RESULTS_DIR}/${SAMPLE_NAME}.${CALLER}.${TYPE}.bed

		done
	fi

done


rm ${RESULTS_DIR}/${SAMPLE_NAME}*txt





# call genotypes with SV2
if [[ "${GENO}" == true ]]
then
	BAM_SAMPLE=$( samtools view -H ${BAM} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq )
	


	# get SV2 reference build
	if [[ "${REF_BUILD}" == "GRCh38"  ]]
	then
		BUILD="hg38"
	else
		BUILD="hg19"
	fi




	# create the pedigree FAM file for the sample
	PED=${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.ped
	
	if [[ ${SEX} == "M" ]]
	then
		SEX_CODE=1

	elif [[ ${SEX} == "F" ]]
	then 
		SEX_CODE=2
	
	else
		log "Error: sex variable (${SEX}) could not be determined" 1
	fi

	echo -e "${BAM_SAMPLE}\t${BAM_SAMPLE}\t0\t0\t${SEX_CODE}\t-9" > ${PED}
	



	if [[ "${FAST}" == true ]]
	then

		# get DEL and DUP in correct format, and restrict to calls under upper limit length
		# NOTE: this will skip if output already exist

		CALLER=cnvpytor
		if [[ ! -f ${RESULTS_DIR}/GENOTYPE/sv2_features/${SAMPLE_NAME}_sv2_features.${CALLER}.txt ]]
		then
			cat ${RESULTS_DIR}/${SAMPLE_NAME}.${CALLER}.DEL.bed ${RESULTS_DIR}/${SAMPLE_NAME}.${CALLER}.DUP.bed | \
				tr '_' '\t' | \
				awk -v OFS="\t" -v UPPER=${UPPER} '$3-$2 <= UPPER  { print $1, $2, $3, $5 }' | \
				bedtools sort -i - > ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.${CALLER}.sv2_input.bed

			sv2 -i ${BAM} \
				-b ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.${CALLER}.sv2_input.bed \
				-snv ${VCF} \
				-ped ${PED} \
				-g ${BUILD} \
				-M \
				-O ${RESULTS_DIR}/GENOTYPE \
				-T ${GT_TEMP_DIR} \
				-L ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.${CALLER}.gt.err \
				-o ${SAMPLE_NAME}.${CALLER}.gt -no-anno \
				| tee ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.${CALLER}.gt.log 

			mv ${RESULTS_DIR}/GENOTYPE/sv2_features/${BAM_SAMPLE}_sv2_features.txt ${RESULTS_DIR}/GENOTYPE/sv2_features/${SAMPLE_NAME}_sv2_features.${CALLER}.txt ;  

		fi


		# re-use the feature file generated from CNVpytor to save time
		for CALLER in erds lumpy manta
		do
			if [[ ! -f ${RESULTS_DIR}/GENOTYPE/sv2_features/${SAMPLE_NAME}_sv2_features.${CALLER}.txt ]]
			then
				cat ${RESULTS_DIR}/${SAMPLE_NAME}.${CALLER}.DEL.bed ${RESULTS_DIR}/${SAMPLE_NAME}.${CALLER}.DUP.bed | \
					tr '_' '\t' | \
					awk -v OFS="\t" -v UPPER=${UPPER} '$3-$2 <= UPPER  { print $1, $2, $3, $5 }' | \
					bedtools sort -i - > ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.${CALLER}.sv2_input.bed

				sv2 -i ${BAM} \
					-b ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.${CALLER}.sv2_input.bed \
					-snv ${VCF} \
					-ped ${PED} \
					-g ${BUILD} \
					-M \
					-O ${RESULTS_DIR}/GENOTYPE \
					-T ${GT_TEMP_DIR} \
					-L ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.${CALLER}.gt.err \
					-o ${SAMPLE_NAME}.${CALLER}.gt \
					-no-anno \
					-pre ${RESULTS_DIR}/GENOTYPE/sv2_preprocessing/ \
					| tee ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.${CALLER}.gt.log 

				mv ${RESULTS_DIR}/GENOTYPE/sv2_features/${BAM_SAMPLE}_sv2_features.txt ${RESULTS_DIR}/GENOTYPE/sv2_features/${SAMPLE_NAME}_sv2_features.${CALLER}.txt ;  
			fi

		done




	# for non-fast mode, genotype all CNVs
	# NOTE: this will skip if output already exist
	else

		# CNVpytor full
		if [[ ! -f ${RESULTS_DIR}/GENOTYPE/sv2_features/${SAMPLE_NAME}_sv2_features.cnvpytor.txt ]]
		then
			sv2 -i ${BAM} \
				-v ${RESULTS_DIR}/CNVPYTOR/${SAMPLE_NAME}.cnvpytor.vcf.gz \
				-snv ${VCF} \
				-ped ${PED} \
				-g ${BUILD} \
				-M \
				-O ${RESULTS_DIR}/GENOTYPE \
				-T ${GT_TEMP_DIR} \
				-L ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.cnvpytor.gt.err \
				-o ${SAMPLE_NAME}.cnvpytor.gt -no-anno \
				| tee ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.cnvpytor.gt.log 

			mv ${RESULTS_DIR}/GENOTYPE/sv2_features/${BAM_SAMPLE}_sv2_features.txt ${RESULTS_DIR}/GENOTYPE/sv2_features/${SAMPLE_NAME}_sv2_features.cnvpytor.txt ;  
		fi


		# ERDS pre 
		if [[ ! -f ${RESULTS_DIR}/GENOTYPE/sv2_features/${SAMPLE_NAME}_sv2_features.erds.txt ]]
		then
			sv2 -i ${BAM} \
				-v ${RESULTS_DIR}/ERDS/${SAMPLE_NAME}.erds.vcf \
				-snv ${VCF} \
				-ped ${PED} \
				-g ${BUILD} \
				-M \
				-O ${RESULTS_DIR}/GENOTYPE \
				-T ${GT_TEMP_DIR} \
				-L ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.erds.gt.err \
				-o ${SAMPLE_NAME}.erds.gt \
				-no-anno \
				-pre ${RESULTS_DIR}/GENOTYPE/sv2_preprocessing/ \
				| tee ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.erds.gt.log 

			mv ${RESULTS_DIR}/GENOTYPE/sv2_features/${BAM_SAMPLE}_sv2_features.txt ${RESULTS_DIR}/GENOTYPE/sv2_features/${SAMPLE_NAME}_sv2_features.erds.txt ;  
		fi


		# LUMPY pre
		if [[ ! -f ${RESULTS_DIR}/GENOTYPE/sv2_features/${SAMPLE_NAME}_sv2_features.lumpy.txt ]]
		then
			sv2 -i ${BAM} \
				-v ${RESULTS_DIR}/LUMPY/${SAMPLE_NAME}.lumpy.vcf.gz \
				-snv ${VCF} \
				-ped ${PED} \
				-g ${BUILD} \
				-M \
				-O ${RESULTS_DIR}/GENOTYPE \
				-T ${GT_TEMP_DIR} \
				-L ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.lumpy.gt.err \
				-o ${SAMPLE_NAME}.lumpy.gt \
				-no-anno \
				-pre ${RESULTS_DIR}/GENOTYPE/sv2_preprocessing/ \
				| tee ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.lumpy.gt.log 

			mv ${RESULTS_DIR}/GENOTYPE/sv2_features/${BAM_SAMPLE}_sv2_features.txt ${RESULTS_DIR}/GENOTYPE/sv2_features/${SAMPLE_NAME}_sv2_features.lumpy.txt ;  
		fi


		# Manta pre
		if [[ ! -f ${RESULTS_DIR}/GENOTYPE/sv2_features/${SAMPLE_NAME}_sv2_features.manta.txt ]]
		then
			sv2 -i ${BAM} \
				-v ${RESULTS_DIR}/MANTA/results/variants/diploidSV.vcf.gz \
				-snv ${VCF} \
				-ped ${PED} \
				-g ${BUILD} \
				-M \
				-O ${RESULTS_DIR}/GENOTYPE \
				-T ${GT_TEMP_DIR} \
				-L ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.manta.gt.err \
				-o ${SAMPLE_NAME}.manta.gt \
				-no-anno \
				-pre ${RESULTS_DIR}/GENOTYPE/sv2_preprocessing/ \
				| tee ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.manta.gt.log 

			mv ${RESULTS_DIR}/GENOTYPE/sv2_features/${BAM_SAMPLE}_sv2_features.txt ${RESULTS_DIR}/GENOTYPE/sv2_features/${SAMPLE_NAME}_sv2_features.manta.txt ;  
		fi

	fi



	# filter calls based on GQ
	for CALLER in cnvpytor erds lumpy manta
	do
		bcftools query -i "QUAL > ${GQ}" -f "%CHROM\t%POS\t%END\t%SVTYPE\t%QUAL\t[%GT]\n" ${RESULTS_DIR}/GENOTYPE/sv2_genotypes/${SAMPLE_NAME}.${CALLER}.gt.vcf 2> /dev/null | awk -v PREFIX=${CALLER} -v OFS="\t" '{ print PREFIX "_" $4 "_S_" $1 ":" $2-1 "-" $3 , $5, $6 }' > ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.${CALLER}.gt.gq${GQ}.txt


		for TYPE in BND DEL DUP INS INV
		do
			grep "${TYPE}"  ${RESULTS_DIR}/GENOTYPE/${SAMPLE_NAME}.${CALLER}.gt.gq${GQ}.txt | \
			cut -f1 | \
			grep -w -f - <( grep "${TYPE}" ${RESULTS_DIR}/${SAMPLE_NAME}.${CALLER}.${TYPE}.bed ) | \
				bedtools sort -i - > ${RESULTS_DIR}/${SAMPLE_NAME}.${CALLER}.gq${GQ}.${TYPE}.bed
		done
	
	done

fi

return 0
