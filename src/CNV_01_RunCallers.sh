#! /bin/bash



##################################################
# call with CNVpytor
##################################################


if [[ ! -f ${RESULTS_DIR}/CNVPYTOR/${SAMPLE_NAME}.cnvpytor.vcf.gz ]]
then 
	CHROMOSOMES=$( for CHROM in `seq 1 22` X Y ; do echo -n "chr$CHROM " ; done  )

	if [[ ! -f ${CNVPYTOR_ROOT} ]]
	then
		log "Creating ROOT file" 4
		cnvpytor -root ${CNVPYTOR_ROOT} -chrom ${CHROMOSOMES} -rd ${BAM} -j ${NPROCS} > ${CNVPYTOR_LOG}
		echo -e "\n\n" >> ${CNVPYTOR_LOG}
	fi

	log "Creating histogram" 4
	cnvpytor -root ${CNVPYTOR_ROOT} -his ${BIN_SIZE} -j ${NPROCS}  >> ${CNVPYTOR_LOG}
	echo -e "\n\n" >> ${CNVPYTOR_LOG}


	log "Partition ROOT file" 4
	cnvpytor -root ${CNVPYTOR_ROOT} -partition ${BIN_SIZE} -j ${NPROCS}  >> ${CNVPYTOR_LOG}
	echo -e "\n\n" >> ${CNVPYTOR_LOG}


	log "Call CNVs from ROOT file" 4
	cnvpytor -root ${CNVPYTOR_ROOT} -call ${BIN_SIZE} -j ${NPROCS}  > ${CNVPYTOR_OUTPUT}


	log "Convert output to VCF, compress and index" 4
	perl ${CURR_SCRIPT_DIR}/cnvnator2VCF.pl -prefix ${SAMPLE_NAME} -reference ${REF_BUILD} ${CNVPYTOR_OUTPUT} > ${CNVPYTOR_VCF}
	bgzip ${CNVPYTOR_VCF}
	tabix ${CNVPYTOR_VCF}.gz


	log "Completed calling with CNVpytor\n" 3

else
	log "Using existing CNVpytor files in ${RESULTS_DIR}/CNVPYTOR\n" 3

fi



##################################################
# call with ERDS 
##################################################


if [[ ! -f ${RESULTS_DIR}/ERDS/${SAMPLE_NAME}.erds.vcf ]]
then
	${TOOL_DIR}/erds_tcag/erds_pipeline.pl -b ${BAM} -v ${VCF} -r ${REF_FASTA} -n ${SAMPLE_NAME} -o ${ERDS_OUT_DIR} 2>&1 > ${ERDS_LOG}


	log "Completed calling with ERDS\n" 3

else
	
	log "Using existing ERDS files in ${RESULTS_DIR}/ERDS\n" 3

fi



##################################################
# call with LUMPY
##################################################


if [[ ! -f ${RESULTS_DIR}/LUMPY/${SAMPLE_NAME}.lumpy.vcf.gz ]]
then 
	#if [[ ! -f ${LUMPY_DISCORDANT} ]]
	#then
	#	log "Extracting Discordant reads, and sorting. " 4
	#	samtools view -h -F 1294 ${BAM} | samtools sort > ${LUMPY_DISCORDANT}
	#fi


	#if [[ ! -f ${LUMPY_SPLITTERS}  ]]
	#then
	#	log "Extracting split reads, and sorting. " 4
	#	samtools view -h ${BAM} | ${TOOL_DIR}/lumpy-sv-0.2.13/scripts/extractSplitReads_BwaMem -i stdin | samtools sort > ${LUMPY_SPLITTERS}
	#fi


	#log "Determining histogram for insert size. " 4


	log "Running lumpy via smoove (and sorting output). " 4
	${TOOL_DIR}/smoove_0.2.8/smoove call --name ${SAMPLE_NAME} \
		--fasta ${REF_FASTA} \
		-e ${EXCLUDE_FILE} \
		-o ${RESULTS_DIR}/LUMPY \
		${BAM} 2> ${LUMPY_LOG}

	vcf-sort -t . -c ${RESULTS_DIR}/LUMPY/${SAMPLE_NAME}-smoove.vcf.gz  > ${LUMPY_VCF}
	rm ${RESULTS_DIR}/LUMPY/${SAMPLE_NAME}-smoove.vcf.gz

	bgzip -f ${LUMPY_VCF}
	tabix -f ${LUMPY_VCF}.gz


	log "Completed calling with LUMPY\n" 3

else

	log "Using existing LUMPY files in ${RESULTS_DIR}/LUMPY\n" 3

fi



##################################################
# icall with Manta
##################################################


if [[ ! -f ${RESULTS_DIR}/MANTA/results/variants/diploidSV.vcf.gz ]]
then
	log "Running configuration file" 4
	python ${MANTA_DIR}/bin/configManta.py --bam ${BAM} --referenceFasta ${REF_FASTA} --callRegions ${INCLUDE_FILE} --runDir ${MANTA_RUN_DIR}


	log "Running workflow file" 4
	${MANTA_RUN_DIR}/runWorkflow.py -m local -j 4 2> ${MANTA_LOG}


	log "Completed calling with Manta" 3

else

	log "Using existing Manta files in ${RESULTS_DIR}/MANTA\n" 3

fi





