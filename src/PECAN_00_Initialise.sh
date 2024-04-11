#!/bin/bash

############################################################
# CNV calling: from BAM to VCF/BED
# 
# This script is to initialise the bash variables
# and check the data and software exists
############################################################






### Text variables
PASS_TEST_LIGHT="[\e[102mPASSED\e[0m]"
PASS_TEST="[\e[42mPASSED\e[0m]"
FAIL_TEST_LIGHT="[\e[101mFAILED\e[0m]"
FAIL_TEST="[\e[41mFAILED\e[0m]"




# Create a logger for messages. 
DEBUG="\e[94mDEBUG\e[0m"
INFO="INFO \e[0m"
NOTE="\e[92mNOTE \e[0m"
WARN="\e[93mWARN \e[0m"
ERROR="\e[91mERROR\e[0m"

log() {
	DATE=`date +"%H:%M:%S %Y/%m/%d"`

	case $2 in
		4)
			MSG=${DEBUG}
			;;

		3)
			MSG=${INFO}
			;;

		2)
			MSG=${WARN}
			;;

		1)
			MSG=${ERROR}
			;;

		0)
			MSG=${NOTE}
			;;

		*)
			MSG=${INFO}
			;;
	esac

	if [[ $2 -le ${VERBOSE} ]]
	then
		MESSAGE="${MSG}\t${DATE}\t$1"
		echo -e "${MESSAGE}" >> ${MAIN_LOG_FILE}
		echo -e "${MESSAGE}"
	fi
}



# Create a function to return pass or fail
testResult(){
	if [ $1 != "0" ]
	then
		log "${FAIL_TEST} - $(awk -v var=$1 -F '\t' '{if($1==var) print $2": "$3}' ${SCRIPT_DIR}/return_errors)" 1
	else
		log "${PASS_TEST}" 3
	fi

	if [ $# -gt 1 ]
	then
		cat ${ERRORS_DIR}/$2.*.err 2> /dev/null | grep -E '^Error|^Warning' 2> /dev/null | while read -r line
		do
			log "$line" 1
		done
	fi

	log " " 3
}



# Error logger for script
err() {
	printf "%s\n" "$*" 1>&2 ;
}


# Print horizontal line for spacing, determined by terminal size
printLine(){
	width=$COLUMNS

	if [[ $width -eq 0 ]]
	then 
		width=$(tput cols)

		if [[ width -eq 0 ]]
		then
			width=80
		fi
	fi

	num=${width-35}

	return $(printf '#%.0s' $(seq 1 ${num}))
}








### Find source data
if [[ ! -f ${BAM} ]]
then
	log "Cannot find BAM input file. " 1
	exit 1
fi


if [[ ! -f ${VCF} ]]
then
	log "Cannot find VCF input file. " 1
	exit 1
fi


log "Found source files" 3
log "BAM - ${BAM}" 4
log "VCF - ${VCF}" 4
log " " 3





### Processors and Memory
TPROCS=$(grep -c ^processor /proc/cpuinfo)

if [[ -z "${NPROCS}" ]]
then
	NPROCS=4
fi

APROCS=$(( ${NPROCS} - 1 ))


log "Using ${NPROCS} processors out of a total ${TPROCS}" 3
log " " 3






### Global Directories and Files
# Test if the source directory exists
if [[ ! -d ${SOURCE_DIR} ]]
then
	log "Problem with source directory: ${SOURCE_DIR}" 1
	exit 8
fi

log  "Current directory: ${SOURCE_DIR}" 4

log " " 4






### Create required files and directories
# Log directory for all log files generated by Plink
if [[ ! -d ${LOG_DIR} ]]
then
	mkdir -p ${LOG_DIR}
	log "Creating log directory: ${LOG_DIR##*/}" 3

else
	log "Log directory already exists at: ${LOG_DIR##*/}" 4
fi



# Results directory for all results files
if [[ ! -d ${RESULTS_DIR} ]]
then
	log "Creating results directory: ${RESULTS_DIR##*/}" 3
	mkdir -p ${RESULTS_DIR}
	mkdir -p ${RESULTS_DIR}/CNVPYTOR
	mkdir -p ${RESULTS_DIR}/ERDS
	mkdir -p ${RESULTS_DIR}/LUMPY
	mkdir -p ${RESULTS_DIR}/MANTA
	mkdir -p ${RESULTS_DIR}/GENOTYPE

else
	log "Results directory already exists at: ${RESULTS_DIR##*/}" 4 
fi




# Temporary directory for any intermediate steps
if [[ ! -d ${TEMP_DIR} ]]
then
	log "Creating temp directory: ${TEMP_DIR##*/}" 3
	mkdir -p ${TEMP_DIR}
fi



log " " 4
log " " 3







### Set the program paths, and load them
# if the variables are empty, throw an error
if [[ -z "$( command -v samtools )" ]]
then
	log "Error: samtools not found." 1
	exit 6
fi
log "samtools was successfully found." 4



if [[ -z "$( command -v bcftools)" ]]
then
	log "Error: bcftools not found." 1
	exit 6
fi
log "bcftools was successfully found." 4






### Find/create the Resource Files, depending on the reference build.
if [[ ! -f ${REF_FASTA} ]]
then
	log "No FASTA file found" 1
else
	log "Reference file: $(basename ${REF_FASTA})" 4
	
	if [[ ! -f ${REF_FASTA}.fai ]]
	then
		log "Indexing FASTA file" 3
		samtools faidx ${REF_FASTA}
	fi


	if [[ -z ${INCLUDE_FILE} ]]
	then
		log "Generating include file" 3
		INCLUDE_FILE=$( echo ${REF_FASTA} | sed -e  's/\.gz$//g ; s/\.fasta$//g ; s/\.fa$//g' ).include.bed
		grep -E "^chr[0-9XY]{1,2}\s" ${REF_FASTA}.fai | awk -v FS="\t" '{ print $1, $2-1, $2}' |  bedtools sort -i - > ${INCLUDE_FILE}
	fi


	if [[ -z ${EXCLUDE_FILE} ]]
	then
		log "Generating exclude file" 3
		EXCLUDE_FILE=$( echo ${REF_FASTA} | sed -e  's/\.gz$//g ; s/\.fasta$//g ; s/\.fa$//g' ).exclude.bed
		grep -Ev "^chr[0-9XY]{1,2}\s" ${REF_FASTA}.fai | awk -v FS="\t" '{ print $1, $2-1, $2}' |  bedtools sort -i - > ${EXCLUDE_FILE}
	fi



fi







### set the input parameters from BAM file
log "Checking/calculating the input parameters" 3

# first, calculate the required bin size for CNVpytor based on the average depth of coverage. 
# This is rounded for convenience. 
if [[ -z "${BIN_SIZE}" && ${S} == 1 ]]
then
	log "Set the bin size by the average depth of coverage (using mosdepth)" 4

	if [[ ! -s ${TEMP_DIR}/${SAMPLE_NAME}.mosdepth.summary.txt ]]
	then
		${TOOL_DIR}/mosdepth_v0.3.2/mosdepth -t ${NPROCS} ${TEMP_DIR}/${SAMPLE_NAME} -f ${REF_FASTA} ${BAM}
	fi
	DP=$( grep 'total' ${TEMP_DIR}/${SAMPLE_NAME}.mosdepth.summary.txt | cut -f4 )
	BIN_SIZE=$( echo ${DP} | awk '{ print int( ((15000/$1)+50)/100 ) * 100 }' )
fi

log "Bin size: ${BIN_SIZE}" 4



# calculate the maximum read length based on the first chromosome in the BAM file
if [[ -z "${READ_LENGTH}" && ${S} == 1 ]]
then
	log "Determine max read length from first chromosome in the BAM" 4

	if [[ ! -f ${TEMP_DIR}/${SAMPLE_NAME}.read_lengths.txt ]]
	then
		CHR=$( samtools view ${BAM} | head -1 | cut -f3 )
		samtools view -T ${REF_FASTA} ${BAM} ${CHR} | awk '{ print length($10) }' | uniq | sort -rn | uniq > ${TEMP_DIR}/${SAMPLE_NAME}.read_lengths.txt
	fi
	READ_LENGTH=$( head -1 ${TEMP_DIR}/${SAMPLE_NAME}.read_lengths.txt )
fi

log "Read length: ${READ_LENGTH}" 4





log " "  4







### Name the output files


# global




# CNVpytor
CNVPYTOR_LOG=${RESULTS_DIR}/CNVPYTOR/${SAMPLE_NAME}.cnvpytor.log
CNVPYTOR_ROOT=${RESULTS_DIR}/CNVPYTOR/${SAMPLE_NAME}.root
CNVPYTOR_STATS=${RESULTS_DIR}/CNVPYTOR/${SAMPLE_NAME}.stats.txt
CNVPYTOR_EVAL=${RESULTS_DIR}/CNVPYTOR/${SAMPLE_NAME}.eval.txt
CNVPYTOR_OUTPUT=${RESULTS_DIR}/CNVPYTOR/${SAMPLE_NAME}.cnvpytor.txt
CNVPYTOR_VCF=${RESULTS_DIR}/CNVPYTOR/${SAMPLE_NAME}.cnvpytor.vcf



# erds
ERDS_LOG=${RESULTS_DIR}/ERDS/${SAMPLE_NAME}.erds.log
ERDS_OUT_DIR=${RESULTS_DIR}/ERDS/


# lumpy 
LUMPY_LOG=${RESULTS_DIR}/LUMPY/${SAMPLE_NAME}.lumpy.log
LUMPY_DISCORDANT=${RESULTS_DIR}/LUMPY/${SAMPLE_NAME}.disc.bam
LUMPY_SPLITTERS=${RESULTS_DIR}/LUMPY/${SAMPLE_NAME}.split.bam
LUMPY_HISTO=${RESULTS_DIR}/LUMPY/${SAMPLE_NAME}.histo.txt
LUMPY_VCF=${RESULTS_DIR}/LUMPY/${SAMPLE_NAME}.lumpy.vcf
LUMPY_GT_VCF=${RESULTS_DIR}/LUMPY/${SAMPLE_NAME}.lumpy.gt.vcf


# manta
MANTA_LOG=${RESULTS_DIR}/MANTA/${SAMPLE_NAME}.manta.log
MANTA_DIR=/home/shared/cathal/tools/manta-1.5.0/
MANTA_RUN_DIR=${RESULTS_DIR}/MANTA



