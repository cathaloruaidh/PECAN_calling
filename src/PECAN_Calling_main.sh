#! /bin/bash



############################################################
# CNV calling: from BAM to VCF/BED
# 
# This script is to initialise the bash variables
# and check the data and software exists
############################################################






### Global directories and files
# The raw data lies here. 
# BED or PED/MAP files must live here. 
SOURCE_DIR=${PWD}

LOG_DIR=${SOURCE_DIR}/logs
RESULTS_DIR=${SOURCE_DIR}/results
TEMP_DIR=${SOURCE_DIR}/temp

MAIN_LOG_FILE=${LOG_DIR}/main.log


### Edit these variables: 
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


# where softare lives
TOOL_DIR=/home/shared/cathal/tools

# where the reference data (FASTA files) live 
REF_DIR=/home/shared/cathal/reference/ReferenceGenome
REFERENCE=GRCh38




# Get general functions defined elsewhere
mkdir ${LOG_DIR}

if [ $? -ne 0 ]
then
	log "Functions script returned an error: $?" 1
	exit 11
fi



set -o pipefail
ulimit -n $(ulimit -Hn) 2> /dev/null


### Argument processing
FAST=false
E=10
SAMPLE_NAME="output_file"
S=1
NPROCS=4  # change to appropriate default
VERBOSE=4
VCF=""
GENO=false
GQ=0
OVLP_FRAC=0.5
COLLAPSE_FRAC=0.25
RLCR_FRAC=0.75
REF_BUILD="GRCh38"
UPPER=30000000
SEX="U"
INCLUDE_FILE=""
EXCLUDE_FILE=""

cmd(){
	echo `basename $0`
}


usage(){
echo -e "\
Usage: `cmd` [OPTIONS ...] \n"

echo -e "\
-b, --bam; <FILE> ; input BAM file; [${BAM}]
-v, --vcf; <FILE>; SNV/indel VCF file; [${VCF}]
-n, --name; <STRING>; output file name; [${SAMPLE_NAME}]
-r, --reference; [GRCh38/GRCh37]; reference genome build; [GRCh38]
-s, --start; <INT>; script to start on (1-6); [${S}]
-e, --end; <INT>; script to finish on (1-6); [${E}]
-u, --upper; <INT>; upper CNV size limit; [${UPPER}]
-x, --sex; [F/M]; the genetic sex of the sample; [${SEX}]
--gq ; <INT> ; SV2 genotype quality score filter ; [${GQ}]
--fracOverlap ; <FLOAT> ; overlap fraction in general ; [${OVLP_FRAC}]
--fracCollapse ; <FLOAT> ; overlap fraction for collapsing ; [${COLLAPSE_FRAC}]
--fracRLCR ; <FLOAT> ; overlap fraction for removing RLCR ; [${RLCR_FRAC}]
  ;
--binSize ; <INT> ; bin size for CNVpytor ; [estimate]
--readLength ; <INT> ; read length ; [estimate]
--include ; <FILE> ; positions to include ; [standard contigs]
--exclude ; <FILE> ; positions to exclude ; [non-standard contigs]
--refDir ; <STRING> ; directory containing reference data ; []
--fasta ; <FILE> ; reference FASTA file ; []
  ;
-f, --faster; ; focus on DEL and DUP only ; [${FAST}]
-g, --genotype ;  ; call genotypes with SV2 ; [${GENO}]
-l, --log; <INT>; set log verbosity level; [${VERBOSE}]
-t, --threads; <INT>; number of threads ; [${NPROCS}]
  ;
-h, --help; ; output this message
" | column -t -s ";"
}


OPTS=`getopt -o b:e:fghl:n:r:s:t:v:u:x: \
	--long bam:,end:,faster,geno,help,log:,name:,reference:,start:,threads:,vcf:,binSize:,readLength:,include:,exclude:,gq:,fracOverlap:,fracCollapse:,fracRLCR:,sex:,upper:,fasta: \
	-n '$(cmd)' -- "$@"`

if [ $? != 0 ]
then 
	echo "Error with arguments. Terminating ..." >&2
	exit 1
fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$OPTS"



while true; do
	case "$1" in
		-b | --bam)
			BAM="$2"
			shift 2
			;;

		-e | --end)
			E="$2"
			shift 2
			;;

		-g | --genotype)
			GENO=true
			shift
			;;

		-f | --faster)
			FAST=true
			shift
			;;

		-h | --help)
			usage
			exit 0
			;;

		-l | --log)
			VERBOSE=$2
			shift 2
			;;

		-n | --name)
			SAMPLE_NAME=$2
			shift 2
			;;

		-r | --reference)
			REF_BUILD="$2"
			shift 2
			;;

		-s | --start)
			S="$2"
			shift 2
			;;

		-t | --threads)
			NPROCS="$2"
			shift 2
			;;

		-u | --upper)
			UPPER="$2"
			shift 2
			;;

		-v | --vcf)
			VCF="$2"
			shift 2
			;;

		-x | --sex)
			SEX="$2"
			shift 2
			;;

		--binSize)
			BIN_SIZE="$2"
			shift 2
			;;

		--readLength)
			READ_LENGTH="$2"
			shift 2
			;;

		--include)
			INCLUDE_FILE="$2"
			shift 2
			;;

		--exclude)
			EXCLUDE_FILE="$2"
			shift 2
			;;

		--fasta)
			REF_FASTA="$2"
			shift 2
			;;

		--gq)
			GQ="$2"
			shift 2
			;;

		--fracOverlap)
			OVLP_FRAC="$2"
			shift 2
			;;

		--fracCollapse)
			COLLAPSE_FRAC="$2"
			shift 2
			;;

		--fracRLCR)
			RLCR_FRAC="$2"
			shift 2
			;;

		--)
			shift
			break
			;;

		\?)
			log "Invalid flags. Exiting ... " 1
			exit 12
			break
			;;

		:)
			log "Flag requires an argument. Exiting ... " 1
			exit 13
			break
			;;

		*)
			log "Error with flags. Exiting ... " 1
			exit 14
			break
			;;
	esac
done



### Script directory
# Test if the scripts directory exists
if [ -z "${SCRIPT_DIR}" ]
then
	log "Problem with script directory. Exiting ... " 1
	exit 3
fi




### Set global variables, perform initial tests and create files

# Create the main log file
if [[ -f ${MAIN_LOG_FILE} ]]
then
	mv ${MAIN_LOG_FILE} ${MAIN_LOG_FILE}.$(date +%Y-%m-%d_%H.%M.%S)
	touch ${MAIN_LOG_FILE}
	log "Renaming previous log file" 3
else
	log "Creating the main log file: ${MAIN_LOG_FILE}" 3
	touch ${MAIN_LOG_FILE}
fi




log "Calling initialisation script" 3 
log " " 3
. ${SCRIPT_DIR}/PECAN_00_Initialise.sh 

if [ $? -ne 0 ]
then
	log "Error - initialisation script returned an error: $?" 1
	exit 11
fi




### Main section of pipeline
log $(printf '#%.0s' $(seq 1 $(($(tput cols)-35)) ) ) 3 2> /dev/null
log " " 3 2>/dev/null
log "Main pipeline section" 3
log " " 3
log "Pre-processing" 3
log " " 3

CURR=0


# Call CNVs with four methods
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	SCRIPT_NAME="PECAN_01_RunCallers"
	log $SCRIPT_NAME 3
	(. ${SCRIPT_DIR}/$SCRIPT_NAME.sh $SCRIPT_NAME) 
	testResult $? $SCRIPT_NAME
fi




# Create the reference dictionary for the reference FASTA file
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	SCRIPT_NAME="PECAN_02_Genotype"
	log $SCRIPT_NAME 3
	(. ${SCRIPT_DIR}/$SCRIPT_NAME.sh $SCRIPT_NAME) 
	testResult $? $SCRIPT_NAME
fi



# Map the reads and convert to a BAM file
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	SCRIPT_NAME="PECAN_03_CollapseWithin"
	log $SCRIPT_NAME 3
	(. ${SCRIPT_DIR}/$SCRIPT_NAME.sh $SCRIPT_NAME) 
	testResult $? $SCRIPT_NAME
fi



# Edit the Read Group
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	SCRIPT_NAME="PECAN_04_WithinCaller_Merge"
	log $SCRIPT_NAME 3
	(. ${SCRIPT_DIR}/$SCRIPT_NAME.sh $SCRIPT_NAME) 
	testResult $? $SCRIPT_NAME
fi



# Reorder the BAM file
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	SCRIPT_NAME="PECAN_05_BetweenCaller_Merge"
	log $SCRIPT_NAME 3
	(. ${SCRIPT_DIR}/$SCRIPT_NAME.sh $SCRIPT_NAME) 
	testResult $? $SCRIPT_NAME
fi



# Sort the BAM file
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	SCRIPT_NAME="PECAN_06_RLCR_Overlap"
	log $SCRIPT_NAME 3
	(. ${SCRIPT_DIR}/$SCRIPT_NAME.sh $SCRIPT_NAME) 
	testResult $? $SCRIPT_NAME
fi




log "Pipeline completed sucessfully." 3
log $(printf '#%.0s' $(seq 1 $(($(tput cols)-35)) ) ) 3 2> /dev/null


