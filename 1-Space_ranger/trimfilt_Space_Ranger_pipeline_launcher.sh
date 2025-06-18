#! /bin/bash

# This script runs trimfilt_Space_Ranger_pipeline.sh for the samples listed in the SAMPLE_LIST_FILE file (-n option).


################
# Requirements #
################

# Docker image: rnaseq_expression:latest
## Conda environment: qc

# Docker image: spaceranger:latest
## Conda environment: spaceranger_env

###########
# Options #
###########

# default values
LOCAL_CORES=4
DRY_RUN=0
LOCAL_MEMORY=16
SPACE_RANGER_NAMES=0

while getopts "a:b:c:de:g:i:m:n:o:pr:s:u:v:w:" param
do
	case $param in
		# suffix to append to EXP_NAME in Docker container names
		a) APPEND_SUFFIX=$OPTARG;;
		# read trimming and filtering configuration file
		b) TRIMMING_FILTERING_CONF=$OPTARG;;
		# max cores Space Ranger may request
		c) LOCAL_CORES=$OPTARG;;
		# dry run
		d) DRY_RUN=1;;
		# date
		e) DATE=$OPTARG;;
		# group ID in hawai
		g) GROUP_ID=$OPTARG;;
		# path to genome index directory
		i) GENOME_INDEX_DIR=$OPTARG;;
		# max RAM Space Ranger may request
		m) LOCAL_MEMORY=$OPTARG;;
		# path to sample list file
		n) SAMPLE_LIST_FILE=$OPTARG;;
		# path to output directory
		o) NFS_OUTPUT_DIR=$OPTARG;; # must be /home/christophe.lepriol/spatial_transcriptomics_nfs...
		# fit intermediary Fastq file basenames with the naming conventions of bcl2fastq and mkfastq used by Space Ranger
		p) SPACE_RANGER_NAMES=1;;
		# path to directory storing raw fastq files
		r) RAW_DATA_DIR=$OPTARG;;
		# path to src directory
		s) SRC_DIR=$OPTARG;;		
		# user in hawai
		u) USER=$OPTARG;;
		# user ID in hawai
		v) USER_ID=$OPTARG;;
		# path to EpiReg working directory
		w) WORK_DIR=$OPTARG;;
		:) echo "Option -$OPTARG expects an argument. Exit." ; exit;;
		\?) echo "Unvalid option -$OPTARG. Exit." ; exit;;
	esac
done

#############
# Functions #
#############

function create_directory {
	local DIR_TO_CREATE=$1
	local DRY_RUN=$2
	
	if [ ! -d ${DIR_TO_CREATE} ]; then
		echo "create directory: ${DIR_TO_CREATE}"
		if [ ${DRY_RUN} -eq 0 ]; then
			mkdir -p ${DIR_TO_CREATE}
		fi
	else
		echo "${DIR_TO_CREATE} already exists."
	fi
}

##########
# Script #
##########

# process options
if [ ! -z ${LOCAL_CORES} ]; then
	C_OPTION="-c ${LOCAL_CORES} "
else
	C_OPTION=""
fi

if [ ! -z ${LOCAL_MEMORY} ]; then
	M_OPTION="-m ${LOCAL_MEMORY} "
else
	M_OPTION=""
fi

if [ ! -z ${APPEND_SUFFIX} ]; then
	P_OPTION="-p ${APPEND_SUFFIX} "
else
	P_OPTION=""
fi

if [ ${SPACE_RANGER_NAMES} -eq 0 ]; then
	Q_OPTION=""
else
	if [ ${SPACE_RANGER_NAMES} -eq 1 ]; then
		Q_OPTION="-q "
	fi
fi

# run trimfilt_Space_Ranger_pipeline.sh for all the samples in the sample list file
while IFS= read -r line;
do
	SAMPLE_ID=$(echo "$line" | cut -f 1)
	SLIDE=$(echo "$line" | cut -f 2)
	AREA=$(echo "$line" | cut -f 3)
	IMAGE=$(echo "$line" | cut -f 4)
	
	# create log directory
	TRIMMING_FILTERING_CONFIG_FILENAME=$(basename ${TRIMMING_FILTERING_CONF})
	GENOME_INDEXES_NAME=$(basename ${GENOME_INDEX_DIR})
	LOG_DIR=${WORK_DIR}/10-ST_analysis/10-Space_Ranger/output/${TRIMMING_FILTERING_CONFIG_FILENAME%.*}/${GENOME_INDEXES_NAME}/00-Samples/${SAMPLE_ID}/logs/${DATE}
	create_directory ${LOG_DIR} ${DRY_RUN}

	# QC and Space Ranger
	GENOME=$(basename ${GENOME_INDEX_DIR})
	LOG_BASENAME=${SAMPLE_ID}_${GENOME}_trimfilt_Space_Ranger_pipeline
	LOG_FILE=${LOG_DIR}/${LOG_BASENAME}.log
	CMD="nohup ${SRC_DIR}/trimfilt_Space_Ranger_pipeline.sh -a ${AREA} -b ${TRIMMING_FILTERING_CONF} ${C_OPTION}-g ${GROUP_ID} -i ${GENOME_INDEX_DIR} -j ${IMAGE} -l ${LOG_DIR} ${M_OPTION}-n ${SAMPLE_ID} -o ${NFS_OUTPUT_DIR} ${P_OPTION}${Q_OPTION}-r ${RAW_DATA_DIR} -s ${SLIDE} -u ${USER} -v ${USER_ID} -w ${WORK_DIR} > ${LOG_FILE}"
	echo ${CMD}
	if [ ${DRY_RUN} -eq 0 ]; then
		eval time ${CMD}
		# if error during trimfilt_Space_Ranger_pipeline.sh, then move on to the next sample
		if [ $? -ne 0 ]; then
			message=${CMD%%*}
			echo "Error while ${message} execution."
			continue
		fi
	fi
done < ${SAMPLE_LIST_FILE}

