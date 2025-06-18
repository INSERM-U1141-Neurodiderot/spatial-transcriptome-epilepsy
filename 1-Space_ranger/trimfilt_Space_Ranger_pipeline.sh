#! /bin/bash

# This script runs:
## read trimming and filtering
## ST Pipeline


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

while getopts "a:b:c:df:g:i:j:l:m:n:o:p:qr:s:t:u:v:w:" param
do
	case $param in
		# visium area identifier
		a) AREA=$OPTARG;;
		# read trimming and filtering configuration file
		b) TRIMMING_FILTERING_CONF=$OPTARG;;
		# max cores Space Ranger may request
		c) LOCAL_CORES=$OPTARG;;
		# dry run
		d) DRY_RUN=1;;
		# group ID in hawai
		g) GROUP_ID=$OPTARG;;
		# path to genome index directory
		i) GENOME_INDEX_DIR=$OPTARG;;
		# path to H&E brightfield image in either TIFF or JPG format
		j) IMAGE=$OPTARG;;
		# path to log directory
		l) LOG_DIR=$OPTARG;;
		# max RAM Space Ranger may request
		m) LOCAL_MEMORY=$OPTARG;;
		# dataset name
		n) EXP_NAME=$OPTARG;;
		# path to output directory
		o) NFS_OUTPUT_DIR=$OPTARG;; # must be /home/christophe.lepriol/spatial_transcriptomics_nfs...
		# suffix to append to EXP_NAME in Docker container names
		p) APPEND_SUFFIX=$OPTARG;;
		# fit output fastq file basenames with the naming conventions of bcl2fastq and mkfastq used by Space Ranger
		q) SPACE_RANGER_NAMES=1;;
		# path to directory storing raw fastq files
		r) RAW_DATA_DIR=$OPTARG;;
		# visium slide serial number
		s) VISIUM_SLIDE=$OPTARG;;
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

function check_file_exists {
	local PATH_TO_FILE=$1
	local DRY_RUN=$2
	
	echo "Check that ${PATH_TO_FILE} exists."
	if [ ! -e ${PATH_TO_FILE} ]; then
		echo "WARNING: the file ${PATH_TO_FILE} does not exist."
		if [ ${DRY_RUN} -eq 0 ]; then
			echo "Exit."
			exit 1
		else
			if [ ${DRY_RUN} -eq 1 ]; then
				echo "Keep on dry run."
			fi
		fi
	fi
}

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

# check files exist
## fastq files
RAW_R1_FASTQ_FILE=${RAW_DATA_DIR}/${EXP_NAME}_R1_001.fastq.gz
check_file_exists ${RAW_R1_FASTQ_FILE} ${DRY_RUN}
RAW_R2_FASTQ_FILE=${RAW_DATA_DIR}/${EXP_NAME}_R2_001.fastq.gz
check_file_exists ${RAW_R2_FASTQ_FILE} ${DRY_RUN}
# read trimming and filtering configuration file
check_file_exists ${TRIMMING_FILTERING_CONF} ${DRY_RUN}

# set sample main output directory
TRIMMING_FILTERING_CONFIG_FILENAME=$(basename ${TRIMMING_FILTERING_CONF})
GENOME_INDEXES_NAME=$(basename ${GENOME_INDEX_DIR})
OUTPUT_DIR=${WORK_DIR}/10-ST_analysis/10-Space_Ranger/output/${TRIMMING_FILTERING_CONFIG_FILENAME%.*}/${GENOME_INDEXES_NAME}/00-Samples/${EXP_NAME}

# set src directories
QC_SRC_DIR=${WORK_DIR}/00-QC/src
SPACE_RANGER_SRC_DIR=${WORK_DIR}/10-ST_analysis/10-Space_Ranger/src

# set Docker container basename
if [ ! -z ${EXP_NAME} ]; then
	if [ ! -z ${APPEND_SUFFIX} ]; then
		CONTAINER_BASENAME="${EXP_NAME}_${APPEND_SUFFIX}"
	else
		CONTAINER_BASENAME="${EXP_NAME}"
	fi
fi

# create and initiate report file
LOG_BASENAME=${EXP_NAME}_${GENOME_INDEXES_NAME}_trimfilt_Space_Ranger_pipeline
REPORT_FILE=${LOG_DIR}/${LOG_BASENAME}.report.txt
echo -e "#################################################\n# trimfilt_Space_Ranger_pipeline.sh report file #\n#################################################\n" > ${REPORT_FILE}
echo -e "# Inputs\nR1 fastq file: ${RAW_R1_FASTQ_FILE}\nR2 fastq file: ${RAW_R2_FASTQ_FILE}" >> ${REPORT_FILE}
echo -e "# Output directory: ${OUTPUT_DIR}" >> ${REPORT_FILE}
echo -e "# Parameters\n## read trimming and filtering configuration file: ${TRIMMING_FILTERING_CONF}" >> ${REPORT_FILE}
echo -e "## Space Ranger\nvisium area identifier: ${AREA}" >> ${REPORT_FILE}
echo "visium slide serial number: ${VISIUM_SLIDE}" >> ${REPORT_FILE}
echo "visium area identifier: ${AREA}" >> ${REPORT_FILE}
echo "H&E brightfield image: ${IMAGE}" >> ${REPORT_FILE}
echo "max cores Space Ranger may request: ${LOCAL_CORES}" >> ${REPORT_FILE}
echo "max RAM Space Ranger may request: ${LOCAL_MEMORY}" >> ${REPORT_FILE}
echo -e "genome index directory: ${GENOME_INDEX_DIR}\n\n" >> ${REPORT_FILE}

# read trimming and filtering
## set variables
### output directory
TRIMMING_FILTERING_OUTPUT_DIR=${OUTPUT_DIR}/00-Trimming_filtering
TRIMMING_FILTERING_TMP_OUTPUT_DIR=${TRIMMING_FILTERING_OUTPUT_DIR}/tmp
create_directory ${TRIMMING_FILTERING_TMP_OUTPUT_DIR} ${DRY_RUN}
### option
if [ ${SPACE_RANGER_NAMES} -eq 0 ]; then
	S_OPTION=""
else
	if [ ${SPACE_RANGER_NAMES} -eq 1 ]; then
		S_OPTION=" -s"
	fi
fi
## Docker run command
### container name
if [ ! -z ${EXP_NAME} ]; then
	NAME_STRING="--name ${CONTAINER_BASENAME}_trim_filt "
else
	NAME_STRING=""
fi
### array of volumes
VOLUMES=("/home/christophe.lepriol/NeuroDev_ADD/:/home/christophe.lepriol/NeuroDev_ADD/" "${TRIMMING_FILTERING_TMP_OUTPUT_DIR}/:/tmp/")
declare -a TRIMMING_FILTERING_VOLUME_STRING
for volume in ${VOLUMES[@]} ; do
  TRIMMING_FILTERING_VOLUME_STRING+=("-v ${volume}")
done
### command
TRIMMING_FILTERING_CMD="\"${QC_SRC_DIR}/trimming_filtering.sh -c ${TRIMMING_FILTERING_CONF} -i ${RAW_R1_FASTQ_FILE} -j ${RAW_R2_FASTQ_FILE} -l ${REPORT_FILE} -o ${TRIMMING_FILTERING_OUTPUT_DIR}${S_OPTION}\""
CMD_STRING="${SPACE_RANGER_SRC_DIR}/hawai_run_cmd_as_user.sh -c ${TRIMMING_FILTERING_CMD} -g ${GROUP_ID} -s ${SPACE_RANGER_SRC_DIR} -u ${USER} -v ${USER_ID}"
DOCKER_CMD="docker run --memory=25G --memory-reservation=20G ${TRIMMING_FILTERING_VOLUME_STRING[@]} ${NAME_STRING}rnaseq_expression:latest ${CMD_STRING}"
echo ${DOCKER_CMD}
DOCKER_LOGS_CMD="docker logs ${CONTAINER_BASENAME}_trim_filt &> ${LOG_DIR}/00-Trimming_filtering.log"
if [ ${DRY_RUN} -eq 0 ]; then
	eval ${DOCKER_CMD}
	eval ${DOCKER_LOGS_CMD}
fi

## cleaning
if [ -d ${TRIMMING_FILTERING_TMP_OUTPUT_DIR} ]; then
	echo "remove tmp directory: ${TRIMMING_FILTERING_TMP_OUTPUT_DIR}"
	rm -Rf ${TRIMMING_FILTERING_TMP_OUTPUT_DIR}
fi

echo -e "\n" >> ${REPORT_FILE}

# run Space Ranger
## check files exist
echo "check trimmed fastq files exist"
### fastq files after trimming and filtering
RAW_R1_FILENAME=$(basename ${RAW_R1_FASTQ_FILE})
if [ ${SPACE_RANGER_NAMES} -eq 0 ]; then
	R1_FASTQ_FILE_BASENAME=${RAW_R1_FILENAME%%.*}.trimmed
else
	if [ ${SPACE_RANGER_NAMES} -eq 1 ]; then
		R1_FASTQ_FILE_BASENAME=${RAW_R1_FILENAME%_R1*}_L001_R1_001
	fi
fi
FASTQ_FILE_1_TRIM=${TRIMMING_FILTERING_OUTPUT_DIR}/${R1_FASTQ_FILE_BASENAME}.fastq
check_file_exists ${FASTQ_FILE_1_TRIM} ${DRY_RUN}
RAW_R2_FILENAME=$(basename ${RAW_R2_FASTQ_FILE})
if [ ${SPACE_RANGER_NAMES} -eq 0 ]; then
	R2_FASTQ_FILE_BASENAME=${RAW_R2_FILENAME%%.*}.trimmed
else
	if [ ${SPACE_RANGER_NAMES} -eq 1 ]; then
		R2_FASTQ_FILE_BASENAME=${RAW_R2_FILENAME%_R2*}_L001_R2_001
	fi
fi
FASTQ_FILE_2_TRIM=${TRIMMING_FILTERING_OUTPUT_DIR}/${R2_FASTQ_FILE_BASENAME}.fastq
check_file_exists ${FASTQ_FILE_2_TRIM} ${DRY_RUN}

## set variables
### output directories: ${OUTPUT_DIR}/10-Pipeline/ must be not be created before running spaceranger count with --id=10-Pipeline, otherwise RuntimeError: ${OUTPUT_DIR}/10-Pipeline is not a pipestance directory
SPACE_RANGER_OUTPUT_DIRNAME=10-Pipeline
SPACE_RANGER_TMP_OUTPUT_DIR=${OUTPUT_DIR}/tmp
create_directory ${SPACE_RANGER_TMP_OUTPUT_DIR} ${DRY_RUN}
SPACE_RANGER_NFS_OUTPUT_DIR=${NFS_OUTPUT_DIR}/${TRIMMING_FILTERING_CONFIG_FILENAME%.*}/${GENOME_INDEXES_NAME}/00-Samples/${EXP_NAME}

### options
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

## Docker run command
### container name
if [ ! -z ${EXP_NAME} ]; then
	NAME_STRING="--name ${CONTAINER_BASENAME}_Space_Ranger "
else
	NAME_STRING=""
fi
### array of volumes
VOLUMES=("/home/christophe.lepriol/NeuroDev_ADD/:/home/christophe.lepriol/NeuroDev_ADD/" "${SPACE_RANGER_TMP_OUTPUT_DIR}/:/tmp/" "/home/christophe.lepriol/NeuroDev_ADD_nfs/spatial_transcriptomics:/home/christophe.lepriol/spatial_transcriptomics_nfs")
declare -a ST_VOLUME_STRING
for volume in ${VOLUMES[@]} ; do
  ST_VOLUME_STRING+=("-v ${volume}")
done
### command
SPACE_RANGER_CMD="\"${SPACE_RANGER_SRC_DIR}/Space_Ranger.sh -a ${AREA} ${C_OPTION}-g ${GENOME_INDEX_DIR} -i ${IMAGE} -j ${SPACE_RANGER_OUTPUT_DIRNAME} -l ${REPORT_FILE} ${M_OPTION}-n ${EXP_NAME} -o ${SPACE_RANGER_NFS_OUTPUT_DIR} -r ${FASTQ_FILE_1_TRIM} -s ${FASTQ_FILE_2_TRIM} -v ${VISIUM_SLIDE}\""
CMD_STRING="${SPACE_RANGER_SRC_DIR}/hawai_run_cmd_as_user.sh -c ${SPACE_RANGER_CMD} -g ${GROUP_ID} -s ${SPACE_RANGER_SRC_DIR} -u ${USER} -v ${USER_ID}"
DOCKER_CMD="docker run --memory=45G --memory-reservation=40G ${ST_VOLUME_STRING[@]} ${NAME_STRING}spaceranger:latest ${CMD_STRING}"
echo ${DOCKER_CMD}
SPACE_RANGER_LOG_FILE=${LOG_DIR}/${SPACE_RANGER_OUTPUT_DIRNAME}.log
DOCKER_LOGS_CMD="docker logs ${CONTAINER_BASENAME}_Space_Ranger &> ${SPACE_RANGER_LOG_FILE}"
if [ ${DRY_RUN} -eq 0 ]; then
	eval ${DOCKER_CMD}
	eval ${DOCKER_LOGS_CMD}
fi

## cleaning
### Space Ranger tmp directory
if [ -d ${SPACE_RANGER_TMP_OUTPUT_DIR} ]; then
	echo "remove tmp directory: ${SPACE_RANGER_TMP_OUTPUT_DIR}"
	rm -Rf ${SPACE_RANGER_TMP_OUTPUT_DIR}
fi
### trimmed and filtered fastq files
FASTQ_FILE_1_TRIM_FASTQC=${TRIMMING_FILTERING_OUTPUT_DIR}/${R1_FASTQ_FILE_BASENAME}_fastqc.html
if [ -e ${FASTQ_FILE_1_TRIM_FASTQC} ]; then
	check_file_exists ${FASTQ_FILE_1_TRIM} ${DRY_RUN}
	echo "remove trimmed and filtered fastq file: ${FASTQ_FILE_1_TRIM}"
	rm -f ${FASTQ_FILE_1_TRIM}
fi
FASTQ_FILE_2_TRIM_FASTQC=${TRIMMING_FILTERING_OUTPUT_DIR}/${R2_FASTQ_FILE_BASENAME}_fastqc.html
if [ -e ${FASTQ_FILE_2_TRIM_FASTQC} ]; then
	check_file_exists ${FASTQ_FILE_2_TRIM} ${DRY_RUN}
	echo "remove trimmed and filtered fastq file: ${FASTQ_FILE_2_TRIM}"
	rm -f ${FASTQ_FILE_2_TRIM}
fi

# check outputs
echo -e "\n# check outputs\n" >> ${REPORT_FILE}
OUTPUT_DIR_BASENAME=$(basename ${OUTPUT_DIR})
echo "## ${OUTPUT_DIR_BASENAME}/:" >> ${REPORT_FILE}
CMD="ls -l ${OUTPUT_DIR}/"
echo ${CMD} >> ${REPORT_FILE}
eval ${CMD} >> ${REPORT_FILE}
CMD="du -sh ${OUTPUT_DIR}/*"
echo ${CMD} >> ${REPORT_FILE}
eval ${CMD} >> ${REPORT_FILE}
echo "### ${OUTPUT_DIR_BASENAME}/$(basename ${TRIMMING_FILTERING_OUTPUT_DIR})/:" >> ${REPORT_FILE}
CMD="ls -l ${TRIMMING_FILTERING_OUTPUT_DIR}/"
echo ${CMD} >> ${REPORT_FILE}
eval ${CMD} >> ${REPORT_FILE}
SPACE_RANGER_OUTPUT_DIR=${OUTPUT_DIR}/${SPACE_RANGER_OUTPUT_DIRNAME}
SPACE_RANGER_OUTPUT_DIR_BASENAME=$(basename ${SPACE_RANGER_OUTPUT_DIR})
echo "### ${OUTPUT_DIR_BASENAME}/${SPACE_RANGER_OUTPUT_DIR_BASENAME}/:" >> ${REPORT_FILE}
CMD="ls -l ${SPACE_RANGER_OUTPUT_DIR}/"
echo ${CMD} >> ${REPORT_FILE}
eval ${CMD} >> ${REPORT_FILE}
echo "#### ${OUTPUT_DIR_BASENAME}/${SPACE_RANGER_OUTPUT_DIR_BASENAME}/outs/:" >> ${REPORT_FILE}
CMD="du -sh ${SPACE_RANGER_OUTPUT_DIR}/outs"
echo ${CMD} >> ${REPORT_FILE}
eval ${CMD} >> ${REPORT_FILE}
CMD="ls -l ${SPACE_RANGER_OUTPUT_DIR}/outs"
echo ${CMD} >> ${REPORT_FILE}
eval ${CMD} >> ${REPORT_FILE}
echo "#### ${OUTPUT_DIR_BASENAME}/${SPACE_RANGER_OUTPUT_DIR_BASENAME}/SPATIAL_RNA_COUNTER_CS/:" >> ${REPORT_FILE}
CMD="du -sh ${SPACE_RANGER_OUTPUT_DIR}/SPATIAL_RNA_COUNTER_CS"
echo ${CMD} >> ${REPORT_FILE}
eval ${CMD} >> ${REPORT_FILE}
echo "### ${OUTPUT_DIR_BASENAME}/logs/$(basename ${LOG_DIR})/:" >> ${REPORT_FILE}
CMD="ls -l ${LOG_DIR}/"
echo ${CMD} >> ${REPORT_FILE}
eval ${CMD} >> ${REPORT_FILE}
CMD="grep -A 37 'Outputs:' ${SPACE_RANGER_LOG_FILE}"
echo ${CMD} >> ${REPORT_FILE}
eval ${CMD} >> ${REPORT_FILE}


