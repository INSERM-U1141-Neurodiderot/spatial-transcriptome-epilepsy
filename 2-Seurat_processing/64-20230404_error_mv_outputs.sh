#! /bin/bash

###########
# Options #
###########

# default values
DRY_RUN=0

while getopts "d:eo:" param
do
	case $param in
		# Sample day
		d) DAY=$OPTARG;;
		# Dry run
		e) DRY_RUN=1;;
		# Path to new output directory
		o) OUTPUT_DIR=$OPTARG;;
		:) echo "Option -$OPTARG expects an argument. Exit." ; exit;;
		\?) echo "Unvalid option -$OPTARG. Exit." ; exit;;
	esac
done



function run_cmd {
	local cmd=$@
	
	echo "${cmd}"
	eval time ${cmd}
	local error_code=$?
	if [ ${error_code} -ne 0 ]; then
		local message=${cmd%%*}
		echo "Error while ${message} execution."
		exit ${error_code}
	fi
}


declare -a FEAT_NB_LIST=(500 1000 1500 2000 2500 3000)
declare -a INTEGRATION_LIST=("no_integration" "Seurat" "Harmony")

for FEAT_NB in ${FEAT_NB_LIST[@]} ; do
	for INTEGRATION in ${INTEGRATION_LIST[@]} ; do
		DEST_DIR=${OUTPUT_DIR}/D${DAY}_samples/Norm_merge/Norm_mito-scale_orig-int_orig/QC-sum1000_det500_mito100_hb100/10-QC_filtering/10-LogNormalize/featnb${FEAT_NB}/${INTEGRATION}/origident_qc_mito_percent_RNA/
		
		DIR2MV=${OUTPUT_DIR}/D${DAY}_samples/Norm_merge/Norm_mito-scale_orig-int_orig/QC-sum1000_det500_mito100_hb100/10-QC_filtering/10-LogNormalize/featnb${FEAT_NB}/${INTEGRATION}/15_dims/
		CMD="cp -Rp ${DIR2MV} ${DEST_DIR}"
		if [ ${DRY_RUN} -eq 0 ]; then
			run_cmd ${CMD}
		else
			echo ${CMD}
		fi
		
		if [ $? -eq 0 ]; then
			CMD="rm -Rf ${DIR2MV}"
			if [ ${DRY_RUN} -eq 0 ]; then
				run_cmd ${CMD}
			else
				echo ${CMD}
			fi
		fi
		
		DIR2MV=${OUTPUT_DIR}/D${DAY}_samples/Norm_merge/Norm_mito-scale_orig-int_orig/QC-sum1000_det500_mito100_hb100/10-QC_filtering/10-LogNormalize/featnb${FEAT_NB}/${INTEGRATION}/origident_qc_mito_percent_RNA/origident_qc_mito_percent_RNA/30_dims/
		CMD="cp -Rp ${DIR2MV} ${DEST_DIR}"
		if [ ${DRY_RUN} -eq 0 ]; then
			run_cmd ${CMD}
		else
			echo ${CMD}
		fi
		
		DIR2RM=${OUTPUT_DIR}/D${DAY}_samples/Norm_merge/Norm_mito-scale_orig-int_orig/QC-sum1000_det500_mito100_hb100/10-QC_filtering/10-LogNormalize/featnb${FEAT_NB}/${INTEGRATION}/origident_qc_mito_percent_RNA/origident_qc_mito_percent_RNA/
		if [ $? -eq 0 ]; then
			CMD="rm -Rf ${DIR2RM}"
			if [ ${DRY_RUN} -eq 0 ]; then
				run_cmd ${CMD}
			else
				echo ${CMD}
			fi
		fi
	done
done

