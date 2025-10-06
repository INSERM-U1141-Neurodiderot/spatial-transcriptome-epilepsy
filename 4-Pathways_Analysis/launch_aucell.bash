#!/bin/bash
#SBATCH -J aucell
#SBATCH -c 15
#SBATCH --mem=550G
#SBATCH --output=/home/adufour/work/logs/auc_epilepsie.log

cd /home/adufour/work
source /home/adufour/.bashrc
source activate scenic
Rscript /home/adufour/work/scripts/R_script/AUCell_epilepsie_spatial.R