#!/bin/bash
#SBATCH -J scpa
#SBATCH -c 1
#SBATCH --mem=60G
#SBATCH --output=/home/adufour/work/logs/SCPA_epilepsie.log

cd /home/adufour/work
source /home/adufour/.bashrc
source activate scpa
Rscript /home/adufour/work/SCPA_epilepsie.R