#!/bin/bash
#SBATCH -J epilepsie_card
#SBATCH -c 3
#SBATCH --mem=250G
#SBATCH --output=/home/adufour/work/logs/epilepsie_card.log

cd /home/adufour/work
source /home/adufour/.bashrc
source activate card_sp
Rscript /home/adufour/work/scripts/R_script/spatial_deconvolution/card_epilepsie.R