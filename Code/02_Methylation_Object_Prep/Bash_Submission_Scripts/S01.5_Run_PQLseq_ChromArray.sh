#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 10
#$ -l h_vmem=1G
#$ -l h_rt=240:00:00
#$ -t 1-29

##############################################################################################

# Created by Charley, 19th Mar 2024
# Run PQLSeq by chromosome in an array

##############################################################################################


### Prep environment ###

module load R/4.2.2

SCRIPT_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/00_Scripts/02_MethylKit_Scripts/Evol_Apps_Project/01_MethylKit_Object_Prep_Scripts


### Call R script ###

# With 2 arguments to pass into R script:

# 1. Chromosome number (which is the job array number i.e. SGE_TASK_ID). This must be the 1st arg.
# Note we have chromosomes 0:28 but 0 can't be an array ID so the array is 1:29, 
# then in the R script the chr number is set to SGE_TASK_ID - 1

# 2. Number of cores available for script

Rscript $SCRIPT_DIR/R01.5_Run_PQLseq_ChromArray.R ${SGE_TASK_ID} ${NSLOTS}
