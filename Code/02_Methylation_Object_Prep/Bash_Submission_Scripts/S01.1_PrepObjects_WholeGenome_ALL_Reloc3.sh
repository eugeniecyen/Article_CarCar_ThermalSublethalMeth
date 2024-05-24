#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 10
#$ -l h_vmem=13G
#$ -l h_rt=240:00:00
#$ -l highmem
#$ -m bea

##############################################################################################


# Create uniteCovALL and MethDiff object


##############################################################################################


### Prep environment ###

module load R/4.2.2

SCRIPT_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/00_Scripts/02_MethylKit_Scripts/Evol_Apps_Project/R_Scripts


### Call R script ###

Rscript $SCRIPT_DIR/R01.1_PrepObjects_WholeGenome_ALL_Reloc3.R ${NSLOTS}
