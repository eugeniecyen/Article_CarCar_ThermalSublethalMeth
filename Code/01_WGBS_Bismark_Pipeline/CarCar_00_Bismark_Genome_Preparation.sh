#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=13G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y
#$ -m bea
#$ -l highmem

module load bismark # v.0.22.1

GENOME_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/02_Bismark_Genome_Index

REPCORES=$((NSLOTS/2))

bismark_genome_preparation \
--verbose --parallel $REPCORES \
$GENOME_DIR
