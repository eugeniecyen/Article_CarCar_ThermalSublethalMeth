#!/bin/bash
#$ -j y
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:00:00
#$ -t 1-40

#######################################################################

# Created by: Charley, 18.08.2023
# Array to merge all destranded meth calls by chromosome into a single genome file

# NB. Listing all chrom manually rather than using * globbing, so that chrom are
# in same order as ref genome

#######################################################################


### Extract sample ID ###

Sample_File=/data/SBCS-EizaguirreLab/Turtle_WGBS/00_General/SampleList_WGBS_2023.08.03.txt
SAMPLE=$(sed -n "${SGE_TASK_ID}p" $Sample_File)


### Set directories ###

DESTRAND_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/04_Bismark_Methylation_Calls/destranded_methylation_calls/$SAMPLE
OUT_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/04_Bismark_Methylation_Calls/destranded_methylation_calls/Whole_Genome

cd $DESTRAND_DIR


### Merge all chrom meth calls ###

# With chr0 at end -> doing it manually rather than with *

cat \
${SAMPLE}_chr1.CpG_merged.cov.gz \
${SAMPLE}_chr2.CpG_merged.cov.gz \
${SAMPLE}_chr3.CpG_merged.cov.gz \
${SAMPLE}_chr4.CpG_merged.cov.gz \
${SAMPLE}_chr5.CpG_merged.cov.gz \
${SAMPLE}_chr6.CpG_merged.cov.gz \
${SAMPLE}_chr7.CpG_merged.cov.gz \
${SAMPLE}_chr8.CpG_merged.cov.gz \
${SAMPLE}_chr9.CpG_merged.cov.gz \
${SAMPLE}_chr10.CpG_merged.cov.gz \
${SAMPLE}_chr11.CpG_merged.cov.gz \
${SAMPLE}_chr12.CpG_merged.cov.gz \
${SAMPLE}_chr13.CpG_merged.cov.gz \
${SAMPLE}_chr14.CpG_merged.cov.gz \
${SAMPLE}_chr15.CpG_merged.cov.gz \
${SAMPLE}_chr16.CpG_merged.cov.gz \
${SAMPLE}_chr17.CpG_merged.cov.gz \
${SAMPLE}_chr18.CpG_merged.cov.gz \
${SAMPLE}_chr19.CpG_merged.cov.gz \
${SAMPLE}_chr20.CpG_merged.cov.gz \
${SAMPLE}_chr21.CpG_merged.cov.gz \
${SAMPLE}_chr22.CpG_merged.cov.gz \
${SAMPLE}_chr23.CpG_merged.cov.gz \
${SAMPLE}_chr24.CpG_merged.cov.gz \
${SAMPLE}_chr25.CpG_merged.cov.gz \
${SAMPLE}_chr26.CpG_merged.cov.gz \
${SAMPLE}_chr27.CpG_merged.cov.gz \
${SAMPLE}_chr28.CpG_merged.cov.gz \
${SAMPLE}_chr0.CpG_merged.cov.gz > $OUT_DIR/${SAMPLE}.CpG_merged.cov.gz


