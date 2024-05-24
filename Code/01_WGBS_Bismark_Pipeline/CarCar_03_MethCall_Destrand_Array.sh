#!/bin/bash
#$ -j y
#$ -cwd
#$ -pe smp 10
#$ -l h_vmem=13G
#$ -l h_rt=240:00:00
#$ -t 1-3
#$ -l highmem

# Created by: Charley
# Date: 05.06.2023

# Array script running Bismark methylation calling per chromosome 
# Runs on whole genome per sample, but splits output into chromosomes
# Outputs in cytosineReport format

# Followed by destranding script on Bismark Cytosine Report files
# Using Robin's script merge_CpG.py: 
# https://github.com/rcristofari/penguin-tools/blob/master/merge_CpG.py 

# Methylation occurs symmetrically at CpG sites in vertebrates, so you can 
# increase coverage by merging adjacent cytosines of both reads per CpG site
# -> look at methylation of CpG sites rather than individual Cs. 
# Also it reduces multiple testing/pseudoreplication issues

#######################################################################

### Prep environment ###

module load bismark # v.0.22.1
module load samtools/1.9 # v.1.10 has incompatible gcc with bismark module

REPCORES=$((NSLOTS/2))


### Extract sample ID ###

Sample_File=/data/SBCS-EizaguirreLab/Turtle_WGBS/00_General/SampleList_New_SVL.txt
SAMPLE=$(sed -n "${SGE_TASK_ID}p" $Sample_File)


### Set directories ###

GENOME_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/02_Bismark_Genome_Index/Whole_Genome
BAM_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/03_Bismark_Alignment/bams/$SAMPLE
METHCALL_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/04_Bismark_Methylation_Calls/methylation_calls/$SAMPLE
CYTOSINE_REPORTS_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/04_Bismark_Methylation_Calls/methylation_calls/$SAMPLE/CytosineReports
DESTRAND_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/04_Bismark_Methylation_Calls/destranded_methylation_calls/$SAMPLE
REPORTS_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/04_Bismark_Methylation_Calls/reports

### Set scripts ###
MERGE_CPG=/data/SBCS-EizaguirreLab/Turtle_WGBS/00_Scripts/01_WGBS_Bismark_Pipeline/merge_CpG.py


# Make output directories, without error if already exists
mkdir -p $METHCALL_DIR


### Run Bismark methylation calling ###

# Separates cytosine report output by chromosome 

bismark_methylation_extractor --parallel ${REPCORES} \
--paired-end --comprehensive \
--split_by_chromosome --no_header --gzip \
--cytosine_report \
--output $METHCALL_DIR \
--genome_folder $GENOME_DIR \
$BAM_DIR/${SAMPLE}_deduplicated.sorted_by_name.bam

echo -e "\n### Methylation calling finished for ${SAMPLE} ###\n"


### Clean up output ###

cd $METHCALL_DIR

echo -e "### Entered directory: $PWD ###\n"

echo -e "\n### Cleaning up output files for ${SAMPLE} ###\n"

# Remove CHG and CHH context files to save space
rm CHG_context_*
rm CHH_context_*

echo -e "\n### Removed CHG and CHH context files ###\n"

# Rename cytosine report files to simple name with just SAMPLE_CHROM.CpG_report.txt.gz
for i in *deduplicated*.CpG_report.txt.gz
do mv $i ${i/"deduplicated.sorted_by_name.CpG_report.txt.chrSLK063_ragtag_"}
done

echo -e "\n### Renamed cytosine report files per chromosome ###\n"

# Create directory for cytosine reports
mkdir -p $CYTOSINE_REPORTS_DIR
mv *.CpG_report.txt.gz $CYTOSINE_REPORTS_DIR

echo -e "\n### Moved cytosine report files to new directory CytosineReports ###\n"

# Move reports to reports folder 
mv *_splitting_report.txt $REPORTS_DIR/splitting_reports

echo -e "\n### Moved splitting reports to splitting_reports directory ###\n"


### Perform destranding ###

# Prep environment 
module unload bismark
module unload samtools/1.9
module load python
source /data/SBCS-EizaguirreLab/Charley/environments/Robin_merge_CpG/bin/activate

echo -e "\n### Running destranding for all chrom of ${SAMPLE} ###\n"

# Navigate to working directory with meth calls
cd $CYTOSINE_REPORTS_DIR

echo -e "### Entered directory: $PWD ###\n"

# Make output directories (without error if already exists)
mkdir -p $DESTRAND_DIR 
mkdir -p $DESTRAND_DIR/both_strands # Put both_strands files into diff directory, as not needed
mkdir -p $REPORTS_DIR/destranding_reports/${SAMPLE}

# Run loop for CytosineReports of all chromosomes
for i in *CpG_report.txt.gz 

do
  # Run destranding script
  python $MERGE_CPG --input $i

  echo -e "\n### Finished destranding: $i ###"
  
  # Tidy output
  gzip *.cov # Zip output files
  chmod a+w *.cov.gz # Add permissions

  mv *_merged.cov.gz $DESTRAND_DIR
  mv *_both_strands.cov.gz $DESTRAND_DIR/both_strands
  mv *_merged.report.txt $REPORTS_DIR/destranding_reports/${SAMPLE}

  echo -e "### Finished tidying up: $i ###\n"
done


### Finished ###
echo -e "\n### Finished destranding for all chrom of ${SAMPLE} ###\n"

echo -e "\n### All done ###\n"

deactivate



