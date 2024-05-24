#!/bin/bash
#$ -j y
#$ -cwd
#$ -pe smp 8
#$ -l h_vmem=16G
#$ -l h_rt=240:00:00
#$ -t 1-192
#$ -l highmem

#######################################################################

# Created by: Charley, adapted from James
# Date: 24.02.2023

# Array script running Bismark whole genome alignment per sample

# 1) Create Bismark alignments per flowcell per sample
# 2) Run Bismark deduplication per flowcell per sample
# 3) Merge flowcells into single bam per sample
# 4) Delete all flowcell bams to save space
# 5) Sort merged bam by name (needed for Bismark methylation calling)

# Output = bam per chromosome per sample ready for methylation calling
# For SNP calling, bams need to be coordinate sorted and indexed after

#######################################################################


### Prep environment ###

module load bismark # v.0.22.1
module load samtools/1.9 # v.1.10 has incompatible gcc with bismark module

REPCORES=$((NSLOTS/2))


### Extract sample ID ###

Sample_File=/data/SBCS-EizaguirreLab/Turtle_WGBS/00_General/SampleList_WGBS_2023.02.20.txt
SAMPLE=$(sed -n "${SGE_TASK_ID}p" $Sample_File)


### Set directories ###

READS_DIR=/data/archive/archive-SBCS-EizaguirreLab/WGBS_trimmed/trimmed_cutadapt/$SAMPLE
GENOME_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/02_Bismark_Genome_Index/Whole_Genome
OUT_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/03_Bismark_Alignment/bams/$SAMPLE

# Make output directories, without error if already exists
mkdir -p $OUT_DIR


### Loop through all flowcells per sample ###

# Check if merged file exists yet -> only run alignment if it doesn't
if [[ -f ${OUT_DIR}/${SAMPLE}_deduplicated.bam ]]; 
 then
   echo -e "\n### Merged bam already exists for ${SAMPLE}. Not re-aligning ###\n"
 else
   echo -e "\n### Merged bam does not exist for ${SAMPLE}. Continue to alignment ###\n"

   # Navigate to READS_DIR
   cd $READS_DIR
  
    # Save unique flowcell info from fastq filenames to variable
    FLOWCELLS=$(ls $READS_DIR | cut -d '_' -f 1,2,3 | uniq )

    # For each flowcell read pair, take read group ID (rgid) from file name
    # Sample name is specified (smid) and bismark is run
  
   for FLOWCELL in $FLOWCELLS

    do
      # Read group ID in the format FlowCell.Lane (replace underscore in filename)
      rgid=$(echo $FLOWCELL | cut -d '_' -f 1,2 | sed 's/_/\./g')
  
      ### Run Bismark alignment ###
      echo -e "\n### Starting alignment for flowcell: ${FLOWCELL}, sample: ${SAMPLE} ###\n"
   
      bismark \
      -p ${REPCORES} \
      --rg_sample $SAMPLE --rg_id $rgid \
      -o $OUT_DIR \
      --basename $FLOWCELL \
      $GENOME_DIR \
      -1 ${FLOWCELL}*_1_trim.fq.gz \
      -2 ${FLOWCELL}*_2_trim.fq.gz

      echo -e "\n### Finished alignment for flowcell: ${FLOWCELL}, sample: ${SAMPLE} ###\n"
  
      ### Run Bismark deduplication ### 
      echo -e "\n### Starting deduplication for flowcell: ${FLOWCELL}, sample: ${SAMPLE} ###\n"
  
      deduplicate_bismark --paired --bam \
      --output_dir $OUT_DIR \
      --outfile ${FLOWCELL}_pe \
      ${OUT_DIR}/${FLOWCELL}_pe.bam
  
      echo -e "\n### Finished deduplication for flowcell: ${FLOWCELL}, sample: ${SAMPLE} ###\n"
  
   done # End flowcell loop
  

   ### Merge deduplicated bams from all flowcells ###
  
   echo -e "\n### Starting merging deduplicated flowcell bams for ${SAMPLE} ###\n"
  
   samtools merge -@ ${NSLOTS} \
   ${OUT_DIR}/${SAMPLE}_deduplicated.bam \
   ${OUT_DIR}/*_pe.deduplicated.bam
  
    echo -e "\n### Finished merging deduplicated flowcell bams for ${SAMPLE} ###\n"
  
fi # End if statement for alignment and deduplication
 

### Delete all flowcell bams, only if merged bam exists ###
  
if [[ -f ${OUT_DIR}/${SAMPLE}_deduplicated.bam ]]; 
    then
      echo -e "\n### Merged bam for has been created for ${SAMPLE} ###"
      echo -e "\n### Deleting flowcell bams ###\n"
     
      # Remove flowcell bams whilst keeping all reports
      rm ${OUT_DIR}/V*.bam
     
       echo -e "### Deleted flowcell bams for ${SAMPLE} ###\n"
       
    else
       echo -e "\n### Merged bam does not exist for ${SAMPLE} ###"
       echo -e "### Will not delete flowcell bams ###\n"
  
fi # End if statement for deleting flowcell bams


### Sort merged deduplicated bam by name ###

# -n flag sorts by read name (required for meth calling)

echo -e "### Sorting merged, deduplicated bam by name for ${SAMPLE} ###"

samtools sort -n -@ ${NSLOTS} \
-o ${OUT_DIR}/${SAMPLE}_deduplicated.sorted_by_name.bam \
${OUT_DIR}/${SAMPLE}_deduplicated.bam

echo -e "### Sorted merged, deduplicated bam by name for ${SAMPLE} ###"


### Delete unsorted bam ###

if [[ -f ${OUT_DIR}/${SAMPLE}_deduplicated.sorted_by_name.bam ]]; 
  then
      echo -e "\n### Sorted bam for has been created for ${SAMPLE} ###"
      echo -e "\n### Deleting unsorted bam ###\n"

      # Delete
      rm ${OUT_DIR}/${SAMPLE}_deduplicated.bam
     
      echo -e "### Deleted unsorted bam for ${SAMPLE} ###\n"
       
   else
       echo -e "\n### Sorted bam does not exist for ${SAMPLE} ###"
       echo -e "### Will not delete unsorted bams ###\n"
fi


echo -e "### All done ###"















