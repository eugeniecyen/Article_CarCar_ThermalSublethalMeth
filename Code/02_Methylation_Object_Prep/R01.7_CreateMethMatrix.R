#### R01.7_CreateMethMatrix.R ###

########################################################################################################################################

# Created by: Charley
# Date: 4th Oct 2023
# For Evol Apps

# Create matrix of methylation values per individual per site for downstream analyses

########################################################################################################################################

###############################
###### Prep environment #######
###############################

####### Load packages ######

library(methylKit)
library(readxl)
library(plyr)
library(dplyr)
library(tidyverse)
library(GenomicRanges)

####### Set directories ######

DIR <- file.path("/data/SBCS-EizaguirreLab/Turtle_WGBS")
DIR_MethMatrix <- file.path(DIR, "05_MethylKit/04_Evol_Apps_Project/02_Output/Tables/Revised_Tables/Candidate_DMS")
metadataPath <- file.path(DIR, "05_MethylKit/04_Evol_Apps_Project/01_Working/01_Metadata")

########################################################################################################################################

###################################
###### Subset candidate DMS #######
###################################

# 29 DMS of interest (>20% meth difference between treatments)
candidate_DMS <- read.csv(file.path(metadataPath, "EvolApps_Candidate_DMS_Revised.csv"))

# Create vector of DMS with just contig name and position in contig
myPos = paste(candidate_DMS$chrom, candidate_DMS$DMS.pos)

# Change the Pos vector into a dataframe): called df
df <- data.frame(chr=sapply(strsplit(myPos, " "), `[`, 1),
                 start=sapply(strsplit(myPos, " "), `[`, 2),
                 end=sapply(strsplit(myPos, " "), `[`, 2))

# Create GRanges object from DMS df
GRanges_candidate_DMS <- makeGRangesFromDataFrame(df)

### Use methylKit select by overlap ###

# Selects records that match list of DMS locations
uniteCov_candidate_DMS <- selectByOverlap(uniteCov, GRanges_candidate_DMS)

# Save
saveRDS(uniteCov_candidate_DMS, file = file.path(DIR_UniteCov, "EvolApps_WholeGenome_uniteCovDM_ALL_CandidateDMS_Revised.RDS"))

# Create percent methylation matrix
perc_methy <- percMethylation(uniteCov_candidate_DMS)
perc_methy_df <- as.data.frame(perc_methy)

# Add 1st column with custom DMS names
perc_methy_df_cols <- add_column(perc_methy_df, DMS.name = candidate_DMS$DMS.name, .before = colnames(perc_methy_df)[1])
# rownames(perc_methy_df_cols) <- paste(uniteCov_candidate_DMS$chr, uniteCov_candidate_DMS$start, sep = ".") # Add row names with chr.pos

# Save 
write.csv(perc_methy_df_cols, file = file.path(DIR_MethMatrix, "EvolApps_WholeGenome_PercMethMatrix_ALL_CandidateDMS_Revised.csv"), row.names = FALSE) 


