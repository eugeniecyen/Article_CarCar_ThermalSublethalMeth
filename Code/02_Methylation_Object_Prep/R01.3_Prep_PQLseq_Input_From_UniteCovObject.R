### R01.3_Prep_PQLseq_Input_From_UniteCovObject.R ###

########################################################################################################################################

# Created by: Charley
# Date: 21st Mar 2024
# For Evol Apps

# Prepare methylKit uniteCov object in correct format for PQLseq requirements
# Split by chromosome, so we can run PQLseq as an array by chromosome -> much faster
# NB. Also works for MACAU2 (same authors so same input files)

########################################################################################################################################

###############################
###### Prep environment #######
###############################

####### Load packages ######

library(methylKit)
library(GenomicRanges)
library(genomation)
library(tidyverse)
library(readxl)

####### Set directories ######

DIR_UniteCov <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/01_UniteCov_Objects"
DIR_PQL <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/02_PQLseq_Objects/Input"
DIR_CHROM <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/02_PQLseq_Objects/Input/By_Chrom"

########################################################################################################################################

########################
###### Functions #######
########################

####### 1) Prep RawCountDataSet input file ######
# Methylated read count file per site -> from numCs columns 
# Dataframe where each row is a genomic site, each column is an individual

# Function that deletes and formats columns
PrepNumCsFun <- function(list_of_dfs) {
  for (i in seq_along(list_of_dfs)) {
    # Delete chr column
    list_of_dfs[[i]] <- list_of_dfs[[i]] %>% dplyr::select(-chr)
    
    # Remove columns containing "coverage" in their name
    col_names <- colnames(list_of_dfs[[i]]) # Delete all coverage columns
    cols_to_delete <- grep("^coverage", col_names)
    list_of_dfs[[i]] <- list_of_dfs[[i]][, -cols_to_delete] # Remove identified column
    
    # Rename coverage column names with just number, so it matches between the two input files
    # Could use sample.id, but easier to spot anything out of order with numbers
    list_of_dfs[[i]] <- list_of_dfs[[i]] %>% rename_with(~str_remove(., "^numCs"))
    colnames(list_of_dfs[[i]]) <- paste0("ind", colnames(list_of_dfs[[i]]))
    
    # Save row names for adding back in, as transformation step removes them (if performing)
    row_names <- rownames(list_of_dfs[[i]])
    
    # # If doing count transformation: add +1 to the no. of methylated reads for every value
    # list_of_dfs[[i]] <- as.data.frame(lapply(list_of_dfs[[i]], function(x) if(is.numeric(x)) x + 1 else x))
    
    # Add row names back in
    rownames(list_of_dfs[[i]]) <- row_names
  }
  return(list_of_dfs)
}

####### 2) Prep LibSize input file ######
# Total read count file (i.e. coverage) per site -> from coverage columns
# Dataframe where each row is a genomic site, each column is an individual

# Function that deletes and formats columns
PrepCoverageFun <- function(list_of_dfs) {
  for (i in seq_along(list_of_dfs)) {
    # Delete chr column
    list_of_dfs[[i]] <- list_of_dfs[[i]] %>% dplyr::select(-chr)
    
    # Delete all numCs columns
    col_names <- colnames(list_of_dfs[[i]]) # Delete all coverage columns
    cols_to_delete <- grep("^numCs", col_names)
    list_of_dfs[[i]] <- list_of_dfs[[i]][, -cols_to_delete] # Remove identified column
    
    # Rename coverage column names with just number, so it matches between the two input files
    # Could use sample.id, but easier to spot anything out of order with numbers
    list_of_dfs[[i]] <- list_of_dfs[[i]] %>% rename_with(~str_remove(., "coverage"))
    colnames(list_of_dfs[[i]]) <- paste0("ind", colnames(list_of_dfs[[i]]))
    
    # Save row names for adding back in, as next step removes them
    row_names <- rownames(list_of_dfs[[i]])
    
    # # If doing count transformation: add +2 to the no. of reads for every value
    # list_of_dfs[[i]] <- as.data.frame(lapply(list_of_dfs[[i]], function(x) if(is.numeric(x)) x + 2 else x))
    
    # Add row names back in
    rownames(list_of_dfs[[i]]) <- row_names
  }
  return(list_of_dfs)
}

########################################################################################################################################

##########################
###### Prep inputs #######
##########################

uniteCov <- readRDS(file.path(DIR_UniteCov, "EvolApps_WholeGenome_uniteCovALL_filt_CtoT_Intersect.RDS"))

####### Prep file columns common to both ######

# Convert uniteCov object to df, as methylBase object doesn't allow column subsetting
uniteCov_df <- as.data.frame(getData(uniteCov))

# Rename rows to site names for record
rownames(uniteCov_df) <- paste(uniteCov_df$chr, uniteCov_df$start, sep = ".")

# Remove site, start, end and strand columns (keep chr column)
uniteCov_df <- uniteCov_df %>% dplyr::select(-start, -end, -strand)

# Delete all numTs columns: not needed by either
col_names <- colnames(uniteCov_df)
cols_to_delete <- grep("^numTs", col_names)
uniteCov_df <- uniteCov_df[, -cols_to_delete]


####### Split into separate dataframes per chromosome ######

# So we can run PQLseq in an array -> faster
uniteCov_df$chr <- as.factor(uniteCov_df$chr)

list_chr_df <- split(uniteCov_df, uniteCov_df$chr)
str(list_chr_df)
names(list_chr_df)


####### 1) Prep RawCountDataSet input file ######

# Methylated read count file per site -> from numCs columns 
# Dataframe where each row is a genomic site, each column is an individual

list_chr_df_numCs <- list_chr_df
list_chr_df_numCs <- PrepNumCsFun(list_chr_df_numCs)

# Save each to separate RDS
lapply(names(list_chr_df_numCs), function(i)
  saveRDS(list_chr_df_numCs[[i]], file.path(DIR_CHROM, paste0("PQLseq_Input_NumCs_CountTransform1_", i , ".RDS"))))


####### 2) Prep LibSize input file ######

# Total read count file (i.e. coverage) per site -> from coverage columns
# Dataframe where each row is a genomic site, each column is an individual

list_chr_df_coverage <- list_chr_df
list_chr_df_coverage <- PrepCoverageFun(list_chr_df_coverage)

# Save each to separate RDS
lapply(names(list_chr_df_coverage), function(i)
  saveRDS(list_chr_df_coverage[[i]], file.path(DIR_CHROM, paste0("PQLseq_Input_Coverage_CountTransform2_", i , ".RDS"))))


