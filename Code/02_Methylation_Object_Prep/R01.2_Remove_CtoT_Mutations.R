#### R01.2_Remove_CtoT_Mutations.R ###

########################################################################################################################################

# Created by: Charley
# Date: 10th Oct 2023
# For Evol Apps paper

# Filter out records that match sites of C -> T mutations, identified via JDG's SNP calling pipeline on WGBS data (Revelio + GATK)
# C->T SNP locations should be removed from the analysis, as they do not represent bisulfite-treatment-associated conversions.

# NB. We called SNPs separately in samples from deep vs shallow treatment, as we were experiencing inflation of DMS being called as
# C->T SNPs, more than expected due to not being masked properly by Revelio algorithm. This was done to reduce interference of DMS, 
# as because siblings in each treatment should have the same genetic background, but different methylation patterns induced by the 
# incubation condition experienced. We further selected for the intersection of SNP sites called in both deep and shallow treated hatchlings
# separately to be extra stringent, as it is more likely to be a true SNP if it is found in both groups, since they are genetically siblings

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
DIR_UniteCov <- file.path(DIR, "05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/01_UniteCov_Objects")
DIR_MethDiff <- file.path(DIR, "05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/02_MethDiff_Objects")

########################################################################################################################################

#########################
###### Load files #######
#########################

####### Prep SNP sites file #######

# Choose SNP sites that intersect in both deep and shallow only
SNP_Locations_Deep <- read.table(file = '/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/01_Metadata/WholeGenome_Hatchlings_Reloc3_2021_CtoT_SNPs_Deep.txt', 
                                 sep = '\t', header=TRUE)

SNP_Locations_Shallow <- read.table(file = '/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/01_Metadata/WholeGenome_Hatchlings_Reloc3_2021_CtoT_SNPs_Shallow.txt', 
                                    sep = '\t', header=TRUE)

myPos_Deep = paste(SNP_Locations_Deep$Chrom, SNP_Locations_Deep$Pos)
myPos_Shallow = paste(SNP_Locations_Shallow$Chrom, SNP_Locations_Shallow$Pos)

# Change the Pos vector into a dataframe): called df
df_Deep <- data.frame(chr=sapply(strsplit(myPos_Deep, " "), `[`, 1),
                 start=sapply(strsplit(myPos_Deep, " "), `[`, 2),
                 end=sapply(strsplit(myPos_Deep, " "), `[`, 2))

df_Shallow <- data.frame(chr=sapply(strsplit(myPos_Shallow, " "), `[`, 1),
                 start=sapply(strsplit(myPos_Shallow, " "), `[`, 2),
                 end=sapply(strsplit(myPos_Shallow, " "), `[`, 2))

# Create GRanges object from DMS df
GRanges_SNP_Locations_Deep <- makeGRangesFromDataFrame(df_Deep)
head(GRanges_SNP_Locations_Deep) # Check

GRanges_SNP_Locations_Shallow <- makeGRangesFromDataFrame(df_Shallow)
head(GRanges_SNP_Locations_Shallow) # Check

# Select intersection by overlaps
GRanges_SNP_Locations_Intersect = GRanges_SNP_Locations_Deep[GRanges_SNP_Locations_Deep  %over% GRanges_SNP_Locations_Shallow,]
SNP_Locations_Intersect <- as.data.frame(GRanges_SNP_Locations_Intersect)

####### Load in methylKit objects to filter #######

uniteCov <- readRDS(file.path(DIR_UniteCov, "EvolApps_WholeGenome_uniteCovALL.RDS")) 
methDiff <- readRDS(file.path(DIR_MethDiff, "EvolApps_WholeGenome_ALL_MethDiff.RDS"))

########################################################################################################################################

###################################
### Filter out C -> T mutations ###
###################################

### Select sites that do not overlap with C->T mutations ###

# For uniteCov objects
sub.uniteCov = uniteCov[! as(uniteCov,"GRanges") %over% GRanges_SNP_Locations_Intersect,]
class(sub.uniteCov)
nrow(sub.uniteCov)
nrow(sub.uniteCov)/nrow(uniteCov)*100 # Percentage remaining
nrow(uniteCov) - nrow(sub.uniteCov) # No. of sites filtered out
100 - (nrow(sub.uniteCov)/nrow(uniteCov)*100) # Percentage filtered out 

# For methDiff objects, should be same no. of rows as sub.uniteCov
sub.methDiff = methDiff[! as(methDiff,"GRanges") %over% GRanges_SNP_Locations_Intersect,]
class(sub.methDiff)

### Save ###
saveRDS(sub.uniteCov, file = file.path(DIR_UniteCov, "EvolApps_WholeGenome_uniteCovALL_filt_CtoT_Intersect.RDS"))
saveRDS(sub.methDiff, file = file.path(DIR_MethDiff, "EvolApps_WholeGenome_ALL_MethDiff_filt_CtoT_Intersect.RDS"))
