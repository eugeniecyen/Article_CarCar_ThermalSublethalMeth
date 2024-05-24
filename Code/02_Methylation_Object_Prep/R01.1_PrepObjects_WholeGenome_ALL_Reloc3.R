### R01.1_PrepObjects_WholeGenome_ALL_Reloc3.R ###

########################################################################################################################################

# Created by Charley, 30th Sept 2023
# To be submitted on Apocrita with corresponding bash script

# Working with reloc 3 data of 2021 hatchery experiment: for Evol Apps paper
# n=10 clutches split into deep and shallow, n=20 per treatment

# (1) Creates methylBase (uniteCov) objects
# Works with destranded, whole genome methylation calls
# 2) Calculate methylation differences for each site across whole genome

########################################################################################################################################

###############################
###### Prep environment #######
###############################

print(paste0("############### Started at ", Sys.time(), " ###############"))
start_time <- Sys.time()

####### Pull arguments set in parent bash submission script ######

args <- commandArgs(trailingOnly = T)

# Set number of cores available (${NSLOTS} in parent bash script)
Num_cores = args[1]

print("Num_cores:")
print(args)

####### Load packages ######

library(methylKit)
library(readxl)
library(plyr)
library(dplyr)
library(ggplot2)

####### Set directories ######

DIR <- file.path("/data/SBCS-EizaguirreLab/Turtle_WGBS")
dataPath <- file.path(DIR, "04_Bismark_Methylation_Calls/destranded_methylation_calls/Whole_Genome")
metadataPath <- file.path(DIR, "05_MethylKit/04_Evol_Apps_Project/01_Working/01_Metadata")
OUTDIR_UniteCov <- file.path(DIR, "05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/01_UniteCov_Objects")
OUTDIR_MethDiff <- file.path(DIR, "05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/02_MethDiff_Objects")
OUTDIR_DMS <- file.path(DIR, "05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/03_DMS_Objects")

########################################################################################################################################

#################################
# Create methylRawList object #
#################################

###### Add metadata ######

print("##### Reading metadata #####")

# Read in Excel file containing metadata
metadata <- readxl::read_xlsx(file.path(metadataPath, "Metadata_Hatchlings_Reloc3.xlsx"))

# MethylKit only allows numeric treatments -> make new variables trt_NUM and Mat_ID
# trt_NUM - where Hatchling = 1 and Mother = 2
# Mat_ID - where SLL is removed
metadata <- mutate(metadata,
                   trt_NUM = case_when(
                     Depth == "Shallow" ~ 1, 
                     Depth == "Deep" ~ 2 
                   )) %>%
  mutate(metadata,
         Mat_ID = as.numeric(substr(metadata$Mother, 4, 6)))

# Check metadata
print("Printing head of metadata:")
head(metadata)


###### Create list of Bismark meth call files to load ######

print("##### Creating list of meth call files to analyse #####")

# Make a list of all destranded methylation call files
temp = list.files(path=dataPath,
                  pattern = ".CpG_merged.cov.gz",
                  full.names = T,
                  recursive = T)

# Subset list for only sample IDs that match those in metadata
temp <- temp[ grep(pattern=paste(metadata$Sample, collapse = "|"), x = temp) ]
# James' explanation for why you need to use collapse="|": 
# When you want to grep multiple things you use | to separate like bash (e.g. "Fogo|BoaVista"). 
# When the things you want to search are in a vector or variable, if you grep df$var, I think it greps the whole thing. 
# But to grep all of the things in the vector one at a time, you collapse (= separate) with | and it treats them separately

# Check correct files are included in list
print("No. of files:")
length(temp) # check number of files

print("Printing head of file list:")
head(temp) # check 1st few files on the list


###### Make a methylRawList object ######

# NB. Can take a while

print("##### Creating methylRawList object with min coverage filter of 5X #####")

# Create object with min coverage filter of 5X
myobj.mincov5=methylKit::methRead(as.list(temp),
                                  pipeline='bismarkCoverage',
                                  mincov=5,
                                  sample.id=as.list(metadata$Sample),
                                  assembly="CarCar_QM_v1_2021_12_Scaff_With_Chr0.fasta",
                                  treatment=metadata$trt_NUM,
                                  context="CpG")

print("##### Finished #####")

########################################################################################################################################

###############################
## Filtering and normalising ##
###############################

####### Filter based on coverage ######

# It's standard practice to filter samples based on coverage.
# Discard bases with very high coverage, as it could be due to PCR bias
# Discard bases with low coverage, as we want high enough to increase statistical power and be confident in meth values (already done)

print("##### Filtering methylRawList object by maximum >99.9% and minimum <5X coverage #####")

filt.myobj.mincov5=filterByCoverage(myobj.mincov5, lo.count=5, lo.perc=NULL,
                                    hi.count=NULL, hi.perc=99.9)

print("##### Finished #####")


####### Normalise coverage #######

# Not necessarily needed if coverage is similar across samples, but good to do in case
# Together with removing extreme coverage, will help reduce the bias in  statistical tests 
# that might occur due to systematic over-sampling of reads in certain samples

print("##### Normalising coverage #####")

normFil.myobj.mincov5=normalizeCoverage(filt.myobj.mincov5)

print("##### Finished #####")


####### Remove files that are no longer needed to save memory #######

rm(myobj.mincov5)
rm(filt.myobj.mincov5)


########################################################################################################################################

#####################
### Merge samples ###
#####################

# In order to do further analysis, we need to get the bases covered by reads in all samples.
# We will merge all samples into a methylBase object for all base-pair locations covered in all samples.

print("##### Creating uniteCovALL object with CpG present in all individuals per treatment #####")

uniteCovALL = methylKit::unite(normFil.myobj.mincov5, destrand=FALSE, mc.cores = Num_cores)
uniteCovALL = as(uniteCovALL,"methylBase")

# Check
print("Printing head of uniteCovALL:")
head(uniteCovALL)

print("No. of rows in uniteCovALL is:")
nrow(uniteCovALL)

# Save
print("##### Saving #####")

saveRDS(uniteCovALL, file = file.path(OUTDIR_UniteCov, "EvolApps_WholeGenome_uniteCovALL.RDS"))

print("##### Finished saving #####")

# Remove to save space
rm(normFil.myobj.mincov5)

########################################################################################################################################

##############################
### Create MethDiff object ###
##############################

# Create a methylKit DiffMeth object containing mean % methylation difference between treatments
# NB. We will be ignoring p/q values generated and use PQLseq instead

MethDiff_Depth <- calculateDiffMeth(uniteCovALL, mc.cores = Num_cores)

print("Printing head of MethDiff_Depth:")
head(MethDiff_Depth)
print("No. of rows in MethDiff_Depth is:")
nrow(MethDiff_Depth)

# Save file
print("##### Saving DiffMeth object #####")
saveRDS(MethDiff_Depth, file = file.path(OUTDIR_MethDiff, "EvolApps_WholeGenome_ALL_MethDiff.RDS"))
print("##### Finished saving DiffMeth object #####")

########################################################################################################################################

end_time <- Sys.time()

runtime <- end_time - start_time

print("Total run time:")
print(runtime)










