### R01.5_Run_PQLseq_ChromArray.R ###

########################################################################################################################################

# Created by: Charley
# Date: 19th Mar 2024
# For submission as an Apocrita job
# For Evol Apps

# Script running PQLSeq as an array by chromosome
# Manual: https://cran.hafro.is/web/packages/PQLseq/PQLseq.pdf 
# Paper: https://academic.oup.com/bioinformatics/article/35/3/487/5055584 

# NB. Changes required to prevent this package from setting off ITS alarms for using too many cores
# Solution provided by Sam Lawson @ ITS, 22/07/2024:
# Edit the PQLseq function: Change to FORKING mode: cl <- makeCluster(numCore, type = "FORK", outfile = "par_log.txt")
# Add to the top of this R script:
# 1) Add Sys.setenv(MKL_NUM_THREADS = 1)
# 2) Install RhpcBLASctl and run library(RhpcBLASctl); blas_set_num_threads(1); omp_set_num_threads(1)

########################################################################################################################################

###############################
###### Prep environment #######
###############################

library(foreach)
library(doParallel)
library(parallel)
library(Rcpp)
library(RcppArmadillo)
library(PQLseq)
library(tidyverse)
library(RhpcBLASctl) # To deal with alarming nodes issue

####### Set multi-threading parameters to avoid alarming nodes ######

Sys.setenv(MKL_NUM_THREADS = 1)
blas_set_num_threads(1)
omp_set_num_threads(1)

####### Set directories ######

DIR_PQL <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/02_PQLseq_Objects/Input"
DIR_CHROM <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/02_PQLseq_Objects/Input/By_Chrom/No_Transform"
DIR_FIT <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/02_PQLseq_Objects/Output_PQLseq/By_Chrom/No_Transform"

####### Pull arguments set in parent bash submission script ######

args <- commandArgs(trailingOnly = T)

# Set desired chromosome/s using task array ID, as specified in parent bash script
# Chrom number (i.e. array task ID) must be 1st arg in parent bash script
# Note this is SGE_TASK_ID - 1 because we have chromosomes 0-28 but the array must start at 1 -> (1:29)

print("Arg1 and Arg2:")
print(args)
chr = as.numeric(args[1])-1

####### Load custom function ######

# Required to prevent over-threading
# Also altered package code to auto-detect no. of cores to use with nslots -> don't need to supply as additional argument
source("/data/SBCS-EizaguirreLab/Turtle_WGBS/00_Scripts/Functions/Custom_PQLseq_Function.R")

# Set up environment for custom function: 
# https://stackoverflow.com/questions/24331690/modify-package-function
environment(custom_pqlseq) <- asNamespace('PQLseq')
assignInNamespace("pqlseq", custom_pqlseq, ns = "PQLseq")


start_time <- Sys.time()
print(paste0("##### Starting chromosome ", chr, " at ", Sys.time(), " #####")) 

########################################################################################################################################

###################
###### Run ########
###################

### Load files ###
print("##### Reading files #####")

# Load methylation count input files
numCs_df <- readRDS(file.path(DIR_CHROM, paste0("PQLseq_Input_NumCs_SLK063_ragtag_chr", chr, ".RDS" )))

print("Head of numCs_df:")
head(numCs_df)
print(paste0("No. of rows: ", nrow(numCs_df)))

coverage_df <- readRDS(file.path(DIR_CHROM, paste0("PQLseq_Input_Coverage_SLK063_ragtag_chr", chr, ".RDS" )))

print("Head of coverage_df:")
head(coverage_df)
print(paste0("No. of rows: ", nrow(coverage_df)))

# Subset rows for testing
# numCs_df <- numCs_df[1:100,]
# coverage_df <- coverage_df[1:100,]

# Load genetic relatedness matrix
relatedness_mat <- read.table(file.path(DIR_PQL, "Evol_Apps_ALL_Theoretical_Relatedness_Matrix.txt"))

print("Head of relatedness matrix:")
head(relatedness_mat)

print("##### Adding Phenotypes vector #####")
treatment.no <- read.table(file.path(DIR_PQL, "Phenotypes_Input_Treatment.txt"), header = F)

print("Phenotypes:")
print(treatment.no)

### Run PQLseq ###
print("##### Running PQLseq in BMM mode #####")

fit = custom_pqlseq(RawCountDataSet = numCs_df , Phenotypes = treatment.no ,
                     RelatednessMatrix = relatedness_mat , LibSize = coverage_df , fit.model = "BMM", verbose=TRUE)

print("Head of fit:")
print(head(fit))
print(paste0("No. of rows: ", nrow(fit)))

### Save ###
print("##### Saving #####")
saveRDS(fit, file.path(DIR_FIT, paste0("PQLseq_fit_chr", chr, ".RDS")))

print(paste0("Saved fit table to: ", DIR_FIT))

########################################################################################################################################

end_time <- Sys.time()

print(paste0("############### All finished for chromosome ", chr, " at ", end_time, " ###############"))

runtime <- end_time - start_time

print("Total run time:")
print(runtime)



