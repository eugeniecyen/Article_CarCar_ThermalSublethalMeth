#### Global methylation stats ###

########################################################################################################################################

# Created by: Charley
# Date: 16th Oct 2023
# For Evol Apps paper

# 1. Prep global methylation objects
# Split uniteCovALL object by samples in deep (n=10) and shallow (n=10)
# Create methylation matrices

# 2. Analyse global methylation patterns
# Filter for sites that are methylated (>0.7 and non-zero) and count
# Lmers and boxplots to compare global methylation level in deep vs shallow

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
library(lmerTest)
library(lme4)

####### Set directories ######

DIR <- file.path("/data/SBCS-EizaguirreLab/Turtle_WGBS")

metadataPath <- file.path(DIR, "05_MethylKit/04_Evol_Apps_Project/01_Working/01_Metadata")
DIR_UniteCov <- file.path(DIR, "05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/01_UniteCov_Objects")
OUTDIR <- file.path(DIR, "05_MethylKit/04_Evol_Apps_Project/02_Output/Tables/Original_Tables/Global_Meth_Per_Ind")
DIR_Figures <- file.path(DIR, "05_MethylKit/04_Evol_Apps_Project/02_Output/Figures/Global_Meth_Per_Ind")

###### Add metadata ######

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

####### Load in uniteCov object #######

uniteCov <- readRDS(file.path(DIR_UniteCov, "EvolApps_WholeGenome_uniteCovALL_filt_CtoT_Intersect.RDS"))

########################################################################################################################################

#############################################
###### Split uniteCovALL by treatment #######
#############################################

####### Prepare deep and shallow metadata #######

# Get sample and treatment vectors for shallow
shallow <- metadata[metadata$Depth == 'Shallow',]
shallow.ids <- as.vector(shallow$Sample)
shallow.treatments <- as.vector(shallow$trt_NUM)

# Get sample and treatment vectors for deep
deep <- metadata[metadata$Depth == 'Deep',]
deep.ids <- as.vector(deep$Sample)
deep.treatments <- as.vector(deep$trt_NUM)

####### Subset uniteCovALL by treatment ####### 

# Shallow
uniteCovShallow <- methylKit::reorganize(uniteCov, sample.ids=shallow.ids , treatment = shallow.treatments)
# saveRDS(uniteCovShallow, file = file.path(DIR_UniteCov, "EvolApps_WholeGenome_uniteCovALL_filt_CtoT_Intersect_Shallow.RDS"))

# Deep
uniteCovDeep <- methylKit::reorganize(uniteCov, sample.ids=deep.ids , treatment = deep.treatments)
# saveRDS(uniteCovDeep, file = file.path(DIR_UniteCov, "EvolApps_WholeGenome_uniteCovALL_filt_CtoT_Intersect_Deep.RDS"))

########################################################################################################################################

##########################################
###### Create methylation matrices #######
##########################################

# Create percent methylation matrices per uniteCov object for deep and shallow

# Create percentage methylation matrix
perc_methy_Shallow=percMethylation(uniteCovShallow)
perc_methy_Deep=percMethylation(uniteCovDeep)

# Convert to df
perc_methy_Shallow_df <- as.data.frame(perc_methy_Shallow)
perc_methy_Deep_df <- as.data.frame(perc_methy_Deep)

# Count number of cells that are 0 in the entire df
sum(perc_methy_Shallow_df == "0") # 1,432,339 sites
sum(perc_methy_Deep_df == "0") # 1,445,430 sites 

# Save
# Shallow
# saveRDS(perc_methy_Shallow, file = file.path(OUTDIR, "EvolApps_WholeGenome_uniteCovALL_filt_CtoT_Intersect_Shallow_PercMethMatrix.RDS"))
# Deep
# saveRDS(perc_methy_Deep, file = file.path(OUTDIR, "EvolApps_WholeGenome_uniteCovALL_filt_CtoT_Intersect_Deep_PercMethMatrix.RDS"))


########################################################################################################################################

################################################################################
###### Filter matrix to count elements with certain meth, per individual #######
################################################################################

####### Create df with sums of elements that meet meth criteria per individual ####### 

### Function to create df with sums of elements that meet meth criteria per individual ###
myMethCountFun <- function(obj.name, treatment){

  # Set up empty df to bind to
  temp_df <-  as.data.frame(matrix(ncol = 4, nrow=0))
  colnames(temp_df) <- c("sample.ID", "count.nonzero.meth", "count.zero.meth", "count.0.7.meth") # Name columns
  
  # Loop through all individuals and bind to temp_df
  # NB. problems with using colname to subset column as it's numerical -> using indices works better
  for (i in 1:length(colnames(obj.name))) {
    print(i)
    print(colnames(obj.name)[i])
  
  # Extract sample ID
  sample.ID <- colnames(obj.name[i])

  # Extract no. of elements that are non-zero
  count.nonzero.meth <- sum(obj.name[,i] > 0)
  print(count.nonzero.meth)
  
  # Extract no. of elements that are zero
  count.zero.meth <- sum(obj.name[,i] == 0)
  print(count.zero.meth)
  
  # Extract no. of elements greater than 0.7
  count.0.7.meth <- sum(obj.name[,i] > 70) # NB. Values in % -> 70 not 0.7
  print(count.0.7.meth)
  
  # Put values them together in a vector
  values <- c(sample.ID, count.nonzero.meth, count.zero.meth, count.0.7.meth)
  
  # Create df of each individual
  temp_df_ind <-  as.data.frame(matrix(ncol = 4, nrow=0)) # Create empty df to bind to
  df_per_ind <- rbind(temp_df_ind, values) # Bind
  colnames(df_per_ind) <- c("sample.ID", "count.nonzero.meth", "count.zero.meth", "count.0.7.meth") # Name columns
  
  # Fill up outer df
  temp_df <- rbind(temp_df, df_per_ind)
}

  # Add column with treatment to temp_df
  if (treatment == "shallow"){
    temp_df$Treatment <- rep("Shallow")
  } else if (treatment == "deep"){
    temp_df$Treatment <- rep("Deep")
  }
  
  return(temp_df)
}

####### Run ####### 

temp_df_shallow <- myMethCountFun(perc_methy_Shallow_df, treatment = "shallow")
temp_df_deep <- myMethCountFun(perc_methy_Deep_df, treatment = "deep")

# Bind together into 1 df
meth_per_ind <- rbind(temp_df_shallow, temp_df_deep)

# Change count columns to numeric
meth_per_ind$count.nonzero.meth <- as.numeric(meth_per_ind$count.nonzero.meth)
meth_per_ind$count.zero.meth <- as.numeric(meth_per_ind$count.zero.meth)
meth_per_ind$count.0.7.meth <- as.numeric(meth_per_ind$count.0.7.meth)
str(meth_per_ind) # Check 

# Add columns with ratios of meth
# NB. not really needed as you divide by the same number, but could be more comparable in future with other datasets
meth_per_ind$ratio.nonzero.meth <- (meth_per_ind$count.nonzero.meth) / nrow(perc_methy_Shallow_df) # Same total sites in deep and shallow
meth_per_ind$ratio.zero.meth <- (meth_per_ind$count.zero.meth) / nrow(perc_methy_Shallow_df) # Same total sites in deep and shallow
meth_per_ind$ratio.0.7.meth <- (meth_per_ind$count.0.7.meth) / nrow(perc_methy_Shallow_df) # Same total sites in deep and shallow

# Add maternal ID to table based on matched sample ID
meth_per_ind$maternal_ID <- metadata$Mother[match(meth_per_ind$sample.ID, metadata$Sample)]

# Save
write.csv(meth_per_ind, file = file.path(OUTDIR, "EvolApps_ALL_filt_CtoT_Intersect_GlobalMeth_Per_Ind.csv"), row.names=FALSE) 


########################################################################################################################################

####################################
###### Final plots and stats #######
####################################

# meth_per_ind <- read.csv(file = file.path(OUTDIR, "EvolApps_ALL_filt_CtoT_Intersect_GlobalMeth_Per_Ind.csv"))

####### Plot ####### 

# Non zero
nonzero.meth.boxplot <- ggplot(meth_per_ind, aes(x=Treatment, y=count.nonzero.meth, fill=Treatment)) +
  geom_boxplot() +
  labs(x="Treatment", y="Count of CpGs with non-zero methylation") +
  scale_fill_manual(values=c("#0D0887","#FA9E3B")) +
  theme_classic() +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
        legend.text = element_text(size=14), legend.title = element_text(size=16)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))

nonzero.meth.boxplot

ggsave(file=file.path(DIR_Figures ,"GlobalMeth_NonZeroMethCount.pdf"), plot=nonzero.meth.boxplot , 
       width=6, height=7, units="in")

# >70%
high.meth.boxplot <- ggplot(meth_per_ind, aes(x=Treatment, y=count.0.7.meth, fill=Treatment)) +
  geom_boxplot() +
  labs(x="Treatment", y="Count of CpGs with methylation > 0.7") +
  scale_fill_manual(values=c("#0D0887","#FA9E3B")) +
  theme_classic() +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
        legend.text = element_text(size=14), legend.title = element_text(size=16)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))

high.meth.boxplot

ggsave(file=file.path(DIR_Figures ,"GlobalMeth_0.7MethCount.pdf"), plot=high.meth.boxplot , 
       width=6, height=7, units="in")


####### Stats ####### 

model <- lmerTest::lmer(count.nonzero.meth ~ Treatment + (1|maternal_ID), meth_per_ind)
summary(model)
anova(model)

model <- lmerTest::lmer(count.0.7.meth ~ Treatment + (1|maternal_ID), meth_per_ind)
summary(model)
anova(model)




