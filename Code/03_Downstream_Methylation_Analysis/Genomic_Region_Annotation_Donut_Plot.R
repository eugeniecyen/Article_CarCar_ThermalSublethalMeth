#### Global methylation stats ###

########################################################################################################################################

# Created by: Charley
# Date: 8th Nov 2023
# For Evol Apps paper

# Functionally annotate genomic region type of all CpGs i.e. global methylation
# Create single and double donut plots with rings showing location of CpGs at DMS vs global
# Chi-squared test comparing if CpGs are similarly distributed in DMS vs global

########################################################################################################################################

###############################
###### Prep environment #######
###############################

###### Load packages ######

library(tidyverse)
library(methylKit)
library(rtracklayer)
library(GenomicRanges)
library(genomation)


####### Set directories ######

# In local

uniteCov_DIR <- file.path("/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/01_UniteCov_Objects")
DMS_DIR <- file.path("/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/03_DMS_Objects")
AnnoPath <- file.path("/data/SBCS-EizaguirreLab/Turtle_Genome/10_Annotation/06_Functional_Annotation/Scaff_With_Chr0/Final_Annotations/HMMER3_Databases/No_Pathways")

OUTDIR_Tables <- file.path("/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/02_Output/Tables/Revised_Tables/Donut_Genomic_Region")
OUTDIR_Figures <- file.path("/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/02_Output/Figures/Revised_Figures/Donut_Genomic_Region")

########################################################################################################################################

###############################
# Load functional annotations #
###############################

####### Load GFF3 file ######

myannotGff3 <- rtracklayer::readGFF(file.path(AnnoPath, "CarCar_QM_v1_2021_12_Scaff_With_Chr0_Functional_Annotation.gff3"))


####### Load bed12 file ######

# NB Promoters are defined by options at genomation::readTranscriptFeatures function. 
# The default option is to take -1000,+1000bp around the TSS and you can change that. 
# following Heckwolf 2020 and Sagonas 2020, we consider 1500bp upstream and 500 bp downstream

# NB. ReadTranscriptFeatures defines the promoters. If you are using annotation with isoforms, you need to use unique.prom=FALSE option, otherwise
# promoter boundaries will not be assigned a gene name!

myannotBed12=readTranscriptFeatures(file.path(AnnoPath,"CarCar_QM_v1_2021_12_Scaff_With_Chr0_Functional_Annotation.bed12"),
                                     remove.unusual = FALSE, up.flank = 1500, down.flank = 500, unique.prom=FALSE)
                                    
head(myannotBed12)

# Recursively change the gene names to keep only ID
getName <- function(x) {sub(";.*", "", sub(".*ID=", "", x))}

for (i in 1:length(myannotBed12)){
  myannotBed12[[i]]$name <- getName(myannotBed12[[i]]$name)
}

########################################################################################################################################

####################################################
## Add feature type annotations to DMS sites: DMS ##
####################################################

####### Create GRanges object from DMS df ######

# Load in file
DMSobject <- readRDS(file.path(DMS_DIR, "EvolApps_WholeGenome_ALL_DMS_PQLseq_10pc_q0.05_no0q.RDS"))

# Create columns with metadata to include in final csv: names can't be those auto detected by makeGRangesFromDataFrame()
# These are included for checking the final functional annotation, due to previous issue where chrom assignment was mixed up
DMSobject$DMS_pos <- DMSobject$start
DMSobject$original_chrom <- DMSobject$chr

# Create vector of DMS with contig name and position in contig
# plus qvalue and meth diff to add as metadata -> can add column in final df
myPos = paste(DMSobject$chr, DMSobject$end, DMSobject$qvalue, DMSobject$meth.diff, DMSobject$DMS_pos, DMSobject$original_chrom)

# Change the Pos vector into a dataframe): called df
df <- data.frame(chr=sapply(strsplit(myPos, " "), `[`, 1),
                 start=sapply(strsplit(myPos, " "), `[`, 2),
                 end=sapply(strsplit(myPos, " "), `[`, 2),
                 qvalue=sapply(strsplit(myPos, " "), `[`, 3),
                 meth.diff=sapply(strsplit(myPos, " "), `[`, 4),
                 DMS_pos=sapply(strsplit(myPos, " "), `[`, 5),
                 original_chrom=sapply(strsplit(myPos, " "), `[`, 6))

# Convert into GRanges object
GRangeOBJ <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)


####### Annotate all DMS by genic parts ######

# Assign feature type to all DMS first, so we can use it to guide finer-tuned annotation in the next step

# Annotate DMS with genic parts
A = annotateWithGeneParts(target = as(GRangeOBJ,"GRanges"), feature = myannotBed12)


# Assign the feature type to GRangeOBJ
# NB. Need to match names in myannotBed12, otherwise you get issues in next step when splitting annotation by genic vs intergenic
names(myannotBed12) # Check

# Using info in A@members, replace rows where there are promoters, exons, introns or intergenic as appropriate
GRangeOBJ$featureType = ifelse(A@members[,1]==1, "promoters",
                               ifelse(A@members[,2]==1, "exons",
                                      ifelse(A@members[,3]==1, "introns", "intergenic")))

# Create df to look at for troubleshooting (not necessary, but useful)
GRangeOBJ_df <- as.data.frame(GRangeOBJ)

####### Create table of proportions for plotting donut chart ######

# Get no. of DMS per architecture type 
gen_arch_DMS <- as.data.frame(table(GRangeOBJ_df$featureType))
colnames(gen_arch_DMS) <- c("Location", "Frequency")
str(gen_arch_DMS)

# Calculate percentages for each location
gen_arch_DMS$Proportion = round((gen_arch_DMS$Frequency / sum(gen_arch_DMS$Frequency))*100, digits = 2)

# Calculate position of labels: halfway through the cumulative position
gen_arch_DMS <- gen_arch_DMS %>%
  arrange(desc(Location)) %>%
  mutate(ypos = cumsum(Proportion) - 0.5*Proportion)

# Save
setwd(OUTDIR_Tables)
write.csv(gen_arch_DMS, file = "EvolApps_GenomicFeatureType_DMS_ALL_PQLSeq_no0q.csv", row.names = FALSE)

########################################################################################################################################

############################################################
## Add feature type annotations to DMS sites: global meth ##
############################################################

####### Prep dataframes for annotation ####### 

# Load in file
uniteCov <- readRDS(file.path(uniteCov_DIR, "EvolApps_WholeGenome_uniteCovALL_filt_CtoT_Intersect.RDS"))
colnames(uniteCov)

# Create columns with metadata to include in final csv: names can't be those auto detected by makeGRangesFromDataFrame()
# These are included for checking the final functional annotation, due to previous issue where chrom assignment was mixed up
uniteCov$original_pos <- uniteCov$start
uniteCov$original_chrom <- uniteCov$chr

# Create vector with contig name and position in contig
myPos = paste(uniteCov$chr, uniteCov$end, uniteCov$original_pos, uniteCov$original_chrom)

# Change the Pos vector into a dataframe): called df
df <- data.frame(chr=sapply(strsplit(myPos, " "), `[`, 1),
                 start=sapply(strsplit(myPos, " "), `[`, 2),
                 end=sapply(strsplit(myPos, " "), `[`, 2),
                 original_pos=sapply(strsplit(myPos, " "), `[`, 3),
                 original_chrom=sapply(strsplit(myPos, " "), `[`, 4))

# Convert into GRanges object
GRangeOBJ <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)


####### Annotate all CpG by genic parts ######

# Assign feature type to all DMS first, so we can use it to guide finer-tuned annotation in the next step

# Annotate DMS with genic parts
A = annotateWithGeneParts(target = as(GRangeOBJ,"GRanges"), feature = myannotBed12)


# Assign the feature type to GRangeOBJ
# NB. Need to match names in myannotBed12, otherwise you get issues in next step when splitting annotation by genic vs intergenic
names(myannotBed12) # Check

# Using info in A@members, replace rows where there are promoters, exons, introns or intergenic as appropriate
GRangeOBJ$featureType = ifelse(A@members[,1]==1, "promoters",
                               ifelse(A@members[,2]==1, "exons",
                                      ifelse(A@members[,3]==1, "introns", "intergenic")))

# Create df
GRangeOBJ_df <- as.data.frame(GRangeOBJ)

####### Create table of proportions for plotting donut chart ######

# Get no. of sites per architecture type 
gen_arch_ALL <- as.data.frame(table(GRangeOBJ_df$featureType))
colnames(gen_arch_ALL) <- c("Location", "Frequency")
str(gen_arch_ALL)

# Calculate percentages for each location
gen_arch_ALL$Proportion = round((gen_arch_ALL$Frequency / sum(gen_arch_ALL$Frequency))*100, digits = 2)

# Calculate position of labels: halfway through the cumulative position
gen_arch_ALL <- gen_arch_ALL %>%
  arrange(desc(Location)) %>%
  mutate(ypos = cumsum(Proportion) - 0.5*Proportion)

sum(gen_arch_ALL$Frequency) # Check total adds up correctly

# Save
setwd(OUTDIR_Tables)
write.csv(gen_arch_ALL, file = "EvolApps_GenomicFeatureType_ALL.csv", row.names = FALSE)

########################################################################################################################################

######################
## Plot donut chart ##
######################

setwd(OUTDIR_Tables)

# Choose colour palette
viridis_custom_palette <- c("#0d0887", "#5b01a5", "#a82296", "#da5b69", "#fa9e3b")
                                       
# Assign labels 
legend.labels=c("Exon", "Intergenic", "Intron", "Promoter")

####### Single plot with legend inside donut #######

donut_ALL <- ggplot(gen_arch_ALL, aes(x = 2, y = Proportion, fill = Location)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) +
  # geom_text(aes(y = ypos, label = Proportion), color = "white") + # Not enough space for labels on graph
  scale_fill_manual(labels= legend.labels, values = c(viridis_custom_palette[5], viridis_custom_palette[1], viridis_custom_palette[2], viridis_custom_palette[3])) +
  theme_void() +
  theme(legend.text = element_text(size=10), legend.title = element_text(size=10, face="bold"),
        legend.position = c(0.5, 0.5)) +
  xlim(0.5, 3)

donut_ALL

ggsave(file=file.path(OUTDIR_Figures ,"Genomic_Location_Donut_ALL.pdf"), plot=donut_ALL , 
       width=14, height=10, units="in")

donut_DMS <- ggplot(gen_arch_DMS, aes(x = 2, y = Proportion, fill = Location)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) +
  # geom_text(aes(y = ypos, label = Proportion), color = "white") + # Not enough space for labels on graph
  scale_fill_manual(labels= legend.labels, values = c(viridis_custom_palette[5], viridis_custom_palette[1], viridis_custom_palette[2], viridis_custom_palette[3])) +
  theme_void() +
  theme(legend.text = element_text(size=10), legend.title = element_text(size=10, face="bold"),
        legend.position = c(0.5, 0.5)) +
  xlim(0.5, 3)

donut_DMS

ggsave(file=file.path(OUTDIR_Figures ,"Genomic_Location_Donut_DMS.pdf"), plot=donut_DMS , 
       width=14, height=10, units="in")

########################################################################################################################################

######################
## Chi squared test ##
######################

# To statistically test whether observed proportion is different between DMS vs global meth 
# in each genomic location category: exon, intergenic, intron, promoter

### Prep df ###

# Merge into same df
gen_arch_DMS$Type = "DMS"
gen_arch_ALL$Type = "Global"
gen_arch_merge = merge(gen_arch_DMS, gen_arch_ALL, all = TRUE)
gen_arch_merge$Type <- as.factor(gen_arch_merge$Type)

setwd(OUTDIR_Tables)
write.csv(gen_arch_merge, file = "EvolApps_GenomicFeatureType_Revised.csv", row.names = FALSE)

# Create contingency table df
df <- data.frame(Location = as.vector(gen_arch_merge$Location[gen_arch_merge$Type == "DMS"]) ,
                 DMS = as.vector(gen_arch_merge$Proportion[gen_arch_merge$Type == "DMS"]), 
                 Global = as.vector(gen_arch_merge$Proportion[gen_arch_merge$Type == "Global"]))

### Perform chi squared test ###

mydata_chi2_matrix <- df %>% 
  dplyr::select(-Location) %>% 
  as.matrix()

chisq.test(mydata_chi2_matrix)
chisq.test(mydata_chi2_matrix, simulate.p.value = TRUE)





