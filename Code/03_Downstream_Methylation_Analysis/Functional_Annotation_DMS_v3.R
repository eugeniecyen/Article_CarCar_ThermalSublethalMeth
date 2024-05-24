### Functional annotation v3 ###

########################################################################################################################################

# Created by: Charley, adapted from Alice
# Date: 20th Oct 2023
# For Evol Apps

# Functionally annotate DMS with genomic region and associated gene 

########################################################################################################################################

########################
### Prep environment ###
########################

###### Load packages ######

library(tidyverse)
library(methylKit)
library(genomation)
library(rtracklayer)
library(GenomicRanges)

####### Set directories ######

DMS_DIR <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/03_DMS_Objects"
AnnoPath <- "/data/SBCS-EizaguirreLab/Turtle_Genome/10_Annotation/06_Functional_Annotation/Scaff_With_Chr0/Final_Annotations/HMMER3_Databases/No_Pathways/longest_isoforms_only"
OUTDIR <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/02_Output/Tables/Functional_Annotation"

########################################################################################################################################

###################################
### Load functional annotations ###
###################################

####### Load GFF3 file ######

myannotGff3 <- rtracklayer::readGFF(file.path(AnnoPath, "CarCar_QM_v1_2021_12_Scaff_With_Chr0_Functional_Annotation_LongestIsoform.gff3"))


####### Load bed12 file ######

# NB Promoters are defined by options at genomation::readTranscriptFeatures function. 
# The default option is to take -1000,+1000bp around the TSS and you can change that. 
# following Heckwolf 2020 and Sagonas 2020, we consider 1500bp upstream and 500 bp downstream

# NB. ReadTranscriptFeatures defines the promoters. For some reason for this annotation, you need to use unique.prom=FALSE option, otherwise 
# promoter boundaries will not be assigned a gene name! This is the case whether annotation has all isoforms or just longest
myannotBed12=readTranscriptFeatures(file.path(AnnoPath, "CarCar_QM_v1_2021_12_Scaff_With_Chr0_Functional_Annotation_LongestIsoform.bed12"),
                                    remove.unusual = FALSE, up.flank = 1500, down.flank = 500, unique.prom=FALSE)

head(myannotBed12)

# Recursively change the gene names to keep only ID
getName <- function(x) {sub(";.*", "", sub(".*ID=", "", x))}

for (i in 1:length(myannotBed12)){
  myannotBed12[[i]]$name <- getName(myannotBed12[[i]]$name)
}

########################################################################################################################################

#################################################
### Add feature type annotations to DMS sites ###
#################################################

####### Create GRanges object from DMS df ######

# Load in file
DMSobject <- readRDS(file.path(DMS_DIR, "EvolApps_WholeGenome_ALL_DMS_PQLseq_10pc_q0.05.RDS"))

# For intersect between old and PQLseq DMS
# DMSobject = new_methdiff[as(new_methdiff,"GRanges") %over% Intersect_q0.05,]

# These are included for checking the final functional annotation, due to previous issue where chrom assignment was mixed up
# NB. column names can't be the same those auto detected by makeGRangesFromDataFrame()
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
# NB. Need to match names(myannotBed12), otherwise you get issues in next step when splitting annotation by genic vs intergenic
# Using info in A@members, replace rows where there are promoters, exons, introns or intergenic as appropriate
GRangeOBJ$featureType = ifelse(A@members[,1]==1, "promoters",
                               ifelse(A@members[,2]==1, "exons",
                                      ifelse(A@members[,3]==1, "introns", "intergenic")))

########################################################################################################################################

###########################################
### Assign gene identities to DMS sites ###
###########################################

####### Case 1: the feature is GENIC -> get the annotation by intersection with bed file #######

# As DMS is on a gene (promoter, intron or exon), we now don't filter by distance from TSS for all DMS beforehand, otherwise 
# you could have a DMS at the end of a large gene that could be filtered out by accident

# Subset GRangesOBJ to retain DMS sites on genes only
GRangeOBJ1=GRangeOBJ[!GRangeOBJ$featureType %in% "intergenic"]

# Add empty column called geneInfo to fill with gene IDs
GRangeOBJ1$feature.name <- NA

add_geneInfo_genic <- function(x, GRangeOBJ, annotBed12){
  ov = GenomicRanges::findOverlaps(
    annotBed12[[x]],
    GRangeOBJ[GRangeOBJ$featureType %in% x,])
  ## Add gene annotation to subject GRanges (i.e. left join)
  mcols(GRangeOBJ[GRangeOBJ$featureType %in% x,])[subjectHits(ov), "feature.name"] = mcols(annotBed12[[x]])[queryHits(ov), "name"]
  return(GRangeOBJ)
}

GRangeOBJ_ex = add_geneInfo_genic("exons", GRangeOBJ1, myannotBed12) # Add gene names for sites on exons
GRangeOBJ_ex_in = add_geneInfo_genic("introns", GRangeOBJ_ex, myannotBed12) # Add gene names for sites on introns
GRangeOBJ_genic = add_geneInfo_genic("promoters", GRangeOBJ_ex_in, myannotBed12) # Add gene names for sites on promoters

### Case 2: the feature is INTERGENIC: get the annotation by proximity to nearest TSS
# Heckwolf 2020: To be associated to a gene, the DMS had to be either inside the gene or,
# if intergenic, not further than 10 kb away from the TSS.

# Subset GRangesOBJ to retain intergenic sites only
GRangeOBJ2=GRangeOBJ[GRangeOBJ$featureType %in% "intergenic"]

a = annotateWithGeneParts(target = as(GRangeOBJ2,"GRanges"), feature = myannotBed12)

# Filter out sites that are further than 10 kb away from the TSS
rows2rm = which((a@dist.to.TSS$dist.to.feature>10000 | a@dist.to.TSS$dist.to.feature< -10000) &
                  rowSums(a@members) %in% 0)
if (is_empty(rows2rm)){  GRangeOBJ2 = GRangeOBJ2
} else { GRangeOBJ2 = GRangeOBJ2[-rows2rm,] }

# Re-annotate the subsetted object
b = annotateWithGeneParts(as(GRangeOBJ2,"GRanges"), myannotBed12)

# Get genes associated with these TSS
c = getAssociationWithTSS(b)

# Add these gene ID associations to GRangeOBJ2
GRangeOBJ2$feature.name=c$feature.name
GRangeOBJ_intergenic <- GRangeOBJ2

####### Merge the 2 cases back together into single GRangesOBJ and clean up #######

# Merge
myGRangeOBJ=c(GRangeOBJ_genic, GRangeOBJ_intergenic)

# Remove columns for original pos/chrom after checking it still matches!
myGRangeOBJ$original_pos <- NULL
myGRangeOBJ$original_chrom <- NULL

########################################################################################################################################

#####################################################
### Collect full annotations for identified genes ###
#####################################################

####### Extract full info for genes from myannotGff3 present in myGRangeOBJ, and merge into 1 df ####### 

annotDF = merge(myGRangeOBJ %>% data.frame() %>% dplyr::rename(chrom = seqnames), # Rename myGRangeOBJ column seqnames to chrom
                data.frame(subset(myannotGff3, ID %in% myGRangeOBJ$feature.name)) %>% # Subset annotGff3 for only genes associated to DMS in myGRangeOBJ
                  dplyr::select(c("ID", "Note","Ontology_term", "start","end","strand", "Parent"))  %>% # Extract only  desired fields then rename
                  dplyr::rename(feature.name=ID, start.gene = start, end.gene = end), by="feature.name")

####### Add column with no. of sites per gene #######
annotDF = merge(annotDF,
                data.frame(table(annotDF$feature.name)) %>% dplyr::rename(feature.name=Var1, nDMSperGene=Freq))

# Add extra info (nbr CpG per gene length, gene length, chrom name)
annotDF = annotDF  %>%
  mutate(geneLengthkb = (end.gene - start.gene)/1000, nDMSperGenekb = round(nDMSperGene/geneLengthkb,2))

# Unlist note (list -> character)
annotDF$Note = unlist(annotDF$Note)

########################################################################################################################################

################################
### Export table to csv file ###
################################

setwd(OUTDIR)

# Delete unneeded columns
anotDMS_minimal <- dplyr::select(annotDF, c('feature.name', 'chrom', 'start.gene', 'end.gene', 'nDMSperGene', 'DMS_pos', 'featureType',
                                            'Note', 'Ontology_term', 'meth.diff', 'qvalue'))

# Homogenise column names
anotDMS_minimal <- anotDMS_minimal %>% 
  dplyr::rename(
    gene.id = feature.name,
    DMS.pos = DMS_pos,
    feature.type = featureType,
    note = Note,
    ontology.term = Ontology_term,
  )

anotDMS_minimal <- apply(anotDMS_minimal,2,as.character)

# Save
write.csv(anotDMS_minimal, file = "EvolApps_ALL_DMS_PQLseq_10pc_q0.05.csv", row.names = FALSE)


