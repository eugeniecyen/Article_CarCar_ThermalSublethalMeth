### NMDS Cluster analysis ###

########################################################################################################################################

# Created by Alice, adapted and added to by Charley
# Date: 3rd May 2024
# For reloc 3 for Evol Apps paper

# Clustering: NMDS and dendrogram plots for global and differential methylation
# Run Adonis tests
# DMS volcano and stats

########################################################################################################################################

###############################
###### Prep environment #######
###############################

###### Load packages ######

library(vegan)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(MASS)
library(plyr)
library(methylKit)
library(dplyr)
library(dendextend)
library(qualpalr)
library(GenomicRanges)
library(viridis)


####### Set directories ######

DMS_DIR <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/03_DMS_Objects"
UniteCov_DIR <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/01_UniteCov_Objects"
metadataPath <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/01_Metadata"


####### Set colour palette ######

# Use viridis magma palette: select discrete colours to use for plotting a few discrete categories
# https://waldyrious.net/viridis-palette-generator/ 
viridis_custom_palette <- c("#0d0887", "#5b01a5", "#a82296", "#da5b69", "#fa9e3b")


####### Load metadata ######

metadata <- readxl::read_xlsx(file.path(metadataPath, "Metadata_Hatchlings_Reloc3.xlsx"))

# MethylKit only allows numeric treatments -> make new variables trt_NUM and Mat_ID
# trt_NUM - where shallow = 1 and deep = 2
# Mat_ID - where SLL is removed
metadata <- mutate(metadata,
                   trt_NUM = case_when(
                     Depth == "Shallow" ~ 1, 
                     Depth == "Deep" ~ 2 
                   )) %>%
  mutate(metadata,
         Mat_ID = as.numeric(substr(metadata$Mother, 4, 6)))

# Check
head(metadata)


####### Load datasets ######

### Load in uniteCov objects ###

# For DMS
uniteCovDM_no0q <- readRDS(file.path(UniteCov_DIR, "EvolApps_WholeGenome_uniteCovDM_ALL_PQLSeq_10pc_q0.05_no0q.RDS")) 

# For global methylation
uniteCovALL <- readRDS(file.path(UniteCov_DIR, "EvolApps_WholeGenome_uniteCovALL_filt_CtoT_Intersect.RDS")) 

### Effect of subsetting out "outliers" 175-2 and 176-5 ###
# subset <- metadata[! metadata$Sample == '175-2',]
# subset <- subset[! subset$Sample == '176-5',]
# subset.ids <- as.vector(subset$Sample)
# subset.treatments <- as.vector(subset$trt_NUM)
# uniteCovALL_NoOutliers <- methylKit::reorganize(uniteCovALL, sample.ids=subset.ids , treatment = subset.treatments) # Create

########################################################################################################################################

###############################
#######    Run NMDS     #######
###############################

OUTDIR <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/02_Output/Figures/Revised_Figures/NMDS"

###### NMDS functions ######

### Make percentage methylation matrix ###
makePercentMetMat <- function(dataset){
  # creates a matrix containing percent methylation values
  perc.meth=percMethylation(dataset)
  # KOSTAS MBE: "Methylated sites and regions with low variation
  # and a standard deviation below 0.3, that is, noninformative
  # sites across individuals, were excluded from the cluster analyses"
  SD=apply(perc.meth,1, sd, na.rm = TRUE)
  if (length(which(SD<0.3)) >0 ) {
    perc.meth <- perc.meth[-which(SD<0.3),]  
  }
  x=t(perc.meth)
  return(x)
}

### Make NMDS function ###
myNMDS <- function(obj.name, mds.to.plot){
  # Set object name and contents to different variables
  # Needed to be able to save figures with name of object used to create it 
  myname <- obj.name
  myobj <- get(myname)
  
  print( paste0( "###### Running NMDS for ", myname , " ######" ) )
  
  ### (1) Run NMDS ###
  # Create percentage methylation matrix
  perc_methy <- makePercentMetMat(myobj)
  
  # Run NMDS based on bray-curtis distances
  # metaMDS finds the most stable NMDS solution by randomly starting from different points in your data
  NMDS <- metaMDS(comm = perc_methy, distance = "bray", maxit=1000, k = 6)
  # Chris: depending on what your methylation value looks like, you may have to homogenize your data 
  # to get a better stress stats, and this would not come with your distance, but with transforming your data
  # when data are bound between 0-1 or 0-100, it is common to either transform everything as log(X+1) or fourth-root(X+1)
  
  # View stressplot (Shepard plot)
  stressplot(NMDS)
  
  ### (2) Plot NMDS ###
  print( paste0( "###### Plotting NMDS with MDS", mds.to.plot , " ######" ) )
  
  # Set plotting variables
  MDS1 = NMDS$points[,1] ; MDS2 = NMDS$points[,2] ; MDS3 = NMDS$points[,3]
  
  NMDS_dt = data.table::data.table(MDS1 = MDS1, MDS2 = MDS2, MDS3 = MDS3,
                                   ID = metadata$Sample,
                                   MAT=as.factor(metadata$Mother), 
                                   DEPTH=as.factor(metadata$Depth))
  
  mycols = c(viridis_custom_palette[1], viridis_custom_palette[5]); myshape = c(21,22) # Colours
  myvar <- "DEPTH" # Treatment variable
  
  # Set polygons, with MDS to plot based on argument mds.to.plot of function 
  if (mds.to.plot == "1v2"){
    dima=1; dimb=2
  } else if (mds.to.plot == "1v3"){
    dima=1; dimb=3
  } else if (mds.to.plot == "2v3"){
    dima=2; dimb=3
  }
  
  hulls <- NMDS_dt[, .SD[chull(get(paste0("MDS", dima)), get(paste0("MDS", dimb)))], by = get(myvar) ]
  
  # Plot 
  myNMDSplot <- ggplot(NMDS_dt, 
                       aes_string(x=paste0("MDS",dima), y=paste0("MDS",dimb))) +
    geom_polygon(data = hulls, aes_string(fill=myvar), alpha=0.4) +
    scale_color_manual(name = "Depth",
                       values = mycols,
                       labels = c("Deep", "Shallow"))+
    scale_fill_manual(name = "Depth",
                      values = mycols,
                      labels = c("Deep", "Shallow"))+
    scale_shape_manual(name = "Depth",
                       values = myshape,
                       labels = c("Deep", "Shallow")) +
    geom_point(aes_string(fill=myvar, shape=myvar), size = 3, alpha = .6) +
    geom_text_repel(aes(label=gsub("SLL", "", NMDS_dt$ID)),
                    size = 3,
                    max.overlaps = Inf,
                    box.padding = 0.4)+
    theme_bw() +
    theme(legend.title=element_blank(), legend.position = "top",
          panel.grid = element_blank(),
          text = element_text(size = 18))
  
  ### (3) Save plots ###
  print( paste0( "###### Saving plot in ", OUTDIR , " ######" ) )

  if (mds.to.plot == "1v2"){
    ggsave(file=file.path(OUTDIR , paste0("NMDS_M1v2_", myname, ".pdf")), plot=myNMDSplot,
           width=9.12, height=7.7, units="in" )
  } else if (mds.to.plot == "1v3"){
    ggsave(file=file.path(OUTDIR , paste0("NMDS_M1v3_", myname, ".pdf")), plot=myNMDSplot,
           width=9.12, height=7.7, units="in" )
  } else if (mds.to.plot == "2v3"){
    ggsave(file=file.path(OUTDIR , paste0("NMDS_M2v3_", myname, ".pdf")), plot=myNMDSplot,
           width=9.12, height=7.7, units="in" )
  }

  # Print needed so it still prints plot when in a loop
  return(print(myNMDSplot))
}

######  Run NMDS function ######

# Argument 1 (my.obj): NAME of uniteCov object to run NMDS analysis on
# Argument 2 (mds.to.plot): which MDS to plot: "1v2" = MDS1 vs MDS2, "1v3" = MDS1 vs MDS3, "2v3" = MDS2 vs MDS3

### Run loop creating all plots for all uniteCovDM objects in environment ###

# NB. Argument 1 in myNMDS function (obj.name) no longer needs to be in quotation marks, as it is 
# already a string after using ls() to find all objects in environment

list_uniteCov_objects <- ls( pattern = "^uniteCovDM")

for (uniteCov_obj_name in list_uniteCov_objects) {
  
  print( paste0( "###### Running NMDS for ", uniteCov_obj_name , " ######" ) )
  
  # Apply myNMDS function for all 3 NMDS plot types
  myNMDS(obj.name = uniteCov_obj_name , mds.to.plot = "1v2")
  myNMDS(obj.name = uniteCov_obj_name , mds.to.plot = "1v3")
  myNMDS(obj.name = uniteCov_obj_name , mds.to.plot = "2v3")
}

########################################################################################################################################

###############################
#######   Dendrogram    #######
###############################

OUTDIR <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/02_Output/Figures/Revised_Figures/Dendrogram"
setwd(OUTDIR)

# To add color bars:
# Generate color palette
depth <- factor(metadata$Depth)
mat <- factor(metadata$Mother)
col_depth <- c(viridis_custom_palette[1], viridis_custom_palette[5])[1:length(unique(metadata$Depth))][depth]
pal <- turbo(length(unique(metadata$Mother)))
col_mat <- pal[mat]

# Run clustering

# For global methylation
mydendro <- clusterSamples(uniteCovALL, dist = "correlation", method = "ward", plot = T)
dend = as.dendrogram(mydendro)

# For DMS
mydendro <- clusterSamples(uniteCovDM_no0q, dist = "correlation", method = "ward", plot = T)
dend = as.dendrogram(mydendro)

# Plot dendrogram
# For DMS (y shift needs to be lower)
dend %>% plot(main=paste("Nest Depth and Maternal ID\nMethylation Clustering\n", 
                         "Distance method: correlation; Clustering method: ward.D"),
              ylab = "Height")

colored_bars(cbind(col_depth, col_mat), dend, y_shift = -0.95,
             rowLabels = c("Depth", "Maternal ID"))


########################################################################################################################################

###############################
###   Statistical Analyses   ##
###############################

###### Create distance matrix ######

makeDatadistFUN <- function(dataset){
  x=makePercentMetMat(dataset)
  # creates a distance matrix. Method: Bray-Curtis, package vegan
  data.dist = as.matrix((vegdist(x, "bray", upper = FALSE))) 
}

# NB. Not possible with uniteCovHALF due to missing values
data.dist = makeDatadistFUN(uniteCovDM_no0q) # DMS
data.dist = makeDatadistFUN(uniteCovALL) # Global

###### Adonis test ######

# Include Maternal ID and relocation in the model
adonis2(data.dist ~ Depth + Mother + Depth:Mother,
        data = metadata)

########################################################################################################################################

##############################
###   Volcano plot of DMS   ##
##############################

OUTDIR <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/02_Output/Figures/Revised_Figures/Volcano"
DMSobject <- readRDS(file.path(DMS_DIR, "EvolApps_WholeGenome_ALL_DMS_PQLseq_10pc_q0.05_no0q.RDS"))

# Vertical line to add to plot
Line_0 <- geom_vline(xintercept=c(0), colour="black")

DMSobject$log10q <- -log10(DMSobject$qvalue)

Volcano_Plot <- ggplot(data=DMSobject, aes(x=meth.diff, y=log10q)) + geom_point(alpha=0.5) +
  scale_x_continuous(limits = c(-35, 35), breaks = seq(-30, 30, by = 10)) +
  scale_y_continuous(limits = c(0, 11), breaks = seq(0,10, by = 5)) +
  labs(x="Methylation difference (%)", y="-log10(q)") +
  Line_0 +
  theme_classic() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

Volcano_Plot

ggsave(file=file.path(OUTDIR ,"Volcano_DMS_no0q.pdf"), plot=Volcano_Plot , 
       width=8, height=9, units="in")


########################################################################################################################################

####################
###   DMS stats   ##
####################

DMSobject$meth.diff.pos <- abs(DMSobject$meth.diff)

mean(DMSobject$meth.diff.pos)
sd(DMSobject$meth.diff.pos)
summary(DMSobject$meth.diff)

nrow(DMSobject[DMSobject$meth.diff > 0 ,]) # 109 hyper
nrow(DMSobject[DMSobject$meth.diff < 0 ,]) # 178 hypo

mean(DMSobject$qvalue)
sd(DMSobject$qvalue)






