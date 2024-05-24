### Create a Manhattan plot highlighting differential methylation and genes of interest ###

########################################################################################################################################

# 8th Nov 2023
# Created by James, edited by Charley
# For Evol Apps

# Manhattan plot of DMS, with 29 DMS of interest labelled

########################################################################################################################################

###############################
###### Prep environment #######
###############################

### Load packages
library( ggplot2 )
library (purrr)
library(genomation)
library( dplyr )
library( ggrepel )
library(methylKit)

### Set directories
AssemDIR <- "/data/SBCS-EizaguirreLab/Turtle_Genome/00_Final_Assemblies"
DIR <- "/data/SBCS-EizaguirreLab/Turtle_WGBS"
DMDIR <- file.path(DIR, "05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/02_MethDiff_Objects")
metadataPath <- file.path(DIR, "05_MethylKit/04_Evol_Apps_Project/01_Working/01_Metadata")
PlotDIR <- file.path(DIR, "05_MethylKit/04_Evol_Apps_Project/02_Output/Figures/Revised_Figures/Manhattan")

##### Load data

# genes of interest to label in Manhattan plot
df_Candidate_DMS <- read.csv(file.path(metadataPath, "EvolApps_Candidate_DMS_Revised.csv"))

# File containing chromosome locations and lengths
df_chrms <- read.table(file.path(AssemDIR, "CarCar_QM_v1_2021_12_Scaff_With_Chr0.fasta.fai") ) %>% dplyr::rename(
  chr = V1,
  Length = V2,
  Start = V3
) %>%
  mutate(chrom_nr = as.numeric(gsub("SLK063_ragtag_chr", "", chr)),
         # chrom_order=factor(chrom_nr) %>% as.numeric()
  )%>% 
  arrange(chrom_nr) %>%
  mutate(gstart=lag(Length,default=0) %>% cumsum(),
         gend=gstart+Length, 
         type=LETTERS[2-(chrom_nr%%2)],
         gmid=(gstart + ((gend - gstart) / 2)  ) )

### Differential methylation data
DiffMeth <- readRDS(file.path(DMDIR, "EvolApps_WholeGenome_ALL_MethDiff_filt_CtoT_Intersect_PQLseq.RDS"))

## DM object needs data on chromosome locations so we can get absolute position of sites, not sites within chromosome
df_DiffMeth = merge( x = getData(DiffMeth), by.x = c("chr"),
                     df_chrms[ , c("chr", "gstart")], by.y = c("chr"))

df_DiffMeth <- df_DiffMeth %>%
  mutate( gpos = start + gstart ) # genome position = base position of site (= start) + base position of start of chromosome ( = gstart)

## Add gene names (only 29) to the diff meth data
df_DiffMeth <- merge(x = df_DiffMeth, by.x = c( "chr", "start" ), all.x = T,
                     y = df_Candidate_DMS[ ,  c( "DMS.name", "chrom", "DMS.pos" ) ], by.y = c("chrom", "DMS.pos"))

####### Create Manhattan Plot

# Alternating colours for chromosomes
mycols=c("grey70","white")

# For smaller chromosomes, x labels overlap. Make custom list of labels to add to x axis
x=12:22
my_labels = c(0:10, x[!x%%2] ) # list is 0-10 and every other number in x
# also want the midpoints of the chromosomes we label to specifyx axis tick location
my_breaks = df_chrms %>% dplyr::filter(chrom_nr %in% my_labels) %>% dplyr::select(gmid) %>% as.vector() # must match label IDs

## We will plot "significant" and "non-significant" points separately
## significant = > 10pc difference and q < 0.01, and remove q=0 
Sig_IDX <- which(abs(df_DiffMeth$meth.diff) > 10 & df_DiffMeth$qvalue < 0.05 & df_DiffMeth$qvalue != 0)

# For plotting from all non-significant points (NB. messy with PQLSeq output for -log10q plot)
NonSig_IDX_All <- which(abs(df_DiffMeth$meth.diff) < 10 | df_DiffMeth$qvalue > 0.05)
length(Sig_IDX) + length(NonSig_IDX_All) # Check it adds up correctly
# Other option? Filter for sites where meth diff is < x% to get rid of some points, but is this misleading?
# NonSig_IDX_All <- which(abs(df_DiffMeth$meth.diff) < 10 & abs(df_DiffMeth$meth.diff) > 1 | df_DiffMeth$qvalue > 0.05)

### For methylation difference of DMS only ###
# Add lines at 10%md and 20%md to show cut offs
Lines_Dash <- geom_hline(yintercept=c(10, -10, 20, -20), colour="grey30", linetype="dotted", linewidth=0.8)

Manhattan_DiffMeth <- ggplot() +
  geom_rect( data = df_chrms, aes(xmin=gstart,xmax=gend,ymin=-Inf,ymax=Inf,fill=type), alpha=.2,
             show.legend = FALSE) +
  geom_point( data =  df_DiffMeth[ Sig_IDX, ], aes( x = gpos, y = meth.diff,
                                                    colour = qvalue)) +
  scale_color_gradient(name = "q-value", low = "#fa9e3b", high = "#0d0887") +
  scale_fill_manual(name = "type", values = mycols) +
  scale_x_continuous(breaks=my_breaks$gmid,
                     labels=my_labels,
                     expand = c(0.01,0)
  ) +
  Lines_Dash +
  geom_label_repel( data = df_DiffMeth[ !is.na(df_DiffMeth$DMS.name), ],
                    aes( x = gpos, y = meth.diff, label = DMS.name),
                    box.padding   = 0.8
  ) +
  theme_bw() +
  theme( panel.grid = element_blank(),
         axis.line = element_blank(),
         axis.title = element_text(size = 14),
         axis.text = element_text(size = 10)
  ) +
  xlab( "Chromosome") + 
  ylab ( "Methylation difference (%)" )

Manhattan_DiffMeth

# save to file
ggsave( file.path(PlotDIR, "Meth_Diff_DMS/Manhattan_WholeGenome_ALL_MethDiff_DMS_PQLseq_no0q.pdf"),
        plot = Manhattan_DiffMeth,
        device = "pdf",
        width=14, height=7, units="in" )
