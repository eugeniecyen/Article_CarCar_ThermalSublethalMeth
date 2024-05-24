### GO annotation of DMS ###

########################################################################################################################################

# Created by: Charley, adapted from Alice
# Date: 13th Oct 2023
# For reloc 3 for Evol Apps paper

# (1) Get gene names from CarCar functional annotation 
# (2) Features annotation of location in gene (e.g. promoter, exon, intron, intergenic)
# Needs functional annotation of genome: GFF3 and bed12 file

########################################################################################################################################

###############################
###### Prep environment #######
###############################

###### Load packages ######

list.of.packages <- c("tidyverse", "methylKit", "genomation", "rtracklayer", "devtools",
                      "forcats", "stringr", "purrr", "gtable") # graphical libraries

# Create function to install packages if they don't exist
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)){
    install.packages(new.pkg, dependencies = TRUE,repos = "http://cran.us.r-project.org")}
  sapply(pkg, require, character.only = TRUE)
}

ipak(list.of.packages)

library(BiocManager)

list.bioc <- c("AnnotationDbi",
               "Category", # for hypergeometric GO test
               "GenomicFeatures",## for annotation
               "GenomicRanges",
               "GOstats", # for GO analysis
               "GSEABase")  # for GO term GeneSetCollection           
ipak2 <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)){
    BiocManager::install(new.pkg)}
  sapply(pkg, require, character.only = TRUE)
}

ipak2(list.bioc)

## For GO enrichment (needs "GOstats" & "GSEABase" installed)
# install_github("asishallab/goEnrichment")
library(goEnrichment)

####### Set directories ######

# In local

DMS_DIR <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/03_DMS_Objects"
uniteCov_DIR <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/01_UniteCov_Objects"
AnnoPath <- "/data/SBCS-EizaguirreLab/Turtle_Genome/10_Annotation/06_Functional_Annotation/Scaff_With_Chr0/Final_Annotations/HMMER3_Databases/No_Pathways/longest_isoforms_only"

OUTDIR_Figures <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/02_Output/Figures/Revised_Figures/GO_Enrichment"
OUTDIR_Tables <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/02_Output/Tables/Revised_Tables/GO_Enrichment"

####### Load GFF3 file ######

annotGff3 <- rtracklayer::readGFF(file.path(AnnoPath, "CarCar_QM_v1_2021_12_Scaff_With_Chr0_Functional_Annotation_LongestIsoform.gff3"))

########################################################################################################################################

##########################
### Prep gene universe ###
##########################

# This is all genes which are present in the full dataset i.e. uniteCovALL
uniteCovALL <- readRDS(file.path(uniteCov_DIR, "EvolApps_WholeGenome_uniteCovALL_filt_CtoT_Intersect.RDS"))

####### Create gene universe #######

gene_universe <- data.frame(
  subsetByOverlaps(GRanges(annotGff3), GRanges(uniteCovALL))) %>% # Subset CpGs in uniteCovALL only
  filter(lengths(Ontology_term)!=0) %>% # rm non existing GO terms
  filter(type %in% "gene")  %>% # keep all the genes with GO terms
  dplyr::select(c("ID", "Ontology_term")) %>%
  mutate(go_linkage_type = "IEA") %>% #NB: IEA but not necessarily true, but needed for further algo
  relocate("Ontology_term","go_linkage_type","ID") %>%
  tidyr::unnest(Ontology_term) %>% # 1 GO per line (was a list before in this column)
  data.frame()

gene_universe$ID %>% unique %>% length # 11,119 genes with GO annotations in uniteCovALL

####### Create gene set collection #######

goFrame <- AnnotationDbi::GOFrame(gene_universe, organism="Caretta caretta")
goAllFrame <- GOAllFrame(goFrame)
gsc_universe <- GeneSetCollection(goAllFrame, setType = GOCollection())

########################################################################################################################################

##############################
### Prep gene sub-universe ###
##############################

# For gene sub universe: all genes on DMS that have GO annotations attached

####### Load in file #######
DMSobject <- readRDS(file.path(DMS_DIR, "EvolApps_WholeGenome_ALL_DMS_PQLseq_10pc_q0.05_no0q.RDS"))

# Split into hypo and hyper meth
hyperDMS <- DMSobject[DMSobject$meth.diff > 0 ,] # 109
hypoDMS <- DMSobject[DMSobject$meth.diff < 0 ,] # 178

####### Extract positions of DMS to GRanges object #######

# Create vector of DMS with just contig name and position in contig
myPos_hyper = paste(hyperDMS$chr, hyperDMS$end)
myPos_hypo = paste(hypoDMS$chr, hypoDMS$end)

# Change the Pos vector into a dataframe): called df
df_hyper <- data.frame(chr=sapply(strsplit(myPos_hyper, " "), `[`, 1),
                 start=sapply(strsplit(myPos_hyper, " "), `[`, 2),
                 end=sapply(strsplit(myPos_hyper, " "), `[`, 2))

df_hypo <- data.frame(chr=sapply(strsplit(myPos_hypo, " "), `[`, 1),
                 start=sapply(strsplit(myPos_hypo, " "), `[`, 2),
                 end=sapply(strsplit(myPos_hypo, " "), `[`, 2))

# Create GRanges object from DMS df
GRanges_DMSobject_hyper <- makeGRangesFromDataFrame(df_hyper)
GRanges_DMSobject_hypo <- makeGRangesFromDataFrame(df_hypo)


####### Create subuniverse on genes with DMS #######

# Hyper
subuniverse_hyper <- data.frame(
  subsetByOverlaps(GRanges(annotGff3), GRanges(GRanges_DMSobject_hyper))) %>% # Subset CpGs in DMS only
  filter(lengths(Ontology_term)!=0) %>% # rm non existing GO terms
  filter(type %in% "gene")  %>% # keep all the genes with GO terms
  dplyr::select(c("ID", "Ontology_term")) %>%
  mutate(go_linkage_type = "IEA") %>% #NB: IEA but not necessarily true, but needed for further algo
  relocate("Ontology_term","go_linkage_type","ID") %>%
  tidyr::unnest(Ontology_term) %>% # 1 GO per line (was a list before in this column)
  data.frame()

subuniverse_hyper$ID %>% unique %>% length # 40 genes with GO annotations in hyper DMS

# Hypo
subuniverse_hypo <- data.frame(
  subsetByOverlaps(GRanges(annotGff3), GRanges(GRanges_DMSobject_hypo))) %>% # Subset CpGs in DMS only
  filter(lengths(Ontology_term)!=0) %>% # rm non existing GO terms
  filter(type %in% "gene")  %>% # keep all the genes with GO terms
  dplyr::select(c("ID", "Ontology_term")) %>%
  mutate(go_linkage_type = "IEA") %>% #NB: IEA but not necessarily true, but needed for further algo
  relocate("Ontology_term","go_linkage_type","ID") %>%
  tidyr::unnest(Ontology_term) %>% # 1 GO per line (was a list before in this column)
  data.frame()

subuniverse_hypo$ID %>% unique %>% length # 54 genes with GO annotations in hypo DMS

########################################################################################################################################

###########################################
### Run conditional hypergeometric test ###
###########################################

## **IMPORTANT NOTE from Mel: why conditional hypergeometric test?**
# The GO ontology is set up as a directed acyclic graph, where a parent term 
# is comprised of all its child terms. If you do a standard hypergeometric, you might e.g., find 'positive regulation of kinase activity' to 
# be significant. If you then test 'positive regulation of catalytic activity', which is a parent term, then it might be significant as well, 
# but only because of the terms coming from positive regulation of kinase activity. The conditional hypergeometric takes this into account, and 
# only uses those terms that were not already significant when testing a higher order (parent) term.

####### Set up Alice's functions ####### 

runTestHypGeom <- function(subuniverse, onto){
  ## Constructing a GOHyperGParams objects or KEGGHyperGParams objects from a GeneSetCollection
  ## Then run hypergeometric test:
  GO_NO_fdr <- hyperGTest(GSEAGOHyperGParams(name="GO_set",
                                             geneSetCollection = gsc_universe,
                                             geneIds = as.vector(unique(subuniverse$ID)), # gene ids for the selected gene set
                                             universeGeneIds = unique(gene_universe$ID),
                                             ontology = onto, # ("BP", "CC", or "MF")
                                             pvalueCutoff = 0.05,
                                             conditional = TRUE, # see note above
                                             testDirection = "over")) # over represented GO terms
  ## Use GOEnrich as a wrapper around GOStat for extra FDR comparison
  ## Does not solve all issues, but better than nothing. See: https://support.bioconductor.org/p/5571/
  GO_fdr <- joinGOEnrichResults(goEnrichTest(gsc=gsc_universe,
                                             gene.ids = as.vector(unique(subuniverse$ID)),# genes in selected gene set
                                             univ.gene.ids = unique(gene_universe$ID),
                                             ontologies = onto, # A string for GO to use ("BP", "CC", or "MF")
                                             pvalue.cutoff = 0.05,
                                             cond = TRUE, # see note above
                                             test.dir = "over"),# over represented GO terms
                                p.adjust.method = "fdr")        
  return(list(GO_NO_fdr=GO_NO_fdr, GO_fdr=GO_fdr))
}  

makeGO <- function(subuniverse){
  GO_MF <- runTestHypGeom(subuniverse = subuniverse, onto = "MF")
  GO_CC <- runTestHypGeom(subuniverse = subuniverse, onto = "CC")
  GO_BP <- runTestHypGeom(subuniverse = subuniverse, onto = "BP")
  ## Get percentage of genes over reppresented in universe
  dfMFperc = GO_MF$GO_NO_fdr %>% summary() %>% mutate(genePercent = Count/Size*100) %>% 
    dplyr::select(c("Term", "genePercent")) %>% dplyr::rename(GO.name=Term)
  dfCCperc = GO_CC$GO_NO_fdr %>% summary() %>% mutate(genePercent = Count/Size*100) %>% 
    dplyr::select(c("Term", "genePercent")) %>% dplyr::rename(GO.name=Term)
  dfBPperc = GO_BP$GO_NO_fdr %>% summary() %>% mutate(genePercent = Count/Size*100) %>% 
    dplyr::select(c("Term", "genePercent")) %>% dplyr::rename(GO.name=Term)
  ## Add this information to FDR corrected table
  GO_MF_all = merge(GO_MF$GO_fdr, dfMFperc)
  GO_CC_all = merge(GO_CC$GO_fdr, dfCCperc)
  GO_BP_all = merge(GO_BP$GO_fdr, dfBPperc)
  ## Merge the df MP and BP
  dfGO = rbind(GO_MF_all, GO_CC_all, GO_BP_all)
  dfGO = dfGO %>% mutate(Term = factor(x = GO.term, levels = GO.term))
  ## Relabel GO group names
  dfGO$GO.category[dfGO$GO.category %in% "CC"]="Cellular components"
  dfGO$GO.category[dfGO$GO.category %in% "BP"]="Biological processes"
  dfGO$GO.category[dfGO$GO.category %in% "MF"]="Molecular functions"
  return(dfGO)
}

###### Split DMS by hyper and hypo ######

# # Run makeGO function on subuniverses
dfGO_DMS_hyper = makeGO(subuniverse_hyper) # 9 hits at p=0.05
dfGO_DMS_hyper$Type = "Hyper"

dfGO_DMS_hypo = makeGO(subuniverse_hypo) # 42 hits at p=0.05
dfGO_DMS_hypo$Type = "Hypo"

dfGO = merge(dfGO_DMS_hyper, dfGO_DMS_hypo, all = TRUE)

# Save
write.csv(dfGO, file.path(OUTDIR_Tables, "GO_Enrichment_EvolApps_DMS_ALL_PQLSeq_p0.05_SplitHyperHypo.csv"), row.names = F)

########################################################################################################################################

#################
### Nice plot ###
#################

Width <- 12
Height <- 12

### Split by hyper and hypo ###
GOplot <- dfGO %>%
  ggplot(aes(x = Type, y = factor(GO.name))) +
  geom_point(aes(color = p.value.adjusted, size = genePercent)) +
  scale_color_gradient(name = "adjusted\np-value", low = "#fa9e3b", high = "#0d0887") +
  scale_size_continuous(name = "% of genes") +
  theme_bw() + ylab("") + xlab("") +
  theme(legend.box.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
        legend.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
        legend.key = element_rect(fill = "#ebebeb", color = "#ebebeb"), 
        legend.position = "left") +  # grey box for legend
  facet_grid(fct_inorder(GO.category) ~ ., scales = "free", space = "free") +
  scale_y_discrete(limits = rev)  # reverse axis to have alphabetical order

print(GOplot)

#Save plot 
ggsave(filename = file.path(OUTDIR_Figures, "GO_Enrichment_ALL_DMS_PQLSeq_SplitHyperHypo_p0.05.pdf"), 
       plot = GOplot,
       height = Height,
       width = Width)

### All DMS together ###
GOplot <- dfGO_DMS %>%
  ggplot(aes(x = Type, y = factor(GO.name))) +
  geom_point(aes(color = p.value.adjusted, size = genePercent)) +
  scale_color_gradient(name = "adjusted\np-value", low = "#fa9e3b", high = "#0d0887") +
  scale_size_continuous(name = "% of genes") +
  theme_bw() + ylab("") + xlab("") +
  theme(legend.box.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
        legend.background = element_rect(fill = "#ebebeb", color = "#ebebeb"),
        legend.key = element_rect(fill = "#ebebeb", color = "#ebebeb"), 
        legend.position = "left") +  # grey box for legend
  facet_grid(fct_inorder(GO.category) ~ ., scales = "free", space = "free") +
  scale_y_discrete(limits = rev)  # reverse axis to have alphabetical order

print(GOplot)

#Save plot 
ggsave(filename = file.path(OUTDIR_Figures, "GO_Enrichment_ALL_DMS_PQLSeq_p0.05.pdf"), 
       plot = GOplot,
       height = Height,
       width = Width)
