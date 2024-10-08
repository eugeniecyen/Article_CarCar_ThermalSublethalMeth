## Summary

This repository contains scripts and data used for the research article: 

**DNA methylation carries signatures of sublethal effects under thermal stress in loggerhead sea turtles**

Authors: Eugenie C. Yen, James D. Gilbert, Alice Balard, Inês O. Afonso, Kirsten Fairweather, Débora Newlands, Artur Lopes, Sandra M. Correia, Albert Taxonera, Stephen J. Rossiter, José M. Martín-Durán, Christophe Eizaguirre

http://dx.doi.org/10.1111/eva.70013 

Raw sequencing reads are available on ENA under study accession PRJEB75968
<br/><br/>
## Contents

### Code

**01_WGBS_Bismark_Pipeline**: Scripts used to generate CpG methylation calls from raw WGBS reads via the Bismark pipeline

**02_Methylation_Object_Prep**: Scripts used to prepare methylation objects across sites present in all individuals and identifying differentially methylated sites (DMS) via PQLseq

**03_Downstream_Methylation_Analysis**: Scripts for downstream methylation analyses and plots (clustering, functional annotation, global methylation counts and GO enrichment analyses)

**04_Phenotypic_Analysis**: Scripts for statistical tests for associations with phenotypic fitness-related traits at the clutch-level and individual, hatchling-level
<br/><br/>
### Data

**Metadata_Hatchlings_Reloc3.xlsx**: Experiment metadata for the 40 hatchlings that were sequenced with WGBS in this study

**EvolApps_ALL_filt_CtoT_Intersect_GlobalMeth_Per_Ind.csv**: Global methylation counts per hatchling

**Evol_Apps_PQLseq_Kinship_Matrix.txt**: Input relatedness matrix between hatchlings for PQLseq

**EvolApps_WholeGenome_ALL_DMS_PQLseq_10pc_q0.05_no0q.RDS**: methylKit methylDiff object of 287 DMS identified via PQLseq from sites present in all hatchlings (q<0.05, meth. difference > 10%)

**EvolApps_WholeGenome_uniteCovDM_ALL_PQLSeq_10pc_q0.05_no0q.RDS**: methylKit methylBase object of coverage and no. of methylated bases per hatchling at 287 DMS

**EvolApps_WholeGenome_PercMethMatrix_ALL_CandidateDMS_Revised.csv**: Matrix of % methylation per hatchling at 29 DMS of interest (meth. difference >20%)

**Reloc3_NestLevelData_WithTemperatureInfo_NERCExp2021.csv**: Clutch-level data of hatching success and incubation temperature

**Reloc3_Hatchling_Fitness_NERC_Exp_2021_WithPatrolDates_WithIDs_T.csv**: Individual-level phenotypic data for all 408 hatchlings, with clutch-level data attached
<br/><br/>
### File structure

```
.
├── Code
│   ├── 01_WGBS_Bismark_Pipeline
│   │   ├── CarCar_00_Bismark_Genome_Preparation.sh
│   │   ├── CarCar_01_CutAdapt_Trimming_Array.sh
│   │   ├── CarCar_02_AlignDedupSort_WholeGenome_Array.sh
│   │   ├── CarCar_03_MethCall_Destrand_Array.sh
│   │   └── CarCar_04_Merge_Destranded_Chrom_MethCalls.sh
│   ├── 02_Methylation_Object_Prep
│   │   ├── Bash_Submission_Scripts
│   │   │   ├── S01.1_PrepObjects_WholeGenome_ALL_Reloc3.sh
│   │   │   └── S01.5_Run_PQLseq_ChromArray.sh
│   │   ├── R01.1_PrepObjects_WholeGenome_ALL_Reloc3.R
│   │   ├── R01.2_Remove_CtoT_Mutations.R
│   │   ├── R01.3_Prep_PQLseq_Input_From_UniteCovObject.R
│   │   ├── R01.4_Prep_PQLseq_Input_Theoretical_Relatedness_Matrix.R
│   │   ├── R01.5_Run_PQLseq_ChromArray.R
│   │   ├── R01.6_MergeChroms_RunSLIM_SubsetDMS.R
│   │   └── R01.7_CreateMethMatrix.R
│   ├── 03_Downstream_Methylation_Analysis
│   │   ├── Clustering_Analyses_Volcano.R
│   │   ├── Functional_Annotation_DMS_v3.R
│   │   ├── Genomic_Region_Annotation_Donut_Plot.R
│   │   ├── Global_Methylation_Stats.R
│   │   ├── GO_Enrichment_Analysis.R
│   │   └── ManhattanPlot.R
│   └── 04_Phenotypic_Analysis
│       ├── Meth_vs_Fitness_DMS.R
│       └── Phenotypic_Stats.R
├── Data
│   ├── EvolApps_ALL_filt_CtoT_Intersect_GlobalMeth_Per_Ind.csv
│   ├── Evol_Apps_PQLseq_Kinship_Matrix.txt
│   ├── EvolApps_WholeGenome_ALL_DMS_PQLseq_10pc_q0.05_no0q.RDS
│   ├── EvolApps_WholeGenome_PercMethMatrix_ALL_CandidateDMS_Revised.csv
│   ├── EvolApps_WholeGenome_uniteCovDM_ALL_PQLSeq_10pc_q0.05_no0q.RDS
│   ├── Metadata_Hatchlings_Reloc3.xlsx
│   ├── Reloc3_Hatchling_Fitness_NERC_Exp_2021_WithPatrolDates_WithIDs_T.csv
│   └── Reloc3_NestLevelData_WithTemperatureInfo_NERCExp2021.csv
└── README.md
```



