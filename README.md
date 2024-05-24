# Article_CarCar_ThermalSublethalMeth

This repository contains all scripts used for the research article: 

**DNA methylation carries signatures of sublethal effects under thermal stress in loggerhead sea turtles**

Authors: Eugenie C. Yen, James D. Gilbert, Alice Balard, Inês O. Afonso, Kirsten Fairweather, Débora Newlands, Artur Lopes, Sandra M. Correia, Albert Taxonera, Stephen J. Rossiter, José M. Martín-Durán, Christophe Eizaguirre

Preprint: https://doi.org/10.1101/2023.11.22.568239

Raw sequencing reads are available on ENA under study accession: PRJEB75968


## Repository structure
```
.
├── 01_WGBS_Bismark_Pipeline
│   ├── CarCar_00_Bismark_Genome_Preparation.sh
│   ├── CarCar_01_CutAdapt_Trimming_Array.sh
│   ├── CarCar_02_AlignDedupSort_WholeGenome_Array.sh
│   ├── CarCar_03_MethCall_Destrand_Array.sh
│   └── CarCar_04_Merge_Destranded_Chrom_MethCalls.sh
├── 02_Methylation_Object_Prep
│   ├── Bash_Submission_Scripts
│   │   ├── S01.1_PrepObjects_WholeGenome_ALL_Reloc3.sh
│   │   └── S01.5_Run_PQLseq_ChromArray.sh
│   ├── R01.1_PrepObjects_WholeGenome_ALL_Reloc3.R
│   ├── R01.2_Remove_CtoT_Mutations.R
│   ├── R01.3_Prep_PQLseq_Input_From_UniteCovObject.R
│   ├── R01.4_Prep_PQLseq_Input_Theoretical_Relatedness_Matrix.R
│   ├── R01.5_Run_PQLseq_ChromArray.R
│   ├── R01.6_MergeChroms_RunSLIM_SubsetDMS.R
│   └── R01.7_CreateMethMatrix.R
├── 03_Downstream_Methylation_Analysis
│   ├── Clustering_Analyses_Volcano.R
│   ├── Functional_Annotation_DMS_v3.R
│   ├── Genomic_Region_Annotation_Donut_Plot.R
│   ├── Global_Methylation_Stats.R
│   ├── GO_Enrichment_Analysis.R
│   └── ManhattanPlot.R
└── README.md
```



