### R01.4_Prep_PQLseq_Input_Theoretical_Relatedness_Matrix.R ###

########################################################################################################################################

# Created by: Charley
# Date: 13th Mar 2024
# OnDemand script
# For Evol Apps

# Create genetic relatedness matrix assuming siblings represented have relatedness of 0.5 and different nests are 0.
# Same method as in von Holdt et al. 2022 MACAU analysis: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.13878?saml_referrer 

# NB. Genetic relatedness matrix must have samples in same order as count data
# Currently no ID checking feature implemented -> need to make sure this is the case

########################################################################################################################################

###############################
###### Prep environment #######
###############################

####### Load packages ######

library(tidyverse)
library(readxl)

####### Set directories ######

DIR_PQL <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/02_MethylKit_Objects/02_PQLseq_Objects/Input"

########################################################################################################################################

############################################################
###### Create theoretical genetic relatedness matrix #######
############################################################

### Load family info ###

metadata <- read.csv(file.path(DIR_PQL , "metadata_uniteCovALL_Slots.csv"))

# Read in Excel file containing metadata
metadata_full <- readxl::read_xlsx("/data/SBCS-EizaguirreLab/Turtle_WGBS/05_MethylKit/04_Evol_Apps_Project/01_Working/01_Metadata/Metadata_Hatchlings_Reloc3.xlsx")

# MethylKit only allows numeric treatments -> make new variables trt_NUM and Mat_ID
# trt_NUM - where Hatchling = 1 and Mother = 2
# Mat_ID - where SLL is removed
metadata_full <- mutate(metadata_full,
                        trt_NUM = case_when(
                          Depth == "Shallow" ~ 1, 
                          Depth == "Deep" ~ 2 
                        )) %>%
  mutate(metadata_full,
         Mat_ID = as.numeric(substr(metadata_full$Mother, 4, 6)))

# Add to macau metadata file
metadata$mat.id <- metadata_full$Mat_ID


### Create matrix ###

# Create matrix with r=0 for all elements
gen_df <- as.data.frame(matrix(c(0), ncol = 40, nrow=40))
rownames(gen_df) <- metadata$sample.id
colnames(gen_df) <- metadata$sample.id
row_names <- rownames(gen_df)
col_names <- colnames(gen_df)

# If maternal ID in metadata matches (siblings), fill cell with 0.5
unique_mat_ids <- unique(metadata$mat.id) # Get unique values of "mat_id" from the metadata data frame

# Run loop:
# For each maternal ID, pull out all sample IDs with this maternal ID
for (id in unique_mat_ids) {
  sample_ids_to_replace <- metadata[metadata$mat.id == id, ] %>% pull(sample.id)
  
  # For each row and each column, if the row name and column name is in sample_ids_to_replace, fill cell with 0.5
  for (row_name in row_names) {
    for (col_name in col_names) {
      if (row_name %in% sample_ids_to_replace && col_name %in% sample_ids_to_replace) {
        gen_df[row_name, col_name] <- 0.5
      }
    }
  }
}

# If row name and column name are the same (i.e. same sample ID), fill cell with 1
for (row_name in row_names) {
  if (row_name %in% col_names) {
    gen_df[row_name, row_name] <- 1 
  }
}

# View
gen_df

# Rename column and rows to ind_name
rownames(gen_df) <- metadata$ind_num
colnames(gen_df) <- metadata$ind_num
rownames(gen_df) <- paste0("ind", rownames(gen_df))
colnames(gen_df) <- paste0("ind", colnames(gen_df))

# Convert to matrix
gen_matrix <- as.matrix(gen_df)

# Save
write.table(gen_matrix, file.path(DIR_PQL, "Evol_Apps_ALL_Theoretical_Relatedness_Matrix.txt"), row.names=TRUE, col.names=TRUE, sep = " ", quote = FALSE)




