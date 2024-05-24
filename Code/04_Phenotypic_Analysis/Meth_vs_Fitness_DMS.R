####### Evol Apps Paper 2023 (Yen et al.) #######
####### Analyse Relationship Between Hatchling #######
####### Fitness Data and DNA Methylation #######

## Explore relationship between hatchling fitness and DNA methylation at key sites
## Methylation sites are identified earlier in the paper/analysis

##### This is v2 of this script, written to analyse new data generated for re-submission of paper #####

# Author: Christophe Eizaguirre, Nov 2023
# Adapted by: James Gilbert, 16 Nov 2023 / May 2024

#################################################################################
####################### Linking DMS with fitness ################################
#################################################################################

#################################
###### Prepare Environment ######
#################################

### Load packages
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(MASS)
library(plyr)
library(dplyr)
library(gclus)
library(corrplot)    
library(tibble)
library(lmerTest)
library( MuMIn )
library( grid )
library(gridExtra)
library(nlme)
library(here)

### Set Directories
DIR <- here() ####### R project file must be in a project dir that contains the following subdirs:
DATADIR <- file.path(DIR, "Data/Fitness_Analysis")
DATAOUT <- file.path(DATADIR, "Revised_Output" )
PLOTOUTDIR <- file.path(DIR, "Figures"); dir.create(PLOTOUTDIR, showWarnings = F)

mycols = c("#0D0887", "#FA9E3B")

### Read and process data
# read in hatchling data
df_hatchling <- read.csv(file.path(DATADIR, "Reloc3_Hatchling_Fitness_NERC_Exp_2021_WithPatrolDates_WithIDs_T.csv") ) %>%
  mutate( Maternal_ID = as.factor(Maternal_ID ) ) %>%
  filter( WGBS == "yes")

# Read in methylation data for candidate DMS
df_DMS <- read.csv(file.path(DATADIR, "EvolApps_WholeGenome_PercMethMatrix_ALL_CandidateDMS_Revised.csv" ),
                   header = T,
                   check.names = FALSE)

# needs processing in order to be merged with df_hatchling
# need Hatchling ID as a column and the rest of the columns are candidate DMS
rownames(df_DMS) <- df_DMS$DMS.name

df_DMS <- df_DMS %>%
  dplyr::select( -DMS.name) %>%
  t() %>% # so hatchlings are rows and DMS are columns
  as.data.frame() %>%
  rownames_to_column(., "Hatchling_ID") # remake variable with hatchling ID
  
df_hatchling_DMS <- merge(x = df_hatchling, all.x = T, 
      y = df_DMS,
      by = "Hatchling_ID")

#################################
###### DMS vs Fitness Data ######
#################################

#### fitness trait to be assessed is specified in the function
# note for Run_average and Flip_average, we use log( + 1) transformation

# store model outputs
List_fitness_DMS_models = list()
List_fitness_DMS_plots = list()
List_mod_stats = list()

# makes plots/stats for all DMS ( all colnames in df_hatchling_DMS that are in the DMS df)
# i.e. all apart from column 1; hatchling ID
COL_IDX <- which(colnames(df_hatchling_DMS[ , ]) %in% colnames(df_DMS))
COL_IDX <- COL_IDX[ ! COL_IDX %in% 1 ]

# Loop through the specified range of columns
for ( col_index in COL_IDX ) { 
  
  # Extract the variable from the data frame
  variable <- df_hatchling_DMS[, col_index]
  variable_name <- names(df_hatchling_DMS)[ col_index] # for naming output
  
  # loop through fitness traits
  for ( fitness_trait in c( "SCL_mm", "Weight_g", "Run_average", "Flip_average" ) ) {
    
    print( paste0(" ####### making model and plot for ", variable_name, " and ", fitness_trait, " #######") )
    
    # make fitness trait label nicer for plot
    if (fitness_trait ==  "Flip_average") {
      fitness_trait_lab = gsub( "Flip_average", "Self-righting time (log(seconds + 1)", fitness_trait )
    }
    if (fitness_trait ==  "Run_average") {
      fitness_trait_lab = gsub( "Run_average", "Crawl time (log(seconds + 1)", fitness_trait )
    }
    if (fitness_trait ==  "SCL_mm") {
      fitness_trait_lab = gsub( "SCL_mm", "Straight carapace length (mm)", fitness_trait )
    }
    if (fitness_trait ==  "Weight_g") {
      fitness_trait_lab = gsub( "Weight_g", "Mass (g)", fitness_trait )
    }
    
    #### Model calculations
    # we want to know whether fitness is related to methylation
    # we  first check if the DMS (variable) is related to fitness while controlling for treatment (additive model)
    # then to check whether any relationship is the same among both deep and shallow
    # we make models where there is an interaction between DMS (variable) and treatment
    # we compare this against a null model containing just maternal ID
    # constructed separately for run and flip so we can add log + 1 transformation
    if (fitness_trait %in% c("Run_average", "Flip_average") == T ) {

      model_int <- lmer(data = df_hatchling_DMS,
                    log( df_hatchling_DMS[ , fitness_trait] + 1 ) ~  variable*Treatment + (1|(Maternal_ID)))
      
      model_add <- lmer(data = df_hatchling_DMS,
                        log( df_hatchling_DMS[ , fitness_trait] + 1 ) ~  variable+Treatment + (1|(Maternal_ID)))
      

      null_model <- lmer(data = df_hatchling_DMS,
                         log( df_hatchling_DMS[ , fitness_trait] + 1 ) ~  (1|(Maternal_ID)))


      } else if (fitness_trait %in% c( "SCL_mm", "Weight_g" ) == T ) {

      model_int <- lmer(data = df_hatchling_DMS,
                    df_hatchling_DMS[ , fitness_trait] ~  variable * Treatment + (1|(Maternal_ID)))
      
      model_add <- lmer(data = df_hatchling_DMS,
                        df_hatchling_DMS[ , fitness_trait] ~  variable + Treatment + (1|(Maternal_ID)))
       
      null_model <- lmer(data = df_hatchling_DMS,
                         df_hatchling_DMS[ , fitness_trait] ~  (1|(Maternal_ID)))

      }
    
    #######################################
    ##### STATS FOR interaction model #####
    #######################################
    # extract p value of the interaction term from the model
    interaction_p = anova(model_int) %>% as.data.frame() %>% select(`Pr(>F)`) %>% filter(rownames(.) == 'variable:Treatment') %>% round(4)
    # get delta AIC vs null model
    delta_AIC_int = round( AIC(model_int) - AIC(null_model), 2)
    
    # prints marginal (fixed effects) and conditional (fixed and random effects) R2
    model_r2_int = r.squaredGLMM(model_int)
    
    R2_total_int = model_r2_int[1, 2] * 100
    R2_fixed_int = model_r2_int[1, 1] * 100
    R2_random_int = ( as.numeric(R2_total_int) - as.numeric(R2_fixed_int) )
    
    List_fitness_DMS_models[[ paste0( variable_name, "_", fitness_trait, "_interaction"  )]] = model_int
    
    List_mod_stats[[ paste0( variable_name, "_", fitness_trait, "_interaction"  )]] = data.frame( DMS = variable_name,
                                                                                  fitness_trait = fitness_trait,
                                                                                  mod_type = "interaction",
                                                                                  AIC_v_mod_matID_only = delta_AIC_int,
                                                                                  DMS_pvalue = interaction_p,
                                                                                  R2_total = R2_total_int,
                                                                                  R2_fixed = R2_fixed_int,
                                                                                  R2_random = R2_random_int)
    
    
    # Create a scatter plot with a smoothed line using ggplot2 for INTERACTION
    plot <- ggplot(df_hatchling_DMS, aes(x = variable, y = df_hatchling_DMS[ , fitness_trait ],
                             color=Treatment
    )) +
      geom_point() +
      geom_smooth(method = "lm")+
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(x=paste0( "DNA Methylation at site ", variable_name, " (%)"),
           y= fitness_trait_lab ) +
      theme(text = element_text(size = 15))+
      scale_color_manual(values=mycols)

    # Print the plot
    # print(plot)

    #### Save list of plots to file
    ggsave(
      filename = file.path(PLOTOUTDIR, "PQLseq_DMS", paste0( "p_", interaction_p, "_", variable_name, "_", fitness_trait, "_interaction.pdf")),
      plot = plot,
      height = 10,
      width = 10,
      units = "in"
    )

    List_fitness_DMS_plots[[ paste0( variable_name, "_", fitness_trait, "_interaction"  )]] = plot
    
    ####################################
    ##### STATS FOR additive model #####
    ####################################
    # extract p value of the DMS when included as a fixed effect in the model
    additive_p = anova(model_add) %>% as.data.frame() %>% select(`Pr(>F)`) %>% filter(rownames(.) == 'variable') %>% round(4)
    
    delta_AIC_add = round( AIC(model_add) - AIC(null_model), 2)
    
    # prints marginal (fixed effects) and conditional (fixed and random effects) R2
    model_r2_add = r.squaredGLMM(model_add)
    
    R2_total_add = model_r2_add[1, 2] * 100
    R2_fixed_add = model_r2_add[1, 1] * 100
    R2_random_add = ( as.numeric(R2_total_add) - as.numeric(R2_fixed_add) )
    
    List_fitness_DMS_models[[ paste0( variable_name, "_", fitness_trait, "_additive"  ) ]] = model_add
    
    List_mod_stats[[ paste0( variable_name, "_", fitness_trait, "_additive"  ) ]] = data.frame( DMS = variable_name,
                                                                                                fitness_trait = fitness_trait,
                                                                                                mod_type = "additive",
                                                                                                AIC_v_mod_matID_only = delta_AIC_add,
                                                                                                DMS_pvalue = additive_p,
                                                                                                R2_total = R2_total_add,
                                                                                                R2_fixed = R2_fixed_add,
                                                                                                R2_random = R2_random_add)
    
    ### stat_smooth plots the interaction but for the additive model we colour points by group
    ### but include a single line representing overall slope
    gradient = model_add@beta[ 2 ] # betas are coefficients and the 2nd is DMS in our model
    # note the points and line are specified differently because we want points by Treatment but line is combined
    plot_add <- ggplot(df_hatchling_DMS) +
      geom_point(aes(x = variable, y = df_hatchling_DMS[ , fitness_trait ],
                     colour = Treatment)) +
      geom_smooth(aes(x = variable, y = df_hatchling_DMS[ , fitness_trait ]),
                  method = "lm", colour = "black") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(x=paste0( "DNA Methylation at site ", variable_name, " (%)"),
           y= fitness_trait_lab ) +
      theme(text = element_text(size = 15))+
      scale_color_manual(values=mycols)

    # Print the plot
    # print(plot)

    #### Save list of plots to file
    ggsave(
      filename = file.path(PLOTOUTDIR, "PQLseq_DMS", paste0( "p_", additive_p, "_", variable_name, "_", fitness_trait, "_additive.pdf")),
      plot = plot_add,
      height = 10,
      width = 10,
      units = "in"
    )

    List_fitness_DMS_plots[[ paste0( variable_name, "_", fitness_trait, "_additive"  )]] = plot_add
    
  }
}

##### FDR calculation
# Number of tests
k = length(List_mod_stats)
# Significance level (alpha)
alpha = 0.05
# Calculate the sum of reciprocals from 1 to k
harmonic_sum <- sum(1 / seq(1, k))
# Calculate the adjusted alpha value using the B-Y method
alpha_adjusted = alpha / harmonic_sum
print( alpha_adjusted )

##### compile model stats
df_model_stats = do.call(rbind, List_mod_stats) %>%
  rename( pvalue = `Pr..F.`) %>%
  mutate( adjusted_alpha = alpha_adjusted )

head(df_model_stats)

write.csv(df_model_stats, file.path(DATAOUT, "df_model_stats_all.csv"))

###### Saving model output
### we will report full model output for the models that we focus on
### i.e. those where DMS term p < adjusted_alpha
### Also show the impact of SCL on model output
IDX_p_lt_adj_alpha = which(df_model_stats$pvalue < alpha_adjusted )

df_model_stats[ IDX_p_lt_adj_alpha, ]

for ( i in IDX_p_lt_adj_alpha) {
  
  # Extract the variable from the data frame
  DMS <- df_model_stats$DMS[ i ]
  fitness_trait <- df_model_stats$fitness_trait[ i ]
  mod_type = df_model_stats$mod_type[ i ]
  
  # re create model
  model = lmer(data = df_hatchling_DMS,
               log( df_hatchling_DMS[ , fitness_trait] + 1 ) ~  df_hatchling_DMS[ , DMS]*Treatment + (1|(Maternal_ID)))
  
  # print model summary to file
  sink( file.path(DATAOUT, paste0(DMS, "_", fitness_trait, "_model_summary.csv") ) )
  print( summary(model) )
  sink(file= NULL)
  
  # show impact from SCL
  if (fitness_trait %in% c("Run_average", "Flip_average") == T ) {
    
    # create new model with SCL included as a fixed effect
    new_model <- lmer(data = df_hatchling_DMS,
                  log( df_hatchling_DMS[ , fitness_trait] + 1 ) ~  df_hatchling_DMS[ , DMS]*Treatment + SCL_mm + (1|(Maternal_ID)))
    
    # print new model summary to file
    sink( file.path(DATAOUT, paste0(DMS, "_", fitness_trait, "_withSCL_model_summary.csv") ) )
    print( summary(new_model) )
    sink(file= NULL)
    print(summary(new_model))
    
    # compare performance of models with and without SCL
    print( anova(new_model, model) )
    
  }
  
}

### end of script


### Post hoc tests for significant interactions ###

# RALYL
model1 <- List_fitness_DMS_models$RALYL_SCL_mm_interaction

# Visualise fitted lines
emmip(model1, variable ~ Treatment, cov.reduce = range)

# Summary of contrasts at max and min of range, by treatment
# cov.reduce required for continuous variable
emm <- emmeans(model1, pairwise ~ variable * Treatment, cov.reduce = range)
pairs(emm, simple = "Treatment")


# TMEM273
model1 <- List_fitness_DMS_models$TMEM273_Run_average_interaction

# Visualise fitted lines
emmip(model1, variable ~ Treatment, cov.reduce = range)

# Summary of contrasts at max and min of range, by treatment
# cov.reduce required for continuous variable
emm <- emmeans(model1, pairwise ~ variable * Treatment, cov.reduce = range)
pairs(emm, simple = "Treatment")






