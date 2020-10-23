
library(tidyverse)
combinations=expand.grid(design_mode = c("case_vs_control_contrast"), continuous_var_aggregation = c("mean","max", "min", "median"), continuous_var_transformation = c("log", "log10", "log2", "zscore", "None"), pheno_codes = c("doid", "hpo", "icd10", "icd9"), post_analysis=c("None", "coloc"), gwas_trait_type=c("binary", "quantitative", "None"))
combinations %>% filter(!((post_analysis != 'None') & (gwas_trait_type == 'None'))) %>%
                 filter(!((post_analysis=='None') & (gwas_trait_type %in% c("binary","quantitative"))))
