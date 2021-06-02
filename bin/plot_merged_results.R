#!/usr/bin/env Rscript

####################
# Import libraries #
####################

suppressPackageStartupMessages({
library(optparse)
library(data.table)
library(tidyverse)
library(PheWAS)
    })

##########################################################
# fix phenotype code function                                    
##########################################################
fix_phenotype = function(x){
    x = as.character(x)
    first_part = strsplit(x, '\\.')[[1]][1] %>% str_pad(3, 'left', '0')
    second_part = strsplit(x, '\\.')[[1]][2]
    if (length(strsplit(x, '\\.')[[1]]) == 1){
        x = first_part
    }
    if (length(strsplit(x, '\\.')[[1]]) > 1){
        x = paste0(first_part,'.', second_part)
    }
    return(x)
}
##########################################################
# Read all files, merge them and saves them back                                     
##########################################################

list_results = list.files('./', pattern="phewas_result.csv", full.=TRUE)

results_d = sapply(list_results, function(x) data.table::fread(x), simplify=FALSE) %>% 
    rbindlist(fill=TRUE) %>%
    select_if(~!all(is.na(.))) %>%
    drop_na() %>%
    select(description, group, snp, everything()) %>% as.data.frame()


results_d$phenotype = sapply(results_d$phenotype, function(x) fix_phenotype(x))


results_d %>%
    write.csv("merged_results.csv", row.names=FALSE)
results_d %>% 
    head(1000) %>%
    write.csv("merged_top_results.csv", row.names=FALSE)

##########################################################
# Plot and saves                                
##########################################################

png("phewas_man.png", width = 1000, height = 750, units = 'px', pointsize=16)
phewasManhattan(results_d %>% select(-description, -group), annotate.angle=0)
dev.off()