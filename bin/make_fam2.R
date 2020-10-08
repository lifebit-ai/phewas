#!/usr/bin/env Rscript

####################
# Import libraries #
####################

suppressPackageStartupMessages({
library(optparse)
library(data.table)
library(tidyverse)

    })

##########################################################
# Parse arguments                                        
##########################################################

option_list = list(
  make_option(c("--pheno_file"), action="store", default='None', type='character',
              help="String containing long format table for phenotypes."),
)

args = parse_args(OptionParser(option_list=option_list))

pheno_file          = args$pheno_file # file


########################################
### Data input
########################################
df = fread(pheno_file)

# remove vcf column
df = df %>% select(-vcf) %>% rename(IID=sampleid)
# make Family ID column
df$FID = df$IID

# make paternal & maternal ID columns
df$PAT = 0
df$MAT = 0

# reorder columns
df = df %>% select(FID,IID,PAT,MAT,SEX, everything())

write.table('sample.phe', sep=' ', quote=F, row.names=F)
