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
# Parse arguments                                        
##########################################################

option_list = list(
  make_option(c("--pheno_file"), action="store", default='None', type='character',
              help="String containing long format table for phenotypes."),
  make_option(c("--geno_file"), action="store", default='None', type='character',
              help="String containing genotypes to be tested."),
  make_option(c("--n_cpus"), action="store", default='1', type='character',
              help="String containing input metadata for columns in Cohort Browser output."),
  make_option(c("--pheno_codes"), action="store", default='None', type='character',
              help="String representing phenotype nomenclature (ie. DOID, ICD9, ICD10, HPO).")
)

args = parse_args(OptionParser(option_list=option_list))

pheno_file          = args$pheno_file # file
geno_file           = args$geno_file
n_cpus              = integer(args$n_cpus) # int
pheno_codes         = args$pheno_codes

########################################
### Data input
########################################
# load genotypes
genotypes=read.table(geno_file,header=TRUE)
genotypes$FID = genotypes$IID
genotypes = genotypes[,c(-2:-6)]
names(genotypes)[1]="id"

# load the pheno data
id.code.count = read.csv(pheno_file,colClasses=c("character","character",'character',"integer"))
# TODO: add option to import something like this:
  # example phenotype.csv:
  # id,T2D,max.a1c
  # 1,T,10
  # 2,F,NA
  # 3,F,6

# if using DOID codes convert to icd9
if (pheno_codes == "doid") {
  mappings = read.csv("assets/mappings.csv")
  # remove prefix
  mappings$curie_id = gsub("^DOID:", "", mappings$curie_id)
  mappings$mapped_curie = gsub("^ICD9CM:", "", mappings$mapped_curie)
  # remove uncessary cols
  keeps = c("curie_id", "mapped_curie")
  mappings = mappings[keeps]
  # rename cols
  names(mappings) = c("doid","icd9")
  # replace DOID codes with icd9
  id.code.count$doid = with(mappings, icd9[match(id.code.count$doid, doid)])
  # rename DOID col
  names(id.code.count)[names(id.code.count) == 'doid'] = 'icd9'
  # remove NA values
  #id.code.count = id.code.count[complete.cases(id.code.count), ]
  # id.code.count$doid =  mappings[match(id.code.count$icd9, mappings$icd9), 1, drop=F]
  phenotypes=createPhewasTable(id.code.count)
}
if (pheno_codes == 'icd10'){
  icd10_data = read.csv("assets/Phecode_map_v1_2_icd10cm_beta.csv", colClasses=rep("character",8)) %>% select(icd10cm, phecode)
  id.code.count = id.code.count %>% 
    left_join(icd10_data, by=c("code" = "icd10cm")) %>%
    select(id, phecode, count)

  phenotypes=createPhewasTable(id.code.count)
}
if (pheno_codes == 'hpo'){
  hpo_data = fread("assets/HPO_to_phecode.txt") %>% select(term_id, phecode)
  id.code.count$code = str_replace_all(id.code.count$code, "HPO:", "")
  id.code.count = id.code.count %>% 
    left_join(hpo_data, by=c("code" = "phecode")) %>%
    select(id, phecode, count)

  phenotypes=createPhewasTable(id.code.count)
}

########################################
### Run the PheWAS
########################################
results=phewas(phenotypes=phenotypes,genotypes=genotypes,cores=as.numeric(n_cpus),significance.threshold=c("bonferroni"))

# Add PheWAS descriptions
results_d=addPhecodeInfo(results)

# Write results to csv ordered by significance
all_res=results_d[order(results_d$p),]
write.csv(all_res, file="phewas_results.csv", row.names=FALSE)
top_res=results_d[order(results_d$p)[1:1000],]
write.csv(top_res, file="top_results.csv", row.names=FALSE)

########################################
### Plotting
########################################
png("phewas_man.png", width = 1000, height = 750, units = 'px', pointsize=16)
phewasManhattan(results, annotate.angle=0)
dev.off()

