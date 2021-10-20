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
  make_option(c("--cov_file"), action="store", type='character',
              help="Covariate file for the corresponding genotype/phenotype data."),            
  make_option(c("--n_cpus"), action="store", default='1', type='character',
              help="Number of cpus used."),
  make_option(c("--pheno_codes"), action="store", default='None', type='character',
              help="String representing phenotype nomenclature (ie. DOID, ICD9, ICD10, HPO)."),
    make_option(c("--min_code_count"), action="store", default='2', type='character',
              help="Minimum code count to define case/control."),
      make_option(c("--add_exclusions"), action="store", default=TRUE, type='logical',
              help="Applying pheWAS exclusions to phecodes."),
  make_option(c("--outprefix"), action="store", default='covid', type='character',
              help="String containing prefix to be added to files.")
)

args = parse_args(OptionParser(option_list=option_list))

pheno_file          = args$pheno_file # file
geno_file           = args$geno_file
covariate_file      = args$cov_file
n_cpus              = as.numeric(args$n_cpus) # int
min_code_count      = as.numeric(args$min_code_count)
add_exclusions      = args$add_exclusions
pheno_codes         = args$pheno_codes
outprefix           = args$outprefix

########################################
### Data input
########################################
# load genotypes
genotypes=read.table(geno_file,header=TRUE)
genotypes$FID = genotypes$IID
genotypes = genotypes[,c(-2:-6)]
names(genotypes)[1]="id"


 if (!is.null(covariate_file)) {
   covariates = read.table(covariate_file,header=TRUE,sep=',')
 }
<<<<<<< HEAD

=======
>>>>>>> d0016ff... {EG} run with covariate file (#52)


# load the pheno data
id.code.count = read.csv(pheno_file,colClasses=c("character",'character',"integer"))
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
  id.code.count$doid = with(mappings, icd9[match(id.code.count$code, doid)])
  # rename DOID col
  names(id.code.count)[names(id.code.count) == 'code'] = 'icd9'

  phenotypes=createPhewasTable(id.code.count, min.code.count = min_code_count, add.exclusions = add_exclusions)
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
if (length(dim(genotypes)) < 1){
results_d = data.frame()
write.csv(results_d, file=paste0(outprefix,"_phewas_results.csv"), row.names=FALSE)
}

if (length(dim(genotypes)) > 1){
  if (is.null(covariate_file)) {
  results=phewas(phenotypes=phenotypes,genotypes=genotypes,cores=as.numeric(n_cpus),significance.threshold=c("bonferroni"))
}
else {
  results=phewas(phenotypes=phenotypes,genotypes=genotypes,cores=as.numeric(n_cpus),covariates=covariates,significance.threshold=c("bonferroni"))
}


# Add PheWAS descriptions
results_d=addPhecodeInfo(results)

# Write results to csv ordered by significance
all_res=results_d[order(results_d$p),]
write.csv(all_res, file=paste0(outprefix,"_phewas_results.csv"), row.names=FALSE)


}
