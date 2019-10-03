#!/usr/bin/env Rscript

###################################
# Perform PheWAS analysis & generate Manhattan plot
args = commandArgs(trailingOnly=TRUE)
# Requires:
pheno_file <- args[1] # file
n_cpus <- args[2] # int
pheno_codes <- args[3] # string eg `doid` to convert DOID codes -> icd9
# also requires `mappings.csv` file in pwd
###################################

library(PheWAS)

########################################
### Data input
########################################
# load genotypes
genotypes=read.table("r_genotypes.raw",header=TRUE)[,c(-2:-6)]
names(genotypes)[1]="id"

# load the pheno data
id.icd9.count <- read.csv(pheno_file,colClasses=c("integer","character","integer"))
# TODO: add option to import something like this:
  # example phenotype.csv:
  # id,T2D,max.a1c
  # 1,T,10
  # 2,F,NA
  # 3,F,6

# if using DOID codes convert to icd9
if (pheno_codes == "doid") {
  mappings <- read.csv("mappings.csv")
  # remove prefix
  mappings$curie_id <- gsub("^DOID:", "", mappings$curie_id)
  mappings$mapped_curie <- gsub("^ICD9CM:", "", mappings$mapped_curie)
  # remove uncessary cols
  keeps <- c("curie_id", "mapped_curie")
  mappings <- mappings[keeps]
  # rename cols
  names(mappings) <- c("doid","icd9")
  # replace DOID codes with icd9
  id.icd9.count$doid <- with(mappings, icd9[match(id.icd9.count$doid, doid)])
  # rename DOID col
  names(id.icd9.count)[names(id.icd9.count) == 'doid'] <- 'icd9'
  # remove NA values
  #id.icd9.count <- id.icd9.count[complete.cases(id.icd9.count), ]
  # id.icd9.count$doid <-  mappings[match(id.icd9.count$icd9, mappings$icd9), 1, drop=F]
}

########################################
### Data transformation
########################################
# Create the PheWAS code table- translates the icd9s, adds exclusions, and reshapes to a wide format
phenotypes=createPhewasTable(id.icd9.count)

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
write.csv(all_res, file="top_results.csv", row.names=FALSE)

########################################
### Plotting
########################################
png("phewas_man.png", width = 1000, height = 750, units = 'px', pointsize=16)
phewasManhattan(results, annotate.angle=0)
dev.off()

