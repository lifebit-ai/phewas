#!/usr/bin/env Rscript

library(PheWAS)

########################################
### Data input
########################################
genotypes=read.table("r_genotypes.raw",header=TRUE)[,c(-2:-6)]
names(genotypes)[1]="id"

# TODO: add option to import something like this:
    # example phenotype.csv:
    # id,T2D,max.a1c
    # 1,T,10
    # 2,F,NA
    # 3,F,6
icd9.codes <- read.csv("id.icd9.count.csv",colClasses=c("integer","character","integer"))


########################################
### Data transformation
########################################
# Create the PheWAS code table- translates the icd9s, adds exclusions, and reshapes to a wide format
phenotypes=createPhewasTable(id.icd9.count)


########################################
### Run the PheWAS
########################################
# TODO: add multithreaded approach, set cores based on cmd line input
results=phewas(phenotypes,genotypes,cores=1,significance.threshold=c("bonferroni"))

#Add PheWAS descriptions
results_d=addPhecodeInfo(results)


########################################
### Plotting
########################################
png("phewas_man.png", width = 1000, height = 750, units = 'px', pointsize=16)
phewasManhattan(results, annotate.angle=0)
dev.off()











