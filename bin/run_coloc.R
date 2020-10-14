#!/usr/bin/env Rscript

####################
# Import libraries #
####################

suppressPackageStartupMessages({
library(optparse)
library(data.table)
library(tidyverse)
library(coloc)
library(ComplexHeatmap)
    })

####################
# Functions        #
####################

run_coloc_analysis = function(phewas_df, gwas_df, phewas_phenotype, trait_type) {
    #get the group in pheWAS you are comparing with the GWAS
    phewas_df = phewas_df %>% filter(description == phewas_phenotype)
    #Run coloc analysis
    if (trait_type == 'binary'){
        coloc_results = coloc.abf(dataset1=list(pvalues=phewas_df$p,N=phewas_df$n_total, s=(phewas_df$n_cases/phewas_df$n_total), type="cc"), #PheWAS
                        dataset2=list(pvalues=gwas_df$p.value, N=(gwas_df$N.Cases + gwas_df$N.Controls), s=(gwas_df$N.Cases/gwas_df$N.Controls), type="cc"), #GWAS
                        MAF=gwas_df$AF_Allele2)
    }
    if (trait_type == 'quantitative'){
        coloc_results = coloc.abf(dataset1=list(pvalues=phewas_df$p,N=phewas_df$n_total,type="cc"), #PheWAS
                        dataset2=list(pvalues=gwas_df$p.value,  N=gwas_df$N, type="quant"), #GWAS
                        MAF=gwas_df$AF_Allele2)
    }
    return(coloc_results$summary)
}

##########################################################
# Parse arguments                                        
##########################################################

option_list = list(
  make_option(c("--phewas_summary"), action="store", default='None', type='character',
              help="String containing long format table for phenotypes."),
  make_option(c("--gwas_summary"), action="store", default='None', type='character',
              help="String containing genotypes to be tested."),
  make_option(c("--outprefix"), action="store", default='covid', type='character',
              help="String containing prefix to be added to files."),
  make_option(c("--gwas_trait_type"), action = "store", default='binary', type='character',
              help="String containing type of trait from GWAS catalogue. Choice between 'quantitative' or 'binary'")
)

args = parse_args(OptionParser(option_list=option_list))

phewas_summary          = args$phewas_summary # file
gwas_summary            = args$gwas_summary
gwas_trait_type         = args$gwas_trait_type
outprefix               = args$outprefix

########################################
### Data input
########################################
# load phewas results
phewas_df = fread(phewas_summary) %>% as_tibble()
phewas_df$snp = sapply(strsplit(phewas_df$snp, '_'), '[', 1)

# load gwas summary table
gwas_df = fread(gwas_summary) %>% as_tibble()

#Find common SNPs

gwas_df = gwas_df %>% filter(SNPID %in% mask_snp) %>% arrange(SNPID)
phewas_df = phewas_df %>% filter(snp %in% mask_snp) %>% arrange(snp)

########################################
### Run coloc
########################################

coloc_results = sapply(unique(phewas_df$description), function(x) run_coloc_analysis(phewas_df, gwas_df, x, gwas_trait_type))

phenotypes = names(coloc_results)

coloc_results = rbindlist(coloc_results) %>% as.data.frame()
rownames(coloc_results) = phenotypes 

########################################
### Save results
########################################
coloc_results %>% write.csv(paste0(outprefix,"_coloc_results.csv"), row.names=FALSE)
########################################
### Prepare to plot
########################################

grouping = phewas_df %>% select(description, group) %>% distinct() %>% arrange(description) %>% select(group)
grouping = grouping$group

########################################
### Plot heatmap with results
########################################

group_col_fun = structure(brewer.pal(length(grouping), "Set3"), 
    names = grouping)
ha = rowAnnotation(
    group = anno_simple(grouping), col = group_col_fun),
    annotation_name_side = "left")
ht = Heatmap(coloc_results, name = "PP", right_annotation = ha, cluster_columns=FALSE, cluster_rows=TRUE, clustering_method_rows='ward.D2')

png(paste0(outprefix, "_coloc_heatmap.png"), width = 1000, height = 750, units = 'px', pointsize=16)
draw(ht)
dev.off()