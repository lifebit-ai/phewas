#!/usr/bin/env Rscript
install.packages("devtools", repos='http://cran.us.r-project.org')
devtools::install_github("PheWAS/PheWAS", ref='v1.0-dev')
library(PheWAS)
devtools::install_github("chr1swallace/coloc")
library(coloc)
install.packages('R.utils', repos='https://cloud.r-project.org/')
devtools::install_github("jokergoo/ComplexHeatmap")
install.packages('RColorBrewer', repos='https://cloud.r-project.org/')