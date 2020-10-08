#!/usr/bin/env Rscript
install.packages("devtools", repos='http://cran.us.r-project.org')
devtools::install_github("PheWAS/PheWAS", ref='v1.0-dev')
library(PheWAS)
