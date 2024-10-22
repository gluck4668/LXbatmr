
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXbatmr")

library(LXbatmr)

rm(list=ls())
gc()

or_files = c("D:/Desktop/LXbatmr_shiny (study ver)/data/OR/GCST90091033_buildGRCh37.tsv.gz___outcome_mr_all_significant (1e5).xlsx",
           "D:/Desktop/LXbatmr_shiny (study ver)/data/OR/GCST90271622.gz___outcome_mr_all_significant (1e5).xlsx")
# devtools::load_all()

#---------------------
venn_analysis(or_files)

