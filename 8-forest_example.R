
library(LXbatmr)


rm(list=ls())

devtools::load_all()

or_file ="D:/Desktop/LXbatmr_shiny (study ver)/data/OR/GCST90091033_buildGRCh37.tsv.gz___outcome_mr_all_significant (1e5).xlsx"

# mr_method=c("Inverse variance weighted","MR Egger")     #选择展示的方法

mr_method=c("Inverse variance weighted")

#---------------------------------
forest_plot(or_file,mr_method)
