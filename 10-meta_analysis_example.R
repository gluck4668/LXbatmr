
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXbatmr")

library(LXbatmr)

rm(list=ls())
gc()


or_files = "D:/Desktop/LXbatmr_shiny (study ver)/data/meta_data/GCST90271622_standard_outcome_mr_all.xlsx"

# devtools::load_all()

#---------------------------------------------------------
meta_analysis (or_files)
