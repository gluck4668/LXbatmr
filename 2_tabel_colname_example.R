
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXbatmr")

library(LXbatmr)

#-------------------------

rm(list=ls())
gc()

# devtools::load_all()

exp_data_dir <- "D:/Desktop/LXbatmr 2024-9-19(v3.16)/exp_gwas"

table_colname(exp_data_dir)

