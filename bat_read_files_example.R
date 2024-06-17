
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXbatmr")

library(LXbatmr)

#-------------------------

rm(list=ls())
gc()

devtools::load_all()

exp_data_dir <- "meta_list"

bat_read_files(exp_data_dir)

