install.packages("devtools")
library(devtools)

install_github("gluck4668/LXbatmr")

library(LXbatmr)


rm(list=ls())
gc()

#devtools::load_all()


data_list = NA # "finngen_R10_NAFLD_mr_significant.xlsx" # NA

data_dir= "1400 metabolites_gwas_data (cleared result)" # "info_test" # NA


#----------------------------------------

bat_get_gwasinfo(data_dir,data_list)

