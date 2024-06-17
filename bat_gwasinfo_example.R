install.packages("devtools")
library(devtools)

install_github("gluck4668/LXbatmr")

library(LXbatmr)


rm(list=ls())
gc()

devtools::load_all()


rm(list = ls())

data_list = NA #"bat_mr_result_significant.xlsx"

data_dir= "NAFLD"  # NA


#----------------------------------------

bat_get_gwasinfo(data_dir,data_list)

