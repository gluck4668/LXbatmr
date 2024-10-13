
library(LXbatmr)


rm(list=ls())

devtools::load_all()

or_file ="D:/Desktop/LXbatmr 2024-10-09(v3.37)/finngen_R10/finngen_R10_NAFLD_outcome_mr_all_significant.xlsx"

# mr_method=c("Inverse variance weighted","MR Egger")     #选择展示的方法

mr_method=c("Inverse variance weighted")

#---------------------------------
forest_plot(or_file,mr_method)
