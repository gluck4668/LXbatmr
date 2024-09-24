
library(LXbatmr)


rm(list=ls())

devtools::load_all()

or_file ="D:/Desktop/LXbatmr 2024-9-23(v3.19)/OR/GCST90271622.gz___outcome_mr_all_significant (1e5).xlsx"

mr_method=c("Inverse variance weighted","MR Egger")     #选择展示的方法

# devtools::load_all()

#---------------------------------
fores_plot(or_file,mr_method)
