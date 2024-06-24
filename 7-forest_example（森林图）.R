
library(LXbatmr)

devtools::load_all()

rm(list=ls())

or_dir ="mr_or"

mr_method=c("Inverse variance weighted")     #选择展示的方法


#---------------------------------
forest_one_by_one (or_dir,mr_method)

forest_all (or_dir,mr_method)
