
library(LXbatmr)

devtools::load_all()

rm(list=ls())


file_list="finngen_R10_NAFLD_mr_significant.xlsx___gwas_info.xlsx"

old_dir="meta_list"

new_dir="meta_new"


#------------

copy_file(file_list,old_dir,new_dir)
