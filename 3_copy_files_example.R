
library(LXbatmr)

devtools::load_all()

rm(list=ls())


file_list="GCST90275041.tsv___outcome_mr_ivw_significant.xlsx"

old_dir="1400 metabolites_gwas_data"

new_dir="meta_GCST90275041"


#------------

copy_file(file_list,old_dir,new_dir)
