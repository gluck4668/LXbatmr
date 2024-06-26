
library(LXbatmr)

#devtools::load_all()

rm(list=ls())


file_list="GCST90275041.tsv___outcome_mr_ivw_significant.xlsx"

from_dir="1400 metabolites_gwas_data (cleared result)"

to_dir="meta_GCST90275041"


#------------

copy_file(file_list,from_dir,to_dir)
