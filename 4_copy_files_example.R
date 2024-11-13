
library(LXbatmr)

#devtools::load_all()

rm(list=ls())


file_list="F:/LXbatmr_shiny (study ver)/drug_MR_results_2024-10-23/significant_drugMR_ivw.xlsx"

from_dir="1400 metabolites_gwas_data (cleared result)"

to_dir="meta_finngen_R10_NAFLD"


#------------

copy_file(file_list,from_dir,to_dir)
