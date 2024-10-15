
library(LXbatmr)

#devtools::load_all()

rm(list=ls())


file_list="D:/Desktop/LXbatmr_shiny 2024-10-14/drug_MR_results_2024-10-13  233/drugMR_ivw_all.xlsx"

from_dir="1400 metabolites_gwas_data (cleared result)"

to_dir="meta_finngen_R10_NAFLD"


#------------

copy_file(file_list,from_dir,to_dir)
