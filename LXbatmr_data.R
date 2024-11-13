
library(openxlsx)

meta_gwas_id <- read.xlsx("1400 metabolites_gwas_data _gwas_info.xlsx")

usethis::use_data(meta_gwas_id,overwrite = T)


data(meta_gwas_id)

