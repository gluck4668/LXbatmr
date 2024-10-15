
library(LXbatmr)

rm(list=ls())
gc()

# devtools::load_all()

data_dir <- "D:/Desktop/LXbatmr 2024-10-08(v3.36)/20241008/exp1008"

# library(data.table)
# head(fread(dir(data_dir,full.names = T)[1]))

#-------exposure data  parameters------
{
snp_exp = "variant_id"
beta_exp = "beta"
se_exp = "standard_error"
effect_allele_exp = "effect_allele"
other_allele_exp = "other_allele"
eaf_exp = "effect_allele_frequency"
pval_exp = "p_value"
#---------------------------------

clum_p=1e-5
clump_kb=10000
clump_r2=0.001

}


#--------------------------------------------------------------
clear_exp_data(data_dir=data_dir,snp_exp=snp_exp,beta_exp=beta_exp,
               se_exp=se_exp,
               effect_allele_exp=effect_allele_exp,
               other_allele_exp=other_allele_exp,
               eaf_exp=eaf_exp,pval_exp=pval_exp,
               clum_p=clum_p,clump_kb=clump_kb,clump_r2=clump_r2)






