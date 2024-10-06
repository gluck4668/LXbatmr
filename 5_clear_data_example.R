
library(LXbatmr)

rm(list=ls())
gc()

# devtools::load_all()

data_dir <- "D:/Desktop/LXbatmr 2024-9-21(v3.17)/data/exp_gwas"

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

samplesize_info= "D:/Desktop/LXbatmr 2024-9-21(v3.17)/exp_gwas_gwas_id_info.xlsx" # "The 233 circulating metabolic biomarkers buildGRCh.xlsx_gwas_info.xlsx" # 先运行：1_gwas_info_example （获取样本数量）

#--------------------------------
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
               samplesize_info=samplesize_info,
               clum_p=clum_p,clump_kb=clump_kb,clump_r2=clump_r2)






