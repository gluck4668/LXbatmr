
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXbatmr")

library(LXbatmr)


rm(list=ls())
gc()

 devtools::load_all()

#------------Batch processing MR analysis---------------------

#----exposure data----------------
{
exp_data_dir="D:/Desktop/LXbatmr 2024-10-12/data/exp_excel"  # 数据所在的文件夹

#---暴露数据主要参数，缺一不可------
snp_exp = "SNP"
beta_exp = "beta.exposure"
se_exp = "se.exposure"
effect_allele_exp = "effect_allele.exposure"
other_allele_exp = "other_allele.exposure"
eaf_exp = "eaf.exposure"
pval_exp = "pval.exposure"


#----暴露数据筛选和去除连锁不平衡的条件------
clum_p=1e-5
clump_kb=10000
clump_r2=0.001


#-----outcome data------------------
out_data_file=c("D:/Desktop/LXbatmr 2024-10-12/data/outcome/GCST90271622_standard.gz",
                "D:/Desktop/LXbatmr 2024-10-12/data/outcome/finngen_R10_NAFLD_standard.gz")

#---结局数据主要参数，缺一不可------
snp_out = "variant_id"
beta_out = "beta"
se_out = "standard_error"
effect_allele_out = "effect_allele"
other_allele_out = "other_allele"
eaf_out = "effect_allele_frequency"
pval_out = "p_value"
}

#----------------------------

bat_tsmr(exp_data_dir,snp_exp,beta_exp,se_exp,effect_allele_exp,
        other_allele_exp,eaf_exp,pval_exp,clum_p,clump_kb,clump_r2,
        out_data_file,snp_out,beta_out,se_out,effect_allele_out,
        other_allele_out,eaf_out,pval_out)

#----------------------------
library(data.table)

head(fread(dir(exp_data_dir,full.names = T)[1]))  # 核对暴露数据“列名称”
head(fread(out_data_file[1]))         #核对结局数据“列名称”



