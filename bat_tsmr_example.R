
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
exp_data_dir="data"  # 数据所在的文件夹,or a xlsx file including the id list of exposrue data

#---如果是ieu网站数据（id或本地文件均可）不需要填下面信息；其他来源的数据要填
snp_exp = "SNP"
beta_exp = "beta.exposure"
se_exp = "se.exposure"
effect_allele_exp = "effect_allele.exposure"
other_allele_exp = "other_allele.exposure"
eaf_exp = "eaf.exposure"
pval_exp = "pval.exposure"

#samplesize_eqtl=265431 #如果表里有，直接用；如果没有，则需要在下载数据时，查看相应信息

clum_p=1e-5
clump_kb=10000
clump_r2=0.001


#-----outcome data------------------
out_data_file="finngen_R10_NAFLD.gz"

snp_out = "rsids"
beta_out = "beta"
se_out = "sebeta"
effect_allele_out = "alt"
other_allele_out = "ref"
eaf_out = "af_alt"
pval_out = "pval"
}


bat_tsmr(exp_data_dir,snp_exp,beta_exp,se_exp,effect_allele_exp,
        other_allele_exp,eaf_exp,pval_exp,clum_p,clump_kb,clump_r2,
        out_data_file,snp_out,beta_out,se_out,effect_allele_out,
        other_allele_out,eaf_out,pval_out)

