
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXbatmr")

library(LXbatmr)


#---Drug target MR -----------------------------------------------------------

rm(list=ls()) #每次都要运行，以清除内存中旧的数据

{

target_gene_data = "D:/Desktop/LXbatmr_shiny 2024-10-12/data/didymin target MR/didymin_target_gene.xlsx"

#---exposure data---
exp_data="D:/Desktop/LXbatmr_shiny 2024-10-12/data/didymin target MR/metabolites" # 暴露，必填

chr_exp="chromosome"
pos_exp="base_pair_location"
snp_exp = "SNP"
beta_exp = "beta" #如果有beta,用beta;如果没有，则用OR或者logOR
se_exp = "standard_error"
effect_allele_exp = "effect_allele"
other_allele_exp = "other_allele"
eaf_exp = "effect_allele_frequency"
pval_exp = "p_value"


#---outcome data---
out_data=c("D:/Desktop/LXbatmr_shiny 2024-10-12/data/didymin target MR/NAFLD_standard/finngen_R10_NAFLD_standard.gz",
          "D:/Desktop/LXbatmr_shiny 2024-10-12/data/didymin target MR/NAFLD_standard/GCST90271622_standard.gz") #* 结局，必填

# outcome_name="NAFLD"  # 结局，一般是疾病名称，选填

chr_out="chromosome"
pos_out="base_pair_location"
snp_out = "variant_id"
beta_out = "beta" #如果有beta,用beta;如果没有，则用OR或者logOR
se_out = "standard_error"
effect_allele_out = "effect_allele"
other_allele_out = "other_allele"
eaf_out = "effect_allele_frequency"
pval_out = "p_value"

#---MR parameters---
clum_p=5e-8        #* 关联性，必填
clump_kb = 100  #* 连锁不平衡kb, 必填
clump_r2 = 0.3  #* 连锁不平衡r2,必填

}
#---run the main function---

devtools::load_all()

#-----------------------------------------------

drugMR (target_gene_data=target_gene_data,
                   exp_data=exp_data,
                   chr_exp=chr_exp,
                   pos_exp=pos_exp,
                   snp_exp=snp_exp,
                   beta_exp=beta_exp,
                   se_exp=se_exp,
                   effect_allele_exp=effect_allele_exp,
                   other_allele_exp=other_allele_exp,
                   eaf_exp=eaf_exp,
                   pval_exp=pval_exp,

                   out_data=out_data,
                   chr_out=chr_out,
                   pos_out=pos_out,
                   snp_out=snp_out,
                   beta_out=beta_out,
                   se_out=se_out,
                   effect_allele_out=effect_allele_out,
                   other_allele_out=other_allele_out,
                   eaf_out=eaf_out,
                   pval_out=pval_out,

                   clum_p=clum_p,
                   clump_kb=clump_kb,
                   clump_r2=clump_r2 )


#---The end--------------------------------------------------------------------




