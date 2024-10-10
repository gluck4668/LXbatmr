
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXbatmr")

library(LXbatmr)



#---Drug target MR -----------------------------------------------------------

rm(list=ls()) #*每次都要运行，以清除内存中旧的数据

{
#---target gene---
# tart_gene ="HMGCR"
# chr_exp= 5
# pos_exp_start = 74632993
# pos_exp_end = 74657941

target_gene_data = "D:/Desktop/LXbatmr 2024-10-10(v3.38)/data/drugMR/drug_target.xlsx"


#---exposure data---
exp_data="D:/Desktop/LXbatmr 2024-10-10(v3.38)/data/drugMR/expouse" # 暴露，必填

# exposure_name="sphingolipid" # 暴露，一般是疾病名称，选填

# 以下6项，如果GWAS data不是来自https://gwas.mrcieu.ac.uk/datasets/，需要提供。如果是，则不需要填
chr_exp="chromosome"
pos_exp="base_pair_location"
snp_exp = "variant_id"
beta_exp = "beta" #如果有beta,用beta;如果没有，则用OR或者logOR
se_exp = "standard_error"
effect_allele_exp = "effect_allele"
other_allele_exp = "other_allele"
eaf_exp = "effect_allele_frequency"
pval_exp = "p_value"


#---outcome data---
out_data=c("D:/Desktop/LXbatmr 2024-10-10(v3.38)/data/drugMR/outcome/finngen_R10_NAFLD_standard.gz",
          "D:/Desktop/LXbatmr 2024-10-10(v3.38)/data/drugMR/outcome/GCST90271622_standard.gz") #* 结局，必填

# outcome_name="NAFLD"  # 结局，一般是疾病名称，选填

# 以下6项，如果GWAS data不是来自https://gwas.mrcieu.ac.uk/datasets/，需要提供。如果是，则不需要填
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
exp_p=1e-6        #* 关联性，必填
clump_kb = 100  #* 连锁不平衡kb, 必填
clump_r2 = 0.3  #* 连锁不平衡r2,必填

}
#---run the main function---
devtools::load_all()

# LXdrugMR()


#---The end--------------------------------------------------------------------




