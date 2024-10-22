
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXbatmr")

library(LXbatmr)


rm(list=ls())
gc()

 devtools::load_all()

#-----eqtl data (exposure data)-------------------------------------
{

#暴露数据：可以是id或本地数据，但如果是vcf.gz格式，还需要下载其index文件（tbi格式）
#eqtlID="THRB_cis_eQTLs_data.xlsx"
eqtlID= "D:/Desktop/LXbatmr_shiny (study ver)/data/coloc_data/expouse"
type1 = "quant" #数据类型，一般有分类变量"cc"和连续变量"quant"

#---查询exposure数据参数
#library(data.table)
# head(fread(dir(eqtlID,full.names = T)[1]))

#---如果是ieu网站数据（id或本地文件均可）不需要填下面信息；其他来源的数据要填
chr_eqtl="chromosome"
pos_eqtl="base_pair_location"
SNP_eqtl="SNP"
beta_eqtl="beta"
pval_eqtl="p_value"
se_eqtl="standard_error"
eaf_eqtl="effect_allele_frequency"

samplesize_eqtl_file= "D:/Desktop/LXbatmr_shiny (study ver)/data/coloc_data/expouse_gwas_id_info.xlsx" # "meta_GCST90091033_622_gwas_id_info.xlsx"
# NA # 如果数据来自GWAS Catalog，可以填NA；如果来自其他数据库，则需要提供含有ID和samplesize的excel表格


#-------gwas data (outcome data) ----------------------------------

#结局数据：可以是id或本地数据，但如果是vcf.gz格式，还需要下载其index文件（tbi格式）
outcomeID =c("D:/Desktop/LXbatmr_shiny (study ver)/data/coloc_data/outcome/finngen_R10_NAFLD_standard.gz",
             " D:/Desktop/LXbatmr_shiny (study ver)/data/coloc_data/outcome/GCST90275041_standard.gz") #结局数据
type2 = "cc"

#---查询outcome数据参数
# head(fread(outcomeID[1]))

#---如果是ieu网站数据（id或本地文件均可）不需要填下面信息；其他来源的数据要填
chr_gwas="chromosome"
pos_gwas="base_pair_location"
SNP_gwas="variant_id"
beta_gwas="beta"
pval_gwas="p_value"
se_gwas="standard_error"
eaf_gwas="effect_allele_frequency"

samplesize_gwas= "D:/Desktop/LXbatmr_shiny (study ver)/data/coloc_data/outcome_gwas_id_info.xlsx" #"outcome_sample_info.xlsx_gwas_info.xlsx"


#-----target gene---------------
# target_gene_file = "D:/Desktop/LXbatmr_shiny (study ver)/data/coloc_data/target.xlsx"  # 如果没有目标基因，可以填NA
 target_gene_file = NA

#-----#共定位阈值-----------------------------------------
SNP_PP_H4=0.75  #共定位阈值

}

#-------------------------------------------------------------------

LXcoloc_gwas (eqtlID=eqtlID,
            type1=type1,
            chr_eqtl=chr_eqtl,
            pos_eqtl=pos_eqtl,
            SNP_eqtl=SNP_eqtl,
            beta_eqtl=beta_eqtl,
            pval_eqtl=pval_eqtl,
            se_eqtl=se_eqtl,
            eaf_eqtl=eaf_eqtl,
            samplesize_eqtl_file=samplesize_eqtl_file,
            outcomeID =outcomeID,
            type2 = type2,
            chr_gwas=chr_gwas,
            pos_gwas=pos_gwas,
            SNP_gwas=SNP_gwas,
            beta_gwas=beta_gwas,
            pval_gwas=pval_gwas,
            se_gwas=se_gwas,
            eaf_gwas=eaf_gwas,
            samplesize_gwas= samplesize_gwas,
            target_gene_file =target_gene_file,
            SNP_PP_H4=SNP_PP_H4 )

#-------------end------------------------------------





