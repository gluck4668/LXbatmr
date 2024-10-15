
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXbatmr")

library(LXbatmr)

#-------------------------

rm(list=ls())
gc()

# devtools::load_all()

gwas_file <- c("D:/Desktop/LXbatmr 2024-10-06(v3.23)/coloc_data/outcome/finngen_R10_NAFLD.gz",
               "D:/Desktop/LXbatmr 2024-10-10(v3.37)/data/drugMR/expouse/ieu-a-300.vcf.gz"
               )

chr="X.chrom"
pos="base_pair_location"
SNP="variant_id"
beta="beta"
pval="p_value"
se="standard_error"
effect_allele="effect_allele"
other_allele="other_allele"
eaf="effect_allele_frequency"


standard_data (gwas_file=gwas_file,
                          chr=chr,
                          pos=pos,
                          SNP=SNP,
                          beta=beta,
                          pval=pval,
                          se=se,
                          effect_allele=effect_allele,
                          other_allele=other_allele,
                          eaf=eaf)





