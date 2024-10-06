
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXbatmr")

library(LXbatmr)

#-------------------------

rm(list=ls())
gc()

# devtools::load_all()

gwas_file <- c("D:/Desktop/LXbatmr 2024-10-01(v3.22)/data/exp_gwas/GCST90199631_buildGRCh38.tsv.gz",
               "D:/Desktop/LXbatmr 2024-10-01(v3.22)/data/exp_gwas/GCST90199634_buildGRCh38.tsv.gz",
               "D:/Desktop/LXbatmr 2024-10-01(v3.22)/data/exp_gwas/GCST90199645_buildGRCh38.tsv.gz")

chr="chromosome"
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

