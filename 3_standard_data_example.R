
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXbatmr")

library(LXbatmr)

#-------------------------

rm(list=ls())
gc()

# devtools::load_all()

gwas_file <- c("E:/NAFLD/finngen_R10_NAFLD.gz"
                   )

chr="#chrom"
pos="pos"
SNP="rsids"
beta="beta"
pval="pval"
se="sebeta"
effect_allele="alt"
other_allele="ref"
eaf="af_alt"

clum_p=5e-8
clum_p=NA

standard_data (gwas_file=gwas_file,
                          chr=chr,
                          pos=pos,
                          SNP=SNP,
                          beta=beta,
                          pval=pval,
                          se=se,
                          effect_allele=effect_allele,
                          other_allele=other_allele,
                          eaf=eaf,
                          clum_p=clum_p)




head(df)
