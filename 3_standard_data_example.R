
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXbatmr")

library(LXbatmr)

#-------------------------

rm(list=ls())
gc()

devtools::load_all()

{
gwas_file <- c("data/ieu/ieu-a-300.vcf.gz"
                   )

chr="chr"
pos="pos"
SNP="SNP"
beta="beta"
pval="pval"
se="se"
effect_allele="effect_allele"
other_allele="other_allele"
eaf="eaf"

clum_p=5e-8
# clum_p=NA

write.gz=TRUE

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
                          clum_p=clum_p,
                          write.gz=write.gz
                           )


}


# library(data.table)
# head(fread(gwas_file[1]))






