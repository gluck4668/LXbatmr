
standard_data <- function(gwas_file=gwas_file,
                          chr=chr,
                          pos=pos,
                          SNP=SNP,
                          beta=beta,
                          pval=pval,
                          se=se,
                          effect_allele=effect_allele,
                          other_allele=other_allele,
                          eaf=eaf){

library(dplyr)
library(stringr)
library(data.table)
library(foreach)


foreach(x=c(1:length(gwas_file)), .errorhandling = "pass") %do% {

  print(paste0("It is number ",x, " of ",length(gwas_file)))

  df <- fread(gwas_file[x])

  colnames(df)[grep("chr",colnames(df),ignore.case = TRUE)] <- "chromosome"

  colnames(df)[which(colnames(df)==chr)] <- "chromosome"

  colnames(df)[which(colnames(df)==pos)] <- "base_pair_location"
  colnames(df)[which(colnames(df)==SNP)] <- "variant_id"
  colnames(df)[which(colnames(df)==beta)] <- "beta"
  colnames(df)[which(colnames(df)==pval)] <- "p_value"
  colnames(df)[which(colnames(df)==se)] <- "standard_error"
  colnames(df)[which(colnames(df)==effect_allele)] <- "effect_allele"
  colnames(df)[which(colnames(df)==other_allele)] <- "other_allele"
  colnames(df)[which(colnames(df)==eaf)] <- "effect_allele_frequency"

  fwrite(df, paste0(str_extract(gwas_file[x],".*(?=\\.)"),"_standard.gz"))

  print (paste0(" The ",x, " of ",length(gwas_file), " was completed "))

} # foreach end

print(paste0("The standarded data was saved in '",str_extract(gwas_file[1],".*(?=\\/)"),"'"))

} # standard_data end



