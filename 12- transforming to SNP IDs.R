
# BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
# BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh38")
# https://bioconductor.org/packages/3.19/data/annotation/src/contrib/SNPlocs.Hsapiens.dbSNP144.GRCh38_0.99.20.tar.gz
# https://bioconductor.org/packages/3.19/data/annotation/src/contrib/SNPlocs.Hsapiens.dbSNP144.GRCh37_0.99.20.tar.gz

library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(data.table)
library(stringr)
library(dplyr)

rm(list=ls())

file_dir = "test"
built="GRCh37"
#-----------------------------------------------------------
df <- head(fread(dir(file_dir,full.names=T)))
df$SNP <- NA

x=1

for (x in nrow(df)){
  df_chr= str_extract(df[x,1],".*?(?=\\:)") %>%  as.character() # 1
  df_pos= str_extract(df[x,1],"(?<=\\:).*?(?=\\:)")  # 930158

  if(grepl("GRCh37",built,ignore.case = T))
    snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37 else
      snps <- SNPlocs.Hsapiens.dbSNP144.GRCh38

  chr_snps<-snpsBySeqname(snps, df_chr) # 提取 chr=1 的全部snps
  pos_idx<-match(df_pos,pos(chr_snps))  # 匹配pos=930158，所在的行
  rsids<-mcols(chr_snps)$RefSNP_id[pos_idx]

  df$SNP[x] <- rsids

}



