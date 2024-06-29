

library(data.table)
library(openxlsx)
library(dplyr)
library(stringr)


pleiotropy = c("Table S3 GCST90275041.tsv___pleiotropy_test.xlsx")  # 多效性

heterogeneity = c("Table S4 GCST90275041.tsv___heterogeneity_test.xlsx") #敏感性

outcome_mr_ivw_significant = c("Table S5 GCST90275041.tsv___outcome_mr_all_significant.xlsx")


#-----------------------------------------------------------------------------------

ple <- read.xlsx(pleiotropy) %>% filter(pval<0.05)  %>% select(id.exposure)
het <- read.xlsx(heterogeneity) %>% filter(Q_pval<0.05) %>% select(id.exposure)

exc <-  unique(c(ple$id.exposure,het$id.exposure)) # 去重


#--------只一个文件-----
#exc <- read.xlsx("finngen_R10_NAFLD.gz___rev_or_ivw .xlsx") %>% filter(pval<0.05)


out_ivw <- read.xlsx(outcome_mr_ivw_significant) 

out_id <-  out_ivw[!(out_ivw$id.exposure %in% exc),]

write.xlsx(out_id,paste0(str_extract(outcome_mr_ivw_significant,".*?(?=___)"),"_ivw_significant ",Sys.Date(),".xlsx"))


