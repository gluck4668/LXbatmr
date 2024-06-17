

bat_tsmr <- function(exp_data_dir,snp_exp,beta_exp,se_exp,effect_allele_exp,
                   other_allele_exp,eaf_exp,pval_exp,clum_p,clump_kb,clump_r2,
                   out_data_file,snp_out,beta_out,se_out,effect_allele_out,
                   other_allele_out,eaf_out,pval_out){

#-----R packages-----------------------

inst_packages()

#-----creating dir--------------------
dir_file=paste0(out_data_file," (results)")
  if(!dir.exists(dir_file))
    dir.create(dir_file)

#------reading exposure data list----
if(grepl("\\.",exp_data_dir))
  exp_list <- read.xlsx(exp_data_dir) else
   exp_list <- dir(exp_data_dir) %>% as.data.frame()
#exp_list <- exp_list[c(1:2),1]  %>% as.data.frame()
colnames(exp_list)[1] <- "id"

#-------reading data ----------------
out_data <- fread(out_data_file) # read outcome data

#------creating empty data frames----
mr_result_all = data.frame() # including all data
mr_result_significant = data.frame() # including data with p<0.05

#-------foreach analysis------------
foreach(x=c(1:nrow(exp_list)), .errorhandling = "pass") %do% {

print(paste0("It is number ",x, " of ",nrow(exp_list)))

file_type <- tolower(str_extract(exp_list[x,],"(?<=\\.)[^\\.]+$"))

exp_data_file <- paste0(exp_data_dir,"/",exp_list[x,1])

# str_extract("exp_list.vcf.gz.txt","(?<=\\.)[^\\.]+$")
if(file_type=="txt"){
  tryCatch(exp_df <-read.table(exp_data_file,header = T,sep = "\t"),
           error=function(e){message("The data cann't be readed by read.table. Try read_table function....")})
           if(!exists("exp_df") | nrow(exp_df)==0){
           exp_df <-readr::read_table(exp_data_file,show_col_types = T)}
  } else {
     if(file_type=="gz")
        exp_df <- fread(exp_data_file) else
        eval(str2expression(exp_df <- paste0("read.",file_type,"('",exp_data_file,"')")))
        }


#--------------format exposure data---------------
exp_df <- eval(str2expression(paste0("subset(exp_df,",pval_exp,"<clum_p)")))

# exp_df <- format_exposure_data(filename = exp_df,
#                                   sep = ",",
#                                   snp_col = snp_exp,
#                                   beta_col = beta_exp,
#                                   se_col = se_exp,
#                                   effect_allele_col = effect_allele_exp,
#                                   other_allele_col = other_allele_exp,
#                                   eaf_col = eaf_exp,
#                                   pval_col = pval_exp)

#-------Linkage Disequilibrium (LD) test--------
# exp_df <- exp_ld_clum(exp_df,clump_kb=clump_kb,clump_r2=clump_r2)
# exp_df <- exp_df$exp_clum

#----------------format outcome data---------------
out_df <- format_outcome_data (filename = out_data,
                               snps = exp_df$SNP,
                               sep = ",",
                               snp_col = snp_out,
                               beta_col = beta_out,
                               se_col = se_out,
                               effect_allele_col = effect_allele_out,
                               other_allele_col = other_allele_out,
                               eaf_col = eaf_out,
                               pval_col = pval_out)

#--------harmonise data--------------------------
harm_df <- harmonise_data(exposure_dat = exp_df,
                          outcome_dat = out_df,
                          action=2)

#-------calculating R2, F,meanF------------------
#harm_df$R2 <- (2 * (harm_df$beta.exposure^2) * harm_df$eaf.exposure * (1 - harm_df$eaf.exposure)) /
#    (2 * (harm_df$beta.exposure^2) * harm_df$eaf.exposure * (1 - harm_df$eaf.exposure) +
#     2 * harm_df$samplesize.exposure*harm_df$eaf.exposure * (1 - harm_df$eaf.exposure) * harm_df$se.exposure^2)

#harm_df$f <- harm_df$R2 * (harm_df$samplesize.exposure - 2) / (1 - harm_df$R2)

#harm_df$meanf<- mean(harm_df$f)

#------screen out the data with F>10------------
#harm_res <- harm_df %>% subset(f>10)

#------MR analysis-----------------------------
#res <- mr(harm_res)
res <- mr(harm_df)

res_or = generate_odds_ratios(res)
res_or <- res_or[3,]

# mr_result=rbind(mr_result,cbind(id=x, method=res_or$method, pvalue= round(res_or$pval,5),
#                                beta=res_or$b,lo_ci=res_or$lo_ci,up_ci=res_or$up_ci,
#                                or=res_or$or,or_lci95=res_or$or_lci95,or_uci95=res_or$or_uci95))

mr_result_all <- rbind(mr_result_all,res_or)

if(res_or$pval<0.05){
  mr_result_significant <- rbind(mr_result_significant,res_or)
 }

#print(paste0("......Processing ",x, "is done......",  ))

p = round(100*x/nrow(exp_list),4)
print(paste0("......",p,"% is completed......"))

#-----foreach end----------
}

write.xlsx(mr_result_all,paste0(dir_file,"/bat_mr_result_all.xlsx"))

if(nrow(mr_result_significant)>0){
write.xlsx(mr_result_significant,paste0(dir_file,"/bat_mr_result_significant.xlsx"))}

#-------bat_tsmr end--------
}










