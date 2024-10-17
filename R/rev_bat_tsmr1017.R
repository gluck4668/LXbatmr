
rev_bat_tsmr <-function(exp_data_file,snp_exp,beta_exp,se_exp,effect_allele_exp,
              other_allele_exp,eaf_exp,pval_exp,clum_p,clump_kb,clump_r2,
              out_data_dir,snp_out,beta_out,se_out,effect_allele_out,
              other_allele_out,eaf_out,pval_out) {

#-------R packages--------------------------
inst_packages()

#-------creating dir-----------------------
#  str_extract(exp_data_file,"[^/]+$") %>% str_extract(".*?(?=\\.)")

dir_file=paste0("Rev_mr_result_",Sys.Date())
  if(!dir.exists(dir_file))
    dir.create(dir_file)

#-------reading exposure data ----------------
exp_data_name <-exp_data_file

x=1

foreach(x=c(1:length(exp_data_name)), .errorhandling = "pass") %do% {

print(paste0("Number: ",x, " of ",length(exp_data_file)," of exposure data is being processed..."))

file_type <- tools::file_ext(exp_data_file[x]) # 获取文件扩展名

exp_file_id <- str_extract(exp_data_name[x],"[^\\/]+$*.(?=\\.)")

if(grepl("_",exp_file_id))
  exp_file_id <- str_extract(exp_file_id,".*(?=_)")

dir_save <- paste0(dir_file,"/",exp_file_id,"_result")
if(!dir.exists(dir_save))
  dir.create(dir_save)

if(file_type=="txt"){
  tryCatch(exp_df <-read.table(exp_data_file[x],header = T,sep = "\t"),
           error=function(e){message("The data cann't be readed by read.table. Try read_table....")})
  exp_df <-readr::read_table(exp_data_file[x])
    } else {
    if(file_type=="gz" | file_type=="tsv")
    exp_df <- fread(exp_data_file[x]) else
      eval(str2expression(exp_df <- paste0("read.",file_type,"('",exp_data_file[x],"')")))
     }

#------exposure beta and OR transformation---------
if(!grepl("log",beta_exp,ignore.case = T) | !grepl("ln",beta_exp,ignore.case = T))
{
  if(grepl("odd",beta_exp,ignore.case = T) | grepl("OR",beta_exp,ignore.case = T)){
    exp_df <- exp_df %>% as.data.frame()
    p <- grep(beta_exp,names(exp_df)) %>% as.numeric()
    exp_df[,p] <- log(exp_df[,p])}
}

#------------------------------------------------------------------
exp_df <- eval(str2expression(paste0("subset(exp_df,",pval_exp,"<clum_p)")))

exp_df <- format_exposure_data(filename = exp_df,
                                   sep = ",",
                                   snp_col = snp_exp,
                                   beta_col = beta_exp,
                                   se_col = se_exp,
                                   effect_allele_col = effect_allele_exp,
                                   other_allele_col = other_allele_exp,
                                   eaf_col = eaf_exp,
                                   pval_col = pval_exp)

#-------Linkage Disequilibrium (LD) test--------
 exp_df <- exp_ld_clum(exp_df,clump_kb=clump_kb,clump_r2=clump_r2)
 exp_df <- exp_df$exp_clum
 exp_df$id.exposure <- exp_file_id



#------creating empty data frames----
 rev_or_all = data.frame() # including all data
 rev_or_ivw = data.frame()
 rev_or_ivw_significant = data.frame() # including data with p<0.05

#-------foreach analysis------------
 out_list <- dir(out_data_dir) %>% as.character()

y=1

foreach(y=c(1:length(out_list)), .errorhandling = "pass") %do% {

print(paste0("Number: ",y, " of ",length(out_list)," of outcome data is being processed..."))

file_type <-  tools::file_ext(out_list[y]) # 获取文件扩展名

out_data_file <- paste0(out_data_dir,"/",out_list[y])

out_file_id <- str_extract(out_data_file,"[^\\/]+$*.(?=\\.)")

if(grepl("_",out_file_id))
  out_file_id <- str_extract(out_file_id,".*(?=_)")

# str_extract("exp_list.vcf.gz.txt","(?<=\\.)[^\\.]+$")
if(file_type=="txt"){
  tryCatch(out_df <-read.table(out_data_file,header = T,sep = "\t"),
           error=function(e){message("The outcome data cann't be readed by read.table. Try read_table function....")})
  if(!exists("out_df") | nrow(out_df)==0){
    out_df <-readr::read_table(out_data_file,show_col_types = F)}
} else {
  if(file_type=="gz" | file_type=="tsv")
    out_df <- fread(out_data_file) else
      eval(str2expression(out_df <- paste0("read.",file_type,"('",out_data_file,"')")))
}

#------outcome beta and OR transformation---------
if(!grepl("log",beta_out,ignore.case = T) | !grepl("ln",beta_out,ignore.case = T))
{
  if(grepl("odd",beta_out,ignore.case = T) | grepl("OR",beta_out,ignore.case = T)){
    out_data <- out_data %>% as.data.frame()
    p <- grep(beta_out,names(out_data)) %>% as.numeric()
    out_data[,p] <- log(out_data[,p])}
}

#-----------------------------------------------
out_df <- format_outcome_data (filename = out_df,
                               snps = exp_df$SNP,
                               sep = ",",
                               snp_col = snp_out,
                               beta_col = beta_out,
                               se_col = se_out,
                               effect_allele_col = effect_allele_out,
                               other_allele_col = other_allele_out,
                               eaf_col = eaf_out,
                               pval_col = pval_out)

out_df$id.outcome <- out_file_id

#--------harmonise data--------------------------
harm_df <- harmonise_data(exposure_dat = exp_df,
                          outcome_dat = out_df,
                          action=2)

#-------calculating R2, F,meanF------------------
# harm_df$R2 <- (2 * (harm_df$beta.exposure^2) * harm_df$eaf.exposure * (1 - harm_df$eaf.exposure)) /
#    (2 * (harm_df$beta.exposure^2) * harm_df$eaf.exposure * (1 - harm_df$eaf.exposure) +
#     2 * harm_df$samplesize.exposure*harm_df$eaf.exposure * (1 - harm_df$eaf.exposure) * harm_df$se.exposure^2)

# harm_df$f <- harm_df$R2 * (harm_df$samplesize.exposure - 2) / (1 - harm_df$R2)

# harm_df$meanf<- mean(harm_df$f)

#------screen out the data with F>10------------
# harm_res <- harm_df %>% subset(f>10)

#------MR analysis-----------------------------
# res <- mr(harm_res)
res <- mr(harm_df)

#res$exposure <- "NAFLD"
#res$outcome <- "Metabolites"
or_all = generate_odds_ratios(res) %>% data.frame()

ivw_sit <- grep("Inverse variance weighted",or_all$method)
or_ivw <- or_all[ivw_sit,]

# mr_result=rbind(mr_result,cbind(id=x, method=res_or$method, pvalue= round(res_or$pval,5),
#                                beta=res_or$b,lo_ci=res_or$lo_ci,up_ci=res_or$up_ci,
#                                or=res_or$or,or_lci95=res_or$or_lci95,or_uci95=res_or$or_uci95))

rev_or_all <- rbind(rev_or_all,or_all)
rev_or_ivw <- rbind(rev_or_ivw,or_ivw)

if(or_ivw$pval<0.05){
  rev_or_ivw_significant  <- rbind(rev_or_ivw_significant,or_ivw) }

#print(paste0("......Processing ",x, "is done......",  ))

print(paste0("Number: ",x, " of ",length(out_list)," of outcome data was done."))


} #-----y foreach outcome end----------

write.xlsx(rev_or_all,paste0(dir_save,"/",exp_file_id,"___rev_or_all .xlsx"))
write.xlsx(rev_or_ivw,paste0(dir_save,"/",exp_file_id,"___rev_or_ivw .xlsx"))

if(nrow(rev_or_ivw_significant)>0){
  write.xlsx(rev_or_ivw_significant,paste0(dir_save,"/",exp_file_id,"___rev_or_ivw_significant.xlsx")) }


print(paste0("Number: ",x, " of ",length(exp_data_file)," of exposure data was completed."))

}# ----x foreach exposure end-----------------------


print(paste0("The MR analysis result can found in the folder of '",dir_file,"' "))
#-------Rev_bat_tsmr end--------
}







