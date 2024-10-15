
bat_tsmr <- function(exp_data_dir,snp_exp,beta_exp,se_exp,effect_allele_exp,
                   other_allele_exp,eaf_exp,pval_exp,clum_p,clump_kb,clump_r2,
                   out_data_file,snp_out,beta_out,se_out,effect_allele_out,
                   other_allele_out,eaf_out,pval_out){

#-----R packages--------------------------------
inst_packages()

#-----mr dir-----------------------------------
mr_dir <- paste0("mr_results","_",Sys.Date())

if(!dir.exists(mr_dir))
    dir.create(mr_dir)

#----exposure and outcome files information----
exp_list <- dir(exp_data_dir,full.names = T) %>% as.character()
out_list <- out_data_file %>% as.character()

#------creating empty data frames----
mr_all <- data.frame()
mr_all_significant <- data.frame()
mr_ivw_significant <- data.frame()
heterogeneity_test <- data.frame()
pleiotropy_test <- data.frame()

#-------x foreach analysis-----------------------------------------------------

foreach(x=c(1:length(exp_list)), .errorhandling = "pass") %do% {

print(paste0("Number ",x, " of ",length(exp_list)," exposure data is being processed..." ))

file_type <- tools::file_ext(exp_list[x]) %>% tolower() # tools::file_ext查看文件扩展名

#-----read exposure data----------------------------
if(file_type=="txt"){
  tryCatch(exp_df <-read.table(exp_list[x],header = T,sep = "\t"),
           error=function(e){message("The exposure data cann't be readed by read.table. Try read_table function....")})
           if(!exists("exp_df") | nrow(exp_df)==0){
           exp_df <-readr::read_table(exp_list[x],show_col_types = F)}
  } else {
     if(file_type=="gz" | file_type=="tsv")
        exp_df <- fread(exp_data_file) else
        exp_df <- eval(str2expression(paste0("read.",file_type,"('",exp_list[x],"')")))
        }

exp_df <- exp_df %>% as.data.frame() # read outcome data
pval_exp_sit <- grep(pval_exp,names(exp_df))
exp_df[,pval_exp_sit] <- as.numeric(exp_df[,pval_exp_sit]) # p_val转换成数字格式

beta_exp_sit <- grep(beta_exp,names(exp_df))
exp_df[,beta_exp_sit] <- as.numeric(exp_df[,beta_exp_sit])

eaf_exp_sit <- grep(eaf_exp,names(exp_df))
exp_df[,eaf_exp_sit] <- as.numeric(exp_df[,eaf_exp_sit])

#------exposure data file name-------------------
exp_data_file_name <- str_extract(exp_list[x],"[^/]+$") %>% str_extract(".*?(?=\\.)")

#------exposure beta and OR transformation---------
if(!grepl("log",beta_exp,ignore.case = T) | !grepl("ln",beta_exp,ignore.case = T))
  {
  if(grepl("odd",beta_exp,ignore.case = T) | grepl("OR",beta_exp,ignore.case = T)){
    exp_df <- exp_df %>% as.data.frame()
    p <- grep(beta_exp,names(exp_df)) %>% as.numeric()
    exp_df[,p] <- log(exp_df[,p])}
   }

#--------------format exposure data---------------
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

 #----outcome ----------------------------------------------------------------

 foreach(y=c(1:length(out_list)), .errorhandling = "pass") %do% {

 message (paste0("Number ",y, " of ",length(out_list)," outcome data is being processed..." ))

 #-------reading data ----------------
 out_data <- fread(out_list[y]) %>% as.data.frame() # read outcome data

 pval_out_sit <- grep(pval_out,names(out_data))
 out_data[,pval_out_sit] <- as.numeric(out_data[,pval_out_sit]) # p_val转换成数字格式

 beta_out_sit <- grep(beta_out,names(out_data))
 out_data[,beta_out_sit] <- as.numeric(out_data[,beta_out_sit])

 eaf_out_sit <- grep(eaf_out,names(out_data))
 out_data[,eaf_out_sit] <- as.numeric(out_data[,eaf_out_sit])

 #------outcome data file name-------------------
 out_data_file_name <- str_extract(out_list[y],"[^/]+$") %>% str_extract(".*?(?=\\.)")

 #------outcome beta and OR transformation---------
 if(!grepl("log",beta_out,ignore.case = T) | !grepl("ln",beta_out,ignore.case = T))
   {
   if(grepl("odd",beta_out,ignore.case = T) | grepl("OR",beta_out,ignore.case = T)){
     out_data <- out_data %>% as.data.frame()
     p <- grep(beta_out,names(out_data))
     out_data[,p] <- log(out_data[,p])}
   }

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

if(exists("harm_df")){
  harm_df$id.exposure <- exp_data_file_name
  harm_df$id.outcome <- out_data_file_name
  }

#------MR analysis-----------------------------
#res <- mr(harm_res)
res <- mr(harm_df)
res$id.exposure <- exp_data_file_name
res$id.outcome <- out_data_file_name

res_or_all <- generate_odds_ratios(res) %>% data.frame()
ivw_sit <- grep("Inverse variance weighted",res_or_all$method,ignore.case = T)
res_or <- res_or_all[ivw_sit,]

# mr_result=rbind(mr_result,cbind(id=x, method=res_or$method, pvalue= round(res_or$pval,5),
#                                beta=res_or$b,lo_ci=res_or$lo_ci,up_ci=res_or$up_ci,
#                                or=res_or$or,or_lci95=res_or$or_lci95,or_uci95=res_or$or_uci95))

#==============================================================================
mr_all <- rbind(mr_all,res_or_all)

if(res_or$pval<0.05){
#-------------------------------------
  mr_all_significant <- rbind(mr_all_significant,res_or_all)
  mr_ivw_significant <- rbind(mr_ivw_significant,res_or)

#-------------------------------------
  #异质性检验，没有异质性的。
  mr_hete <- mr_heterogeneity(harm_df,method_list=c("mr_ivw"))
  heterogeneity_test <- rbind(heterogeneity_test,mr_hete)


  #多水平校验，这里是没有多水平效应的
  mr_pleio <- mr_pleiotropy_test(harm_df)
  pleiotropy_test <- rbind(pleiotropy_test,mr_pleio)

#-------------------------------------
 is.het = all(!is.na(mr_hete$Q_pval),mr_hete$Q_pval<0.05)
 is.ple = all(!is.na(mr_pleio$pval),mr_pleio$pval<0.05)

 #-----creating outcome dir-----------------------------
 if(is.het | is.ple)
 out_dir=paste0(mr_dir,"/",out_data_file_name,"__(",exp_data_file_name,")__excluded") else
     out_dir=paste0(mr_dir,"/",out_data_file_name,"__(",exp_data_file_name,")")

 if(!dir.exists(out_dir))
   dir.create(out_dir)

 mr_save_file <-paste0(out_dir,"/",out_data_file_name,"__(",exp_data_file_name,")")

#-----------save files---------------
  write.xlsx(res_or, paste0(mr_save_file,"_OR_data.xlsx"))
  write.xlsx(harm_df, paste0(mr_save_file,"_harmonise_data.xlsx"))
  write.xlsx(mr_hete,  paste0(mr_save_file,"__heterogeneity_test.xlsx"))
  write.xlsx(mr_pleio, paste0(mr_save_file,"_pleiotropy_test.xlsx"))

#进行了一个（MR-PRESSO）检验，这个也是多水平效应检验，P值应该要大于0.05
  tryCatch(
  pres <- mr_presso(BetaOutcome="beta.outcome",
                    BetaExposure ="beta.exposure",
                    SdOutcome ="se.outcome",
                    SdExposure = "se.exposure",
                    OUTLIERtest =TRUE,DISTORTIONtest = TRUE,
                    data =harm_df, NbDistribution = 1000,
                    SignifThreshold = 0.05),
  error=function(e){message("There is no enough intrumental variables for PRESSO analysis")}
       )

  if(exists("pres")){
    capture.output(pres,file = paste0(mr_save_file,"_presso_test.txt"))
    write.xlsx(pres, paste0(mr_save_file,"_presso_test.xlsx"))
    }

  # 查看离群值
  #Outlier <- pres$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`

  #if(!is.null(Outlier) & is.numeric(Outlier) )
  #  harm_df <- mr_data[-Outlier,] # 剔除离群值61和70


  #Leave-one-out analysis是指逐步剔除SNP后观察剩余的稳定性，理想的是剔除后变化不大
  while (!is.null(dev.list()))  dev.off()#关闭Plots
  par(cex = 0.6)
  par(mar = c(0,4,2,0))
  png(filename = paste0(mr_save_file,"_leaveoneout_plot.png"),
      width=800, height=600,units = "px",res = 100)
  par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

  single<- mr_leaveoneout(harm_df)
  p1<- mr_leaveoneout_plot(single)
  print(p1)

  dev.off()

  #散点图
  while (!is.null(dev.list()))  dev.off()#关闭Plots
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  png(filename = paste0(mr_save_file,"_scatter_plot.png"),
      width=800, height=600,units = "px",res = 100)
  par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

  p2 <- mr_scatter_plot(res,harm_df)
  print(p2)

  dev.off()


  #绘制森林图
  while (!is.null(dev.list()))  dev.off()#关闭Plots
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  png(filename = paste0(mr_save_file,"_forest_plot.png"),
      width=800, height=600,units = "px",res = 100)
  par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

  res_single<- mr_singlesnp(harm_df)
  p3 <- mr_forest_plot(res_single)
  print(p3)

  dev.off()

  #绘制漏斗图，主要是看蓝线周围的散点是否对称
  while (!is.null(dev.list()))  dev.off()#关闭Plots
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  png(filename = paste0(mr_save_file,"_funnel_plot.png"),
      width=800, height=600,units = "px",res = 100)
  par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

  p4 <- mr_funnel_plot(res_single)
  print(p4)

  dev.off()#关闭Plots

} #---if(res_or$pval<0.05)  end ------------------------------

message(paste0("Number ",y, " of ",length(out_list)," outcome data was done." ))

} #----- y foreach(y=c(1:length(out_list)), .errorhandling = "pass") end-------


#print(paste0("......Processing ",x, "is done......",  ))

print(paste0("Number ",x, " of ",length(exp_list)," exposure data was completed." ))


}#-------x foreach analysis end------------------------------------------------

write.xlsx(mr_all,paste0(mr_dir,"/",out_data_file_name,"_outcome_mr_all.xlsx"))
# write.xlsx(heterogeneity_test,paste0(mr_dir,"/",out_data_file_name,"_heterogeneity_test.xlsx"))
# write.xlsx(pleiotropy_test,paste0(mr_dir,"/",out_data_file_name,"_pleiotropy_test.xlsx"))

if(nrow(mr_all_significant)>0){

  is.hete <- all(!is.na(heterogeneity_test$Q_pval),any(heterogeneity_test$Q_pval<0.05))
  if(is.hete)
     het_id <- heterogeneity_test %>% filter(Q_pval<0.05) %>% .$id.exposure else
       het_id <- NA

  is.plei <- all(!is.na(pleiotropy_test$pval),any(pleiotropy_test$pval<0.05))
  if(is.plei)
    ple_id <- pleiotropy_test %>% filter(pval<0.05) %>% .$id.exposure else
      ple_id <- NA

  exclud_id <- c(het_id,ple_id) %>% unique() %>% na.omit()

  if(length(exclud_id)>0){
    mr_all_significant <- mr_all_significant[!(mr_all_significant$id.exposure %in% exclud_id),]
    mr_ivw_significant <- mr_ivw_significant[!(mr_ivw_significant$id.exposure %in% exclud_id),]
  }


  write.xlsx(mr_all_significant,paste0(mr_dir,"/",out_data_file_name,"_outcome_mr_all_significant.xlsx"))
  write.xlsx(mr_ivw_significant,paste0(mr_dir,"/",out_data_file_name,"_outcome_mr_ivw_significant.xlsx"))


  }



print(paste0("The MR analysis result can be found in the folder of '",mr_dir,"' "))
#-------bat_tsmr end--------
}










