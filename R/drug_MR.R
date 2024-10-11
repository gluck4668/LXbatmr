
drugMR <- function(target_gene_data=target_gene_data,
                     exp_data=exp_data,
                     chr_exp=chr_exp,
                     pos_exp=pos_exp,
                     snp_exp=snp_exp,
                     beta_exp=beta_exp,
                     se_exp=se_exp,
                     effect_allele_exp=effect_allele_exp,
                     other_allele_exp=other_allele_exp,
                     eaf_exp=eaf_exp,

                     out_data=out_data,
                     snp_out=snp_out,
                     beta_out=beta_out,
                     se_out=se_out,
                     effect_allele_out=effect_allele_out,
                     other_allele_out=other_allele_out,
                     eaf_out=eaf_out,
                     pval_out=pval_out,

                     pval_exp=pval_exp,
                     clump_kb=clump_kb,
                     clump_r2=clump_r2 ){

#-------R packages ----------------------
inst_packages ()

#---create dir----------------------------
dir_save <- paste0("drug_MR_results_",Sys.Date())
if(!dir.exists(dir_save))
    dir.create(dir_save)

#---data information----------------------
exp_files <- dir(exp_data,full.names = T)
target_gene_df <- read.xlsx(target_gene_data)
out_files <- out_data

#-----建立一个空白的dataframe，用于存储结果------------------------------------
drugMR_mr_all <- data.frame()
drugMR_ivw_all <- data.frame()

drugMR_mr_significant <- data.frame()
drugMR_ivw_significant <- data.frame()

#-----read exposure data file--------------------------------------------------
foreach(x=1:length(exp_files),.errorhandling = "pass") %do% {

print(paste0("Number ",x, " of ",length(exp_files)," exposure data is being processed."))

exp_file_x <- exp_files[x]

exp_file_type <- tools::file_ext(exp_file_x)  # 文件扩展名

if(tolower(exp_file_type)=="gz")
  exp_df <- fread(exp_file_x) else
  {if(tolower(exp_file_type)=="txt")
    exp_df <- read_table(exp_file_x) else
      exp_df <- eval(str2expression( paste0("read.",exp_file_type,"('",exp_file_x,"')")
                ))
  }

# head(exp_df)

#-----筛选exp data-----
exp_df <- eval(str2expression(paste0("distinct(exp_df,",snp_exp,",.keep_all = T)"))) #snp去重
exp_df <-eval(str2expression(paste0("subset(exp_df,",pval_exp,"<",exp_p,")"))) #筛选p值

if(nrow(exp_df)>0)
  print(paste0("Note: At the condition of exp_p<",exp_p,", there were ",nrow(exp_df), " SNPs in '",
               str_extract(exp_file_x,"(?<=\\/)[^/]+$"),"'.")) else
    stop(paste0("Note: At the condition of exp_p<",exp_p,", there was no any SNP in '",
                str_extract(exp_file_x,"(?<=\\/)[^/]+$"),"'."))

#---save the exp data---
exp_file_name <- str_extract(exp_file_x,"(?<=\\/)[^/]+$") %>% str_extract(".*?(?=[.])")

exp_file_save <- paste0(dir_save,"/",exp_file_name," (cleared).csv")
write.csv(exp_df,exp_file_save,row.names = F)

#----format the exposure data-----
exp_df <- read_exposure_data(filename = exp_file_save,
                             sep = ",",
                             chr_col = chr_exp,
                             pos_col = pos_exp,
                             snp_col = snp_exp,
                             beta_col = beta_exp,
                             se_col = se_exp,
                             effect_allele_col = effect_allele_exp,
                             other_allele_col = other_allele_exp,
                             eaf_col = eaf_exp,
                             pval_col = pval_exp,
                             ncase_col = NA,
                             ncontrol_col = NA
                             )

exp_df$id.exposure <- exp_file_name

if(file.exists(exp_file_save))
  file.remove(exp_file_save)

#------exposure beta and OR transformation---------
if(!grepl("log",beta_exp,ignore.case = T) | !grepl("ln",beta_exp,ignore.case = T)){
  if(grepl("odd",beta_exp,ignore.case = T) | grepl("OR",beta_exp,ignore.case = T)){
    exp_df <- exp_df %>% as.data.frame()
    exp_df$beta.exposure <- log(as.numeric(exp_df$beta.exposure))
    }
  }

#-------Linkage Disequilibrium (LD) test--------
ld_clum <- exp_ld_clum(exp_df,clump_kb,clump_r2)
exp_df <- ld_clum$exp_clum

#----target gene SNP---------------------------
foreach(y=1:nrow(target_gene_df),.errorhandling = "pass") %do% {

message(paste0("Number ",y, " of ",nrow(target_gene_df)," target genes is being processed."))

target_y <- target_gene_df[y,]

target_gene <- target_y[1,1]
chr_target <- target_y[1,2]
pos_start <- target_y[1,3]
pos_end <- target_y[1,4]

target_gene_snp <- subset(exp_df,chr.exposure==chr_target & pos.exposure>pos_start-100000 & pos.exposure<pos_end+100000)

if(nrow(target_gene_snp)==0)
  stop(paste0("Note: The '",exp_file_name, "' dataset didn't have any SNPs for the target gene (",target_gene,")"))

target_gene_name <- paste0(dir_save,"/",target_gene," snp.csv")
# write.csv(target_gene_snp,target_gene_name)

exp_data <- target_gene_snp
exp_snp <- data.frame(target_gene_snp$SNP)
names(exp_snp) <- "SNP"

#----read outcome data---------------------------------------------------------
foreach(z=1:length(out_files),.errorhandling = "pass") %do% {

message(paste0("Number ",z, " of ",length(out_files)," outcome data is being processed."))

out_df <- fread(out_files[z])
head(out_df)

colnames(out_df)[which(colnames(out_df)==snp_out)] <- "SNP"

out_df <- merge(exp_snp,out_df,by="SNP")
out_df <- distinct(out_df,SNP,.keep_all = T)

out_file_save <- str_extract(out_files[z],"(?<=\\/)[^/]+$") %>% str_extract(".*(?=\\.)")
out_file <- paste0(dir_save,"/",out_file_save," out_data (cleared).csv")
write.csv(out_df,out_file)

out_data <- read_outcome_data(
            snps=exp_snp$SNP,
            filename=out_file,
            sep = ",",
            snp_col = "SNP",
            beta_col = beta_out,
            se_col = se_out,
            effect_allele_col = effect_allele_out,
            other_allele_col = other_allele_out,
            eaf_col = eaf_out,
            pval_col = pval_out,
            ncase_col = NA,
            ncontrol_col = NA)

if(file.exists(out_file))
  file.remove(out_file)

#----harmonise data------------------------------------------------------------
harmonise_df <- harmonise_data(exposure_dat=data.frame(exp_data),
                              outcome_dat=data.frame(out_data),
                              action= 1)

gene_name <- target_gene_df[y,1]
exp_file_name <- str_extract(exp_file_x,"(?<=\\/)[^/]+$") %>% str_extract(".*?(?=[.])")
out_file_name <- str_extract(out_files[z],"(?<=\\/)[^/]+$") %>% str_extract(".*?(?=[.])")

harmonise_df$exposure <- exp_file_name
harmonise_df$id.exposure <- exp_file_name
harmonise_df$outcome <- out_file_name
harmonise_df$id.outcome <- out_file_name
harmonise_df$target.gene <- gene_name

#---MR analysis----------------------------------------------------------------
mr_analysis <- mr(harmonise_df) %>% as.data.frame()
mr_analysis$target.gene <- gene_name

p_ivw <- subset(mr_analysis,method=="Inverse variance weighted")$pval

if(p_ivw < 0.05){
result_file <- paste0(gene_name,"_(",exp_file_name,")_(",out_file_name,")___significant")
result_dir <- paste0(dir_save,"/",result_file)
  } else{
    result_file <- paste0(gene_name,"_(",exp_file_name,")_(",out_file_name,")")
    result_dir <- paste0(dir_save,"/",result_file)
    }

if(!dir.exists(result_dir))
  dir.create(result_dir)

write.xlsx(harmonise_df, paste0(result_dir,"/","harmonise_result (",result_file,").xlsx"))

write.xlsx(mr_analysis, paste0(result_dir,"/","MR_result (",result_file,").xlsx"))

#-----生成OR和可信区间---------------------------------------------------------
odds_ratios <- function (or_mr) {
               or_mr$lo_ci <- or_mr$b - 1.96 * or_mr$se
               or_mr$up_ci <- or_mr$b + 1.96 * or_mr$se
               or_mr$OR <- exp(or_mr$b)
               or_mr$OR_lci95 <- exp(or_mr$lo_ci)
               or_mr$OR_uci95 <- exp(or_mr$up_ci)
               or_mr$OR_inhibitor <- 1/exp(or_mr$b)
               or_mr$OR_lci95_inhibitor <- or_mr$OR_inhibitor - 1.96 * or_mr$se
               or_mr$OR_uci95_inhibitor <- or_mr$OR_inhibitor + 1.96 * or_mr$se
               return(or_mr)
               }

mr_OR <- odds_ratios(mr_analysis)
mr_OR$target.gene <- gene_name

write.xlsx(mr_OR, paste0(result_dir,"/","MR_OR_result (",result_file,").xlsx"))

#---MR-PRESSO检验，这个也是多水平效应检验，P值应该要大于0.05-------------------
presso_analysis <- mr_presso(BetaOutcome="beta.outcome",
                  BetaExposure ="beta.exposure",
                  SdOutcome ="se.outcome",
                  SdExposure = "se.exposure",
                  OUTLIERtest =TRUE,DISTORTIONtest = TRUE,
                  data =harmonise_df, NbDistribution = 1000,SignifThreshold = 0.05)

presso_analysis$target.gene <- gene_name

write.xlsx(presso_analysis, paste0(result_dir,"/","PRESSO_result (",result_file,").xlsx"))

#-----heterogeneity异质性检验，没有异质性的------------------------------------
heterogeneity_analysis <- mr_heterogeneity(harmonise_df,
                                           method_list=c("mr_ivw","mr_egger_regression"))

heterogeneity_analysis$target.gene <- gene_name

write.xlsx(heterogeneity_analysis, paste0(result_dir,"/","Heterogeneity_result (",result_file,").xlsx"))

#----pleiotropy多水平校验，这里是没有多水平效应的------------------------------
pleiotropy_analysis <- mr_pleiotropy_test(harmonise_df)

pleiotropy_analysis$target.gene <- gene_name

write.xlsx(pleiotropy_analysis, paste0(result_dir,"/","Pleiotropy_result (",result_file,").xlsx"))

#----Leave-one-out analysis是指逐步剔除SNP后观察剩余的稳定性，理想的是剔除后变化不大
while (!is.null(dev.list()))  dev.off()#关闭Plots
par(cex = 0.6)
par(mar = c(0,4,2,0))
png(filename = paste0(result_dir,"/","mr_leaveoneout_plot (",result_file,").png"),
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

single<- mr_leaveoneout(harmonise_df)
p1<- mr_leaveoneout_plot(single)
print(p1)

dev.off()

#------散点图------------------------------------------------------------------
while (!is.null(dev.list()))  dev.off()#关闭Plots
par(cex = 0.6);
par(mar = c(0,4,2,0))
png(filename = paste0(result_dir,"/","mr_smessageter_plot (",result_file,").png"),
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

p2 <- mr_scatter_plot(mr_analysis,harmonise_df)
print(p2)

dev.off()

#------绘制森林图--------------------------------------------------------------
while (!is.null(dev.list()))  dev.off()#关闭Plots
par(cex = 0.6);
par(mar = c(0,4,2,0))
png(filename = paste0(result_dir,"/","mr_forest_plot (",result_file,").png"),
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

res_single<- mr_singlesnp(harmonise_df)
p3 <- mr_forest_plot(res_single)
print(p3)

dev.off()

#-------绘制漏斗图，主要是看蓝线周围的散点是否对称-----------------------------
while (!is.null(dev.list()))  dev.off()#关闭Plots
par(cex = 0.6);
par(mar = c(0,4,2,0))
png(filename = paste0(result_dir,"/","funnel_plot (",result_file,").png"),
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

p4 <- mr_funnel_plot(res_single)
print(p4)

dev.off()#关闭Plots

#---save the results-----------------------------------------------------------
drugMR_mr_all <- rbind(drugMR_ivw_all,mr_OR)

ivw_set <- subset(mr_OR,method=="Inverse variance weighted")
drugMR_ivw_all <- rbind(drugMR_mr_all,ivw_set)

if(ivw_set$pval < 0.05){
  drugMR_mr_significant <- rbind(drugMR_mr_significant,mr_OR)
  drugMR_ivw_significant <- rbind(drugMR_ivw_significant,ivw_set)
  }

message(paste0("Number ",z, " of ",length(out_files)," outcome data was doned."))

} # foreach(z=1:length(out_files),.errorhandling = "pass") end

message(paste0("Number ",y, " of ",nrow(target_gene_df)," target genes was doned."))

} #foreach(y=1:nrow(target_gene_df),.errorhandling = "pass") end

print(paste0("Number ",x, " of ",length(out_files)," exposure data was doned."))
} #foreach(x=1:length(out_files),.errorhandling = "pass") end


#---所有循环end----------------------------------------------------------------

#----save the all results------------------------------------------------------
write.xlsx(drugMR_mr_all, paste0(dir_save,"/","drugMR_mr_all.xlsx"))
write.xlsx(drugMR_ivw_all, paste0(dir_save,"/","drugMR_ivw_all.xlsx"))

if(nrow(drugMR_ivw_significant)>0){
write.xlsx(drugMR_mr_significant, paste0(dir_save,"/","significant_drugMR_mr.xlsx"))
write.xlsx(drugMR_ivw_significant, paste0(dir_save,"/","significant_drugMR_ivw.xlsx"))
}

#------分析结束----------------------------------------------------------------
print(paste0("The results of drug MR analysis have been saved in the folder '",dir_save,"'."))

}
