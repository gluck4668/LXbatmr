
#--------------------------------------
clear_exp_data <- function(data_dir=data_dir,snp_exp=snp_exp,beta_exp=beta_exp,
               se_exp=se_exp,
               effect_allele_exp=effect_allele_exp,
               other_allele_exp=other_allele_exp,
               eaf_exp=eaf_exp,pval_exp=pval_exp,
               clum_p=clum_p,clump_kb=clump_kb,clump_r2=clump_r2) {

#--------
if(!dir.exists(data_dir) | length(dir(data_dir))==0)
stop(paste0("The folder '",data_dir,"' was not existed or it was empty, please check it."))

#-----R packages-----------------------
 inst_packages()

#-----creating dir--------------------
dir_file=paste0(data_dir," (cleared result)")
 if(!dir.exists(dir_file))
  dir.create(dir_file)

#-------reading data files------------
exp_list <- dir(data_dir)
exp_list <- paste0(data_dir,"/",exp_list) %>% data.frame()
names(exp_list) <- "file_id"

#-------foreach analysis-----------------------------------------
foreach(x=c(1:nrow(exp_list)),.errorhandling = "pass") %do% {

print(paste0("It is number ",x," of ", nrow(exp_list)))

file_type <-  tools::file_ext(exp_list[x,1])  # 文件扩展名

if(file_type=="txt"){
  tryCatch(exp_df <-read.table(exp_list[x,1],header = T,sep = "\t"),
             error=function(e){message("The data cann't be readed by read.table. Try read_table....")})
    exp_df <-readr::read_table(exp_list[x,1])
  } else {
    if(file_type=="gz" | file_type=="tsv")
      exp_df <- fread(exp_list[x,1]) else
        eval(str2expression(exp_df <- paste0("read.",file_type,"('",exp_list[x,1],"')")))
  }

#------exposure beta and OR transformation---------
  if(!grepl("log",beta_exp,ignore.case = T) | !grepl("ln",beta_exp,ignore.case = T))
  {
    if(grepl("odd",beta_exp,ignore.case = T) | grepl("OR",beta_exp,ignore.case = T)){
      exp_df <- exp_df %>% as.data.frame()
      p <- grep(beta_exp,names(exp_df)) %>% as.numeric()
      exp_df[,p] <- log(exp_df[,p])}
  }

#-----format data----------------------------------
exp_df <- data.frame(exp_df)
p_sit <- grep(pval_exp,names(exp_df))
exp_df[,p_sit] <- as.numeric(exp_df[,p_sit])

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

idx01 <- str_extract(exp_list[x,1],"(?<=/)[^/]+$")
if(grepl("_",idx01))
  data_id <- str_extract(idx01,".*(?=_)") else
    data_id <- idx01

exp_df$id.exposure <- data_id

#----calculating F value-------------
# b <- exp_df$beta.exposure
# eaf <- exp_df$eaf.exposure
# se <- exp_df$se.exposure
# sam <- exp_df$samplesize.exposure

# exp_df$R2<-(2*(b^2)*eaf*(1-eaf)/(2*(b^2)*eaf*(1-eaf)+2*sam*(se^2)*eaf*(1-eaf)))
# exp_df$F<-exp_df$R2*(exp_df$sam-2)/(1-exp_df$R2)

file_name <- gsub("[.]","_",idx01)
write.csv(exp_df,paste0(dir_file,"/",file_name,".csv"))

p = round(100*x/nrow(exp_list),4)
message(paste0("......",p,"% was completed......"))
#----the end of foreach--------------
}

print(paste0("The cleared data can be found in the folder of '", dir_file,"', good luck!"))
#------the end-----
}





