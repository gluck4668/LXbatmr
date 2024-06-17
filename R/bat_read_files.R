

bat_read_files <- function(exp_data_dir){

inst_packages()

#dir_file="analysis results"
#  if(!dir.exists(dir_file))
#    dir.create(dir_file)

exp_list <- dir(exp_data_dir) %>% as.data.frame()
#exp_list <- exp_list[c(1:2),1]  %>% as.data.frame()
colnames(exp_list) <- "id"

data_list <- list()

files_list <- data.frame(vid=rep(NA,50))

#---------------------------------
foreach(x=c(1:nrow(exp_list)), .errorhandling = "pass") %do% {

print (paste0("Reading number ",x," of ", nrow(exp_list)))

file_type <- tolower(str_extract(exp_list[x,],"(?<=\\.)[^\\.]+$"))

exp_data_file <- paste0(exp_data_dir,"/",exp_list[x,1])

# str_extract("exp_list.vcf.gz.txt","(?<=\\.)[^\\.]+$")
if(file_type=="txt"){
  exp_df <-read_table2(exp_data_file)} else {
    if(file_type=="gz" | file_type=="tsv")
      exp_df <- fread(exp_data_file) else
        eval(str2expression(exp_df <- paste0("read.",file_type,"('",exp_data_file,"')")))
  }

#data_list[[exp_list[x,]]] <- exp_df[1,]

#id= str_extract(exp_list[x,],".*(?=_)")
val <- c(exp_list[x,],names(exp_df))
name_list <- c(val,rep(NA,50-length(val)))
files_list <- cbind(files_list,id=name_list)

message(paste0("Reading number ",x," of ", nrow(exp_list)," is done"))

#-------foreach end---------
}

files_title <- files_list[,-1]
colnames(files_title) <- files_title[1,]
files_title <- files_title[-1,]

write.xlsx(files_title,paste0(exp_data_dir,"/exposure_file_col_names.xlsx"))


result <- list(data_list=data_list,files_title=files_title)

return(result)

#-----bat_read_files end-----
}





