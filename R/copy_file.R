
#-----------------------------------

copy_file <- function(file_list,old_dir,new_dir){

#------------
inst_packages()

#------------

file.list <- read.xlsx(file_list) %>% .[,1]
old.list <- dir(old_dir)

if(!dir.exists(new_dir))
   dir.create(new_dir)

nonfile <-data.frame()

foreach(x=file.list,.errorhandling = "pass") %do% {

print(paste0(x, " is being copied......"))
if(!any(grepl(x,old.list,ignore.case = T))){
  message (paste0(x, " can not be found in the '",old_dir,"' folder"))
  nonfile <- rbind(nonfile,x)
  }


copy.file <- grep(x,old.list,ignore.case = T,value = T)

from.file <- paste0(old_dir,"/",copy.file)
to.file <- paste0(new_dir,"/",copy.file)

file.copy(from = from.file,to = to.file, overwrite = T,copy.mode = TRUE )

#-----foreach end--------
}

colnames(nonfile)[1] <- "mis_files"
write.xlsx(nonfile,paste0(new_dir,"/mis_files.xlsx"))

print(paste0("The result can be found in '",new_dir,"' folder"))

}
