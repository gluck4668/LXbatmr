
venn_analysis <- function(or_dir){

#--------
if(!dir.exists(or_dir) | length(dir(or_dir))==0)
  stop(paste0("The folder '",or_dir,"' was not existed or it was empty, please check it."))

#---R packages--------
com_packages()

#----------------------
dir_file <- paste0("venn_result ",Sys.Date())
if(!dir.exists(dir_file))
  dir.create(dir_file)

#------reading data----------------
or_file_list <- dir(or_dir,full.names = T)

or_list <- str_extract(or_file_list,"(?<=/)[^/]+$") %>% str_extract(".*(?=___)") %>% str_extract(".*(?=\\.)")


#-----merge data------------------
venn_id_list <- list()
venn_files <- list()

for(x in or_file_list){
  df <- read.xlsx(x)
  venn_id_list[[x]] <- df$id.exposure
  venn_files[[x]] <- df
}


#-----venn analysis-----------------------------------------
names(venn_id_list) <- or_list

gg_ven <- ggvenn(venn_id_list,columns = or_list,
                 stroke_linetype = 1,text_size = 6,stroke_size = 0.5,
                 show_percentage = T,set_name_size = 4,
                 fill_color = brewer.pal(n=length(names(venn_id_list)),name = "Set3"))

ggsave(filename = paste0(dir_file,"/venn_plot ",Sys.Date(),".png"),plot = gg_ven,
       width = 800,height = 600,units = "px",dpi = 150)

# int <- get.venn.partitions(venn_id_list) # 取交集
#  inter_id <- int$..values..[[1]]

#------------------
inter_id <-reduce(venn_id_list,intersect) # 取交集，用更简单的purrr包的reduce()
venn_df <- data.frame(inter_id)
write.xlsx(venn_df,paste0(dir_file,"/venn_data ",Sys.Date(),".xlsx"))


#---------------------------------------
for(x in names(venn_files)){
  df <- venn_files[[x]] %>% data.frame()
  df <- df[(df$id.exposure %in% inter_id),]

  filename <- str_extract(x,"(?<=/)[^/]+$") %>% str_extract(".*(?=___)") %>% str_extract(".*(?=\\.)")

  write.xlsx(df,paste0(dir_file,"/",filename,"_venn_intersect ",Sys.Date(),".xlsx"))
  }


print(paste0("The venn analysis results can be found in the folder of '",dir_file,"' "))

gg_ven

}#----the end of venn functioin----
