
#----------------------------------------------------------#
cir_pic <- function(or_file,mr_method,cir_txt_sit,
                    gap.degree,
                    start.degree,
                    track.height){
#-----R packages-----------------------
com_packages()

#-------------------------------------
dir_or <- paste0("circos_plot_result_",Sys.Date())
if(!dir.exists(dir_or))
   dir.create(dir_or)

#-----reading mr or data---------------
file_type <- tools::file_ext(or_file) #文件扩展名

if(file_type=="txt"){
  tryCatch(or_data <-read.table(or_file,header = T,sep = "\t"),
           error=function(e){message("The OR data file cann't be readed by read.table. Try read_table function....")})
  if(!exists("or_df") | nrow(or_df)==0){
    or_data <-readr::read_table(or_file,show_col_types = F)}
} else {
  if(file_type=="gz" | file_type=="tsv")
    or_data <- fread(or_file) else
      or_data <- eval(str2expression(paste0("read.",file_type,"('",or_file,"')")))
}

# or_data$method[or_data$method=="Inverse variance weighted"]="IVW"
or_df <- or_data[(or_data$method %in% mr_method),] %>% as.data.frame()

for(a in 1:nrow(or_df)){
if(grepl("_",or_df[a,1]))
  or_df[a,1] <- str_extract(or_df[a,1],".*?(?=_)")

if(grepl("\\.",or_df[a,1]))
  or_df[a,1] <- str_extract(or_df[a,1],".*?(?=\\.)")
}

for(b in 1:nrow(or_df)){
  if(grepl("_",or_df[b,2]))
    or_df[b,2] <- str_extract(or_df[b,2],".*?(?=_)")

  if(grepl("\\.",or_df[b,2]))
    or_df[b,2] <- str_extract(or_df[b,2],".*?(?=\\.)")
}

outcome_num <- or_df$id.outcome[!duplicated(or_df$id.outcome)]
out_list <- list()

library(foreach)
foreach(x=outcome_num,.errorhandling = "pass") %do% {
out_list[[x]] <- subset(or_df,id.outcome==x)
write.csv(subset(or_df,id.outcome==x),paste0(dir_or,"/",x,"_circos_outcome_data.csv"))
}

#-----foreach(y=outcome_num--------------------------
foreach(y=outcome_num,.errorhandling = "pass") %do% {
mat=acast(out_list[[y]], id.exposure~method, value.var="pval")

#mat <- na.omit(mat)

col_color = colorRamp2(c(0.001,0.01, 0.02, 0.05,0.09,
                         0.1, 0.2,0.5,0.7, 0.9),
                       c("#ff0000","#ff3300","#ff6600", "#ff9900","#ffff00",
                         "white", "#00ccff","#0066ff","#0033ff","#0000ff" ))


#while (!is.null(dev.list()))   dev.off()
#circos.clear()
#circos.par(gap.after=c(20))
#pdf(file =paste0(out_data_file,"_circos_plot.pdf"),
#    width=8, height=8)


#--------
while (!is.null(dev.list()))   dev.off()
{
circos.clear()
circos.par(gap.after=c(20))

file_save <-str_extract(or_file,"(?<=/)[^/]+$") %>% str_extract(".*?(?=\\.)")
png(filename = paste0(dir_or,"/",y,"_circos_plot.png"),
    width=1200, height=1200,units = "px",res = 150)

par(mai = c(0.5, 0.5, 0.8, 0.5)) # 边距 c(bottom, left, top, right)

# 设置布局参数
circos.par(
  gap.degree = gap.degree,  # 圆环之间的角度间隔
 # cell.padding = c("cell", 0.02, "mm"),  # 单元填充
  start.degree = start.degree  # 起始角度
)

circos.heatmap(mat, # 数据矩阵
               col = col_color,
               rownames.cex = 0.8,
               dend.side = "inside",
               rownames.side = "outside",
               track.height = track.height,  # 扇形高度
               bg.border = "black")
# 添加标题
title( paste0("Circos.heatmap (outcome data: ",y,")"))

}

# while (!is.null(dev.list()))   dev.off()

cn = colnames(mat)

n = length(cn)
#cir_txt_sit = 0.35+(n:1)*0.2 # methods 的位置和行距

circos.text(rep(CELL_META$cell.xlim[2], n) +
              convert_x(1, "mm"),
              cir_txt_sit,
             # 0.35+(n:1)*0.2, # methods 的位置和行距
              cn,
              cex =0.8,  # 字体大小
              adj = c(0, 0.5),
              facing = "inside")

circos.clear()

lgd = Legend(title="Pvalue", col_fun=col_color,at = c(0.01, 0.05, 0.1,0.2, 0.3),size = 12)
grid.draw(lgd)


while (!is.null(dev.list()))   dev.off()

} # foreach (y) end


print(paste0("The circos plot can be found in the folder of '",dir_or,"' "))


#---------the end-------------
}





