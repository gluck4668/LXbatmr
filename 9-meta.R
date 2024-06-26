
devtools::load_all()

rm(list=ls())
#---------------------------------------------------------
# inst_packages()
library(meta)
library(openxlsx)

#--------------
or_dir <- dir("mr_or",full.names = T)

i=1
or_df <- read.xlsx(or_dir[i])

log_or <- log(or_df$or)
log_uci<-log(or_df$or_uci95)
log_lci<-log(or_df$or_lci95)
se_log_or <-(log_uci-log_lci)/(2*1.96)


meta_df <- data.frame(id= paste0 (or_df$id.exposure,"-",or_df$id.outcome),
                      or= or_df$or,
                      or_uci95=or_df$or_uci95,
                      or_lci95=or_df$or_lci95
                      )

meta_analysis <- meta:: metagen(log_or,se_log_or,sm="OR",data=meta_df, studlab=meta_df$id)

summary(meta_analysis)


while (!is.null(dev.list()))   dev.off()
png(filename = paste0("meta_plot.png"),
    width=1600, height=1400,units = "px",res = 150)
par(mai = c(0.5, 0.5, 0.5, 0.5)) # 边距 c(bottom, left, top, right)

meta::forest(meta_analysis,plotwidth = "10cm",fontsize = 8)

while (!is.null(dev.list()))   dev.off()

