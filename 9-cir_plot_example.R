
library(LXbatmr)

# devtools::load_all()

rm(list=ls())

or_file = "D:/Desktop/LXbatmr 2024-10-24/cicos_data/significant_drugMR_mr.xlsx" # OR 表

mr_method=c("Inverse variance weighted","MR Egger","Weighted mode","Weighted median","Simple mode")   #选择展示的方法

# mr_method=c("Inverse variance weighted")   #选择展示的方法

n_meth = length(mr_method)

cir_txt_sit = 0.8+(n_meth:1)*0.40# methods 的位置和行距

#----------------------------------------------------------#
cir_pic (or_file,mr_method,cir_txt_sit,
         gap.degree=35,   # 环形热图开口大小
         start.degree=350, # 开口位置
         track.height=0.2) # 扇形高度


