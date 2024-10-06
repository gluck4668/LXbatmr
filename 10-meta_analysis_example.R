
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXbatmr")

library(LXbatmr)

rm(list=ls())
gc()


or_dir = "D:/Desktop/LXbatmr 2024-9-29(v3.21)/OR"

# devtools::load_all()

#---------------------------------------------------------
meta_analysis (or_dir)
