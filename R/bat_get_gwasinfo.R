
#-------------------------------------------
bat_get_gwasinfo <- function(data_dir,data_list){

inst_packages()

if(!is.na(data_list)) {
  id_list <- read.xlsx(data_list) %>% data.frame()
}


if(!is.na(data_dir)){
gwas_data_info <- data.frame()

file_list <- dir(data_dir) %>%
            str_extract(".*?(?=\\.)")

file_table <- data.frame()

for(i in file_list){
  if(grepl("_",i))
    x=str_extract(i,".*?(?=_)") else
      x=i
  file_table <- rbind(file_table,cbind(id=x))
}

id_list <- file_table %>% unique()
}


#-------foreach analysis----------------------
foreach(x=c(1:nrow(id_list)), .errorhandling = "pass") %do% {

gwas_info <- get_studies(study_id =id_list[x,1])

if(grepl("level",gwas_info@studies$reported_trait,ignore.case = T))
 info <- str_extract(gwas_info@studies$reported_trait,".*(?= )") else
   info <- gwas_info@studies$reported_trait

id_list[x,2] <-info

#samplesize <- str_extract(gwas_info@studies$initial_sample_size,"\\d.*\\d")
samplesize <- gwas_info@ancestries$number_of_individuals
#pop <- str_extract(gwas_info@studies$initial_sample_size,"(?<=\\ ).*?(?=\\ )")
pop <- gwas_info@ancestral_groups$ancestral_group
nsnp=gwas_info@studies$snp_count
pubmed_id <- gwas_info@publications$pubmed_id

m_info <- data.frame(id=id_list[x,1],
                     name=info,
                     samplesize=samplesize,
                     pop=pop,
                     nsnp=nsnp,
                     pubmed_id)

gwas_data_info <- rbind(gwas_data_info,m_info)


#-----foreach end-------------
}

if(!is.na(data_dir)){
  write.xlsx(gwas_data_info,paste0(data_dir,"/gwas_data_info.xlsx"))
} else{
write.xlsx(id_list,"bat_mr_result_significant (metabolite_info).xlsx")
write.xlsx(gwas_data_info,"gwas_data_info.xlsx")
}

#----the end----------
}


#bat_get_gwasinfo(data_dir,data_list)


