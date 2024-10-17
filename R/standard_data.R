
standard_data <- function(gwas_file=gwas_file,
                          chr=chr,
                          pos=pos,
                          SNP=SNP,
                          beta=beta,
                          pval=pval,
                          se=se,
                          effect_allele=effect_allele,
                          other_allele=other_allele,
                          eaf=eaf,
                          clum_p=clum_p){

#------R package------
inst_packages()

#------read data------
x=1
foreach(x=c(1:length(gwas_file)), .errorhandling = "pass") %do% {

  print(paste0("It is number ",x, " of ",length(gwas_file)))

  file_path <- gwas_file[x]

  file_name <- str_extract(file_path,"(?<=\\/)[^/]+$")

  if(grepl("vcf.gz",file_name,ignore.case = TRUE)){
    df <- VariantAnnotation::readVcf(file_path) %>%  # 读取vcf数据
      gwasglue::gwasvcf_to_TwoSampleMR(type = "exposure")

    } else {
          if(grepl("vcf",file_name,ignore.case = TRUE)){
             if(!file.exists("variants.tsv.bgz")){
               ukbb_url <- "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz"
               download.file(ukbb_url,destfile = "variants.tsv.bgz",method = "libcurl" )}
               ukbb_anno <- fread("variants.tsv.bgz")
               df <- fread(file_path)
               df <- inner_join(ukbb_anno,exp_df,by = "variant")
               rm(ukbb_anno)
                  } else
                    df <- fread(file_path)
       }

  #------standard data------
  colnames(df)[grep("chr",colnames(df),ignore.case = TRUE)] <- "chromosome"

  colnames(df)[which(colnames(df)==chr)] <- "chromosome"

  colnames(df)[which(colnames(df)==pos)] <- "base_pair_location"
  colnames(df)[which(colnames(df)==SNP)] <- "variant_id"

  colnames(df)[which(colnames(df)==beta)] <- "beta"

  # beta = log(OR), 如果只提供OR，需要转换成log(OR)，才是beta
  if(!grepl("Log", beta,ignore.case = T)) {
      if(grepl("odd",beta,ignore.case = T) | grepl("OR",beta,ignore.case = T)){
        df$beta <- log(as.numeric(df$beta))
      }
    }

  colnames(df)[which(colnames(df)==pval)] <- "p_value"
  colnames(df)[which(colnames(df)==se)] <- "standard_error"

  colnames(df)[which(colnames(df)==effect_allele)] <- "effect_allele"
  df$effect_allele <- toupper(df$effect_allele)

  colnames(df)[which(colnames(df)==other_allele)] <- "other_allele"
  df$other_allele <- toupper(df$other_allele)

  colnames(df)[which(colnames(df)==eaf)] <- "effect_allele_frequency"


  #--------
  clump_df <- subset(df,p_value<clum_p) #筛选p值

  fwrite(clump_df, paste0(str_extract(file_path,".*(?=\\.)"),"_ (pval less than ",clum_p,").csv"))

  fwrite(df, paste0(str_extract(file_path,".*(?=\\.)"),"_standard.gz"))

  print (paste0(" The ",x, " of ",length(gwas_file), " was completed "))

} # foreach end

print(paste0("The standarded data was saved in '",str_extract(gwas_file[1],".*(?=\\/)"),"'"))

} # standard_data end



