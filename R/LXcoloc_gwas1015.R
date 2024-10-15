
LXcoloc_gwas <- function(eqtlID=eqtlID,
                       type1=type1,
                       chr_eqtl=chr_eqtl,
                       pos_eqtl=pos_eqtl,
                       SNP_eqtl=SNP_eqtl,
                       beta_eqtl=beta_eqtl,
                       pval_eqtl=pval_eqtl,
                       se_eqtl=se_eqtl,
                       eaf_eqtl=eaf_eqtl,
                       samplesize_eqtl_file=samplesize_eqtl_file,
                       outcomeID =outcomeID,
                       type2 = type2,
                       chr_gwas=chr_gwas,
                       pos_gwas=pos_gwas,
                       SNP_gwas=SNP_gwas,
                       beta_gwas=beta_gwas,
                       pval_gwas=pval_gwas,
                       se_gwas=se_gwas,
                       eaf_gwas=eaf_gwas,
                       samplesize_gwas= samplesize_gwas,
                       target_gene_file =target_gene_file,
                       SNP_PP_H4=SNP_PP_H4){

#----R packages---------------
  LXcoloc_inst_packages()

#---create dir----------------
  dir_file= paste0("coloc analysis results_",Sys.Date())

  if(!dir.exists(dir_file))
      dir.create(dir_file)

#------data information-----------------------------------
exp_list <- dir(eqtlID,full.names = T) %>% as.character()
out_list <- outcomeID %>% as.character()

if(!is.na(target_gene_file)){
  target_list <- read.xlsx(target_gene_file)
}

#------sample information---------------------------------
sam_eqtl <- read.xlsx(samplesize_eqtl_file) # 暴露样本数量
sam_gwas <- read.xlsx(samplesize_gwas) # 结局样本数量

#------foreach 01--------------------------------------------

x=1

foreach(x=c(1:length(exp_list)), .errorhandling = "pass") %do% {

print(paste0("Number ",x, " of ",length(exp_list)," exposure data is being processed..."))

file_type <- tools::file_ext(exp_list[x]) # 文件扩展名

eqtl_df <- fread(exp_list[x]) %>%
           dplyr::rename(chr.exposure=chr_eqtl,
                         pos.exposure=pos_eqtl,
                         SNP=SNP_eqtl,
                         beta.exposure=beta_eqtl,
                         se.exposure=se_eqtl,
                         pval.exposure=pval_eqtl,
                         eaf.exposure=eaf_eqtl
                         )

# head(eqtl_df)

eqtl_df$pval.exposure <- as.numeric(eqtl_df$pval.exposure)  # 转换成数据格式
eqtl_df$beta.exposure <- as.numeric(eqtl_df$beta.exposure)
eqtl_df$eaf.exposure <- as.numeric(eqtl_df$eaf.exposure)
eqtl_df$pos.exposure <- as.numeric(eqtl_df$pos.exposure)

id01 <- str_extract(exp_list[x],"(?<=\\/)[^/]+$") %>% str_extract(".*?(?=\\.)")

if(grepl("_",id01))
  id_exp <- str_extract(id01,".*?(?=_)") else
    id_exp <- id01

idx <- grep(id_exp,sam_eqtl[,1])

eqtl_df$samplesize.exposure <- sam_eqtl$samplesize[idx] #添加例数

eqtl_df <- eqtl_df[!is.na(eqtl_df$beta.exposure),]%>% .[!duplicated(.$SNP),]

#------target gene -------------------

if(!is.na(target_gene_file))
target_n <- nrow(target_list) else
  target_n <- 1

y=1

foreach(y=c(1:target_n), .errorhandling = "pass") %do% {

  message(paste0("Number ",y, " of ",target_n," target genes is being processed..."))

if(!is.na(target_gene_file)){
   target_gene <- target_list[y,1]
   chr_target = target_list[y,2]
   pos_target_start = target_list[y,3]
   pos_target_end = target_list[y,4]

#---共定位区域------------------------
   chrpos <- paste0(chr_target,":",pos_target_start-100000,"-",pos_target_end+100000)
   geneChr <- chr_target
   geneStart <- pos_target_start
   geneEnd <- pos_target_end

   } else {
      top_snp <- eqtl_df %>% dplyr::arrange(pval.exposure) %>% .[1,]
      chrpos <- paste0(top_snp$chr.exposure,":",top_snp$pos.exposure-100000,"-",top_snp$pos.exposure+100000)
      geneChr <- top_snp$chr.exposure
      geneStart <- top_snp$pos.exposure
      geneEnd <- top_snp$pos.exposure
  } #--if(!is.na(target_gene_file)) end-----

#---eqtl data 取共定位区域----
 eqtl_coloc <- eqtl_df %>%
    subset(chr.exposure==geneChr & pos.exposure>geneStart-100000 & pos.exposure<geneEnd+100000)

#----清除内存中的eqtl_df------
  rm(eqtl_df)

#------foreach outcome gwas --------------------------------------------
z=1

foreach(z=c(1:length(out_list)), .errorhandling = "pass") %do% {

   message(paste0("Number ",z, " of ",length(out_list)," outcome data is being processed..."))

   gwas_df <- fread(out_list[z])

   gwas_df <- gwas_df %>% dplyr::rename(chr.outcome=chr_gwas,
                                        pos.outcome=pos_gwas,
                                        SNP=SNP_gwas,
                                        beta.outcome=beta_gwas,
                                        pval.outcome=pval_gwas,
                                        se.outcome=se_gwas,
                                        eaf.outcome=eaf_gwas
                                       )

   # head(gwas_df)

   gwas_df$se.outcome <- as.numeric(gwas_df$se.outcome)
   gwas_df$pval.outcome <- as.numeric(gwas_df$pval.outcome)
   gwas_df$beta.outcome <- as.numeric(gwas_df$beta.outcome)
   gwas_df$eaf.outcome <- as.numeric(gwas_df$eaf.outcome)
   gwas_df$pos.outcome <- as.numeric(gwas_df$pos.outcome)

#-----------添加结局样本数量----------------------------------------#
   id02 <- str_extract(out_list[z],"(?<=\\/)[^/]+$") %>% str_extract(".*?(?=\\.)")

   if(grepl("_",id02))
     id_gwas <- str_extract(id02,".*?(?=_)") else
       id_gwas <- id02

   idx_gwas <- grep(id_gwas,sam_gwas[,1])

   gwas_df$samplesize.outcome <- sam_gwas$samplesize[idx_gwas] #添加例数

#------去掉beta NA和重复的SNP------------------------------------------------#
   gwas_df <- gwas_df[!is.na(gwas_df$beta.outcome),]%>% .[!duplicated(.$SNP),]

#---outcome data 取共定位区域----
   gwas_coloc <- gwas_df %>%
     subset(chr.outcome==geneChr & pos.outcome>geneStart-100000 & pos.outcome<geneEnd+100000)

#----清空内存中的gwas_df--------
   rm(gwas_df)


#---整合eqtl和gwas数据-------
   overlapping_snps <- intersect(eqtl_coloc$SNP,gwas_coloc$SNP)

   merg_data <- merge(eqtl_coloc,gwas_coloc,by="SNP")
   if(nrow(merg_data)==0)
     if(nontarget)
       stop("There is no any overlapping SNP between the eqtl and gwas data.") else
         stop("There is no any overlapping SNP between the eqtl and gwas data. Thus, the target gene may be unsuitable.")

   dataset1 <- list(pvalues=merg_data$pval.exposure,
                    N=merg_data$samplesize.exposure,
                    MAF= ifelse(merg_data$eaf.exposure>0.5,1-merg_data$eaf.exposure,merg_data$eaf.exposure),
                    beta=merg_data$beta.exposure,
                    varbeta= merg_data$se.exposure^2,
                    type=type1,
                    snp=merg_data$SNP,
                    z= merg_data$beta.exposure/merg_data$se.exposure,
                    chr=merg_data$chr.exposure,
                    pos=merg_data$pos.exposure,
                    id=str_extract(dir(eqtlID)[x],".*(?=_)"))


   dataset2 <- list(pvalues=merg_data$pval.outcome,
                    N=merg_data$samplesize.outcome,
                    MAF= ifelse(merg_data$eaf.outcome>0.5,1-merg_data$eaf.outcome,merg_data$eaf.outcome),
                    beta=merg_data$beta.outcome,
                    varbeta= merg_data$se.outcome^2,
                    type=type2,
                    snp=merg_data$SNP,
                    z= merg_data$beta.outcome/merg_data$se.outcome,
                    chr=merg_data$chr.outcome,
                    pos=merg_data$pos.outcome,
                    id=str_extract(outcomeID,".*(?=\\.)"))

   #------共定位分析--------------------------------------------------------
   col_result <- coloc::coloc.abf(dataset1, dataset2)
   col_data <- col_result$results %>% dplyr::arrange(-SNP.PP.H4)  # -SNP.PP.H4降序

   #筛选共定位的位点
   col_snp <- col_data %>% filter(SNP.PP.H4 > SNP_PP_H4)
   lead_snp <- col_snp$snp[1]

   if(is.na(lead_snp))
     print(paste0("The 'SNP.PP.H4' is ",round(col_data$SNP.PP.H4[1],3),
                  ", suggesting that there is no any colocalization snp.")) else
                    print(paste0("The colocalization snp may be '",lead_snp,
                                 "' with a 'SNP.PP.H4' value of ",round(col_snp$SNP.PP.H4,4)))


   id_out <- id_gwas

   if(round(col_data$SNP.PP.H4[1],3)>SNP_PP_H4){
      if(!is.na(target_gene_file))
         snp_dir=paste0(dir_file,"/",target_gene,"-",id_exp,"-",id_out," (SNP.PP.H4 ",round(col_data$SNP.PP.H4[1],3),")") else
           snp_dir=paste0(dir_file,"/",id_exp,"-",id_out," (SNP.PP.H4 ",round(col_data$SNP.PP.H4[1],3),")")

      } else{
         if(!is.na(target_gene_file))
          snp_dir=paste0(dir_file,"/",target_gene,"-",id_exp,"-",id_out) else
            snp_dir=paste0(dir_file,"/",id_exp,"-",id_out)
          }


   if(!dir.exists(snp_dir))
     dir.create(snp_dir)

   if(!is.na(target_gene_file))
     file=paste0(snp_dir,"/",target_gene, "-(",id_exp,"-",id_out,") ",Sys.Date(),".csv") else
        file=paste0(snp_dir,"/",id_exp,"-",id_out," ",Sys.Date(),".csv")

   write.csv(col_data, file=file, row.names=F)

   #------locuscompare可视化-------------------
   fn1_gwas <- data.frame(rsid=merg_data$SNP,pval=merg_data$pval.outcome)
   fn2_eqtl <- data.frame(rsid=merg_data$SNP,pval=merg_data$pval.exposure)

   coloc_plot <- locuscompare(in_fn1 = fn1_gwas, in_fn2 = fn2_eqtl,
                              title1 = id_out, title2 = id_exp)

   coloc_plot

   if(!is.na(target_gene_file))
     plot01_file=paste0(snp_dir,"/locuscompare_plot (",target_gene,"-",id_exp,"-",id_out,") ",Sys.Date(),".png") else
       plot01_file=paste0(snp_dir,"/locuscompare_plot (",id_exp,"-",id_out,") ",Sys.Date(),".png")

   ggsave(plot01_file,coloc_plot, width = 1000,height = 800,units = "px",dpi = 150)

   #---------geni.plots可视化--------------------------------------------------
   snp_list=merg_data$SNP
   table(duplicated(snp_list))

   print("Begin to calculate the correlation coefficient,it may take a long time.......")

   corr_list <-LXcoloc_ld_snplist (variants=snp_list, with_alleles = TRUE, pop = "EUR",
                                 opengwas_jwt = get_opengwas_jwt()) %>% data.frame()

   if(ncol(corr_list)<2 | nrow(corr_list)<2)
     stop(paste0("There is no enough corr number between the eqtl and gwas data.
       Maybe the chrompos is inappropriate or the target gene is inappropriate."))

   print("Calculating the correlation coefficient is done.")
   corr_snp <- rownames(corr_list) %>% str_extract(.,".*?(?=_)")
   colnames(corr_list) <- corr_snp
   rownames(corr_list) <- corr_snp
   corr_list$snp <- corr_snp
   corr_list <- corr_list %>% dplyr::arrange(snp)
   corr_list <- corr_list[,-ncol(corr_list)]

   set_01 <- data.frame(marker=dataset1$snp,pvalue_1=dataset1$pvalues) %>%
     subset(marker %in% corr_snp)

   set_02 <- data.frame(marker=dataset2$snp,pvalue_2=dataset2$pvalues) %>%
     subset(marker %in% corr_snp)

   mark=data.frame(marker = dataset1$snp,chr = dataset1$chr, pos = dataset1$pos)

   assoc <- Reduce(function(x,y)inner_join(x,y,by="marker"),list(mark,set_01,set_02)) %>%
     data.frame() %>% dplyr::arrange(marker)

   traits=c(dataset1$id[1], str_extract(dataset2$id[1],"(?<=/)[^/]+$"))

   if(is.na(lead_snp))
     highlights= col_data$snp[1] else
       highlights = lead_snp

   library(geni.plots)
   assoc$pos <-as.numeric(assoc$pos) %>% as.integer()

   plot_geni <- geni.plots::fig_region_stack(data = assoc, corr = corr_list,
                                             traits = traits,highlights=highlights,
                                             build = 37,title_center = TRUE)
   plot_geni


   if(!is.na(target_gene_file))
     plot02_file=paste0(snp_dir,"/geni_plot (",target_gene,"-",id_exp,"-",id_out,") ",Sys.Date(),".png") else
       plot02_file=paste0(snp_dir,"/geni_plot (",id_exp,"-",id_out,") ",Sys.Date(),".png")

   ggsave(plot02_file,plot_geni,width = 1000,height = 1000,units = "px",dpi = 150)


#------------------------------------------

   message(paste0("Number ",z, " of ",length(outcomeID),"  outcome data was completed"))

    } # z foreach outcome gwas end

#------------------------------------------

  message(paste0("Number ",y, " of ",nrow(target_list)," target genes was completed."))

 }#---- y foreach target gene end---

  print(paste0("Number ",x, " of ",length(exp_list)," exposure data was done."))

#------------------------------------------

}#---- x foreach exposrue end---

print( paste0("The colocalization analysis  result can be found in the folder of '",dir_file,"'."))

}#----------the end---------------------------------------------
