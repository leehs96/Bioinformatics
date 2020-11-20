library(data.table)
library(utils)
library(vcfR)

getwd()
setwd("~/cancer/WGS_hg37/variant_calling/varscan_results")


path <- "/users/hslee/cancer/WGS_hg37/variant_calling/filtered/ADB_LL//"
vcf <-list.files(path = path, pattern = '*.vcf', full.names = T)

#####################################################################
for (i in seq_along(vcf)){
  a <- gsub(path,"",vcf[i])
  b <- gsub('/','',a)
  name <- (gsub(".pass.vcf","",b))
  assign(name, fread(vcf[i], skip = "CHROM"))
}#######################################
#####################################################################
strelka2_snp_filter <- function(mylist_1){

  #mylist_1 <- mylist_1[!grepl('##', mylist_1$V1)]
  #mylist_1$V1[1] <- 'CHROM'
  #header <- as.character(mylist_1[1])
  #names(mylist_1) <- header
  #mylist_1 <- mylist_1[-1,]
  
  df <- data.frame()

  for (i in 1:nrow(mylist_1)){
      ref = mylist_1$REF[i]
      alt = mylist_1$ALT[i]
      normal = strsplit(mylist_1$NORMAL[i],':')
      tumor = strsplit(mylist_1$TUMOR[i],':')
      ifelse(alt == 'A', germ_var <- strsplit(normal[[1]][5],',')[[1]][1],
             ifelse(alt == 'C', germ_var <- strsplit(normal[[1]][6],',')[[1]][1],
                    ifelse(alt == 'G', germ_var <- strsplit(normal[[1]][7],',')[[1]][1],
                           ifelse(alt == 'T', germ_var <- strsplit(normal[[1]][8],',')[[1]][1],break
                                  ))))
      ifelse(alt == 'A', soma_var <- strsplit(tumor[[1]][5],',')[[1]][1],
             ifelse(alt == 'C', soma_var <- strsplit(tumor[[1]][6],',')[[1]][1],
                    ifelse(alt == 'G', soma_var <- strsplit(tumor[[1]][7],',')[[1]][1],
                           ifelse(alt == 'T', soma_var <- strsplit(tumor[[1]][8],',')[[1]][1],break
                                  ))))
      ifelse(ref == 'A', germ_ref_read <- strsplit(normal[[1]][5],',')[[1]][1],
             ifelse(ref == 'C', germ_ref_read <- strsplit(normal[[1]][6],',')[[1]][1],
                    ifelse(ref == 'G', germ_ref_read <- strsplit(normal[[1]][7],',')[[1]][1],
                           ifelse(ref == 'T', germ_ref_read <- strsplit(normal[[1]][8],',')[[1]][1],break
                           ))))
      
      ifelse(ref == 'A', soma_ref_read <- strsplit(tumor[[1]][5],',')[[1]][1],
             ifelse(ref == 'C', soma_ref_read <- strsplit(tumor[[1]][6],',')[[1]][1],
                    ifelse(ref == 'G', soma_ref_read <- strsplit(tumor[[1]][7],',')[[1]][1],
                           ifelse(ref == 'T', soma_ref_read <- strsplit(tumor[[1]][8],',')[[1]][1],break
                           ))))
      
      germ_vaf = as.numeric(germ_var)/as.numeric(germ_var)+as.numeric(germ_ref_read)
      if(is.nan(germ_vaf)){
        germ_vaf = 0
      }
      soma_vaf = as.numeric(soma_var)/as.numeric(soma_var)+as.numeric(soma_ref_read)
      if(is.nan(soma_vaf)){
        soma_vaf = 0
      }
      
      
      if(as.numeric(soma_var) > 2 && 
         as.numeric(germ_var) <2 && 
         as.numeric(soma_vaf) > 20){
       df <- rbind(df, mylist_1[i,])
      } 
      
  }
  return(df)
}#######################
#####################################################################
varscan2_snp_filter <- function(mylist_2){
  #mylist_2 <- mylist_2[!grepl('##', mylist_2$V1)]
  #mylist_2$V1[1] <- 'CHROM'
  #header <- as.character(mylist_2[1])
  #names(mylist_2) <- header
  #mylist_2 <- mylist_2[-1,]
  
  df <- data.frame()
  for (i in 1:nrow(mylist_2)){
    
    input_germ <- strsplit(mylist_2$NORMAL[i],':')[[1]][1]
    input_somatic <- strsplit(mylist_2$TUMOR[i],':')[[1]][1]
    
    germ_var <- strsplit(mylist_2$NORMAL[i],':')[[1]][5]
    soma_var <- strsplit(mylist_2$TUMOR[i],':')[[1]][5]
    
    germ_vaf <- gsub('%','',strsplit(mylist_2$NORMAL[i],':')[[1]][6])
    soma_vaf <- gsub('%','',strsplit(mylist_2$TUMOR[i],':')[[1]][6])

    ifelse(input_germ == '0/0' &&
           input_somatic == '0/1' &&
           as.numeric(germ_var) < 2 && 
           as.numeric(soma_var) >= 3,
           df <- rbind(df, mylist_2[i,]),
           ifelse(input_germ == '0/0' &&
                    input_somatic == '1/1' &&
                    as.numeric(germ_var) < 2 && 
                    as.numeric(soma_var) >= 3,
                  df <- rbind(df, mylist_2[i,]),
                  ifelse(input_germ == '0/1' &&
                           input_somatic == '0/0' &&
                           as.numeric(germ_vaf) >= 25 && 
                           as.numeric(soma_var) < 3,
                         df <- rbind(df, mylist_2[i,]),
                         ifelse(input_germ == '0/1' &&
                                  input_somatic == '1/1' && 
                                  as.numeric(germ_vaf) >= 25,
                                df <- rbind(df, mylist_2[i,]),
                                ifelse(input_germ == '1/1' &&
                                         input_somatic == '0/1' && 
                                         as.numeric(soma_vaf) >= 25,
                                       df <- rbind(df, mylist_2[i,]),
                                       ifelse(input_germ == '1/1' &&
                                                input_somatic == '0/0' &&
                                                as.numeric(soma_var) < 2,
                                              df <- rbind(df, mylist_2[i,]),NA
                                              )
                                       )
                                )
                         )
                   )
              )
  } 
  return(df)
}#######################
#####################################################################
strelka2_indel_filter <- function(df){
  out <- data.frame()
  
  for (i in 1:nrow(df)){
    ref = df$REF[i]
    alt = 
  }
  
}



strelka_snp_out <- strelka2_snp_filter(ADB_LL.strelka.somatic.snvs)
varscan_snp_out <- varscan2_snp_filter(ADB_LL.varscan.snp)



read.table()