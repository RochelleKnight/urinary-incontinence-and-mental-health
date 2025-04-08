library(data.table)
library(stringr)
library(dplyr)

# Read in allele frequencies
allele_freq <- fread(file="data/raw_data/allele_freq_1000G_10092023.txt", fill = T)
allele_freq$chr_position <- paste0("chr",allele_freq$CHR,":",allele_freq$position)

for(outcome_name in c("moderate_ANY_UI","moderate_UUI","moderate_SUI","severe_ANY_UI","severe_UUI","severe_SUI")){
  print(paste0("Working on ",outcome_name ))
  
  #Read in UI GWAS data
  ui_gwas <- fread(file=paste0("data/raw_data/original_UI_gwas_results/",outcome_name,"_GWAS.tbl"), select=c("MarkerName","Allele1","Allele2","Effect","StdErr","P-value"), fill = T)
  
  print(unique(substr(ui_gwas$MarkerName,1,2)))
  
  ui_gwas_rsid <- ui_gwas %>% filter(str_detect(MarkerName, "^rs")) %>% 
    left_join(allele_freq, by = c("MarkerName"="RSID"))
  
  ui_gwas_chr_pos <- ui_gwas %>% filter(str_detect(MarkerName, "^chr")) %>% 
    left_join(allele_freq, by = c("MarkerName"="chr_position"))
  
  ui_gwas_chr_pos$MarkerName[which(str_detect(ui_gwas_chr_pos$RSID,"^rs"))] <- ui_gwas_chr_pos$RSID[which(str_detect(ui_gwas_chr_pos$RSID,"^rs"))]
  
  ui_gwas_chr_pos <- ui_gwas_chr_pos %>% filter(!MarkerName %in% ui_gwas_rsid$MarkerName)
  
  ui_gwas_rsid <- ui_gwas_rsid %>% select(intersect(colnames(ui_gwas_chr_pos),colnames(ui_gwas_rsid)))
  ui_gwas_chr_pos <- ui_gwas_chr_pos %>% select(intersect(colnames(ui_gwas_chr_pos),colnames(ui_gwas_rsid)))
  
  ui_gwas <- rbind(ui_gwas_rsid,ui_gwas_chr_pos)
  rm(ui_gwas_rsid,ui_gwas_chr_pos)
  
  ui_gwas$Allele1 <- str_to_upper(ui_gwas$Allele1)
  ui_gwas$Allele2 <- str_to_upper(ui_gwas$Allele2)
  ui_gwas$issue_with_alleles <- NA
  
  ui_gwas$MAF <- as.numeric(ui_gwas$MAF)
  ui_gwas$eaf <- ifelse(ui_gwas$Allele1 == ui_gwas$minor_allele, ui_gwas$MAF, 1-ui_gwas$MAF)
  ui_gwas$issue_with_alleles <- ifelse((ui_gwas$Allele1 == ui_gwas$minor_allele & ui_gwas$Allele2 == ui_gwas$major_allele) |
                                         (ui_gwas$Allele2 == ui_gwas$minor_allele & ui_gwas$Allele1 == ui_gwas$major_allele),
                                       "no","yes")
  
  ui_gwas <- ui_gwas %>% rename("A_allele_from_ALSPAC" = "A_allele",
                                "B_allele_from_ALSPAC" = "B_allele",
                                "minor_allele_from_ALSPAC" = "minor_allele",
                                "major_allele_from_ALSPAC" = "major_allele") %>%
    select("MarkerName","Allele1","Allele2","Effect","StdErr","P-value","MAF","eaf","A_allele_from_ALSPAC","B_allele_from_ALSPAC","minor_allele_from_ALSPAC","major_allele_from_ALSPAC","issue_with_alleles")
  fwrite(ui_gwas, paste0("data/raw_data/",outcome_name,".csv.gz"),compress = "gzip")
}
  
  
  

