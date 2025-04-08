library(dplyr)

snps_used <- tidyr::crossing(c("urinary_incontinence","stress_urinary_incontinence","urgency_urinary_incontinence","mixed_urinary_incontinence"),c("MDD","broad_depression_phenotype","neuroticism","anxiety"))
colnames(snps_used) <- c("exposure","outcome")
tmp1 <- tidyr::crossing(c("anxiety"),c("urinary_incontinence","stress_urinary_incontinence","urgency_urinary_incontinence","mixed_urinary_incontinence"))
colnames(tmp1) <- c("exposure","outcome")
snps_used <- rbind(snps_used,tmp1)

# Number of SNPs in exposure GWAS
exposure_snps <- as.data.frame(matrix(nrow = 0, ncol = 3))
for (exposure in c("urinary_incontinence","stress_urinary_incontinence","urgency_urinary_incontinence","mixed_urinary_incontinence",
                   "anxiety")){
  if(file.exists(paste0("data/exposure/",exposure,"_p_value_5e-06_clumped.csv"))){
    df <- read.csv(paste0("data/exposure/",exposure,"_p_value_5e-06_clumped.csv"))
    df <- df[,c("rsid","phenotype","proxy")]
    exposure_snps <- rbind(exposure_snps,df)
  }
}

exposure_snps <- exposure_snps %>% 
  group_by(phenotype,proxy) %>%
  summarise(count = n())

exposure_snps$phenotype <- gsub("_p_value_5e-06","",exposure_snps$phenotype)
exposure_snps$exposure <- exposure_snps$phenotype
exposure_snps[,c("proxy","phenotype")] <- NULL

snps_used <- snps_used %>% left_join(exposure_snps)
snps_used <- snps_used %>% rename("Number of SNPs in exposure GWAS"="count")

#Number of exposure SNPs in outcome
exposure_snps <- as.data.frame(matrix(nrow = 0, ncol = 3))
for (exposure in c("urinary_incontinence","stress_urinary_incontinence","urgency_urinary_incontinence","mixed_urinary_incontinence",
                   "anxiety")){
  if(file.exists(paste0("data/exposure/",exposure,"_p_value_5e-06_clumped.csv"))){
    df <- read.csv(paste0("data/exposure/",exposure,"_p_value_5e-06_clumped.csv"))
    df <- df[,c("rsid","phenotype","proxy")]
    exposure_snps <- rbind(exposure_snps,df)
  }
}

exposure_snps$phenotype <- gsub("_p_value_5e-06","",exposure_snps$phenotype)
exposure_snps <- exposure_snps %>% rename("exposure"="phenotype",
                                          "SNP"="rsid")
exposure_snps$proxy <- NULL

outcome_snps <- as.data.frame(matrix(nrow = 0, ncol = 3))
for (exposure in c("urinary_incontinence","mixed_urinary_incontinence","stress_urinary_incontinence","urgency_urinary_incontinence",
                   "MDD","broad_depression_phenotype","neuroticism","anxiety")) {
  
  if(file.exists(paste0("data/outcome/outcome_data_for_",exposure,"_p_value_5e-06.csv"))){
    df <- read.csv(paste0("data/outcome/outcome_data_for_",exposure,"_p_value_5e-06.csv"))
    df <- df[,c("SNP","phenotype")]
    df$exposure <- exposure
    df$outcome <- df$phenotype
    df$phenotype <- NULL
    outcome_snps <- rbind(outcome_snps,df)
  }
}

tmp<-outcome_snps %>% filter(exposure == "mixed_urinary_incontinence" & outcome == "broad_depression_phenotype")
outcome_snps[duplicated(outcome_snps),]
tmp <- exposure_snps %>% inner_join(outcome_snps)
tmp <- tmp %>% 
  group_by(exposure,outcome) %>%
  summarise(count = n())

snps_used <- snps_used %>% left_join(tmp)
snps_used <- snps_used %>% rename("Number of exposure SNPs in outcome GWAS"="count")

# Number of exposure SNPs missing from outcome GWAS 
snps_used$"Number of exposure SNPs missing from outcome GWAS" <- snps_used$`Number of SNPs in exposure GWAS` - snps_used$`Number of exposure SNPs in outcome GWAS`

#Number of missing SNPs there are proxies for 
proxies <- as.data.frame(matrix(nrow = 0, ncol = 4))

for (exposure in c("urinary_incontinence","mixed_urinary_incontinence","stress_urinary_incontinence","urgency_urinary_incontinence")) {
  for (outcome in c("MDD","broad_depression_phenotype","neuroticism","anxiety")) {
    tmp <- read.csv(paste0("data/exposure/exposure_",exposure,"_p_value_5e-06_outcome_",outcome,"_with_proxies.csv"))
    tmp <- tmp[,c("rsid","phenotype","outcome","proxy")]
    tmp <- tmp %>% rename("exposure"="phenotype",
                          "SNP"="rsid")
    proxies <- rbind(proxies,tmp)
  }
  
}

for (exposure in c("anxiety")) {
  for (outcome in c("urinary_incontinence","mixed_urinary_incontinence","stress_urinary_incontinence","urgency_urinary_incontinence")) {
    tmp <- read.csv(paste0("data/exposure/exposure_",exposure,"_p_value_5e-06_outcome_",outcome,"_with_proxies.csv"))
    tmp <- tmp[,c("rsid","phenotype","outcome","proxy")]
    tmp <- tmp %>% rename("exposure"="phenotype",
                          "SNP"="rsid")
    proxies <- rbind(proxies,tmp)
  }
}

proxies$exposure <- gsub("_p_value_5e-06","",proxies$exposure)
proxies <- proxies %>% filter(proxy == "proxy_snp") %>% select(-proxy)

tmp <- outcome_snps %>% inner_join(proxies) 
tmp <- tmp %>% 
  group_by(exposure,outcome) %>%
  summarise(count = n())

snps_used <- snps_used %>% left_join(tmp)
snps_used$count <- ifelse(is.na(snps_used$count) & snps_used$`Number of exposure SNPs missing from outcome GWAS` !=0, 0,snps_used$count)
snps_used <- snps_used %>% rename("Number of missing SNPs there are proxies for"="count")                      


snps_used$"Total non-proxy and proxy SNPs used as exposure" <- snps_used$`Number of exposure SNPs in outcome GWAS` + snps_used$`Number of missing SNPs there are proxies for`

# Total SNPs remaining following harmonisation

# Load UI MR results
mr_results_urinary_incontinence <- read.csv(paste0("results/UI_MR/urinary_incontinence_p_value_5e-06.csv"))
mr_results_mixed_urinary_incontinence <- read.csv(paste0("results/UI_MR/mixed_urinary_incontinence_p_value_5e-06.csv"))
mr_results_stress_urinary_incontinence <- read.csv(paste0("results/UI_MR/stress_urinary_incontinence_p_value_5e-06.csv"))
mr_results_urgency_urinary_incontinence <- read.csv(paste0("results/UI_MR/urgency_urinary_incontinence_p_value_5e-06.csv"))

# Load MH MR results
mr_results_anxiety <- read.csv(paste0("results/MH_MR/anxiety_p_value_5e-06.csv"))

mr_results <- rbind(mr_results_urinary_incontinence,mr_results_mixed_urinary_incontinence,
                    mr_results_stress_urinary_incontinence,mr_results_urgency_urinary_incontinence,mr_results_anxiety)

rm(mr_results_urinary_incontinence,mr_results_mixed_urinary_incontinence,
   mr_results_stress_urinary_incontinence,mr_results_urgency_urinary_incontinence,mr_results_anxiety)

mr_results <- mr_results %>% filter(method == "Debiased IVW") %>%
  select(exposure,outcome,nsnp)
mr_results$exposure <- gsub("_p_value_5e-06","",mr_results$exposure)
mr_results <- mr_results[!duplicated(mr_results),]

snps_used <- snps_used %>% left_join(mr_results)
snps_used <- snps_used %>% rename("Total SNPs remaining following harmonisation"="nsnp")

# Tidy exposure and outcome variable names
snps_used <- snps_used %>% dplyr::mutate(exposure = dplyr::case_when(exposure == "urinary_incontinence" ~ "Urinary incontinence",
                                                                       exposure == "stress_urinary_incontinence" ~ "Stress urinary incontinence",
                                                                       exposure == "urgency_urinary_incontinence" ~ "Urgency urinary incontinence",
                                                                       exposure == "mixed_urinary_incontinence" ~ "Mixed urinary incontinence",
                                                                       exposure == "anxiety" ~ "Anxiety" ,
                                                                       TRUE ~ exposure))

snps_used <- snps_used %>% dplyr::mutate(outcome = dplyr::case_when(outcome == "urinary_incontinence" ~ "Urinary incontinence" ,
                                                                      outcome == "stress_urinary_incontinence" ~ "Stress urinary incontinence",
                                                                      outcome == "urgency_urinary_incontinence" ~ "Urgency urinary incontinence",
                                                                      outcome == "mixed_urinary_incontinence" ~ "Mixed urinary incontinence",
                                                                      outcome == "MDD" ~ "Depression" ,
                                                                      outcome == "broad_depression_phenotype" ~ "Broad depression phenotype",
                                                                      outcome == "anxiety" ~ "Anxiety" ,
                                                                      outcome == "neuroticism" ~ "Neuroticism" ,
                                                                      TRUE ~ outcome))


# Order rows
snps_used$exposure <- factor(snps_used$exposure, levels = c("Urinary incontinence","Stress urinary incontinence","Urgency urinary incontinence","Mixed urinary incontinence", "Anxiety"))
snps_used$outcome <- factor(snps_used$outcome, levels = c("Anxiety","Depression","Broad depression phenotype","Neuroticism",
                                                          "Urinary incontinence","Stress urinary incontinence","Urgency urinary incontinence","Mixed urinary incontinence"))

snps_used <-  snps_used[order(snps_used$exposure,snps_used$outcome),]

write.csv(snps_used,"results/formatted_tables/supp_table_3.csv",row.names = F)
