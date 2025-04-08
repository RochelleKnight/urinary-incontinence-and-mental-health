df <- as.data.frame(matrix(nrow = 0, ncol = 7))
colnames(df) <- c("rsid","effect_allele","other_allele","beta","se","pval","phenotype")
df[nrow(df)+1,] <- c("rs34998271","a","g","0.488","0.0835","1.7e-09","UUI_overall")
df[nrow(df)+1,] <- c("rs138724718","a","g","0.598","0.0835","3.39e-09","SUI_overall")
write.csv(df, "data/exposure/exposure_UI_overall.csv")

results <- NULL

# MDD
outcome <- fread(file="data/raw_data/MDD2018_ex23andMe.txt", fill = T)
outcome$beta <- log(outcome$OR)
outcome$eaf <- (outcome$FRQ_A_59851*59851 + outcome$FRQ_U_113154*113154)/(59851+113154)
outcome <- outcome %>% select(SNP,eaf, A1,A2,beta, SE, P) %>%
  rename(effect_allele = A1,
         other_allele = A2,
         se = SE)

outcome_snp <- outcome[outcome$SNP %in% df$rsid]
outcome_snp$phenotype <- "MDD"
results <- rbind(results,outcome_snp)


# Broad depression phenotype
outcome <- fread(file="data/raw_data/PGC_UKB_depression_genome-wide.txt", fill = T)
outcome <- outcome %>% select(MarkerName,Freq, A1,A2,LogOR, StdErrLogOR, P) %>%
  rename(SNP = MarkerName,
         effect_allele = A1,
         other_allele = A2,
         se = StdErrLogOR,
         eaf = Freq,
         beta = LogOR)

outcome_snp <- outcome[outcome$SNP %in% df$rsid]
outcome_snp$phenotype <- "broad_depression_phenotype"
results <- rbind(results,outcome_snp)

# Neuroticism
outcome <- fread(file="data/raw_data/SummaryStats.txt", fill = T)
outcome_eaf <- fread(file="data/raw_data/UKB_Neu_AF.txt", fill = T)
outcome <- outcome %>% left_join(outcome_eaf)
rm(outcome_eaf)

outcome <- outcome %>% select(rsid,af, a_1,a_0,N_res_beta, N_res_se, p_value) %>%
  rename(SNP = rsid,
         effect_allele = a_1,
         other_allele = a_0,
         se = N_res_se,
         eaf = af,
         beta = N_res_beta,
         P = p_value)

outcome_snp <- outcome[outcome$SNP %in% df$rsid]
outcome_snp$phenotype <- "neuroticism"
results <- rbind(results,outcome_snp)

# Anxiety
outcome <- fread(file="data/raw_data/TotAnx_effect_sumstats.txt", fill = T)
outcome <- outcome %>% select(SNP,af,A2,A1,effect, SE, P) %>%
  rename(effect_allele = A2,
         other_allele = A1,
         se = SE,
         eaf = af,
         beta = effect)

outcome_snp <- outcome[outcome$SNP %in% df$rsid]
outcome_snp$phenotype <- "anxiety"
results <- rbind(results,outcome_snp)

write.csv(results, "data/outcome/exposure_UI_overall.csv")
