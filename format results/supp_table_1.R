library(dplyr)
library(stringr)
rm(list = ls())

dir.create("results/formatted_tables", showWarnings = FALSE)

# Load UI MR results
f_stat_urinary_incontinence_pval_5e08 <- read.csv(paste0("results/UI_MR/results_f_stat_urinary_incontinence_p_value_5e-08.csv"))

f_stat_urinary_incontinence_pval_5e06 <- read.csv(paste0("results/UI_MR/results_f_stat_urinary_incontinence_p_value_5e-06.csv"))
f_stat_mixed_urinary_incontinence_pval_5e06 <- read.csv(paste0("results/UI_MR/results_f_stat_mixed_urinary_incontinence_p_value_5e-06.csv"))
f_stat_stress_urinary_incontinence_pval_5e06 <- read.csv(paste0("results/UI_MR/results_f_stat_stress_urinary_incontinence_p_value_5e-06.csv"))
f_stat_urgency_urinary_incontinence_pval_5e06 <- read.csv(paste0("results/UI_MR/results_f_stat_urgency_urinary_incontinence_p_value_5e-06.csv"))

# Load MH MR results
f_stat_MDD_pval_5e08 <- read.csv(paste0("results/MH_MR/results_f_stat_MDD_p_value_5e-08.csv"))
f_stat_broad_depression_phenotype_pval_5e08 <- read.csv(paste0("results/MH_MR/results_f_stat_broad_depression_phenotype_p_value_5e-08.csv"))
f_stat_neuroticism_pval_5e08 <- read.csv(paste0("results/MH_MR/results_f_stat_neuroticism_p_value_5e-08.csv"))
f_stat_anxiety_pval_5e08 <- read.csv(paste0("results/MH_MR/results_f_stat_anxiety_p_value_5e-08.csv"))

f_stat_MDD_pval_5e06 <- read.csv(paste0("results/MH_MR/results_f_stat_MDD_p_value_5e-06.csv"))
f_stat_broad_depression_phenotype_pval_5e06 <- read.csv(paste0("results/MH_MR/results_f_stat_broad_depression_phenotype_p_value_5e-06.csv"))
f_stat_neuroticism_pval_5e06 <- read.csv(paste0("results/MH_MR/results_f_stat_neuroticism_p_value_5e-06.csv"))
f_stat_anxiety_pval_5e06 <- read.csv(paste0("results/MH_MR/results_f_stat_anxiety_p_value_5e-06.csv"))

f_stat_results <- rbind(f_stat_urinary_incontinence_pval_5e08,f_stat_MDD_pval_5e08,f_stat_broad_depression_phenotype_pval_5e08,
                        f_stat_neuroticism_pval_5e08,f_stat_anxiety_pval_5e08,
                        f_stat_urinary_incontinence_pval_5e06,f_stat_mixed_urinary_incontinence_pval_5e06,f_stat_stress_urinary_incontinence_pval_5e06,f_stat_urgency_urinary_incontinence_pval_5e06,
                        f_stat_MDD_pval_5e06,f_stat_broad_depression_phenotype_pval_5e06,f_stat_neuroticism_pval_5e06,f_stat_anxiety_pval_5e06)

rm(f_stat_urinary_incontinence_pval_5e08,f_stat_MDD_pval_5e08,f_stat_broad_depression_phenotype_pval_5e08,
f_stat_neuroticism_pval_5e08,f_stat_anxiety_pval_5e08,
f_stat_urinary_incontinence_pval_5e06,f_stat_mixed_urinary_incontinence_pval_5e06,f_stat_stress_urinary_incontinence_pval_5e06,f_stat_urgency_urinary_incontinence_pval_5e06,
f_stat_MDD_pval_5e06,f_stat_broad_depression_phenotype_pval_5e06,f_stat_neuroticism_pval_5e06,f_stat_anxiety_pval_5e06)


# Format names
f_stat_results$analysis <- as.character(NA)
f_stat_results$analysis[which(str_detect(f_stat_results$exposure,"5e-08"))] <- "main"
f_stat_results$analysis[which(str_detect(f_stat_results$exposure,"5e-06"))] <- "sensitivity"

f_stat_results <- f_stat_results %>% dplyr::mutate(exposure = dplyr::case_when(exposure %in% c("urinary_incontinence_p_value_5e-08","urinary_incontinence_p_value_5e-06") ~ "Urinary incontinence",
                                                                               exposure == "mixed_urinary_incontinence_p_value_5e-06" ~ "Mixed urinary incontinence",
                                                                               exposure == "stress_urinary_incontinence_p_value_5e-06" ~ "Stress urinary incontinence",
                                                                               exposure == "urgency_urinary_incontinence_p_value_5e-06" ~ "Urgency urinary incontinence",
                                                                               exposure %in% c("MDD_p_value_5e-08","MDD_p_value_5e-06") ~ "Depression",
                                                                               exposure %in% c("broad_depression_phenotype_p_value_5e-08","broad_depression_phenotype_p_value_5e-06") ~ "Broad depression phenotype",
                                                                               exposure %in% c("anxiety_p_value_5e-08","anxiety_p_value_5e-06") ~ "Anxiety",
                                                                               exposure %in% c("neuroticism_p_value_5e-08","neuroticism_p_value_5e-06") ~ "Neuroticism",
                                                                               TRUE ~ exposure))

f_stat_results <- f_stat_results %>% dplyr::mutate(outcome = dplyr::case_when(outcome == "urinary_incontinence" ~ "Urinary incontinence" ,
                                                                      outcome == "stress_urinary_incontinence" ~ "Stress urinary incontinence",
                                                                      outcome == "urgency_urinary_incontinence" ~ "Urgency urinary incontinence",
                                                                      outcome == "mixed_urinary_incontinence" ~ "Mixed urinary incontinence",
                                                                      outcome == "MDD" ~ "Depression" ,
                                                                      outcome == "broad_depression_phenotype" ~ "Broad depression phenotype",
                                                                      outcome == "anxiety" ~ "Anxiety" ,
                                                                      outcome == "neuroticism" ~ "Neuroticism" ,
                                                                      TRUE ~ outcome))

# Order results
f_stat_results$exposure <- factor(f_stat_results$exposure, levels=c("Urinary incontinence","Stress urinary incontinence","Urgency urinary incontinence","Mixed urinary incontinence",
                                                                    "Anxiety","Depression","Broad depression phenotype","Neuroticism")) 

f_stat_results$outcome <- factor(f_stat_results$outcome, levels=c("Urinary incontinence","Stress urinary incontinence","Urgency urinary incontinence","Mixed urinary incontinence",
                                                                  "Anxiety","Depression","Broad depression phenotype","Neuroticism")) 

f_stat_results$analysis <- factor(f_stat_results$analysis, levels=c("main","sensitivity")) 
f_stat_results <- f_stat_results[order(f_stat_results$exposure,f_stat_results$outcome,f_stat_results$analysis),]

f_stat_results <- f_stat_results %>% select(SNP,exposure,outcome,analysis,effect_allele.exposure,other_allele.exposure,
                                            eaf.exposure,beta.exposure,se.exposure,pval.exposure,f_stat)

write.csv(f_stat_results, "results/formatted_tables/supp_table_1.csv",row.names = F)

