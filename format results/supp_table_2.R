library(dplyr)
rm(list = ls())

# Load UI MR results
mr_results_urinary_incontinence <- read.csv(paste0("results/UI_MR/urinary_incontinence_p_value_5e-08.csv"))

# Load MH MR results
mr_results_MDD <- read.csv(paste0("results/MH_MR/MDD_p_value_5e-08.csv"))
mr_results_broad_depression_phenotype <- read.csv(paste0("results/MH_MR/broad_depression_phenotype_p_value_5e-08.csv"))
mr_results_neuroticism <- read.csv(paste0("results/MH_MR/neuroticism_p_value_5e-08.csv"))
mr_results_anxiety <- read.csv(paste0("results/MH_MR/anxiety_p_value_5e-08.csv"))

mr_results <- rbind(mr_results_urinary_incontinence,mr_results_MDD,mr_results_broad_depression_phenotype,
                    mr_results_neuroticism,mr_results_anxiety)

rm(mr_results_urinary_incontinence,mr_results_MDD,mr_results_broad_depression_phenotype,
   mr_results_neuroticism,mr_results_anxiety)

# Select IVW and Wald ratio results
mr_results <- mr_results %>% filter(method %in% c("Inverse variance weighted","MR Egger","Weighted median","Weighted mode")) %>%
  select(exposure,outcome,nsnp,method,b,se,pval)

# Load heterogeneity results
heterogeneity_urinary_incontinence <- read.csv("results/UI_MR/heterogeneity_urinary_incontinence_p_value_5e-08.csv")
heterogeneity_MDD <- read.csv("results/MH_MR/heterogeneity_MDD_p_value_5e-08.csv")
heterogeneity_broad_depression_phenotype <- read.csv("results/MH_MR/heterogeneity_broad_depression_phenotype_p_value_5e-08.csv")
heterogeneity_neuroticism <- read.csv("results/MH_MR/heterogeneity_neuroticism_p_value_5e-08.csv")
heterogeneity_anxiety <- read.csv("results/MH_MR/heterogeneity_anxiety_p_value_5e-08.csv")

heterogeneity <- rbind(heterogeneity_urinary_incontinence,heterogeneity_MDD,heterogeneity_broad_depression_phenotype,heterogeneity_neuroticism,heterogeneity_anxiety)
rm(heterogeneity_urinary_incontinence,heterogeneity_MDD,heterogeneity_broad_depression_phenotype,heterogeneity_neuroticism,heterogeneity_anxiety)

heterogeneity$id.exposure <- NULL
heterogeneity$id.outcome <- NULL

mr_results <- mr_results %>% left_join(heterogeneity)
rm(heterogeneity)

#Format estimate
mr_results$causal_estimate <- ifelse(mr_results$outcome != "neuroticism", exp(mr_results$b),mr_results$b)
mr_results$lower_CI <- ifelse(mr_results$outcome != "neuroticism", exp(mr_results$b - qnorm(0.975)*mr_results$se),mr_results$b - qnorm(0.975)*mr_results$se)
mr_results$upper_CI <- ifelse(mr_results$outcome != "neuroticism", exp(mr_results$b + qnorm(0.975)*mr_results$se),mr_results$b + qnorm(0.975)*mr_results$se)
mr_results$est <- ifelse(mr_results$causal_estimate > 0.1, paste0(sprintf("%.2f",mr_results$causal_estimate)," (",sprintf("%.2f",mr_results$lower_CI),"-",sprintf("%.2f",mr_results$upper_CI),"); ",sprintf("%.2f",mr_results$pval)),
                                        ifelse(mr_results$causal_estimate < -0.0009,paste0(sprintf("%.3f",mr_results$causal_estimate)," (",sprintf("%.3f",mr_results$lower_CI),"-",sprintf("%.3f",mr_results$upper_CI),"); ",sprintf("%.2f",mr_results$pval)),
                                               paste0(sprintf("%.4f",mr_results$causal_estimate)," (",sprintf("%.4f",mr_results$lower_CI),"-",sprintf("%.4f",mr_results$upper_CI),"); ",sprintf("%.2f",mr_results$pval))))

mr_results[c("b","se","lower_CI","upper_CI","causal_estimate","pval")] <- NULL

mr_results <- mr_results %>% 
  mutate(Q = case_when(
    Q < 1 ~ as.character(sprintf("%.3f", Q)),
    Q > 1 & Q < 10 ~ as.character(sprintf("%.2f", Q)),
    Q >10 & Q < 100 ~ as.character(sprintf("%.1f", Q)),
    Q > 100 ~ as.character(sprintf("%.0f", Q)),
    TRUE ~ as.character(Q)
  ),
  Q_pval = case_when(
    Q_pval > 0.1 ~ as.character(sprintf("%.2f", Q_pval)),
    Q_pval > 0.01 & Q_pval < 0.1 ~ as.character(sprintf("%.2f", Q_pval)),
    Q_pval < 0.01 ~ as.character(sprintf("%.3f", Q_pval)),
    TRUE ~ as.character(Q_pval)
  ),
  heterogeneity = paste0(Q, "; ",Q_pval)) %>%
  select(-c(Q,Q_pval,Q_df))

# Format names
mr_results <- mr_results %>% dplyr::mutate(exposure = dplyr::case_when(exposure == "urinary_incontinence_p_value_5e-08" ~ "Urinary incontinence" ,
                                                                       exposure == "MDD_p_value_5e-08" ~ "Depression" ,
                                                                       exposure == "broad_depression_phenotype_p_value_5e-08" ~ "Broad depression phenotype",
                                                                       exposure == "anxiety_p_value_5e-08" ~ "Anxiety" ,
                                                                       exposure == "neuroticism_p_value_5e-08" ~ "Neuroticism",
                                                                       TRUE ~ exposure))

mr_results <- mr_results %>% dplyr::mutate(outcome = dplyr::case_when(outcome == "urinary_incontinence" ~ "Urinary incontinence" ,
                                                                      outcome == "stress_urinary_incontinence" ~ "Stress urinary incontinence",
                                                                      outcome == "urgency_urinary_incontinence" ~ "Urgency urinary incontinence",
                                                                      outcome == "mixed_urinary_incontinence" ~ "Mixed urinary incontinence",
                                                                      outcome == "MDD" ~ "Depression" ,
                                                                      outcome == "broad_depression_phenotype" ~ "Broad depression phenotype",
                                                                      outcome == "anxiety" ~ "Anxiety" ,
                                                                      outcome == "neuroticism" ~ "Neuroticism" ,
                                                                      TRUE ~ outcome))

# Order results
mr_results$order <- paste0(mr_results$exposure, " , ", mr_results$outcome)

mr_results$order <- factor(mr_results$order, levels = c(paste0("Urinary incontinence , ",c("Anxiety","Depression","Broad depression phenotype","Neuroticism")),
                                                        paste0("Anxiety , ",c("Urinary incontinence","Stress urinary incontinence","Urgency urinary incontinence","Mixed urinary incontinence")),
                                                        paste0("Depression , ",c("Urinary incontinence","Stress urinary incontinence","Urgency urinary incontinence","Mixed urinary incontinence")),
                                                        paste0("Broad depression phenotype , ",c("Urinary incontinence","Stress urinary incontinence","Urgency urinary incontinence","Mixed urinary incontinence")),
                                                        paste0("Neuroticism , ",c("Urinary incontinence","Stress urinary incontinence","Urgency urinary incontinence","Mixed urinary incontinence"))))
mr_results <- mr_results[order(mr_results$order),]
mr_results$order <- NULL

mr_results <- tidyr::pivot_wider(mr_results, names_from = method, values_from = c(est,heterogeneity))
colnames(mr_results)

mr_results <- mr_results %>% select("exposure","outcome","nsnp","est_Inverse variance weighted","heterogeneity_Inverse variance weighted",
                                    "est_MR Egger","heterogeneity_MR Egger",
                                    "est_Weighted median","heterogeneity_Weighted median",
                                    "est_Weighted mode","heterogeneity_Weighted mode")
                                    
                                    
# Save results
write.csv(mr_results, "results/formatted_tables/supp_table_2.csv", row.names = F)

