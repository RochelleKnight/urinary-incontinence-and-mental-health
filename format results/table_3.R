library(dplyr)

rm(list = ls())
dir.create("results/formatted_tables", showWarnings = FALSE)

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
   mr_results_stress_urinary_incontinence,mr_results_urgency_urinary_incontinence,
   mr_results_anxiety)

# Select IVW and Wald ratio results
mr_results <- mr_results %>% filter(method %in% c("Debiased IVW")) %>%
  select(exposure,outcome,nsnp,method,b,se,pval)

#Format estimate
mr_results$causal_estimate <- ifelse(mr_results$outcome != "neuroticism", exp(mr_results$b),mr_results$b)
mr_results$lower_CI <- ifelse(mr_results$outcome != "neuroticism", exp(mr_results$b - qnorm(0.975)*mr_results$se),mr_results$b - qnorm(0.975)*mr_results$se)
mr_results$upper_CI <- ifelse(mr_results$outcome != "neuroticism", exp(mr_results$b + qnorm(0.975)*mr_results$se),mr_results$b + qnorm(0.975)*mr_results$se)
mr_results$causal_estimate_CI <- ifelse(mr_results$causal_estimate > 0.1, paste0(sprintf("%.2f",mr_results$causal_estimate)," (",sprintf("%.2f",mr_results$lower_CI),"-",sprintf("%.2f",mr_results$upper_CI),")"),
                                        ifelse(mr_results$causal_estimate < -0.0009,paste0(sprintf("%.3f",mr_results$causal_estimate)," (",sprintf("%.3f",mr_results$lower_CI),"-",sprintf("%.3f",mr_results$upper_CI),")"),
                                               paste0(sprintf("%.4f",mr_results$causal_estimate)," (",sprintf("%.4f",mr_results$lower_CI),"-",sprintf("%.4f",mr_results$upper_CI),")")))

mr_results$causal_estimate_pval <- sprintf("%.2f",mr_results$pval)

mr_results[c("b","se","lower_CI","upper_CI","causal_estimate","method","pval")] <- NULL

# Format names
mr_results <- mr_results %>% dplyr::mutate(exposure = dplyr::case_when(exposure == "urinary_incontinence_p_value_5e-06" ~ "Urinary incontinence",
                                                                       exposure == "stress_urinary_incontinence_p_value_5e-06" ~ "Stress urinary incontinence",
                                                                       exposure == "urgency_urinary_incontinence_p_value_5e-06" ~ "Urgency urinary incontinence",
                                                                       exposure == "mixed_urinary_incontinence_p_value_5e-06" ~ "Mixed urinary incontinence",
                                                                       exposure == "anxiety_p_value_5e-06" ~ "Anxiety" ,
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
                                                        paste0("Stress urinary incontinence , ",c("Anxiety","Depression","Broad depression phenotype","Neuroticism")),
                                                        paste0("Urgency urinary incontinence , ",c("Anxiety","Depression","Broad depression phenotype","Neuroticism")),
                                                        paste0("Mixed urinary incontinence , ",c("Anxiety","Depression","Broad depression phenotype","Neuroticism")),
                                                        paste0("Anxiety , ",c("Urinary incontinence","Stress urinary incontinence","Urgency urinary incontinence","Mixed urinary incontinence"))))
mr_results <- mr_results[order(mr_results$order),]
mr_results$order <- NULL


colnames(mr_results) <- c("Exposure","Outcome","nsnp","Causal estimate (95% CI)","P-value")
write.csv(mr_results, "results/formatted_tables/table_3.csv",row.names = F)
