library(TwoSampleMR)
library(MRInstruments)
library(LDlinkR)
library(dplyr)
library(data.table)
library(MendelianRandomization)

source("MR/mr_function.R")
#exposure_name = "MDD"
#outcome_name = "urinary_incontinence"

for(exposure_name in c("MDD","broad_depression_phenotype","neuroticism","anxiety")){
  
  compiled_mr_results <- NULL
  compiled_heterogeneity <- NULL
  compiled_pleiotropy <- NULL
  compiled_results_single <- NULL
  compiled_results_loo <- NULL
  compiled_f_stat <- NULL
  
  for (outcome_name in c("urinary_incontinence","mixed_urinary_incontinence","stress_urinary_incontinence","urgency_urinary_incontinence")) {
    
    print(paste0("working on ",exposure_name," ",outcome_name))
    
    exposure <- read.csv(paste0("data/exposure/exposure_", exposure_name,"_p_value_5e-08_outcome_",outcome_name, "_with_proxies.csv"))
    outcome <- read.csv(paste0("data/outcome/exposure_", exposure_name,"_p_value_5e-08_outcome_",outcome_name, ".csv"))
    
    
    if(nrow(exposure)>0 & nrow(outcome)>0){
      print("Running MR analyses")
      #Run main analyses using p-value of 5e-08 to select exposure SNPs 
      #(SNPs have been selected and clumped in pre-process script - see preprocess/format_UI_exposures.R
      # & preprocess/format_MH_exposures.R)
      
      # Load exposure data
      exposure_main <- read_exposure_data(file = paste0("data/exposure/exposure_", exposure_name,"_p_value_5e-08_outcome_",outcome_name, "_with_proxies.csv"),
                                          clump = FALSE,
                                          sep = ",",
                                          snp_col = "rsid",
                                          beta_col = "beta",
                                          se_col = "se",
                                          effect_allele_col = "effect_allele",
                                          other_allele_col = "other_allele",
                                          pval_col = "pval",
                                          eaf_col = "eaf",
                                          phenotype_col = "phenotype")
      
      # Load outcome data
      
      outcome_main <- read_outcome_data(filename = paste0("data/outcome/exposure_", exposure_name,"_p_value_5e-08_outcome_",outcome_name, ".csv"),
                                        snps = exposure_main$rsid,
                                        sep = ",",
                                        snp_col = "SNP",
                                        beta_col = "beta",
                                        se_col = "se",
                                        effect_allele_col = "effect_allele",
                                        other_allele_col = "other_allele",
                                        pval_col = "P",
                                        eaf_col = "effect_allele_frequency",
                                        phenotype_col = "phenotype")
      
      
      results <- run_mr(exposure_main, outcome_main)
    
      mr_results <- results[[1]]
      heterogeneity <- results[[2]]
      pleiotropy <-  results[[3]]
      results_single <- results[[4]]
      results_loo <- results[[5]]
      f_stat <- results[[6]]
      
      compiled_mr_results <- rbind(compiled_mr_results, mr_results)
      compiled_heterogeneity <- rbind(compiled_heterogeneity,heterogeneity)
      compiled_pleiotropy <- rbind(compiled_pleiotropy,pleiotropy)
      compiled_results_single <- rbind(compiled_results_single,results_single)
      compiled_results_loo <- rbind(compiled_results_loo,results_loo)
      compiled_f_stat <- rbind(compiled_f_stat,f_stat)
    }
  }
  
  write.csv(compiled_mr_results,paste0("results/MH_MR/",exposure_name,"_p_value_5e-08.csv"), row.names = F)
  write.csv(compiled_heterogeneity,paste0("results/MH_MR/heterogeneity_",exposure_name,"_p_value_5e-08.csv"), row.names = F)
  write.csv(compiled_pleiotropy,paste0("results/MH_MR/pleiotropy_",exposure_name,"_p_value_5e-08.csv"), row.names = F)
  write.csv(compiled_results_single,paste0("results/MH_MR/results_single_",exposure_name,"_p_value_5e-08.csv"), row.names = F)
  write.csv(compiled_results_loo,paste0("results/MH_MR/results_loo_",exposure_name,"_p_value_5e-08.csv"), row.names = F)
  write.csv(compiled_f_stat,paste0("results/MH_MR/results_f_stat_",exposure_name,"_p_value_5e-08.csv"), row.names = F)
  
}


for(exposure_name in c("MDD","broad_depression_phenotype","neuroticism","anxiety")){
  
  compiled_mr_results <- NULL
  compiled_heterogeneity <- NULL
  compiled_pleiotropy <- NULL
  compiled_results_single <- NULL
  compiled_results_loo <- NULL
  compiled_f_stat <- NULL
  
  for (outcome_name in c("urinary_incontinence","mixed_urinary_incontinence","stress_urinary_incontinence","urgency_urinary_incontinence")) {
    
    print(paste0("working on ",exposure_name," ",outcome_name))
    
    exposure <- read.csv(paste0("data/exposure/exposure_", exposure_name,"_p_value_5e-06_outcome_",outcome_name, "_with_proxies.csv"))
    outcome <- read.csv(paste0("data/outcome/exposure_", exposure_name,"_p_value_5e-06_outcome_",outcome_name, ".csv"))
    
    if(nrow(exposure)>0 & nrow(outcome)>0){
      print("Running MR analyses")
      #Run main analyses using p-value of 5e-06 to select exposure SNPs 
      #(SNPs have been selected and clumped in pre-process script - see preprocess/format_UI_exposures.R
      # & preprocess/format_MH_exposures.R)
      
      # Load exposure data
      exposure_main <- read_exposure_data(file = paste0("data/exposure/exposure_", exposure_name,"_p_value_5e-06_outcome_",outcome_name, "_with_proxies.csv"),
                                          clump = FALSE,
                                          sep = ",",
                                          snp_col = "rsid",
                                          beta_col = "beta",
                                          se_col = "se",
                                          effect_allele_col = "effect_allele",
                                          other_allele_col = "other_allele",
                                          pval_col = "pval",
                                          phenotype_col = "phenotype",
                                          eaf_col = "eaf")
      
      # Load outcome data
      
      outcome_main <- read_outcome_data(filename = paste0("data/outcome/exposure_", exposure_name,"_p_value_5e-06_outcome_",outcome_name, ".csv"),
                                        snps = exposure_main$rsid,
                                        sep = ",",
                                        snp_col = "SNP",
                                        beta_col = "beta",
                                        se_col = "se",
                                        effect_allele_col = "effect_allele",
                                        other_allele_col = "other_allele",
                                        pval_col = "P",
                                        eaf_col = "effect_allele_frequency",
                                        phenotype_col = "phenotype")
      
      results <- run_mr(exposure_main, outcome_main)
      
      mr_results <- results[[1]]
      heterogeneity <- results[[2]]
      pleiotropy <-  results[[3]]
      results_single <- results[[4]]
      results_loo <- results[[5]]
      f_stat <- results[[6]]
      
      compiled_mr_results <- rbind(compiled_mr_results, mr_results)
      compiled_heterogeneity <- rbind(compiled_heterogeneity,heterogeneity)
      compiled_pleiotropy <- rbind(compiled_pleiotropy,pleiotropy)
      compiled_results_single <- rbind(compiled_results_single,results_single)
      compiled_results_loo <- rbind(compiled_results_loo,results_loo)
      compiled_f_stat <- rbind(compiled_f_stat,f_stat)
    }
  }
  
  write.csv(compiled_mr_results,paste0("results/MH_MR/",exposure_name,"_p_value_5e-06.csv"), row.names = F)
  write.csv(compiled_heterogeneity,paste0("results/MH_MR/heterogeneity_",exposure_name,"_p_value_5e-06.csv"), row.names = F)
  write.csv(compiled_pleiotropy,paste0("results/MH_MR/pleiotropy_",exposure_name,"_p_value_5e-06.csv"), row.names = F)
  write.csv(compiled_results_single,paste0("results/MH_MR/results_single_",exposure_name,"_p_value_5e-06.csv"), row.names = F)
  write.csv(compiled_results_loo,paste0("results/MH_MR/results_loo_",exposure_name,"_p_value_5e-06.csv"), row.names = F)
  write.csv(compiled_f_stat,paste0("results/MH_MR/results_f_stat_",exposure_name,"_p_value_5e-06.csv"), row.names = F)
}
