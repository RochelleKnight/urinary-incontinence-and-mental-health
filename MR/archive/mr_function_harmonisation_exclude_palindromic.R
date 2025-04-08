#Function to run and save MR
#exposure_data <- exposure_main
#outcome_data <- outcome_main
#include_palindromic_snps <- "1"
run_mr_including_excluding_palindromic <- function(exposure_data, outcome_data){
  
  exposure_data[exposure_data$beta.exposure < 0, c("effect_allele.exposure","other_allele.exposure")] <- exposure_data[exposure_data$beta.exposure < 0, c("other_allele.exposure","effect_allele.exposure")]
  exposure_data$beta.exposure <- abs(exposure_data$beta.exposure)
  
  results <- NULL
  compiled_het <- NULL
  compiled_pleio <- NULL
  compiled_res_single <- NULL
  compiled_res_loo <- NULL
  
  for(include_palindromic_snps in c("1","3")){
    # Harmonise data
    # If action = 1 then palindromic SNPs are kept in the analyses and if action = 3 then they are dropped
    df <- harmonise_data(exposure_data, outcome_data, action = include_palindromic_snps)
    
    # Input for debiased IVW
    input <- mr_input(bx = df$beta.exposure,
                      bxse = df$se.exposure,
                      by = df$beta.outcome,
                      byse = df$se.outcome,
                      exposure = unique(df$exposure),
                      outcome = unique(df$outcome),
                      snps = df$SNP)
    
    # Run MR
    mr_results <- mr(df)
    
    mr_results <- mr_results %>% dplyr::select(exposure,outcome,method,nsnp,b,se,pval)
    
    if(include_palindromic_snps == "1"){
      mr_results$palindromic_snps <- "included"
    }else{
      mr_results$palindromic_snps <- "excluded"
    }
    
    results <- rbind(results,mr_results)
    
    # Run debiased IVW
    divw <- mr_divw(input)
    results[nrow(results)+1,] <- c(unique(df$exposure),unique(df$outcome),"Debiased IVW",divw@SNPs, divw@Estimate,divw@StdError,divw@Pvalue,ifelse(include_palindromic_snps == "1","included","excluded"))
    
    # Run sensitivity analyses
    # Is there evidence of heterogeneity in the genetic effects?
    
    if(nrow(df)>1){
      het <- mr_heterogeneity(df)
      if(include_palindromic_snps == "1"){
        het$palindromic_snps <- "included"
      }else{
        het$palindromic_snps <- "excluded"
      }
    }else{
      het <- as.data.frame(matrix(nrow = 0, ncol = 8))
      colnames(het) <- c("id.exposure","id.outcome","outcome","exposure","method","Q","Q_df","Q_pval")
    }
    
    compiled_het <- rbind(compiled_het,het)
    
    #Is there evidence of horizontal pleiotropy?
    if(nrow(df)>1){
      pleio <- mr_pleiotropy_test(df)
      if(include_palindromic_snps == "1"){
        pleio$palindromic_snps <- "included"
      }else{
        pleio$palindromic_snps <- "excluded"
      }
    }else{
      pleio <- as.data.frame(matrix(nrow = 0, ncol = 8))
      colnames(pleio) <- c("id.exposure","id.outcome","outcome","exposure","egger_intercept","se","pval","palindromic_snps")
    }
    
    compiled_pleio <- rbind(compiled_pleio,pleio)
    
    #Perform 2 sample MR on each SNP individually
    res_single <- mr_singlesnp(df)
    res_single <- res_single %>% dplyr::select(exposure, outcome,SNP,b,se,p)
    if(include_palindromic_snps == "1"){
      res_single$palindromic_snps <- "included"
    }else{
      res_single$palindromic_snps <- "excluded"
    }
    compiled_res_single <- rbind(compiled_res_single,res_single)
    
    # Leave one out sensitivity analysis
    res_loo <- mr_leaveoneout(df)
    res_loo <- res_loo %>% dplyr::select(exposure, outcome,SNP,b,se,p)
    if(include_palindromic_snps == "1"){
      res_loo$palindromic_snps <- "included"
    }else{
      res_loo$palindromic_snps <- "excluded"
    }
    compiled_res_loo <- rbind(compiled_res_loo,res_loo)
    
  }
  
  # Save results
  return(list(results,compiled_het,compiled_pleio,compiled_res_single,compiled_res_loo))
}
