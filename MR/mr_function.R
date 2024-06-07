#Function to run and save MR
#exposure_data <- exposure_main
#outcome_data <- outcome_main

run_mr <- function(exposure_data, outcome_data){
  
  exposure_data[exposure_data$beta.exposure < 0, c("effect_allele.exposure","other_allele.exposure")] <- exposure_data[exposure_data$beta.exposure < 0, c("other_allele.exposure","effect_allele.exposure")]
  exposure_data$beta.exposure <- abs(exposure_data$beta.exposure)
  
  results <- NULL
  compiled_het <- NULL
  compiled_pleio <- NULL
  compiled_res_single <- NULL
  compiled_res_loo <- NULL
  
  # Harmonise data
  # If action = 1 then palindromic SNPs are kept in the analyses and if action = 3 then they are dropped
  df <- harmonise_data(exposure_data, outcome_data, action = 2)
  df$f_stat <- (df$beta.exposure^2)/(df$se.exposure^2)
  
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
  
  results <- rbind(results,mr_results)
  
  # Run debiased IVW
  divw <- mr_divw(input)
  results[nrow(results)+1,] <- c(unique(df$exposure),unique(df$outcome),"Debiased IVW",divw@SNPs, divw@Estimate,divw@StdError,divw@Pvalue)
  
  # Run sensitivity analyses
  # Is there evidence of heterogeneity in the genetic effects?
  
  if(nrow(df)>1){
    het <- mr_heterogeneity(df)
  }else{
    het <- as.data.frame(matrix(nrow = 0, ncol = 9))
    colnames(het) <- c("id.exposure","id.outcome","outcome","exposure","method","Q","Q_df","Q_pval")
  }
  
  compiled_het <- rbind(compiled_het,het)
  
  #Is there evidence of horizontal pleiotropy?
  if(nrow(df)>1){
    pleio <- mr_pleiotropy_test(df)
  }else{
    pleio <- as.data.frame(matrix(nrow = 0, ncol = 8))
    colnames(pleio) <- c("id.exposure","id.outcome","outcome","exposure","egger_intercept","se","pval")
  }
  
  compiled_pleio <- rbind(compiled_pleio,pleio)
  
  #Perform 2 sample MR on each SNP individually
  res_single <- mr_singlesnp(df)
  res_single <- res_single %>% dplyr::select(exposure, outcome,SNP,b,se,p)

  compiled_res_single <- rbind(compiled_res_single,res_single)
  
  # Leave one out sensitivity analysis
  res_loo <- mr_leaveoneout(df)
  res_loo <- res_loo %>% dplyr::select(exposure, outcome,SNP,b,se,p)

  compiled_res_loo <- rbind(compiled_res_loo,res_loo)
  
  if(nrow(df)>0){
    df <- df %>% select(SNP,exposure,outcome,effect_allele.exposure,other_allele.exposure,eaf.exposure,
                        beta.exposure,se.exposure,pval.exposure,f_stat)
  }else{
    df <- as.data.frame(matrix(nrow = 0, ncol = 10))
    colnames(df) <- c("SNP","exposure","outcome","effect_allele.exposure","other_allele.exposure","eaf.exposure",
                      "beta.exposure","se.exposure","pval.exposure","f_stat")
  }
  
  
  # Save results
  return(list(results,compiled_het,compiled_pleio,compiled_res_single,compiled_res_loo,df))
}
