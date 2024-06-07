# Function to format all the outcome data sets
#exposure_data_set <- exposure
#exposure_name <- "MDD"
format_UI_outcomes <- function(exposure_data_set,exposure_name){
  
  main <- exposure_data_set %>% filter(P < 5e-08)
  main$phenotype <- paste0(exposure_name, "_p_value_5e-08")
  main$proxy <- "non_proxy_snp"
  
  sensitivity <- exposure %>% filter(P < 5e-06)
  sensitivity$phenotype <- paste0(exposure_name, "_p_value_5e-06")
  sensitivity$proxy <- "non_proxy_snp"
  
  main <- main %>% filter(!grepl("chr",SNP))
  sensitivity <- sensitivity %>% filter(!grepl("chr",SNP))
  
  # Save analyses exposure datasets
  write.csv(main, paste0("data/exposure/", exposure_name, "_p_value_5e-08_raw.csv"), row.names = F)
  write.csv(sensitivity, paste0("data/exposure/", exposure_name, "_p_value_5e-06_raw.csv"), row.names = F)
  
  # Clump SNPs
  if(nrow(main)>0){
    main <- rename(main, pval = P, rsid = SNP) #rename for ld_clump
    # To run ld_clump through the ieugwas API use the below line however the server can be busy
    # so this command won't always run
    #main <-  ld_clump(main,clump_kb=10000, clump_r2=0.001, pop = "EUR")
    # To run ld_clump locally, use the below line. Referenece panel has been downloaded so this will run
    main <-  ld_clump(main,clump_kb=10000, clump_r2=0.001,
                      plink_bin = get_plink_exe(),
                      bfile = "data/reference_panels/EUR")
    main$id <- NULL
    write.csv(main, paste0("data/exposure/", exposure_name, "_p_value_5e-08_clumped.csv"), row.names = F)
  }
  
  if(nrow(sensitivity)>0){
    sensitivity <- rename(sensitivity, pval = P, rsid = SNP) #rename for ld_clump
    # To run ld_clump through the ieugwas API use the below line however the server can be busy
    # so this command won't always run
    #sensitivity <-  ld_clump(sensitivity,clump_kb=10000, clump_r2=0.001, pop = "EUR") 
    
    # To run ld_clump locally, use the below line. Referenece panel has been downloaded so this will run
    sensitivity <-  ld_clump(sensitivity,clump_kb=10000, clump_r2=0.001,
                             plink_bin = get_plink_exe(),
                             bfile = "data/reference_panels/EUR")
    sensitivity$id <- NULL
    write.csv(sensitivity, paste0("data/exposure/", exposure_name, "_p_value_5e-06_clumped.csv"), row.names = F)
  }
  
  analyses_to_run <- NULL
  
  if(nrow(main)>0){
    analyses_to_run <- append(analyses_to_run,"main")
  }
  if(nrow(sensitivity)>0){
    analyses_to_run <- append(analyses_to_run,"sensitivity")
  }
  
  for(analyses in analyses_to_run){
    print(paste0("Working on ", analyses))
    results <- NULL
    
    # Load and format outcome data
    for(outcome_name in c("urinary_incontinence","mixed_urinary_incontinence","stress_urinary_incontinence","urgency_urinary_incontinence")){
      print(paste0("Working on ",outcome_name ))
      
      if(analyses == "main"){
        exposure_df <- main
        save_name <- "_p_value_5e-08"
      }else if(analyses == "sensitivity"){
        exposure_df <- sensitivity
        save_name <- "_p_value_5e-06"
      }
      
      outcome <- read.table("data/raw_data/",outcome_name,".tsv", header = T)
      
      outcome <- outcome %>% rename("SNP" = "variant_id",
                                      "se" = "standard_error",
                                      "P" = "p_value") 
                                          
      outcome$phenotype <- outcome_name
      outcome_snp <- outcome[outcome$SNP %in% exposure_df$rsid,]
      
      # Find proxy SNPs
      missing_snps <- exposure_df$rsid[which(!exposure_df$rsid %in% outcome_snp$SNP)]
      
      if(length(missing_snps)>0){
        print(paste0(length(missing_snps), " missing SNPs - finding proxies"))
        
        if(file.exists("combined_query_snp_list_grch37.txt")){
          file.remove("combined_query_snp_list_grch37.txt")
        }
        
        LDproxy_batch(missing_snps, pop="EUR", r2d = "r2", token = "5cb32cd8137f", append=TRUE) #using missing snps extracts proxies and saves into txt file (combined_query_snp_list)
        proxySNPs <- read.table("combined_query_snp_list_grch37.txt", header=TRUE, row.names=NULL)
        file.remove("combined_query_snp_list_grch37.txt")
        write.table(proxySNPs,paste0("lib/proxy_snps_exposure_", exposure_name, "_outcome_", outcome_name, save_name, ".txt"), row.names=F,quote=F)
        #proxySNPs <- read.table(paste0("lib/proxy_snps_exposure_", exposure_name, "_outcome_", outcome_name, save_name, ".txt"), header = T)
        
        # Remove SNPs where query snp == RS_number (this is just the same as the query SNP) and without an RS number
        proxySNPs <- proxySNPs %>% filter(query_snp != RS_Number &
                                            RS_Number != ".") 
       
        for(snp in unique(proxySNPs$query_snp)){
          tmp_proxy <- proxySNPs %>% filter(query_snp == snp 
                                            & R2 > 0.8) #keep only ld snps with r2 > 0.8
          
          # check which proxy snps are present in both the exposure and outcome gwas
          tmp_proxy <- tmp_proxy[which(tmp_proxy$RS_Number %in% exposure_data_set$SNP &
                                         tmp_proxy$RS_Number %in% outcome$SNP),]
        
          tmp_proxy <- tmp_proxy %>% 
            slice(c(which.min(row.names))) # select rsid with lowest row.names (i.e., higher R2)
          
          if(nrow(tmp_proxy)>0){
            #select the proxy snp from the exposure gwas data set
            exposure_proxy_snps <- exposure_data_set[exposure_data_set$SNP %in% tmp_proxy$RS_Number,]
            exposure_proxy_snps <- rename(exposure_proxy_snps, pval = P, rsid = SNP)
            exposure_proxy_snps$proxy <- "proxy_snp"
            exposure_df <- plyr::rbind.fill(exposure_df,exposure_proxy_snps)
          }
        }
        
        exposure_df$phenotype <- paste0(exposure_name, save_name)
        exposure_df$outcome <- outcome_name
        exposure_df <- exposure_df[!duplicated(exposure_df),]
        print(paste0("Exposure duplicated SNPs: ",anyDuplicated(exposure_df$rsid))) 
        
        # Get outcome data
        print(paste0("Number of outcome SNPs before adding proxies: ", nrow(outcome_snp)))
        outcome_snp <- outcome[outcome$SNP %in% exposure_df$rsid,]
        print(paste0("Number of outcome SNPs after ading proxies: ", nrow(outcome_snp)))
      }else{
        print("No missing SNPs")
        exposure_df$outcome <- outcome_name
      }
      
      results <- rbind(results, outcome_snp)
      write.csv(exposure_df,paste0("data/exposure/exposure_", exposure_name,save_name,"_outcome_",outcome_name, "_with_proxies.csv"), row.names = F)
      write.csv(outcome_snp, paste0("data/outcome/exposure_", exposure_name,save_name,"_outcome_",outcome_name, ".csv"), row.names = F )
      
    }
    
    write.csv(results, paste0("data/outcome/outcome_data_for_", exposure_name,save_name,".csv") )
  }
}

