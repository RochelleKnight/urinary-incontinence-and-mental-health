
check_gwas_data <- function(df){
  print(paste0("Max beta value (log scale) ", max(df$beta)))
  print(paste0("Max beta value (exponentiated) ", exp(max(df$beta))))
  
  print(paste0("Min beta value (log scale) ", min(df$beta)))
  print(paste0("Min beta value (exponentiated) ", exp(min(df$beta))))
  
  print(paste0("Max standard error ", max(df$se)))
  
  print(paste0("Any issues: ", exp(max(df$beta)) > 10 | exp(min(df$beta)) >10 | max(df$se)>10))
  
}