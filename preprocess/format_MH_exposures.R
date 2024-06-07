library(data.table)
library(dplyr)
library(LDlinkR)
library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)
library(R.utils)
library(plinkbinr)
library(gwasvcf)
# Prep all instruments
# Find all proxy SNPs -

setwd("~/GitHub/urinary-incontinence-and-mental-health")

source("preprocess/function_format_MH_exposures.R")
source("preprocess/check_gwas_data.R")

########
# MDD
########
exposure <- fread(file="data/raw_data/MDD.gz", fill = T)
exposure$beta <- log(exposure$OR)
exposure$eaf <- (exposure$FRQ_A_294322*294322 + exposure$FRQ_U_741438*741438)/(294322+741438)
exposure <- exposure %>% select(SNP,eaf, A1,A2,beta, SE, P) %>%
  rename(effect_allele = A1,
         other_allele = A2,
         se = SE)

check_gwas_data(exposure)
format_UI_outcomes(exposure, "MDD")


############################
# Broad depression phenotype
############################

exposure <- fread(file="data/raw_data/PGC_UKB_depression_genome-wide.txt", fill = T)
exposure <- exposure %>% select(MarkerName,Freq, A1,A2,LogOR, StdErrLogOR, P) %>%
  rename(SNP = MarkerName,
         effect_allele = A1,
         other_allele = A2,
         se = StdErrLogOR,
         eaf = Freq,
         beta = LogOR)

check_gwas_data(exposure)
format_UI_outcomes(exposure, "broad_depression_phenotype")





##############
# Neuroticism
##############

exposure <- fread(file="data/raw_data/SummaryStats.txt", fill = T)
exposure_eaf <- fread(file=paste0(data_file_path,"/data/raw_data/UKB_Neu_AF.txt"), fill = T)
exposure <- exposure %>% left_join(exposure_eaf)
rm(exposure_eaf)

exposure <- exposure %>% select(rsid,af, a_1,a_0,N_res_beta, N_res_se, p_value) %>%
  rename(SNP = rsid,
         effect_allele = a_1,
         other_allele = a_0,
         se = N_res_se,
         eaf = af,
         beta = N_res_beta,
         P = p_value)

check_gwas_data(exposure)
format_UI_outcomes(exposure, "neuroticism")

  
##########
# Anxiety
#########
  
exposure <- fread(file="data/raw_data/TotAnx_effect_sumstats.txt", fill = T)
exposure <- exposure %>% select(SNP,af,A2,A1,effect, SE, P) %>%
  rename(effect_allele = A2,
         other_allele = A1,
         se = SE,
         eaf = af,
         beta = effect)

check_gwas_data(exposure)
format_UI_outcomes(exposure, "anxiety")

