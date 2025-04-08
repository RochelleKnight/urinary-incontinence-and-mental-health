#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

 
library(devtools)
#devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(MRInstruments)
library(purrr)
library(openxlsx)
library(meta)
library(metafor)
library(cowplot)
library(gridGraphics)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(readxl)
library(data.table)

if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "purrr", "remotes", "ggplot2", "ggforce", "data.table", "ggforestplot", "gtools", "openxlsx")
pacman::p_load_gh("Spiller/MVMR")

ao <- available_outcomes()

# Set working directory - folder in my computer
setwd(" ")
#---------------------------------------------------------------------#
#                        Smoking Exposure                            #----
#---------------------------------------------------------------------#

#currently used https://gwas.mrcieu.ac.uk/datasets/ieu-a-1239/
education <- extract_instruments("ieu-a-1239") #23 gwas significant snps


#---------------------------------------------------------------------#
#                         Combined Exposures                          #----
#---------------------------------------------------------------------#
#using previous Locke et al
exp_dat <- extract_instruments("ieu-b-40") #79 gwas significant snps

snp_list <- smartbind(exp_dat, education)

write.table(unique(snp_list$SNP), "./MVMR_rsids.txt", row.names = F, col.names = F, quote = F)
x <- read.table("./MVMR_rsids.txt", header=F)

#---------------------------------------------------------------------#
#                   rsIDs for MVMR analysis                           #----
#---------------------------------------------------------------------#
#read exposure

rsid_func <- function(bmi_var) {
  exp_dat <- extract_instruments("ieu-b-40")
  bmi_exp_dat <- exp_dat[which(exp_dat$id.exposure==bmi_var),]
  
  education_exp_dat <- extract_instruments("ieu-a-1239")
  
  x <- smartbind(bmi_exp_dat, education_exp_dat)
  
  exp_snps <- extract_outcome_data(
    snps = x$SNP,
    outcomes = c(bmi_var)
  )
  
  clumped <- clump_data(exp_snps, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1)
}

bmi_snps <- rsid_func("ieu-b-40")

#---------------------------------------------------------------------#
#                   Read in Covid outcome datasets                    #----
#---------------------------------------------------------------------#

## Extract Covid data: Covid vs lab or self-reported negative - release 5
C1 = fread("./COVID19_HGI_C1_ALL_leave_23andme_20201020.b37.txt.gz")


#---------------------------------------------------------------------#
#                   MVMR exposures and outcomes                       #----
#---------------------------------------------------------------------#

mvmr_prep_func <- function(snps_data, bmi_var, myoutcome, outcome_name) {
#BMI exposure
bmi_MVMR <- extract_outcome_data(
  snps = snps_data,
  outcomes = c(bmi_var)
)


bmi_MVMR <- split_outcome(bmi_MVMR)
bmi_MVMR <- format_data(bmi_MVMR, type = "exposure", 
                        phenotype_col = "outcome",
                        snp_col = "SNP",
                        beta_col = "beta.outcome", 
                        se_col = "se.outcome", 
                        eaf_col = "eaf.outcome", 
                        effect_allele_col = "effect_allele.outcome", 
                        other_allele_col = "other_allele.outcome", 
                        pval_col = "pval.outcome", 
                        samplesize_col = "samplesize.outcome", 
                        chr_col = "chr", 
                        pos_col = "pos")  

#Smoking exposure
education_MVMR <- extract_outcome_data(
  snps = snps_data,
  outcomes = c("ieu-a-1239")
)
education_MVMR <- split_outcome(education_MVMR)
education_MVMR <- format_data(education_MVMR, type = "exposure", 
                            phenotype_col = "outcome",
                            snp_col = "SNP",
                            beta_col = "beta.outcome", 
                            se_col = "se.outcome", 
                            eaf_col = "eaf.outcome", 
                            effect_allele_col = "effect_allele.outcome", 
                            other_allele_col = "other_allele.outcome", 
                            pval_col = "pval.outcome", 
                            samplesize_col = "samplesize.outcome", 
                            chr_col = "chr", 
                            pos_col = "pos")

#covid outcome
out_dat<-myoutcome
out_dat$outcome <- outcome_name
out_dat<- rename(out_dat, c("SNP"="rsid",
                               "chr.outcome"="#CHR", 
                               "effect_allele.outcome"="ALT",
                               "other_allele.outcome"="REF",
                               "beta.outcome"="all_inv_var_meta_beta",
                               "se.outcome"="all_inv_var_meta_sebeta",
                               "pval.outcome"="all_inv_var_meta_p",                              
                               "eaf.outcome" = "all_meta_AF",
                               "samplesize.outcome" = "all_meta_sample_N",
                               "SNPid"="SNP",
                               "pval.het"="all_inv_var_het_p"))
out_dat$id.outcome <- rep(outcome_name, nrow(out_dat))


#harmonise
MVMR_data1 <- harmonise_data(
  exposure_dat = bmi_MVMR, 
  outcome_dat = out_dat
)     
MVMR_data2 <- harmonise_data(
  exposure_dat = education_MVMR, 
  outcome_dat = out_dat
)      

data_func <- function(MVMR_DATA) {
  x <- subset(MVMR_DATA, select = c("SNP", "outcome", "exposure","beta.exposure","beta.outcome","se.exposure",
                                    "se.outcome", "pval.exposure", "pval.outcome")) 
}

MVMR_data <- merge(data_func(MVMR_data1), data_func(MVMR_data2), by = c("SNP", "outcome", "beta.outcome", "se.outcome", "pval.outcome"))
}

#Written as a function, assuming each covid outcome dataset is already read in

bmi_education_C1_data <- mvmr_prep_func(bmi_snps$SNP, "ieu-b-40", C1, "C1 Covid vs lab or self-reported negative")


#---------------------------------------------------------------------#
#                           MVMR analysis                             #----
#---------------------------------------------------------------------#

#MVMR function
mvmr_func <- function(data, outcome) {
  x <- format_mvmr(BXGs = data[,c("beta.exposure.x", "beta.exposure.y")], BYG = data[,"beta.outcome"], seBXGs = data[,c("se.exposure.x", "se.exposure.y")], seBYG = data[,"se.outcome"], RSID = data[,"SNP"])
  a <- as.data.frame(t(strength_mvmr(x, gencov = 0)), row.names = c(paste0(data$exposure.x[1], "_MVMR", sep = ""), paste0(data$exposure.y[1], "_MVMR", sep = "")))
  y <- as.data.frame(ivw_mvmr(x, gencov = 0), row.names = c(paste0(data$exposure.x[1], "_MVMR", sep = ""), paste0(data$exposure.y[1], "_MVMR", sep = "")))
  y$outcome <- outcome
  y$exposure.1 <- data$exposure.x[1]
  y$exposure.2 <- data$exposure.y[1]
  m <- cbind(y, a)
  return(m)
}

#C1
mvmr_bmi_education_C1_result <- mvmr_func(bmi_education_C1_data, "C1 Covid vs lab or self-reported negative")


#put all outcomes together
results <- rbind(mvmr_bmi_education_C1_result)

write.xlsx(results, "./results/MVMR_EDUCATION_C1.xlsx", row.names= T)
