#install.packages("BiocManager")
#BiocManager::install("BiocGenerics")
#BiocManager::install("Biostrings")
#BiocManager::install("GenomicRanges")
#BiocManager::install("Rsamtools")
#BiocManager::install("SummarizedExperiment")
#BiocManager::install("VariantAnnotation")
#remotes::install_github("mrcieu/gwasvcf")
#install.packages("LDlinkR")

library(data.table)
library(dplyr)
library(LDlinkR)
library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)
library(plinkbinr)
#library(gwasvcf)
# Prep all instruments
# Find all proxy SNPs - 

setwd("~/GitHub/urinary-incontinence-and-mental-health")

source("preprocess/function_format_UI_exposures.R")
source("preprocess/check_gwas_data.R")

############################
# Urinary Incontinence
############################
exposure <- read.table("data/raw_data/urinary_incontinence.tsv",header = T)
exposure <- exposure %>% rename("SNP" = "variant_id",
                                "se" = "standard_error",
                                "P" = "p_value") 

check_gwas_data(exposure)
format_UI_exposure_data(exposure, "urinary_incontinence")

############################
# Mixed Urinary Incontinence
############################
exposure <- read.table("data/raw_data/mixed_urinary_incontinence.tsv",header = T)
exposure <- exposure %>% rename("SNP" = "variant_id",
                                "se" = "standard_error",
                                "P" = "p_value") 

check_gwas_data(exposure)
format_UI_exposure_data(exposure, "mixed_urinary_incontinence")

 ############################
# Stress Urinary Incontinence
############################
exposure <- read.table("data/raw_data/stress_urinary_incontinence.tsv",header = T)
exposure <- exposure %>% rename("SNP" = "variant_id",
                                "se" = "standard_error",
                                "P" = "p_value") 

check_gwas_data(exposure)
format_UI_exposure_data(exposure, "stress_urinary_incontinence")

############################
# Urgency Urinary Incontinence
############################
exposure <- read.table("data/raw_data/urgency_urinary_incontinence.tsv",header = T)
exposure <- exposure %>% rename("SNP" = "variant_id",
                                "se" = "standard_error",
                                "P" = "p_value") 

check_gwas_data(exposure)
format_UI_exposure_data(exposure, "urgency_urinary_incontinence")

