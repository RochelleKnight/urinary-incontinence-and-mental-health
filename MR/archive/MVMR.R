library(TwoSampleMR)
library(MRInstruments)
library(LDlinkR)
library(dplyr)
library(ieugwasr)
library(data.table)

# Extract BMI exposure data
bmi_exposure <- extract_instruments("ieu-a-974") #37 gwas significant snps
bmi_exposure <- bmi_exposure %>% rename("rsid" = "SNP", "pval"="pval.exposure")

# Load MDD exposure data
MDD_exposure <- read.csv(paste0("data/exposure/exposure_MDD_p_value_5e-08_outcome_severe_UUI_with_proxies.csv"))
MDD_exposure <- MDD_exposure %>% filter(rsid != "rs7531118")

# Combine SNPs and clump
exposure_snps <- rbind(MDD_exposure[c("rsid","pval")],bmi_exposure[c("rsid","pval")])
exposure_snps %>% filter(rsid == "rs7531118")

exposure_snps <-  ld_clump(exposure_snps,clump_kb=10000, clump_r2=0.001, pop = "EUR") # 1 SNP removed due to high LD with other SNPs

######################### MVMR exposure data ###################################

# Extract BMI data for MVMR
bmi_MVMR_exposure <- extract_outcome_data(
  snps = exposure_snps$rsid,
  outcomes = "ieu-a-974")

# Format BMI data ready for harmonisation
bmi_MVMR_exposure <- format_data(bmi_MVMR_exposure, type = "exposure",
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
                              

# Extract MDD data for MVMR
MDD_MVMR_exposure <- fread(file="data/raw_data/MDD2018_ex23andMe.txt",select=c("SNP","A1","A2","OR","SE","P","FRQ_A_59851","FRQ_U_113154"), fill = T)
MDD_MVMR_exposure$beta <- log(MDD_MVMR_exposure$OR)
MDD_MVMR_exposure$eaf <- (MDD_MVMR_exposure$FRQ_A_59851*59851 + MDD_MVMR_exposure$FRQ_U_113154*113154)/(59851+113154)
MDD_MVMR_exposure <- MDD_MVMR_exposure %>% select(SNP,eaf, A1,A2,beta, SE, P) %>%
  rename(effect_allele = A1,
         other_allele = A2,
         se = SE)
MDD_MVMR_exposure <- MDD_MVMR_exposure %>% filter(SNP %in% exposure_snps$rsid)
MDD_MVMR_exposure$phenotype <- "MDD"
length(unique(MDD_MVMR_exposure$SNP))

# Format data ready for harmonisation
MDD_MVMR_exposure <- format_data(MDD_MVMR_exposure,
                        type = "exposure",
                        phenotype_col = "phenotype",
                        snp_col = "SNP",
                        beta_col = "beta", 
                        se_col = "se", 
                        eaf_col = "eaf", 
                        effect_allele_col = "effect_allele", 
                        other_allele_col = "other_allele", 
                        pval_col = "P")

######################### MVMR outcome data ####################################

# Extract UUI data for MVMR
UUI_MVMR_outcome <- fread(file=paste0("data/raw_data/severe_UUI.csv.gz"), select=c("MarkerName","Allele1","Allele2","Effect","StdErr","P-value","MAF","issue_with_alleles"), fill = T)
UUI_MVMR_outcome <- UUI_MVMR_outcome %>% rename("SNP" = "MarkerName",
                              "effect_allele" = "Allele1",
                              "other_allele" = "Allele2",
                              "beta" = "Effect",
                              "se" = "StdErr",
                              "P" = "P-value",
                              "eaf" = "MAF")

UUI_MVMR_outcome$phenotype <- "severe_UUI"
UUI_MVMR_outcome <- UUI_MVMR_outcome %>% filter(SNP %in% exposure_snps$rsid)
UUI_MVMR_outcome <- UUI_MVMR_outcome %>% filter(issue_with_alleles != "yes")

UUI_MVMR_outcome <- format_data(UUI_MVMR_outcome,
                                 type = "outcome",
                                 phenotype_col = "phenotype",
                                 snp_col = "SNP",
                                 beta_col = "beta", 
                                 se_col = "se", 
                                 eaf_col = "eaf", 
                                 effect_allele_col = "effect_allele", 
                                 other_allele_col = "other_allele", 
                                 pval_col = "P")

################################# Harmonise data ###############################

MVMR_data1 <- harmonise_data(
  exposure_dat = MDD_MVMR_exposure, 
  outcome_dat = UUI_MVMR_outcome
)
MVMR_data1 <- MVMR_data1 %>% select("SNP", "outcome", "exposure","beta.exposure","beta.outcome","se.exposure",
                                     "se.outcome", "pval.exposure", "pval.outcome")

MVMR_data2 <- harmonise_data(
  exposure_dat = bmi_MVMR_exposure, 
  outcome_dat = UUI_MVMR_outcome
)
MVMR_data2 <- MVMR_data2 %>% select("SNP", "outcome", "exposure","beta.exposure","beta.outcome","se.exposure",
                                    "se.outcome", "pval.exposure", "pval.outcome")

MVMR_data <- merge(MVMR_data1, MVMR_data2, by = c("SNP", "outcome", "beta.outcome", "se.outcome", "pval.outcome"))

############################### Format for MVMR ################################

x <- format_mvmr(BXGs = MVMR_data[,c("beta.exposure.x", "beta.exposure.y")], BYG = MVMR_data[,"beta.outcome"], seBXGs = MVMR_data[,c("se.exposure.x", "se.exposure.y")], seBYG = MVMR_data[,"se.outcome"], RSID = MVMR_data[,"SNP"])

####################### Test for weak instruments ##############################

sres <- strength_mvmr(r_input = x, gencov = 0)

# Weak instrument strength for MDD - could need to use less BMI SNPs?

pres <- pleiotropy_mvmr(r_input = x, gencov = 0)

########################## Estimate IVW estimates ##############################

# exposure1 is MDD
# exposure2 is BMI

ivw_mvmr(r_input = x, gencov = 0)
exp(-1.4554330)
exp(-0.4216059)
