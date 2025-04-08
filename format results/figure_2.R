library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)

rm(list = ls())
dir.create("results/forest_plots", showWarnings = FALSE)

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
   mr_results_stress_urinary_incontinence,mr_results_urgency_urinary_incontinence,mr_results_anxiety)

mr_results_neuroticism <- mr_results %>% filter(method %in% c("Debiased IVW")
                                                & outcome == "neuroticism") %>%
  select(exposure,outcome,method,b,se)

mr_results <- mr_results %>% filter(method %in% c("Debiased IVW")
                                    & !outcome %in% c("neuroticism")) %>%
  select(exposure,outcome,method,b,se)

mr_results$odds_ratio <- exp(mr_results$b)
mr_results$lower_CI <- exp(mr_results$b - qnorm(0.975)*mr_results$se)
mr_results$upper_CI <- exp(mr_results$b + qnorm(0.975)*mr_results$se)

mr_results[c("b","se")] <- NULL

mr_results$est <- paste0(sprintf("%.2f",mr_results$odds_ratio)," (",sprintf("%.2f",mr_results$lower_CI),"-",sprintf("%.2f",mr_results$upper_CI),")")


# Format exposure name
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


mr_results$outcome <- ifelse(mr_results$exposure == "Urinary incontinence", paste0(strrep(" ",2), mr_results$outcome,strrep(" ",1)),mr_results$outcome)
mr_results$outcome <- ifelse(mr_results$exposure == "Stress urinary incontinence", paste0(strrep(" ",2),mr_results$outcome,strrep(" ",2)),mr_results$outcome)
mr_results$outcome <- ifelse(mr_results$exposure == "Urgency urinary incontinence", paste0(strrep(" ",2),mr_results$outcome,strrep(" ",3)),mr_results$outcome)
mr_results$outcome <- ifelse(mr_results$exposure == "Mixed urinary incontinence", paste0(strrep(" ",2),mr_results$outcome,strrep(" ",4)),mr_results$outcome)
mr_results$outcome <- ifelse(mr_results$exposure == "Anxiety", paste0(strrep(" ",2),mr_results$outcome,strrep(" ",5)),mr_results$outcome)

mr_results[nrow(mr_results)+1,] <- c("Urinary incontinence →","Urinary incontinence →", NA,NA,NA,NA,NA)
mr_results[nrow(mr_results)+1,] <- c("Stress urinary incontinence →","Stress urinary incontinence →", NA,NA,NA,NA,NA)
mr_results[nrow(mr_results)+1,] <- c("Urgency urinary incontinence →","Urgency urinary incontinence →", NA,NA,NA,NA,NA)
mr_results[nrow(mr_results)+1,] <- c("Mixed urinary incontinence →","Mixed urinary incontinence →", NA,NA,NA,NA,NA)
mr_results[nrow(mr_results)+1,] <- c("Anxiety →","Anxiety →",NA,NA,NA,NA,NA)

mr_results <- mr_results %>% mutate(across(c(odds_ratio,upper_CI,lower_CI), as.numeric))

mr_results$outcome <- factor(mr_results$outcome, levels = c("Urinary incontinence →",
                                                            paste0(strrep(" ",2),"Anxiety", strrep(" ",1)),
                                                            paste0(strrep(" ",2),"Depression",strrep(" ",1)),
                                                            paste0(strrep(" ",2),"Broad depression phenotype",strrep(" ",1)),
                                                            "Stress urinary incontinence →",
                                                            paste0(strrep(" ",2),"Anxiety", strrep(" ",2)),
                                                            paste0(strrep(" ",2),"Depression",strrep(" ",2)),
                                                            paste0(strrep(" ",2),"Broad depression phenotype",strrep(" ",2)),
                                                            "Urgency urinary incontinence →",
                                                            paste0(strrep(" ",2),"Anxiety", strrep(" ",3)),
                                                            paste0(strrep(" ",2),"Depression",strrep(" ",3)),
                                                            paste0(strrep(" ",2),"Broad depression phenotype",strrep(" ",3)),
                                                            "Mixed urinary incontinence →",
                                                            paste0(strrep(" ",2),"Anxiety", strrep(" ",4)),
                                                            paste0(strrep(" ",2),"Depression",strrep(" ",4)),
                                                            paste0(strrep(" ",2),"Broad depression phenotype",strrep(" ",4)),
                                                            "Anxiety →",
                                                            paste0(strrep(" ",2),"Urinary incontinence",strrep(" ",5)),
                                                            paste0(strrep(" ",2),"Stress urinary incontinence",strrep(" ",5)),
                                                            paste0(strrep(" ",2),"Urgency urinary incontinence",strrep(" ",5)),
                                                            paste0(strrep(" ",2),"Mixed urinary incontinence",strrep(" ",5))))
                                                            
                                                            
mr_results <- mr_results[order(mr_results$outcome),]

mr_results$colour <- c(rep(c("gray95","white"), 10,),"gray95")

p <- ggplot(mr_results, aes(x = odds_ratio, y = outcome, xmin = lower_CI, xmax = upper_CI, fill = outcome, shape = method)) +
  #geom_hline(aes(yintercept = outcome, colour = colour), linewidth = 12) +
  # geom_hline(aes(yintercept = order, colour = colour)) +
  geom_pointrange(shape = c(rep(15,21)), fill = "black") +
  geom_vline(xintercept = 1, linetype = 3) +
  xlab("Odds ratio") +
  ylab("Adjusted Relative Risk with 95% Confidence Interval") +
  theme_classic() +
  scale_colour_identity() +
  scale_y_discrete(limits = rev(mr_results$outcome)) +
  scale_x_continuous(lim = c(0,2.5), breaks = c(0, 0.5,1,1.5,2,2.5))+
  # scale_x_log10(limits = c(0.25, 4),
  #              breaks = c(0.25, 0.5, 1, 2, 4),
  #              labels = c("0.25", "0.5", "1", "2", "4"), expand = c(0,0)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())
p


table_results <- mr_results %>% select(outcome,est,colour)
table_results$outcome <- as.character(table_results$outcome)

table_results$outcome <- factor(table_results$outcome, levels = rev(c("Urinary incontinence →",
                                                                      paste0(strrep(" ",2),"Anxiety", strrep(" ",1)),
                                                                      paste0(strrep(" ",2),"Depression",strrep(" ",1)),
                                                                      paste0(strrep(" ",2),"Broad depression phenotype",strrep(" ",1)),
                                                                      "Stress urinary incontinence →",
                                                                      paste0(strrep(" ",2),"Anxiety", strrep(" ",2)),
                                                                      paste0(strrep(" ",2),"Depression",strrep(" ",2)),
                                                                      paste0(strrep(" ",2),"Broad depression phenotype",strrep(" ",2)),
                                                                      "Urgency urinary incontinence →",
                                                                      paste0(strrep(" ",2),"Anxiety", strrep(" ",3)),
                                                                      paste0(strrep(" ",2),"Depression",strrep(" ",3)),
                                                                      paste0(strrep(" ",2),"Broad depression phenotype",strrep(" ",3)),
                                                                      "Mixed urinary incontinence →",
                                                                      paste0(strrep(" ",2),"Anxiety", strrep(" ",4)),
                                                                      paste0(strrep(" ",2),"Depression",strrep(" ",4)),
                                                                      paste0(strrep(" ",2),"Broad depression phenotype",strrep(" ",4)),
                                                                      "Anxiety →",
                                                                      paste0(strrep(" ",2),"Urinary incontinence",strrep(" ",5)),
                                                                      paste0(strrep(" ",2),"Stress urinary incontinence",strrep(" ",5)),
                                                                      paste0(strrep(" ",2),"Urgency urinary incontinence",strrep(" ",5)),
                                                                      paste0(strrep(" ",2),"Mixed urinary incontinence",strrep(" ",5)))))

table_results <- table_results[order(table_results$outcome),]

data_table <- ggplot(data = table_results, aes(y = outcome, fill = outcome)) +
  #geom_hline(aes(yintercept = outcome, colour = colour), linewidth =12) +
  geom_text(aes(x = 0, label = outcome), hjust = 0) +
  #geom_text(aes(x = 1, label = method), hjust = 0) +
  geom_text(aes(x = 3, label = est),hjust = 1) +
  scale_colour_identity() +
  theme_void()+
  theme(plot.margin = margin(5, 0, 35, 0))
data_table

top <- grid.arrange(data_table,p, ncol = 2, top = textGrob(paste0("(A)",strrep(" ",60),"OR (95% CI)",strrep(" ",90)),gp=gpar(fontsize=12,font=3)))


# Plot figure for neuroticism outcome
mr_results_neuroticism$lower_CI <- mr_results_neuroticism$b - qnorm(0.975)*mr_results_neuroticism$se
mr_results_neuroticism$upper_CI <- mr_results_neuroticism$b + qnorm(0.975)*mr_results_neuroticism$se

mr_results_neuroticism$est <- NA
mr_results_neuroticism$est[1] <- paste0(sprintf("%.3f",mr_results_neuroticism$b[1])," (",sprintf("%.4f",mr_results_neuroticism$lower_CI[1]),"-",sprintf("%.2f",mr_results_neuroticism$upper_CI[1]),")")
mr_results_neuroticism$est[2:4] <- paste0(sprintf("%.3f",mr_results_neuroticism$b[2:4])," (",sprintf("%.2f",mr_results_neuroticism$lower_CI[2:4]),"-",sprintf("%.2f",mr_results_neuroticism$upper_CI[2:4]),")")

# Format exposure name
mr_results_neuroticism <- mr_results_neuroticism %>% dplyr::mutate(exposure = dplyr::case_when(exposure == "urinary_incontinence_p_value_5e-06" ~ "Urinary incontinence",
                                                                       exposure == "stress_urinary_incontinence_p_value_5e-06" ~ "Stress urinary incontinence",
                                                                       exposure == "urgency_urinary_incontinence_p_value_5e-06" ~ "Urgency urinary incontinence",
                                                                       exposure == "mixed_urinary_incontinence_p_value_5e-06" ~ "Mixed urinary incontinence",
                                                                       TRUE ~ exposure))
mr_results_neuroticism$outcome <- "Neuroticism"


mr_results_neuroticism$outcome <- ifelse(mr_results_neuroticism$exposure == "Urinary incontinence", paste0(strrep(" ",2), mr_results_neuroticism$outcome,strrep(" ",1)),mr_results_neuroticism$outcome)
mr_results_neuroticism$outcome <- ifelse(mr_results_neuroticism$exposure == "Stress urinary incontinence", paste0(strrep(" ",2),mr_results_neuroticism$outcome,strrep(" ",2)),mr_results_neuroticism$outcome)
mr_results_neuroticism$outcome <- ifelse(mr_results_neuroticism$exposure == "Urgency urinary incontinence", paste0(strrep(" ",2),mr_results_neuroticism$outcome,strrep(" ",3)),mr_results_neuroticism$outcome)
mr_results_neuroticism$outcome <- ifelse(mr_results_neuroticism$exposure == "Mixed urinary incontinence", paste0(strrep(" ",2),mr_results_neuroticism$outcome,strrep(" ",4)),mr_results_neuroticism$outcome)

mr_results_neuroticism[nrow(mr_results_neuroticism)+1,] <- c("Urinary incontinence →","Urinary incontinence →", NA,NA,NA,NA,NA,NA)
mr_results_neuroticism[nrow(mr_results_neuroticism)+1,] <- c("Stress urinary incontinence →","Stress urinary incontinence →", NA,NA,NA,NA,NA,NA)
mr_results_neuroticism[nrow(mr_results_neuroticism)+1,] <- c("Urgency urinary incontinence →","Urgency urinary incontinence →", NA,NA,NA,NA,NA,NA)
mr_results_neuroticism[nrow(mr_results_neuroticism)+1,] <- c("Mixed urinary incontinence →","Mixed urinary incontinence →", NA,NA,NA,NA,NA,NA)

mr_results_neuroticism <- mr_results_neuroticism %>% mutate(across(c(b,upper_CI,lower_CI), as.numeric))

mr_results_neuroticism$outcome <- factor(mr_results_neuroticism$outcome, levels = c("Urinary incontinence →",
                                                                                    paste0(strrep(" ",2),"Neuroticism",strrep(" ",1)),
                                                                                    "Stress urinary incontinence →",
                                                                                    paste0(strrep(" ",2),"Neuroticism",strrep(" ",2)),
                                                                                    "Urgency urinary incontinence →",
                                                                                    paste0(strrep(" ",2),"Neuroticism",strrep(" ",3)),
                                                                                    "Mixed urinary incontinence →",
                                                                                    paste0(strrep(" ",2),"Neuroticism",strrep(" ",4))))

mr_results_neuroticism <- mr_results_neuroticism[order(mr_results_neuroticism$outcome),]

mr_results_neuroticism$colour <- rep(c("gray95","white"), 4)

p_neuroticism <- ggplot(mr_results_neuroticism, aes(x = b, y = outcome, xmin = lower_CI, xmax = upper_CI, fill = outcome)) +
  #geom_hline(aes(yintercept = outcome, colour = colour), linewidth = 15) +
  # geom_hline(aes(yintercept = order, colour = colour)) +
  geom_pointrange(shape = rep(15,8), fill = "black") +
  geom_vline(xintercept = 0, linetype = 3) +
  xlab("\u03b2") +
  ylab("Adjusted Relative Risk with 95% Confidence Interval") +
  theme_classic() +
  scale_colour_identity() +
  scale_y_discrete(limits = rev(mr_results_neuroticism$outcome)) +
  #scale_x_continuous(lim = c(-0.01,0.01), breaks = c(-0.01,0,0.01))+
  # scale_x_log10(limits = c(0.25, 4),
  #              breaks = c(0.25, 0.5, 1, 2, 4),
  #              labels = c("0.25", "0.5", "1", "2", "4"), expand = c(0,0)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())
p_neuroticism


table_results_neuroticism <- mr_results_neuroticism %>% select(outcome,est,colour)

table_results_neuroticism$outcome <- factor(table_results_neuroticism$outcome, levels = rev(c("Urinary incontinence →",
                                                                                              paste0(strrep(" ",2),"Neuroticism",strrep(" ",1)),
                                                                                              "Stress urinary incontinence →",
                                                                                              paste0(strrep(" ",2),"Neuroticism",strrep(" ",2)),
                                                                                              "Urgency urinary incontinence →",
                                                                                              paste0(strrep(" ",2),"Neuroticism",strrep(" ",3)),
                                                                                              "Mixed urinary incontinence →",
                                                                                              paste0(strrep(" ",2),"Neuroticism",strrep(" ",4)))))
table_results_neuroticism <- table_results_neuroticism[order(table_results_neuroticism$outcome),]

data_table_neuroticism <- ggplot(data = table_results_neuroticism, aes(y = outcome, fill = outcome)) +
  #geom_hline(aes(yintercept = outcome, colour = colour), linewidth =15) +
  geom_text(aes(x = 0, label = outcome), hjust = 0) +
  #geom_text(aes(x = 1, label = method), hjust = 0) +
  geom_text(aes(x = 3, label = est),hjust = 1) +
  scale_colour_identity() +
  theme_void()+
  theme(plot.margin = margin(5, 0, 35, 0))
data_table_neuroticism

bottom = grid.arrange(data_table_neuroticism,p_neuroticism, ncol = 2,top = textGrob(paste0("(B)",strrep(" ",60),"\u03b2 (95% CI)",strrep(" ",90)),gp=gpar(fontsize=12,font=3)))

grid.arrange(top,bottom, ncol = 1)

png("results/forest_plots/figure_2.png", units = "mm", width=210, height=297,res = 1000)
print(plot_grid(top, bottom,align = "v", nrow = 2, rel_heights = c(0.75, 0.25)))
dev.off() 

