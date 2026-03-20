################################################################################
# Figure S3.6 and Table S3.7 in the chapter supplement 
# Cut-off identification for hemoglobin and overall survival [OS]

# This script creates a four-panel figure looking at  (A) the distribution of 
# hemoglobin, (B) its association of hemoglobin with OS (HR), 
# (C) maximal log-rank test statistic, and (D) Kaplan-Meier curves and HR 
# by a selected cut-off for hemoglobin. It also 
# identifies an optimal cutoff using multiple statistical methods (LS1,
# LS2, CO, CPE, PLRT, YI, ER, and CZ) - see more details in the main chapter and the
# chapter supplement (Table S3.6)

# Note: All analyses were conducted using R version 4.5.1. The packages and functions 
# described in this chapter and its supplement were available and verified at the 
# time of writing. Future updates to these packages may cause code examples to behave differently.
################################################################################

#--------------------------------------------------------
# Load packages
#--------------------------------------------------------

library(tidyverse) #a collection of core R packages for data science
library(survival) #for survival analysis
library(boot) # for bootstrapping and 95% bootstrap CIs
library(rms) #Regression Modeling Strategies: using rcs() to calculate basis terms for restricted cubic splines
library(CPE) #phcpe to calculate CPE
library(timeROC) #time dependent ROC, AUC
library(SurvEval) #calculate CPE_projection 
library(maxstat)# Cut-point identification using maximal statistics
library(rms) # for Cox models
library(rmsMD) # for HR curves
library(ggplot2)
library(survminer) # Survival analysis
library(survMisc) # Contal and O'Quigley method
library(gridExtra)
library(boot) # boostrap resampling
library(survivalROC) # time-dependent ROCs

#--------------------------------------------------------
#Load data
#--------------------------------------------------------
# Read the simulated data set (.csv) from GitHub repository (n=300)
file_url <- "https://raw.githubusercontent.com/chapman-biomarker-book-desachy-dangelo/chapter3_prognostic-biomarkers/main/simulated_hnca_300.csv"

d1 = read_csv(file_url) %>%
  # Convert binary variables of interest to factors for modeling and plotting
  mutate(across(c("zubrod", "t_stage_bin", "n_stage_bin", "p16_status",  "egfr_level_80", "ajcc_stage"),
                ~as.factor(.x), .names="{.col}_F")) %>%
  # Set p16_status reference group to positive (i.e., compares neg vs pos)
  mutate(p16_status_F=relevel(p16_status_F, ref="1"))
         
str(d1)


#--------------------------------------------------------
# Panel A in Figure S3.6: Histogram of Hemoglobin values
#--------------------------------------------------------

hemo_histo <- ggplot(d1) + aes(x=hemoglobin) + geom_histogram() +
  labs(x="Hemoglobin (g/dl)", y="Number of Patients") + ggtitle("A")
hemo_histo

#--------------------------------------------------------
# Panel B in Figure S3.6: Smooth curve for OS HR vs. Hemoglobin
#                         based on rcs with 4 knots
#--------------------------------------------------------

# Define the reference value
reference_hemoglobin <- median(d1$hemoglobin) # reference hemoglobin for HR curve
dd <- datadist(d1)
options(datadist="dd")
dd$limits$hemoglobin[2] <- reference_hemoglobin

# Fit the Cox model with a restricted cubic spline
fit.hemog <- cph(Surv(survival_years, survival) ~ rcs(hemoglobin, 4), data=d1, x=TRUE, y=TRUE)

# Create the plot (Figure S3.6, Panel B)
hemo_HR_curve <- ggrmsMD(fit.hemog, d1,
                   xlabs=list(hemoglobin="Hemoglobin (g/dl)"), titles=list(hemoglobin="B") 
)
hemo_HR_curve

#--------------------------------------------------------
# Table S3.7: LS1--Lau92 (Lausen and Schumacher, 1992)
# Statistic (p-value), used in Panel C figure
#--------------------------------------------------------

cutp.maxstat.lau92 <- maxstat.test(Surv(survival_years, survival) ~ hemoglobin,
                                   smethod="LogRank",
                                   data = d1, pmethod = "Lau92",
                                   alpha = 0.05
)
cutp.maxstat.lau92

#--------------------------------------------------------
# Panel C in Figure S3.6: Standardized log-rank test for
#                         potential hemoglobin cutoffs (LS1/LS2)
#--------------------------------------------------------

cutp.maxstat.lau92.df <- as.data.frame(cbind(cutp.maxstat.lau92$stats, cutp.maxstat.lau92$cuts))
colnames(cutp.maxstat.lau92.df) <- c("stats", "cuts")
plot_lau92 <- ggplot(cutp.maxstat.lau92.df) + geom_point(aes(y=stats, x=cuts)) + geom_line(aes(y=stats, x=cuts))+
  geom_vline(xintercept=13.1, col="blue", linetype="dashed") +
  labs(x="Hemoglobin (g/dl)", y="Standardized Log-Rank Statistic")+
  ggtitle("C")
plot_lau92

#--------------------------------------------------------
# Panel D in Figure S3.6: Kaplan-Meier for low and high-risk subgroups
#--------------------------------------------------------

d1$hemoglobin_ <- cut(d1$hemoglobin,
                      breaks = c(0, 13.1, max(d1$hemoglobin)),
                      labels = c("Low (<=13.1)", "High (>13.1)"),
                      right = TRUE, include.lowest = TRUE)
d1$hemoglobin_ <- factor(d1$hemoglobin_, levels=c("High (>13.1)", "Low (<=13.1)"))
table(d1$hemoglobin_) 
length(which(d1$hemoglobin <= 13.1)) # double-checking 

fit.hemoglobin.cat <- survfit(Surv(survival_years, survival) ~ hemoglobin_, data = d1)
fit.hemoglobin.cat

# Fit a Cox proportional hazards model to get the HR and 95% CI
# This model is required to calculate the Hazard Ratio and its confidence interval.
cox_fit <- coxph(Surv(survival_years, survival) ~ hemoglobin_, data = d1)
# Extract HR and CI
hr_data <- broom::tidy(cox_fit, conf.int = TRUE, exponentiate = TRUE)
hr_text <- sprintf("HR: %.2f (95%% CI: %.2f - %.2f)",
                   hr_data$estimate, hr_data$conf.low, hr_data$conf.high)

# Generate the plot object
p <- ggsurvplot(
  fit.hemoglobin.cat,
  data = d1,
  pval = TRUE,             # Add p-value
  conf.int = FALSE,         # Add confidence intervals for the curves
  risk.table = TRUE,       # Add risk table
  risk.table.col = "strata", # Color risk table by groups
  linetype = "strata",     # Change line type by groups
  ggtheme = theme_bw(),    # Use a clean theme
  palette = c("black", "gray"), # Custom colors
  title = "D"
)

# Add the Hazard Ratio and 95% CI as an annotation to the plot
# specific dataset's time scale and survival probability.
# Here, we place it in the bottom left. Time unit is in years (0-6yrs).
# e.g., to place it on the top right, change x=4, y=0.95
p$plot <- p$plot +
  ggplot2::annotate("text", x = 0, y = 0.10, label = hr_text, size = 4, hjust = 0)
print(p)

#--------------------------------------------------------
# Final Figure S3.6: Arrange Panels A/B/C/D and save
#--------------------------------------------------------

# Multipanel Figure S3.6
grid.arrange(hemo_histo, hemo_HR_curve,
             plot_lau92, p$plot,
             nrow = 2)
# Create a grob from the arranged plots
grob <- arrangeGrob(hemo_histo, hemo_HR_curve,
                    plot_lau92, p$plot,
                    nrow = 2)

# Use ggsave to save the grob
ggsave("Figure S3.6.png", grob, width = 10, height = 10, path=output)



#############################################################
### Comparisons of methods based on maximal test statistics
#############################################################

#--------------------------------------------------------
# Table S3.7: LS2--Lau94 (Lausen, Sauerbrei and Schumacher, 1994) 
# Statistic (p-value)
#--------------------------------------------------------

cutp.maxstat.lau94 <- maxstat.test(Surv(survival_years, survival) ~ hemoglobin,
                                   smethod="LogRank",
                                   data = d1, pmethod = "Lau94", 
                                   alpha = 0.05
)
cutp.maxstat.lau94


#--------------------------------------------------------
# Table S3.7: CO--Contal C, O'Quigley J, 1999
# Statistic (p-value)
#--------------------------------------------------------

cutp2 <- survMisc::cutp(coxph(Surv(survival_years, survival) ~ hemoglobin, data=d1), defCont = 3, plot=TRUE)
c1 <- cutp2$hemoglobin
data.table::setorder(c1, "p")
c1[1,]     

#--------------------------------------------------------
# Table S3.7: LS1/LS2/CO summary table
# Statistic (p-value)
#--------------------------------------------------------
all.results <- rbind(
  c(cutp.maxstat.lau92$method, round(cutp.maxstat.lau92$estimate, 2), round(cutp.maxstat.lau92$statistic, 4), round(cutp.maxstat.lau92$p.value, 8)),
  c(cutp.maxstat.lau94$method, round(cutp.maxstat.lau94$estimate, 2), round(cutp.maxstat.lau94$statistic, 4), round(cutp.maxstat.lau94$p.value, 8)),
  c("Contal-OQuingley", round(c1[1,"hemoglobin"] , 2), round(c1[1,"Q"] , 4), round(c1[1,"p"], 8))
)
colnames(all.results) <- c("Method", "Cut Point", "Statistic", "p-value")
all.results

#--------------------------------------------------------
# Table S3.7: Bootstrap CI for cut-off based on LS1
#--------------------------------------------------------

cutpoint.function <- function(x, index) {
  d <- x[index, ]
  cutp.maxstat.lau92 <- maxstat.test(Surv(survival_years, survival) ~ hemoglobin,
                                     smethod="LogRank",
                                     data = d, pmethod = "Lau92",
                                     alpha = 0.05
  )
  return(cutp.maxstat.lau92$estimate)  
}
BootDist_LS1 <- boot(data = d1, statistic = cutpoint.function, R=999) 
BootDist_LS1 
boot.ci( BootDist_LS1, conf = .95, type='all')
 
#--------------------------------------------------------
# Table S3.7: Bootstrap CI for cut-off based on LS2
#--------------------------------------------------------

cutpoint.function <- function(x, index) {
  d <- x[index, ]
  cutp.maxstat.lau94 <- maxstat.test(Surv(survival_years, survival) ~ hemoglobin,
                                     smethod="LogRank",
                                     data = d, pmethod = "Lau94",
                                     alpha = 0.05
  )
  return(cutp.maxstat.lau94$estimate)  
}
BootDist_LS2 <- boot(data = d1, statistic = cutpoint.function, R=999) 
BootDist_LS2
boot.ci( BootDist_LS2, conf = .95, type='all')
 

#--------------------------------------------------------
# Table S3.7: Bootstrap CI for cut-off based on CO
#--------------------------------------------------------

 cutpoint.function <- function(x, index) {
  d <- x[index, ]
  cutp2 <- survMisc::cutp(coxph(Surv(survival_years, survival) ~ hemoglobin, data=d), plot=FALSE)
  c1 <- cutp2$hemoglobin
  data.table::setorder(c1, "p")
  return(c1[1,][[1]])  
 }
set.seed(3)
BootDist_CO <- boot(data = d1, statistic = cutpoint.function, R=999) 
BootDist_CO
boot.ci( BootDist_CO, conf = .95, type='all')
 

#--------------------------------------------------------
# Table S3.7: Cut-off based on PLRT
#--------------------------------------------------------

cutp <- seq(quantile(d1$hemoglobin, 0.10), quantile(d1$hemoglobin, 0.90), by=0.1)
PLRT <- NULL
HRs <- NULL
for(i in 1:length(cutp)){
  d1$hemoglobin_cutp <- as.numeric(d1$hemoglobin > cutp[i])
  fit.cutp <- coxph(Surv(survival_years, survival) ~ hemoglobin_cutp, data=d1)
  PLRT[i] <- anova(fit.cutp)$Chisq[2]
  HRs[i] <- exp(-coefficients(fit.cutp))
}
plot(cutp, PLRT, type="b")

#PLRT statistic (p-val)
PLRT[which(PLRT==max(PLRT))]

#cutoff
opt_cutp_PLRT <- cutp[which(PLRT==max(PLRT))]
opt_cutp_PLRT
abline(v=opt_cutp_PLRT)

plot(cutp, HRs, type="b")
abline(v=opt_cutp_PLRT)

#--------------------------------------------------------
# Table S3.7: Boostrap CI for cut-off based on PLRT
#--------------------------------------------------------

set.seed(3)
opt_cutp_plrt <- NULL
for(b in 1:999){
  d <- d1[sample(1:nrow(d1), replace = TRUE), ]
  cutp <- seq(quantile(d$hemoglobin, 0.10), quantile(d$hemoglobin, 0.90), by=0.1)
  PLRT <- NULL
  for(i in 1:length(cutp)){
    d$hemoglobin_cutp <- as.numeric(d$hemoglobin > cutp[i])
    fit.cutp <- coxph(Surv(survival_years, survival) ~ hemoglobin_cutp, data=d)
    PLRT[i] <- anova(fit.cutp)$Chisq[2]
  }
  opt.cutp <- cutp[which.max(PLRT)] 
  opt_cutp_plrt <- c(opt_cutp_plrt, opt.cutp)  
}
hist(opt_cutp_plrt)
quantile(opt_cutp_plrt, c(0.025, 0.975))

#--------------------------------------------------------
# Table S3.7: Cut-off based on CPE
#--------------------------------------------------------

cutp <- seq(quantile(d1$hemoglobin, 0.10), quantile(d1$hemoglobin, 0.90), by=0.1)
CPE_index <- NULL
for(i in 1:length(cutp)){
  d1$hemoglobin_cutp <- as.numeric(d1$hemoglobin >= cutp[i])
  fit.cutp <- coxph(Surv(survival_years, survival) ~ hemoglobin_cutp, data=d1)
  CPE_index[i] <- phcpe(fit.cutp)$CPE
  #CPE_index[i] <- fit.cutp$concordance[6]
}
plot(cutp, CPE_index, type="b")

#CPE statistic
CPE_index[which(CPE_index==max(CPE_index))]

#CPE cutoff
opt_cutp_CPE <- cutp[which.max(CPE_index)]
opt_cutp_CPE
abline(v=opt_cutp_CPE)

#--------------------------------------------------------
# Table S3.7: Bootstrap CI for cut-off based on CPE
#--------------------------------------------------------

set.seed(3)
opt_cutp_CPE <- NULL
for(b in 1:999){
  d <- d1[sample(1:nrow(d1), replace = TRUE), ]
  cutp <- seq(quantile(d$hemoglobin, 0.10), quantile(d$hemoglobin, 0.90), by=0.1)
  CPE_index <- NULL
  for(i in 1:length(cutp)){
    d$hemoglobin_cutp <- as.numeric(d$hemoglobin > cutp[i])
    fit.cutp <- coxph(Surv(survival_years, survival) ~ hemoglobin_cutp, data=d)
    CPE_index[i] <- phcpe(fit.cutp)$CPE
  }
  opt_cutp <- cutp[which.max(CPE_index)]  
  opt_cutp_CPE <- c(opt_cutp_CPE, opt_cutp)  
}
length(opt_cutp_CPE)
quantile(opt_cutp_CPE, c(0.025, 0.975))

#--------------------------------------------------------
# Table S3.7: Bootstrap CI for cut-off based on CPE—a different approach
#--------------------------------------------------------

cutpoint.function <- function(x, index) {
  d <- x[index, ]
  cutp <- seq(quantile(d$hemoglobin, 0.10), quantile(d$hemoglobin, 0.90), by=0.1)
  CPE_index <- NULL
  for(i in 1:length(cutp)){
    d$hemoglobin_cutp <- as.numeric(d$hemoglobin > cutp[i])
    fit.cutp <- coxph(Surv(survival_years, survival) ~ hemoglobin_cutp, data=d)
    CPE_index[i] <- phcpe(fit.cutp)$CPE
  }
  opt_cutp <- cutp[which.max(CPE_index)][1]  
  return(opt_cutp)
 }
set.seed(3)
BootDist_CPE_boot <- boot(data = d1, statistic = cutpoint.function, R=999) 
BootDist_CPE_boot
boot.ci( BootDist_CPE_boot, conf = .95, type='all')
plot(BootDist_CPE_boot)


#--------------------------------------------------------
# Table S3.7: Cut-offs based on Youden (YI), 
#             Euclidean distance (ER), and 
#             Concordance Probability (CZ)
# for t=3 years (time-specific)
#--------------------------------------------------------

td_roc <- survivalROC.C(
  Stime = d1$survival_years,
  status = d1$survival,
  marker = -d1$hemoglobin, # since low hemoglobin is associated with worse OS
  predict.time = c(3), # Time point of interest in days
)
Sen_t3 <- td_roc$TP # sensitivity
Spe_t3 <- 1-td_roc$FP # specificity
YI <- Sen_t3 + Spe_t3 - 1 # Youden
ER <- sqrt((0-td_roc$FP)^2+(1-td_roc$TP)^2) # Euclidean distance
CZ <- Sen_t3*Spe_t3
plot(td_roc$cut.values, YI, ylim=c(0,1), type="l", col="blue")
lines(td_roc$cut.values, ER)
lines(td_roc$cut.values, CZ)

# Performance measures at the optimal cut-offs
YI[which.max(YI)]
ER[which.min(ER)]
CZ[which.max(CZ)]

# Optimal cut-offs based on youden (YI), ER, and CZ
td_roc$cut.values[which.max(YI)]
td_roc$cut.values[which.min(ER)]
td_roc$cut.values[which.max(CZ)]

#--------------------------------------------------------
# Table S3.7: Bootstrap CI for cut-off based on YI
#--------------------------------------------------------

cutpoint.function <- function(x, index) {
  d <- x[index, ]
  td_roc_cutp <- survivalROC.C(
      Stime = d$survival_years,
      status = d$survival,
      marker = -d$hemoglobin, # since low hemoglobin is associated with worse OS
      predict.time = c(3), # Time point of interest in days
    )
  Sen_t3 <- td_roc_cutp$TP
  Spe_t3 <- 1 - td_roc_cutp$FP
  YI <- Sen_t3 + Spe_t3 - 1
  opt_cutp <- td_roc_cutp$cut.values[which.max(YI)]
  return(opt_cutp)
}
set.seed(3)
BootDist_YI_boot <- boot(data = d1, statistic = cutpoint.function, R=999) 
BootDist_YI_boot
boot.ci( BootDist_YI_boot, conf = .95, type='all')
# hist(BootDist_YI_boot)
hist(as.numeric(BootDist_YI_boot$t)) #there are 9 -inf in the bootstrapped results, for demonstration purposes, plotting removing those 
plot(BootDist_YI_boot)

#--------------------------------------------------------
# Table S3.7: Bootstrap CI for cut-off based on ER
#--------------------------------------------------------

cutpoint.function <- function(x, index) {
  d <- x[index, ]
  td_roc_cutp <- survivalROC.C(
    Stime = d$survival_years,
    status = d$survival,
    marker = -d$hemoglobin, # since low hemoglobin is associated with worse OS
    predict.time = c(3), # Time point of interest in days
  )
  ER <- sqrt((0-td_roc_cutp$FP)^2+(1-td_roc_cutp$TP)^2)
  opt_cutp <- td_roc_cutp$cut.values[which.min(ER)]
  return(opt_cutp)
}
set.seed(3)
BootDist_ER_boot <- boot(data = d1, statistic = cutpoint.function, R=999) 
BootDist_ER_boot
boot.ci( BootDist_ER_boot, conf = .95, type='all')
hist(BootDist_ER_boot)
plot(BootDist_ER_boot)

#--------------------------------------------------------
# Table S3.7: Bootstrap CI for cut-off based on CZ
#--------------------------------------------------------

cutpoint.function <- function(x, index) {
  d <- x[index, ]
  td_roc_cutp <- survivalROC.C(
    Stime = d$survival_years,
    status = d$survival,
    marker = -d$hemoglobin, # since low hemoglobin is associated with worse OS
    predict.time = c(3), # Time point of interest in days
  )
  Sen_t3 <- td_roc_cutp$TP
  Spe_t3 <- 1 - td_roc_cutp$FP
  CZ <- Sen_t3*Spe_t3
  opt_cutp <- td_roc_cutp$cut.values[which.max(CZ)]
  return(opt_cutp)
}
set.seed(3)
BootDist_CZ_boot <- boot(data = d1, statistic = cutpoint.function, R=999) 
BootDist_CZ_boot
boot.ci( BootDist_CZ_boot, conf = .95, type='all')
hist(BootDist_CZ_boot)
plot(BootDist_CZ_boot)

#--------------------------------------------------------
# Table S3.7: Corrected biomarker effect using the heuristic shrinkage factor
#             and Uncorrected biomarker effect
#--------------------------------------------------------

cutp_optimal <- 14.4 #change the cutoff identified by each method to assess its effect on HR
d1$hemoglobin_ <- cut(d1$hemoglobin,
                      breaks = c(0, cutp_optimal, max(d1$hemoglobin)),
                      labels = c("Low", "High"), # Low =< cut-off
                      right = TRUE, include.lowest = TRUE)
d1$hemoglobin_ <- factor(d1$hemoglobin_, levels=c("High", "Low"))
cox_fit <- coxph(Surv(survival_years, survival) ~ hemoglobin_, data = d1)

# uncorrected HR
summary(cox_fit)

# HR adjusted by shrinkage factor c
c <- 1 - cox_fit$var/cox_fit$coefficients^2
c
HR.adj <- exp(c*cox_fit$coefficients) 
round(HR.adj, 2)

# CI based on the adjusted HR and the model variance
round(c(exp(c*cox_fit$coefficients-qnorm(0.975)*sqrt(cox_fit$var)), exp (c*cox_fit$coefficients+qnorm(0.975)*sqrt(cox_fit$var))), 2)

#--------------------------------------------------------
# Table S3.7: Uncorrected absolute risk difference at 3 years for a given cut-off
#--------------------------------------------------------

cutp_optimal <- 13.4 #change the cutoff identified by each method to assess its effect on HR
d1$hemoglobin_ <- cut(d1$hemoglobin,
                      breaks = c(0, cutp_optimal, max(d1$hemoglobin)),
                      labels = c("Low", "High"), # Low =< cut-off
                      right = TRUE, include.lowest = TRUE)
d1$hemoglobin_ <- factor(d1$hemoglobin_, levels=c("High", "Low"))
KM <- survfit(Surv(survival_years, survival) ~ hemoglobin_, data = d1)
rates3yr <- summary(KM, times = 3)
rates3yr
# Absolute risk difference (low - high)
risk_diff <- rates3yr$surv[1] - rates3yr$surv[2]
# 95% CI 
std_err_diff <- sqrt( sum(rates3yr$std.err^2 ) )
lb_risk_diff <- risk_diff - qnorm(0.975)*std_err_diff  
ub_risk_diff <- risk_diff + qnorm(0.975)*std_err_diff  
round(100*risk_diff, 1)
round(100*lb_risk_diff, 1)
round(100*ub_risk_diff, 1)

#--------------------------------------------------------
# Other options in finding cutoffs
#--------------------------------------------------------
### Package timeROC 

l_cut_points <- seq(min(-d1$hemoglobin), max(-d1$hemoglobin), by=0.1)
TP <- NULL; FP <- NULL
for(i in 1:length(l_cut_points)){
  td_measu <- SeSpPPVNPV(
    cutpoint = l_cut_points[i],
    T = d1$survival_years, 
    marker = -d1$hemoglobin, 
    cause = 1,
    delta = d1$survival, 
    times = 3,
    weighting ="marginal")
  TP <- c(TP, td_measu$TP[2])
  FP <- c(FP, td_measu$FP[2])
}
TP
FP

Sen_t3 <- TP
Spe_t3 <- 1-FP
YI <- Sen_t3 + Spe_t3 - 1
ER <- sqrt((0-FP)^2+(1-TP)^2)
CZ <- Sen_t3*Spe_t3
plot(l_cut_points, YI, ylim=c(0,1))
lines(l_cut_points, ER)
lines(l_cut_points, CZ)

# Optimal cut-offs based on youden (YI), ER, and CZ
l_cut_points[which.max(YI)]
l_cut_points[which.min(ER)]
l_cut_points[which.max(CZ)]


library(cenROC)
ROC_cenROC <- cenROC(Y = d1$survival_years, 
       M = -d1$hemoglobin, 
       censor = d1$survival, 
       t = 3, 
       method = "emp", 
       B=999) # not smoothing
ROC_cenROC$AUC
