################################################################################
# Table S3.3 (in the Chapter Supplement - see link below)
# https://github.com/chapman-biomarker-book-desachy-dangelo/chapter3_prognostic-biomarkers
# Performance measures for biomarkers with a binary outcome (3-year OS)

# This script assesses the incremental prognostic value of p16 status and EGFR 
# expression for 3-year OS (binary outcome) in the presence of known prognostic
# factors (see more details in the chapter supplement) 

# Models evaluated:
#   1) Base (reduced) model: Age, Zubrod performance status, 
#      smoking pack-years (via restricted cubic splines),
#      T stage, and N stage
#   2) Full model: Base model + p16 status 
#   3) Full model: Base model + EGFR level (binary cutoff at 80%)

# Evaluates binary biomarkers using performance measures (PMs): 
#   - AIC/BIC, Log likelihood, AUC (with 95% CI, Delong method), 
#     Brier score, and scaled Brier score. 
# Quantifies the incremental prognostic value using performance 
# improvement measures (PIMs): 
#   - Delta AUC (Delong test), delta Brier, LRT, continuous NRI(0), 
#     modified NRI, and IDI.

# Confidence intervals for Brier, Scaled Brier, Delta Brier, NRI(0), 
# mNRI, and IDI are based on the non-parametric percentile bootstrap method 

# Note: All analyses were conducted using R version 4.5.1. The packages and functions 
# described in this chapter and its supplement were available and verified at the 
# time of writing. Future updates to these packages may cause code examples to behave differently.
################################################################################

#--------------------------------------------------------
# Load packages
#--------------------------------------------------------

library(tidyverse) # for data management
library(pROC) # for ROC, AUC
library(nricens) # for Net reclassification index (NRI)
library(lmtest) # for likelihood ratio test LRT
library(boot) # for bootstrapping and 95% bootstrap CIs
library(rms) # Regression Modeling Strategies: using rcs() to calculate basis terms for restricted cubic splines
library(DescTools) # can use BrierScoreCI() to calculate Brier Score
library(broom) # used in boot helper function; augment() adds columns to the original dataset such as predictions
library(readr) # to read .csv file from GitHub

#--------------------------------------------------------
#Load data
#--------------------------------------------------------
# Read the simulated data set (.csv) from GitHub repository (n=300)
file_url <- "https://raw.githubusercontent.com/chapman-biomarker-book-desachy-dangelo/chapter3_prognostic-biomarkers/main/simulated_hnca_300.csv"

d1 <- read_csv(file_url) %>%
  # Convert binary variables of interest to factors for modeling and plotting
  mutate(across(c("zubrod", "t_stage_bin", "n_stage_bin", "p16_status", "egfr_level_80"), as.factor)) %>%
  # Set p16_status reference group to positive (i.e., compares neg vs pos)
  mutate(p16_status_F=relevel(p16_status, ref="1"),
  # define a binary variable for 3-year overall survival
  OS_3yrs=ifelse(survival==1 & survival_years<=3, 1, 0))

str(d1)

#--------------------------------------------------------
# Define and fit logistic regression models 
# for 3-year OS (binary)
#--------------------------------------------------------

# Base model
os3yr_base=glm(OS_3yrs~age+zubrod+rcs(pack_years,4)+t_stage_bin+n_stage_bin, data=d1, family = "binomial", x=TRUE)
# Base model+p16
os3yr_base_p16=glm(OS_3yrs~age+zubrod+rcs(pack_years,4)+t_stage_bin+n_stage_bin+p16_status, data=d1, family = "binomial", x=TRUE)
# Base model +EGFR
os3yr_base_egfr=glm(OS_3yrs~age+zubrod+rcs(pack_years,4)+t_stage_bin+n_stage_bin+egfr_level_80, data=d1, family = "binomial", x=TRUE)

# Group models in a list for iteration
os3yr_models=list(os3yr_base, os3yr_base_p16, os3yr_base_egfr)
names(os3yr_models)=c("os3yr_base", "os3yr_base_p16", "os3yr_base_egfr")

#--------------------------------------------------------
# Calculate performance measures: AIC, BIC, loglik, AUC, Brier's score
#--------------------------------------------------------

# Define lists to store results
OS_roc_results=list() # to store the ROC objects, can use for plotting the ROC curve
OS_measure_results=NULL # to store AIC, BIC, log likelihood, AUC, and other measures, will be a data frame/matrix

# Loop through models and calculate performance measures
for (i in 1:length(os3yr_models)){
  testing_model=os3yr_models[[i]]
  
  model_aic=round(AIC(testing_model),2) #AIC
  model_bic=round(BIC(testing_model),2) #BIC
  model_loglik=round(logLik(testing_model),2) #Log likelihood
  
  # Save the ROC objects (roc(actual outcome values, predicted probs))
  OS_roc_results[[i]]=roc(testing_model$y, testing_model$fitted.values) 
  
  # AUC with 95%CI
  # ci.auc() returns lower (index 1), mid, upper (index 3)
  model_auc=round(auc(OS_roc_results[[i]]),2) 
  model_auc_ci=paste0("(", round(ci.auc(OS_roc_results[[i]])[1],2), ", ", 
                      round(ci.auc(OS_roc_results[[i]])[3],2), ")") 
  
  # Brier Score with bootstrap CI using DescTools::BrierScoreCI
  set.seed(3695)
  brier=BrierScoreCI(testing_model, R=1000) #calculate the brier score and 95% bootstrap percentile CI
  brier_score=round(brier[1], 2)  #Brier score
  brier_95CI=paste0("(", round(brier[2], 2), ", ", round(brier[3], 2), ")") #brier 95% CI 
  brier_score_95CI=paste(brier_score, " ", brier_95CI)
  
  
  # Combine all results (rows) into a data table format in OS_measure_results
  OS_measure_results=rbind(OS_measure_results, cbind(names(os3yr_models)[i], model_aic, model_bic, model_loglik, 
                                                     model_auc, model_auc_ci, brier_score_95CI))
  colnames(OS_measure_results)=c("Model", "AIC", "BIC", "logLik", "AUC", "AUC_CI", "Brier_95CI")
  
}

OS_measure_results <- as.data.frame(OS_measure_results, stringsAsFactors = FALSE)
View(OS_measure_results)


#--------------------------------------------------------
# Model Comparison Function: computes delta AUC, LRT, NRI,
# these functions also output the corresponding 95% CI
#--------------------------------------------------------

# List to store the model comparison results
model_comparison_results=list()

model_comparison=function(Outcome, Comparison_pair, roc_new, roc_std, reg_base, reg_new){
  
  ## Outcome=text label, e.g "3-year OS"
  ## Comparison_pair=text: description of which models are being compared,
  ## roc_new=for deltaAUC: OS_roc_results[[j]] enter integer j, referring to the model with additional biomarker saved in the prior "roc_results" list,
  ## roc_std=for deltaAUC: OS_roc_results[[i]]; enter integer i, index referring to the base model saved in the prior "roc_results" list,
  ## reg_base=for LRT and NRI : base model regression object,
  ## reg_new=for LRT and NRI: new model with additional biomarker regression object,

  Outcome=Outcome
  Comparison=Comparison_pair
  
  # Delta AUC: difference in AUCs and Delong test
  # roc_new, roc_std are roc objects from list OS_roc_results
  delta.auc=round(auc(roc_new)-auc(roc_std),2) # delta AUC
  model.roctest=roc.test(roc_new, roc_std) # Delong test
  CI.delta.auc=paste0("(", round(model.roctest$conf.int[1],2), ", ", round(model.roctest$conf.int[2],2), ")") # delta AUC 95% CI
  pval.delta.auc=round(model.roctest$p.value,4) #Delong test p-value
  
  # Likelihood ratio test 
  LRT=lmtest::lrtest(reg_base, reg_new)
  LRT.test.stat=round(LRT[2, "Chisq"], 2) #likelihood ratio test statistic
  LRT.pval=round(LRT[2, "Pr(>Chisq)"], 4)  # LRT p-value
  
  ## Continuous NRI(0)
  NRI0_results=nribin(mdl.std=reg_base, mdl.new=reg_new, cut=0, updown = 'diff', alpha=0.05) 
  NRI0_CI=paste0(round(NRI0_results$nri$Estimate[1],2), " (", round(NRI0_results$nri$Lower[1],2), ", ", round(NRI0_results$nri$Upper[1],2), ")") #overall NRI(0) and 95%CI
  NRI0_event_CI=paste0(round(NRI0_results$nri$Estimate[2],2), " (", round(NRI0_results$nri$Lower[2],2), ", ", round(NRI0_results$nri$Upper[2],2), ")") #event NRI(0) and 95%CI
  NRI0_NonEvent_CI=paste0(round(NRI0_results$nri$Estimate[3],2), " (", round(NRI0_results$nri$Lower[3],2), ", ", round(NRI0_results$nri$Upper[3],2), ")") #non-event NRI(0) and 95%CI

  # Collect and append results to comparison_results, assign model_comparison_results to the environment
  comparison_results=cbind(Outcome, Comparison, 
                           delta.auc, CI.delta.auc, pval.delta.auc,
                           LRT.test.stat, LRT.pval,
                           NRI0_CI, NRI0_event_CI, NRI0_NonEvent_CI)
  model_comparison_results[[length(model_comparison_results)+1]]<<-comparison_results
}

#--------------------------------------------------------
#  OS model comparisons
#--------------------------------------------------------
# Indices for the OS models
# index=1: os3yr_base
# index=2: os3yr_base_p16
# index=3: os3yr_base_egfr

# Example calls
#### 3-year OS: Model 2 (base+p-16) vs. 1 (base)
model_comparison(Outcome="3-year OS", Comparison_pair="base+p-16 vs. base", 
                 roc_new=OS_roc_results[[2]], roc_std=OS_roc_results[[1]], 
                 reg_base=os3yr_base, reg_new=os3yr_base_p16)

#### 3-year OS: Model 3 (base+EGFR) vs. 1 (base)
model_comparison(Outcome="3-year OS", Comparison_pair="base+EGFR vs. base", 
                 roc_new=OS_roc_results[[3]], roc_std=OS_roc_results[[1]], 
                 reg_base=os3yr_base, reg_new=os3yr_base_egfr)

###############
## Convert Model comparison list to data frame
model_comparison_results=as.data.frame(do.call(rbind, model_comparison_results))

#--------------------------------------------------------
# IDI: Manual Calculations and its 95% percentile bootstrap CI
#--------------------------------------------------------

###function for calculating the IDI
calculate_idi=function(data, indices, formula_new, formula_base, Outcome){
  
  # data=data frame
  # indices=integer vector: bootstrap row indices supplied automatically by boot::boot (do NOT pass manually).
  # formula_new: formula for the new model that includes the biomarker of interest, e.g., OS_3yrs~age+...+zubrod+p16
  # formula_base: formula for the base_model, e.g., OS_3yrs~age+...+zubrod
  # Outcome: name of the outcome column, e.g, "OS_3yrs"
  
  d_boot = data[indices, ]  
  # Fit the models to the sample
  reg_new = glm(formula_new, data = d_boot, family = "binomial")
  reg_base = glm(formula_base, data = d_boot, family = "binomial")
  
  # Calculate predicted probabilities using the models and add to dataset
  d_boot =  d_boot %>%
    mutate(
      new =
        broom::augment(reg_new, type.predict = "response") %>%
        pull(".fitted"),
      base =
        broom::augment(reg_base, type.predict = "response") %>%
        pull(".fitted")
    )
  
  avg_event=mean(d_boot[d_boot[[Outcome]]==1,]$new)-mean(d_boot[d_boot[[Outcome]]==1,]$base)
  avg_nonevent=mean(d_boot[d_boot[[Outcome]]==0,]$new)-mean(d_boot[d_boot[[Outcome]]==0,]$base)
  IDI=avg_event-avg_nonevent
  return(IDI)
}

## Function to perform bootstrapping and return CI
set.seed(456)
IDI_results=list()

IDI_boot_results=function(formula_new, formula_base, Outcome, Comparison){
  # formula_new: formula for the new model that includes the biomarker of interest, e.g., OS_3yrs~age+...+zubrod+p16
  # formula_base: formula for the base_model, e.g., OS_3yrs~age+...+zubrod
  # Outcome: name of the outcome column, e.g, "OS_3yrs"
  # Comparison: text, label for the output, e.g., "base+p16 vs base"
  
  bootstrap_output=boot(data=d1,
                        statistic=calculate_idi, 
                        R=1000,
                        formula_new=formula_new,
                        formula_base=formula_base,
                        Outcome=Outcome)
  IDI=bootstrap_output$t0
  IDI_percCI=boot.ci(bootstrap_output, type = "perc") 
  IDI_CI_lower=IDI_percCI$percent[4]
  IDI_CI_upper=IDI_percCI$percent[5]
  IDI_CI=paste0(round(IDI,4), " (",  round(IDI_CI_lower,4), ", ", round(IDI_CI_upper,4), ")") 
  
  IDI_bootresults=cbind(Outcome, Comparison,  IDI_CI)
  IDI_results[[length(IDI_results)+1]]<<-IDI_bootresults
}

# Example calls
#### 3-year OS: Model 2 (base+p-16) vs. 1 (base)
IDI_boot_results(formula_new=OS_3yrs~age+zubrod+rcs(pack_years,4)+t_stage_bin+n_stage_bin+p16_status,
                 formula_base=OS_3yrs~age+zubrod+rcs(pack_years,4)+t_stage_bin+n_stage_bin, 
                 Outcome="OS_3yrs",
                 Comparison="Model 2 (base+p-16) vs. 1 (base)")

#### 3-year OS: Model 3 (base+EGFR) vs. 1 (base)
IDI_boot_results(formula_new=OS_3yrs~age+zubrod+rcs(pack_years,4)+t_stage_bin+n_stage_bin+egfr_level_80,
                 formula_base=OS_3yrs~age+zubrod+rcs(pack_years,4)+t_stage_bin+n_stage_bin, 
                 Outcome="OS_3yrs",
                 Comparison="Model 3 (base+EGFR) vs. 1 (base)")

#################
## Convert IDI results list to data frame
IDI_results=as.data.frame(do.call(rbind, IDI_results))

#--------------------------------------------------------
# Modified NRI: Manual Calculations and its 95% percentile bootstrap CI
# As described in Section 2 in Heller (2023), "A modified net reclassification 
# improvement statistic", Journal of the Statistical Planning and Inference 227, 18-33.
#--------------------------------------------------------

###function for calculating the mNRI
calculate_mNRI=function(data, indices, formula_new, formula_base, Outcome){
  
  # data=data frame
  # indices=integer vector: bootstrap row indices supplied automatically by boot::boot (do NOT pass manually).
  # formula_new: formula for the new model that includes the biomarker of interest, e.g., OS_3yrs~age+...+zubrod+p16
  # formula_base: formula for the base_model, e.g., OS_3yrs~age+...+zubrod
  # Outcome: name of the outcome column, e.g, "OS_3yrs"
  
  d_boot = data[indices, ]   
  
  # Fit the models to the sample
  reg_new = glm(formula_new, data = d_boot, family = "binomial")
  reg_base = glm(formula_base, data = d_boot, family = "binomial")
  
  # Calculate predicted probabilities using the models and add to dataset
  d_boot =  d_boot %>%
    mutate(
      new =
        broom::augment(reg_new, type.predict = "response") %>%
        pull(".fitted"),
      base =
        broom::augment(reg_base, type.predict = "response") %>%
        pull(".fitted")
    )
  
  #calculate the absolute difference in predicted probabilities between the 2 models
  abs_diff_probs = abs(d_boot$new - d_boot$base)
  #calculate the expected value of the absolute differences
  mean_abs_diff = mean(abs_diff_probs)
  #calculate the observed event rate in the data
  pi_0 = mean(d_boot[[Outcome]]) 
  #calculate the proportionality constant
  prop_constant = 1 / (2 * pi_0 * (1 - pi_0))
  #calculate the mNRI
  mNRI = prop_constant * mean_abs_diff
  
  return(mNRI)
}


## Function to perform bootstrapping and return results
set.seed(789)
mNRI_results=list()

mNRI_boot_results=function(formula_new, formula_base, Outcome, Comparison){
  
  # formula_new: formula for the new model that includes the biomarker of interest, e.g., OS_3yrs~age+...+zubrod+p16
  # formula_base: formula for the base_model, e.g., OS_3yrs~age+...+zubrod
  # Outcome: name of the outcome column, e.g, "OS_3yrs"
  # Comparison: text, label for the output, e.g., "base+p16 vs base"
  
  bootstrap_output=boot(data=d1,
                        statistic=calculate_mNRI, 
                        R=1000,
                        formula_new=formula_new,
                        formula_base=formula_base,
                        Outcome=Outcome)
  mNRI=bootstrap_output$t0
  mNRI_percCI=boot.ci(bootstrap_output, type = "perc") 
  mNRI_CI_lower=mNRI_percCI$percent[4]
  mNRI_CI_upper=mNRI_percCI$percent[5]
  mNRI_CI=paste0(round(mNRI,2), " (",  round(mNRI_CI_lower,2), ", ", round(mNRI_CI_upper,2), ")") 
  
  mNRI_bootresults=cbind(Outcome, Comparison,  mNRI_CI)
  mNRI_results[[length(mNRI_results)+1]]<<-mNRI_bootresults
}

# Example calls
#### 3-year OS: Model 2 (base+p-16) vs. 1 (base)
mNRI_boot_results(formula_new=OS_3yrs~age+zubrod+rcs(pack_years,4)+t_stage_bin+n_stage_bin+p16_status,
                 formula_base=OS_3yrs~age+zubrod+rcs(pack_years,4)+t_stage_bin+n_stage_bin, 
                 Outcome="OS_3yrs",
                 Comparison="Model 2 (base+p-16) vs. 1 (base)")

#### 3-year OS: Model 3 (base+EGFR) vs. 1 (base)
mNRI_boot_results(formula_new=OS_3yrs~age+zubrod+rcs(pack_years,4)+t_stage_bin+n_stage_bin+egfr_level_80,
                 formula_base=OS_3yrs~age+zubrod+rcs(pack_years,4)+t_stage_bin+n_stage_bin, 
                 Outcome="OS_3yrs",
                 Comparison="Model 3 (base+EGFR) vs. 1 (base)")

#################
## Convert mNRI results list to data frame
mNRI_results=as.data.frame(do.call(rbind, mNRI_results))

#--------------------------------------------------------
# Scaled Brier: manual calculations and its 95% percentile bootstrap CI
#--------------------------------------------------------

###function for calculating the scaled Brier
calculate_scaledBrier=function(data, indices, formula, model){
  
  # data=data frame
  # indices=integer vector: bootstrap row indices supplied automatically by boot::boot (do NOT pass manually).
  # formula: for each model, base, base+p16, base+EGFR
  # model: name of the model, e.g, "base"
  
  d_boot = data[indices, ]   
  
  # Fit the model to the sample
  reg_fit = glm(formula, data = d_boot, family = "binomial")
  
  # Calculate predicted probabilities using the models and add to dataset
  d_boot =  d_boot %>%
    mutate(
      predicted.p =
        broom::augment(reg_fit, type.predict = "response") %>%
        pull(".fitted")) %>%
    mutate(difference_sq=(OS_3yrs-predicted.p)^2)
  
  #calculate the Brier
  Brier=mean(d_boot$difference_sq)
  #calculate the Brier max
  Brier_max=mean(d_boot$predicted.p)*(1-mean(d_boot$predicted.p))^2+(1-mean(d_boot$predicted.p))*(mean(d_boot$predicted.p))^2
  #calculate the Scaled brier
  scaled_Brier=1-(Brier/Brier_max)
  
  return(scaled_Brier)
}


## Function to perform bootstrapping and return results
set.seed(9454)
scaled_Brier_results=list()

scaled_Brier_boot_results=function(formula, model){
  
  # formula: base, base+p16, or base+EGFR
  # model: text, name of the model, e.g, "base".
  bootstrap_output=boot(data=d1,
                        statistic=calculate_scaledBrier, 
                        R=1000,
                        formula=formula, 
                        model=model)
  scaled_Brier=bootstrap_output$t0
  scaled_Brier_percCI=boot.ci(bootstrap_output, type = "perc") 
  scaled_Brier_CI_lower=scaled_Brier_percCI$percent[4]
  scaled_Brier_CI_upper=scaled_Brier_percCI$percent[5]
  scaled_Brier_CI=paste0(round(scaled_Brier,2), " (",  round(scaled_Brier_CI_lower,2), ", ", round(scaled_Brier_CI_upper,2), ")") 
  
  scaled_Brier_bootresults=cbind(model, scaled_Brier_CI)
  scaled_Brier_results[[length(scaled_Brier_results)+1]]<<-scaled_Brier_bootresults
}

#### 3-year OS: Model 1 (base) 
scaled_Brier_boot_results(formula=OS_3yrs~age+zubrod+rcs(pack_years,4)+t_stage_bin+n_stage_bin, model="base")

#### 3-year OS: Model 2 (base+p-16) 
scaled_Brier_boot_results(formula=OS_3yrs~age+zubrod+rcs(pack_years,4)+t_stage_bin+n_stage_bin+p16_status, model="base+p16")

#### 3-year OS: Model 3 (base+EGFR) 
scaled_Brier_boot_results(formula=OS_3yrs~age+zubrod+rcs(pack_years,4)+t_stage_bin+n_stage_bin+egfr_level_80, model="base+egfr")

#################
## Convert scaled Brier results list to data frame
scaled_Brier_results=as.data.frame(do.call(rbind, scaled_Brier_results))

#--------------------------------------------------------
# Delta Brier: Manual calculation and its 95% percentile bootstrap CI
# DescTools package has a function that can calculate BrierScore() for each model
#--------------------------------------------------------
###function for calculating the delta Brier
calculate_deltaBrier=function(data, indices, formula_base, formula_new, model_comp){
  
  # data=data frame
  # indices=integer vector: bootstrap row indices supplied automatically by boot::boot (do NOT pass manually).
  # formula_base: formula for base model
  # formula_new: formula for base model + biomarker of interest
  # model_comp: text of the comparison pair, e.g, "base+p16 vs. base"
  
  d_boot = data[indices, ]   
  
  model_base = glm(formula_base, data=d_boot, family = "binomial") #base model
  model_new = glm(formula_new,, data=d_boot, family = "binomial") #new model, with the biomarker
  
  d_boot =  d_boot %>%
    mutate(
      new =
        broom::augment(model_new, type.predict = "response") %>%
        pull(".fitted"),
      base =
        broom::augment(model_base, type.predict = "response") %>%
        pull(".fitted")
    )
  
  # Manual Calculations
  # brier_base=mean((d_boot$OS_3yrs-d_boot$base)^2)
  # brier_new=mean((d_boot$OS_3yrs-d_boot$new)^2)
  
  # Use DescTools::BrierScore()
  brier_base=BrierScore(d_boot$OS_3yrs, d_boot$base)
  brier_new=BrierScore(d_boot$OS_3yrs, d_boot$new)

  delta_Brier=as.numeric(brier_new-brier_base)
  
  return(delta_Brier)
}

## Function to perform bootstrapping and return results
set.seed(9454)
delta_Brier_results=list()

delta_Brier_boot_results=function(formula_base, formula_new, model_comp){
  
  # formula_base: base model
  # formula_new: base model + biomarker of interest
  # model_comp: text of the comparison pair, e.g, "base+p16 vs. base"

    bootstrap_output=boot(data=d1,
                        statistic=calculate_deltaBrier, 
                        R=1000,
                        formula_base=formula_base,
                        formula_new=formula_new, 
                        model_comp=model_comp)
  delta_Brier=bootstrap_output$t0
  delta_Brier_percCI=boot.ci(bootstrap_output, type = "perc") 
  delta_Brier_CI_lower=delta_Brier_percCI$percent[4]
  delta_Brier_CI_upper=delta_Brier_percCI$percent[5]
  delta_Brier_CI=paste0(round(delta_Brier,2), " (",  round(delta_Brier_CI_lower,2), ", ", round(delta_Brier_CI_upper,2), ")") 
  
  delta_Brier_bootresults=cbind(model_comp, delta_Brier_CI)
  delta_Brier_results[[length(delta_Brier_results)+1]]<<-delta_Brier_bootresults
}

# Example calls
#### 3-year OS: Model 2 (base+p-16) vs. 1 (base)
delta_Brier_boot_results(formula_base=OS_3yrs ~ age + zubrod + rcs(pack_years, 4) + t_stage_bin + n_stage_bin,
                         formula_new=OS_3yrs ~ age + zubrod + rcs(pack_years, 4) + t_stage_bin + n_stage_bin+p16_status,
                 model_comp="Model 2 (base+p-16) vs. 1 (base)")

#### 3-year OS: Model 3 (base+EGFR) vs. 1 (base)
delta_Brier_boot_results(formula_base=OS_3yrs ~ age + zubrod + rcs(pack_years, 4) + t_stage_bin + n_stage_bin,
                         formula_new=OS_3yrs ~ age + zubrod + rcs(pack_years, 4) + t_stage_bin + n_stage_bin+egfr_level_80,
                 model_comp="Model 3 (base+EGFR) vs. 1 (base)")

#################
## Convert delta Brier results list to data frame
delta_Brier_results=as.data.frame(do.call(rbind, delta_Brier_results))

