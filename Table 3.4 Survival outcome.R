################################################################################
# Table 3.4 (in the main chapter): Performance measures for biomarkers
# with a time-to-event outcome (overall survival [OS])

# This script assesses the incremental prognostic value of p16 status and EGFR 
# expression for OS (survival outcome) in the presence of known prognostic
# factors (see more details in the chapter) 

# Models evaluated:
#   1) Base (reduced) model: Age, Zubrod performance status, 
#      smoking pack-years (via restricted cubic splines),
#      T stage, and N stage
#   2) Full model: Base model + p16 status 
#   3) Full model: Base model + EGFR level (binary cutoff at 80%)

# Evaluates binary biomarkers using performance measures (PMs): 
#   - AIC/BIC, Log likelihood, Harrell's c, Uno's c, CPE, 
#     AUC (with 95% CI, Delong method), Brier score and scaled Brier score. 
# Quantifies the incremental prognostic value using 
# performance improvement measures (PIMs): 
#   - LRT, Delta Harrell's c-index, Delta Uno's c-index, Delta CPE, 
#     Delta UAC, Delta Brier, continuous NRI(0), modified NRI, and IDI

# Confidence intervals (CIs) for Brier, Scaled Brier, delta Harrell's c, 
# Uno's c, delta Uno's c, delta CPE, delta AUC, NRI(0) are based on the 
# non-parametric percentile bootstrap method. 
# CIs for IDI were based on the pertubation-resampling method

# Note: All analyses were conducted using R version 4.5.1. The packages and functions 
# described in this chapter and its supplement were available and verified at the 
# time of writing. Future updates to these packages may cause code examples to behave differently.
################################################################################


#--------------------------------------------------------
# Load packages
#--------------------------------------------------------

library(tidyverse) #a collection of core R packages for data science
library(survival) #for survival analysis
library(nricens) #for calculating NRI 
library(boot) # for bootstrapping and 95% bootstrap CIs
library(rms) # Regression Modeling Strategies: using rcs() to calculate basis terms for restricted cubic splines
library(Hmisc) # rcorr.cens() to calculate Harrell's c
library(survAUC) # UnoC to calculate Uno's c
library(CPE) # phcpe to calculate CPE
library(timeROC) # time dependent ROC, AUC
library(survIDINRI) # IDI using censored data
library(SurvEval) # calculate CPE_projection 

#--------------------------------------------------------
#Load data
#--------------------------------------------------------
# Read the simulated data set (.csv) from GitHub repository (n=300)
file_url <- "https://raw.githubusercontent.com/chapman-biomarker-book-desachy-dangelo/chapter3_prognostic-biomarkers/main/simulated_hnca_300.csv"

d1 <- read_csv(file_url) %>%
  # Convert binary variables of interest to factors for modeling and plotting
  mutate(across(c("zubrod", "t_stage_bin", "n_stage_bin", "p16_status",  "egfr_level_80"),
                ~as.factor(.x), .names="{.col}_F")) %>%
  # Set p16_status reference group to positive (i.e., compares neg vs pos)
  mutate(p16_status_F=relevel(p16_status_F, ref="1"))

str(d1)


# Adding the basis function columns of spline components for pack years to dataset, 
# since some functions used later (e.g. IDI.INF(), ProjectionCPE()) for performance measures 
# do not take 'rcs()'
# Note: rcs from rms returns object; extract from column 2 for basis
b=rcs(d1$pack_years,4) # basis functions evaluated at pack-years
d1$pack_years1=b[,2]
d1$pack_years2=b[,3]

#--------------------------------------------------------
# define and fit 3-year OS cox regression models
#--------------------------------------------------------

os_base=coxph(Surv(survival_years, survival)~age+zubrod_F+rcs(pack_years,4)+t_stage_bin_F+n_stage_bin_F, data=d1, x=TRUE)
os_base_p16=coxph(Surv(survival_years, survival)~age+zubrod_F+rcs(pack_years,4)+t_stage_bin_F+n_stage_bin_F+p16_status_F, data=d1, x=TRUE)
os_base_egfr=coxph(Surv(survival_years, survival)~age+zubrod_F+rcs(pack_years,4)+t_stage_bin_F+n_stage_bin_F+egfr_level_80_F, data=d1, x=TRUE)

# group them into a list for iteration
os_models=list(os_base, os_base_p16, os_base_egfr)
names(os_models)=c("os_base", "os_base_p16", "os_base_egfr")

#--------------------------------------------------------
# Calculate model-level performance measures: 
# (AIC/BIC, logLik, Harrell's c, Uno's c, CPE, 
# time-specific AUC, Brier, scaled Brier)
#--------------------------------------------------------

# initialize lists
OS_roc_results=list() # to store the OS ROC results, can use for plotting the ROC curve
OS_measure_results=NULL # to store AIC, BIC, log likelihood, AUC, and other measures

for (i in 1:length(os_models)){
  testing_model=os_models[[i]]
  
  model_aic=round(AIC(testing_model),2) #AIC
  model_bic=round(BIC(testing_model),2) #BIC
  model_aic_bic=paste0(model_aic, "/", model_bic)
  model_logplik=round(logLik(testing_model),2) #Log partial likelihood
  
  #For Harrell's c-index, Uno's c-index, AUC
  riskscore=predict(testing_model, type="lp") # linear predictor of cox model
  surv_obj=Surv(d1$survival_years, d1$survival) #surv object
  
  #calculating Harrell's c-index and its 95% CI, rcorr.cens(x, S), the "direction" of the index
  #is opposite from the concordance() from survival package, simply change the 
  #sign of x in rcorr.cens to get the index
  Harrells_c=rcorr.cens(I(-1*riskscore), surv_obj)[1]
  Harrells_c_SE=rcorr.cens(I(-1*riskscore), surv_obj)[3]/2
  Harrells_c_lower=round(Harrells_c+qnorm(0.025)*Harrells_c_SE,2)
  Harrells_c_upper=round(Harrells_c+qnorm(0.975)*Harrells_c_SE,2)
  
  Harrells_text=paste0(round(Harrells_c,2), " (", Harrells_c_lower, ", ", Harrells_c_upper,")")
  
  #calculating Uno's c-index
  Unos_text=round(UnoC(surv_obj, surv_obj, riskscore),2)
  
  #CPE and its 95%CI 
  CPE_index=phcpe(testing_model, CPE.SE = T)[[1]][1]
  CPE_SE=phcpe(testing_model, CPE.SE = T)[[2]][1]
  CPE_lower=round(CPE_index+qnorm(0.025)*CPE_SE,2)
  CPE_upper=round(CPE_index+qnorm(0.975)*CPE_SE,2)
  CPEs_text=paste0(round(CPE_index,2), " (", CPE_lower, ", ", CPE_upper,")")
  
  #Time specific AUC at t=3 years
  OS_roc_results[[i]]=timeROC(T=d1$survival_years,
                  delta=d1$survival,
                  marker=riskscore,
                  cause=1,
                  times=3, 
                  ROC=T, 
                  iid=T) #required TRUE for computation of all inference procedures, e.g. CI 
  AUC_t3_value=OS_roc_results[[i]]$AUC[2]
  AUC_t3_CI=paste0(round(AUC_t3_value,2), "(", round(confint(OS_roc_results[[i]])$CI_AUC[1]/100,2),
                   ", ", round(confint(OS_roc_results[[i]])$CI_AUC[2]/100,2), ")") 
  
  #Time-specific Brier score and scaled Brier at t=3 years: brier() from package survival
  brier_score=brier(testing_model, times=3)
  brier_t3=round(brier_score$brier,2)
  brier_scaled_t3=round(brier_score$rsquared,4)
  
  # Combine all results into a data table format
  OS_measure_results=rbind(OS_measure_results, 
                           cbind(names(os_models)[i], model_aic_bic,  model_logplik, Harrells_text,
                                 Unos_text, CPEs_text, AUC_t3_CI, brier_t3, brier_scaled_t3))
  colnames(OS_measure_results)=c("Model", "AIC/BIC", "logLik", "Harrell's c", "Uno's c",
                                 "CPE", "3yrAUC", "3yrBrier", "3yrScaledBrier")
  
}

View(OS_measure_results)


#--------------------------------------------------------
# Model Comparison Function (LRT, NRI)
#--------------------------------------------------------

model_comparison_results=list()

model_comparison=function(Comparison_pair, reg_base, reg_new){
  
  # Inputs
  ## Comparison_pair=text: description of which models are being compared,
  ## reg_base=for LRT and NRI : base model regression object,
  ## reg_new=for LRT and NRI: new model with additional biomarker regression object,
  
  Comparison=Comparison_pair
  
  #Likelihood ratio test
  LRT=anova(reg_base, reg_new)
  LRT.test.stat=round(LRT$Chisq[2], 2) #likelihood ratio test statstic
  LRT.pval=round(LRT$`Pr(>|Chi|)`[2],4)  # LRT p-value
  
  ##Continuous NRI(0)
  NRI0_results=nricens(mdl.std=reg_base, mdl.new=reg_new, t0=3, updown = "diff", cut=0, alpha=0.05) 
  NRI0_CI=paste0(round(NRI0_results$nri$Estimate[1],2), " (", round(NRI0_results$nri$Lower[1],2), ", ", round(NRI0_results$nri$Upper[1],2), ")") #overall NRI(0) and 95%CI
  NRI0_event_CI=paste0(round(NRI0_results$nri$Estimate[2],2), " (", round(NRI0_results$nri$Lower[2],2), ", ", round(NRI0_results$nri$Upper[2],2), ")") #event NRI(0) and 95%CI
  NRI0_NonEvent_CI=paste0(round(NRI0_results$nri$Estimate[3],2), " (", round(NRI0_results$nri$Lower[3],2), ", ", round(NRI0_results$nri$Upper[3],2), ")") #non-event NRI(0) and 95%CI
  
  #Append results to performance_results
  comparison_results=cbind(Comparison, 
                           LRT.test.stat, LRT.pval, 
                           NRI0_CI, NRI0_event_CI, NRI0_NonEvent_CI)
  model_comparison_results[[length(model_comparison_results)+1]]<<-comparison_results
}


# Example calls

# OS: Model 2 (base+p-16) vs. 1 (base)
model_comparison(Comparison_pair="base+p-16 vs. base",
                 reg_base=os_base, reg_new=os_base_p16)

# OS: Model 3 (base+EGFR) vs. 1 (base)
model_comparison(Comparison_pair="base+EGFR vs. base", 
                 reg_base=os_base, reg_new=os_base_egfr)

# convert list to data frame
model_comparison_results=as.data.frame(do.call(rbind, model_comparison_results))


#--------------------------------------------------------
# Percentile bootstrap confidence intervals for Delta Harrell's c, Uno's c and
# Delta Uno's c, Delta CPE, Delta CPE_proj, Brier, Scaled Brier, Delta AUC, Delta Brier
#--------------------------------------------------------

set.seed(325)
measures=NULL

calculate_measures=function(data, indices, formula_base, formula_new, base_covs, new_cov){
  
  # Inputs:
  ## data=datafile,
  ## indices=integer vector: bootstrap row indices supplied automatically by boot::boot (do NOT pass manually).
  ## formula_new: formula for the new model that includes the biomarker of interest, e.g., OS_3yrs~age+...+zubrod+p16
  ## formula_base: formula for the base_model, e.g., OS_3yrs~age+...+zubrod
  ## base_covs=base model variable columns, only for ProjectionCPE()
  ## new_cov=new biomarker variable column name, only for ProjectionCPE()
  
  d_boot=data[indices,]
  
  model_base = coxph(as.formula(formula_base), data=d_boot, x=TRUE) #cox base model
  model_new = coxph(as.formula(formula_new), data=d_boot, x=TRUE) #cox new model, with the biomarker
  
  surv_obj=Surv(d_boot$survival_years, d_boot$survival) #surv object, using in rcorr.cens(), UnoC()
  riskscore_base=predict(model_base, type="lp", newdata=d_boot) # linear predictor of cox base model
  riskscore_new=predict(model_new, type="lp", newdata=d_boot) # linear predictor of cox new model
  
  Harrells_c_base=rcorr.cens(I(-1*riskscore_base), surv_obj)[1] #Harrell's c for base model
  Harrells_c_new=rcorr.cens(I(-1*riskscore_new), surv_obj)[1] #Harrell's c for new model
  delta_Harrell=Harrells_c_new-Harrells_c_base

  Unos_c_base=UnoC(surv_obj, surv_obj, riskscore_base) #Uno's c for base model
  Unos_c_new=UnoC(surv_obj, surv_obj, riskscore_new) #Uno's c for new model
  delta_Uno=Unos_c_new-Unos_c_base
  
  CPE_base=phcpe(model_base, CPE.SE = T)[[1]][1] #CPE for base model
  CPE_new=phcpe(model_new, CPE.SE = T)[[1]][1] #CPE for new model
  delta_CPE=CPE_new-CPE_base
  
  #projection based CPE
  covs0=as.matrix(d_boot[, base_covs]) #base model variable columns
  covs1=as.matrix(d_boot[, new_cov]) #additional biomarker
  CPE_proj_base=ProjectionCPE(d_boot$survival_years, d_boot$survival,
                              StandardMarkers = covs0, #base model variable columns
                              NewMarkers = covs1, #additional biomarker
                              tau = max(d_boot$survival_years[d_boot$survival==1]), #max observed failure time
                              Block = T) #default Block=T, recommended for dataset>150 obs
  CPE_proj_base_value=CPE_proj_base[[1]]
  delta_CPEproj=CPE_new-CPE_proj_base_value
  
  AUC_base_ROC=timeROC(T=d_boot$survival_years, delta=d_boot$survival, marker=riskscore_base,
                       cause=1, times=3, ROC=T, iid=T) #AUC for base model
  AUC_base=AUC_base_ROC$AUC[2] #extract the AUC value
  
  AUC_new_ROC=timeROC(T=d_boot$survival_years, delta=d_boot$survival, marker=riskscore_new,
                      cause=1, times=3, ROC=T, iid=T) #AUC for new model
  AUC_new=AUC_new_ROC$AUC[2] #extract the AUC value
  delta_AUC=AUC_new-AUC_base
  
  brier_base_output=brier(model_base, times=3, newdata=d_boot)
  brier_base=brier_base_output$brier #brier score for base model
  brier_scaled_base=brier_base_output$rsquared #scaled brier score for base model
  
  brier_new_output=brier(model_new, times=3, newdata=d_boot)
  brier_new=brier_new_output$brier #brier score for new model
  brier_scaled_new=brier_new_output$rsquared #scaled brier score for new model
  
  delta_Brier=brier_new-brier_base

  #Combine all measures together
  measures=rbind(measures, cbind(Harrells_c_base, Harrells_c_new, delta_Harrell,
                             Unos_c_base, Unos_c_new, delta_Uno,
                             CPE_base, CPE_new, delta_CPE, CPE_proj_base_value, delta_CPEproj, 
                             AUC_base, AUC_new, delta_AUC,
                             brier_base, brier_scaled_base,
                             brier_new, brier_scaled_new, delta_Brier))

  return(measures)
  
}

results=list() #store all outputs from the bootstrap
all_CIs=list()

boot_results=function(formula_base, formula_new, base_covs, new_cov, Comparison){
  bootstrap_output=boot(data=d1,
                        statistic=calculate_measures, 
                        R=1000, # bootstrap samples
                        formula_new=formula_new,
                        formula_base=formula_base,
                        base_covs=base_covs,
                        new_cov=new_cov)
  
  stat_names=c("Harrells_c_base", "Harrells_c_new", "delta_Harrell",
                       "Unos_c_base", "Unos_c_new", "delta_Uno",
                       "CPE_base", "CPE_new", "delta_CPE", "CPE_proj_base_value", "delta_CPEproj",
                       "AUC_base", "AUC_new", "delta_AUC",
                       "brier_base", "brier_scaled_base",
                       "brier_new", "brier_scaled_new", "delta_Brier")
  
  results[[length(results)+1]]<<-bootstrap_output ##save all bootstrap outputs
  
  original_stat=bootstrap_output$t0
  
  
  for (i in 1:length(stat_names)) {
    
    original_val = round(original_stat[i],2)
    ci_result = boot.ci(bootstrap_output, index = i, type = "perc") #get the percentile CI for each stat
    lower_bound = ci_result$percent[4] # Access the 4th element (lower bound)
    upper_bound = ci_result$percent[5] # Access the 5th element (upper bound)
    CI=paste0("(",  round(lower_bound,2), ", ", round(upper_bound,2), ")")
    # Store the results
    all_CIs[[length(all_CIs)+1]] <<- c(Comparison, stat_names[i], original_val, CI)
  }
  
}

# Example calls

# Comparison="Model 2 (base+p-16) vs. 1 (base)"
boot_results(Comparison="Model 2 (base+p-16) vs. 1 (base)",
             formula_base=Surv(survival_years, survival)~age+zubrod_F+rcs(pack_years,4)+t_stage_bin_F+n_stage_bin_F,
             formula_new=Surv(survival_years, survival)~age+zubrod_F+rcs(pack_years,4)+t_stage_bin_F+n_stage_bin_F+p16_status_F,
             base_covs=c("age", "zubrod", "pack_years", "pack_years1", "pack_years2", "t_stage_bin", "n_stage_bin"),
             new_cov="p16_status")

# Comparison="Model 3 (base+egfr) vs. 1 (base)"
boot_results(Comparison="Model 3 (base+EGFR) vs. 1 (base)",
             formula_base=Surv(survival_years, survival)~age+zubrod_F+rcs(pack_years,4)+t_stage_bin_F+n_stage_bin_F,
             formula_new=Surv(survival_years, survival)~age+zubrod_F+rcs(pack_years,4)+t_stage_bin_F+n_stage_bin_F+egfr_level_80_F,
             base_covs=c("age", "zubrod", "pack_years", "pack_years1", "pack_years2", "t_stage_bin", "n_stage_bin"),
             new_cov="egfr_level_80")

# Convert to a dataframe, all est (95%CI) for both comparisons
all_CIs=as.data.frame(do.call(rbind, all_CIs))

#Original estimates and bootstrap est for each sample
original=as.data.frame(t(rbind(results[[1]]$t0, results[[2]]$t0)))
colnames(original)=c("new=+p16", "new=+EGFR")
bootstrap_est=as.data.frame(results[[1]]$t)
colnames(bootstrap_est)=c("Harrells_c_base", "Harrells_c_new", "delta_Harrell",
  "Unos_c_base", "Unos_c_new", "delta_Uno",
  "CPE_base", "CPE_new", "delta_CPE", "CPE_proj_base_value", "delta_CPEproj",
  "AUC_base", "AUC_new", "delta_AUC",
  "brier_base", "brier_scaled_base",
  "brier_new", "brier_scaled_new", "delta_Brier")

#--------------------------------------------------------
# IDI and its 95% CI based on perturbation resampling using IDI.INF() 
# from package survIDINRI 
#--------------------------------------------------------

##IDI.INF()- If any factor variable is involved in the set of predictors, 
##use dummy coding 0/1 indicator varaibles (the variable names without the _F. 

set.seed(564) #for reproducible results, function below uses perturbation resampling

# function for calculating the IDI
calculate_idi=function(variable_to_incl, var_to_remove){
  d_boot=d1[,variable_to_incl]
  tp=3 #specify time point for analysis
  indata1=d_boot #data for the new model
  indata0=d_boot[,-which(names(d_boot) == var_to_remove)] #data for the base model, specify the column to remove
  covs1=as.matrix(indata1[,c(-1, -2)]) #need to be in matrix form, 1 and 2 refers to the survival outcome columns
  covs0=as.matrix(indata0[,c(-1, -2)])
  d_boot2= apply(d_boot[,1:2], 2, as.numeric) ### make sure the survival time and status are numeric variables
  
  x=IDI.INF(d_boot2, covs0, covs1, tp, npert=300)  #npert=300 is the default
  IDI=round(x$m1[1], 4)
  IDI_CI_lower=round(x$m1[2],4)
  IDI_CI_upper=round(x$m1[3],4)
  IDI_CI=paste0(IDI, " (",  IDI_CI_lower, ", ", IDI_CI_upper, ")")
  
  return(IDI_CI)
}

# Example calls

# Comparison="Model 2 (base+p-16) vs. 1 (base)"
calculate_idi(variable_to_incl=c("survival_years", "survival","age", "zubrod", "t_stage_bin", "n_stage_bin", "pack_years", 
                                 "pack_years1", "pack_years2", "p16_status"),
              var_to_remove="p16_status")

# Comparison="Model 3 (base+EGFR) vs. 1 (base)"
calculate_idi(variable_to_incl=c("survival_years", "survival","age", "zubrod", "t_stage_bin", "n_stage_bin", "pack_years", 
                                    "pack_years1", "pack_years2", "egfr_level_80"),
                 var_to_remove="egfr_level_80")

