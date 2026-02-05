# -----------------------------------------------------------------------------
# SCRIPT DESCRIPTION
# -----------------------------------------------------------------------------
# This script investigates effect modification, aiming to determine if the
# effect of heatwave exposure on the cardiometabolic disorders (including type 2 
# diabetes, dyslipidaemie, hypertension) differs across various subgroups of 
# the population.
#
# The key workflow is as follows:
# 1.  Data Preparation for Stratification: Based on the potential modifiers,
#     it creates clean, stratified subsets of the data (e.g., a northern, central, 
#     or southern Taiwan group). A careful approach is taken to ensure
#     individuals remain in a consistent stratum over time.
# 2.  Stratified Analysis: It runs the regression model separately within
#     each stratum. This allows for a direct comparison of the exposure effect
#     size between groups (e.g., comparing the heatwave effect between males
#     and females).
# 3.  Interaction Testing: It systematically tests for statistical interaction 
#     between each heatwave exposure and potential effect modifiers
#     (e.g., age, sex , education level, and geographic region) using cox 
#     regression model (coxph).
# -----------------------------------------------------------------------------
# DATA STRUCTURES
# -----------------------------------------------------------------------------
#
# OUTCOME VARIABLE:
#   - event
#
# EFFECT MODIFIERS (Stratification Variables):
#   - age2_base (<60 or ≥60 years), 
#   - sex (male or female), 
#   - edu2_base (≤ high school or > high school), 
#   - region3_base (northern, central, or southern Taiwan)
#
# PREDICTOR VARIABLES:
#   - Heatwave Exposures: 
#     HWE:hmean90_2d_sum_p1y, 
#         hmean95_2d_sum_p1y, 
#         hmean90_3d_sum_p1y, 
#         hmean95_3d_sum_p1y;
#
#     HWD:hmean90_2d_sumdays_p1y, 
#         hmean95_2d_sumdays_p1y, 
#         hmean90_3d_sumdays_p1y, 
#         hmean95_3d_sumdays_p1y;
#     
#   - Covariates: age, sex, edu2, yr, smoke3, drink3, vegefruit3, pa2, bmi2, 
#                 income2, region3, pm25, ndvi
# -----------------------------------------------------------------------------

# --- SCRIPT START ---

###### Load Required Libraries
library("dplyr")
library(survival)
library("survminer")

# Helper Function ----
HRfun <- function(model) {
  hr<-round(summary(model)$conf.int, digits = 2)
  hr<-ifelse(hr==1, '1.00', hr)
  hr_result<-paste0(hr[1,1], ' (' , hr[1,3], ', ', hr[1,4], ')'  )
  return(hr_result)
}

intfun <- function(model, model1) {
  re<-anova(model, model1)
  pint<-round(re[2,4], digits = 3)
  return(pint)
}


###### import data for analysis 
setwd('D:/hw')
listn <- list.files(getwd(), pattern=".csv", full.names=FALSE)

###### 1. Loop for stratified analysis for three outcomes separately using each dataset
df_result_all<-data.frame()# Initialize results data frame for all outcomes
outcome_name<-c('type 2 diabetes', 'dyslipidaemia', 'hypertension')

for (j in 1:length(listn)){
   df <- read.csv(listn[j])
   
   df$smoke3<-as.factor(df$smoke3)
   df$drink3<-as.factor(df$drink3)
   df$vegefruit3<-as.factor(df$vegefruit3)
   df$region3<-as.factor(df$region3)

   ### Calculate IQR and scale exposure variable by IQR
   #HWE
   df$hmean90_2d_sum_p1y_piqr<-df$hmean90_2d_sum_p1y/IQR(df$hmean90_2d_sum_p1y)
   df$hmean95_2d_sum_p1y_piqr<-df$hmean95_2d_sum_p1y/IQR(df$hmean95_2d_sum_p1y)
   df$hmean90_3d_sum_p1y_piqr<-df$hmean90_3d_sum_p1y/IQR(df$hmean90_3d_sum_p1y)
   df$hmean95_3d_sum_p1y_piqr<-df$hmean95_3d_sum_p1y/IQR(df$hmean95_3d_sum_p1y)
   #HWD
   df$hmean90_2d_sumdays_p1y_piqr<-df$hmean90_2d_sumdays_p1y/IQR(df$hmean90_2d_sumdays_p1y)
   df$hmean95_2d_sumdays_p1y_piqr<-df$hmean95_2d_sumdays_p1y/IQR(df$hmean95_2d_sumdays_p1y)
   df$hmean90_3d_sumdays_p1y_piqr<-df$hmean90_3d_sumdays_p1y/IQR(df$hmean90_3d_sumdays_p1y)
   df$hmean95_3d_sumdays_p1y_piqr<-df$hmean95_3d_sumdays_p1y/IQR(df$hmean95_3d_sumdays_p1y)
   
   ### Create Baseline Stratification Variables 
   # Selects baseline (first) observation for each person and extracts stratification variables.
   base <- df %>% arrange(pid, n) %>% group_by(pid) %>% slice(1) %>% select(pid, age, edu2, region3)
   colnames(base) <- c("pid","age_base","edu2_base","region3_base")
   base$age2_base<-ifelse(base$age_base<60, 0, 1)
   base$region3_base<-as.factor(base$region3_base)
   
   ### Merge Baseline Strata into Main Data 
   # Joins baseline stratification variables back to the full dataset by 'pid'.
   df <- merge(df, base, by="pid")
   
   ###select exposure variable for loop
   expo_list<-df[, c('hmean90_2d_sum_p1y_piqr', 'hmean95_2d_sum_p1y_piqr', 'hmean90_3d_sum_p1y_piqr', 'hmean95_3d_sum_p1y_piqr',
                     'hmean90_2d_sumdays_p1y_piqr',  'hmean95_2d_sumdays_p1y_piqr',  'hmean90_3d_sumdays_p1y_piqr',  'hmean95_3d_sumdays_p1y_piqr')]
   expo_list_name<- c('hmean90_2d_sum_p1y_piqr', 'hmean95_2d_sum_p1y_piqr', 'hmean90_3d_sum_p1y_piqr', 'hmean95_3d_sum_p1y_piqr',
                      'hmean90_2d_sumdays_p1y_piqr',  'hmean95_2d_sumdays_p1y_piqr',  'hmean90_3d_sumdays_p1y_piqr',  'hmean95_3d_sumdays_p1y_piqr')
   
   df_result_expo<-data.frame() ## initial results data frame for each exposure variable
   
   ### loop for stratified analysis for each exposure matrix
   for (i in 1:ncol(expo_list)){
   
   df$expo_hmean<-expo_list[, i]
   
   ## Define Covariates for Stratified Models 
   # Re-defines covariate strings, excluding the stratification variable from the 'cov' string itself.
   cov <- "+ yr+ smoke3 + drink3 + vegefruit3 + pa2 + bmi2 + income2 + pm25 + ndvi"

   # Define model formulas for each stratified analysis
   formula.age <- as.formula(paste0("Surv(tstart, tstop, event) ~ expo_hmean",cov,"+ sex + edu2 + region3"))
   formula.sex <- as.formula(paste0("Surv(tstart, tstop, event) ~ expo_hmean",cov,"+ age + edu2+ region3"))
   formula.edu <- as.formula(paste0("Surv(tstart, tstop, event) ~ expo_hmean",cov,"+ age + sex + region3"))
   formula.region <- as.formula(paste0("Surv(tstart, tstop, event) ~ expo_hmean",cov,"+ age + sex+ edu2"))
   
   # Fit models for each stratum
   model_age0<- coxph(as.formula(formula.age), data = df, subset=(age2_base==0), ties = "efron")
   model_age1<- coxph(as.formula(formula.age), data = df, subset=(age2_base==1), ties = "efron")
   model_sex1<- coxph(as.formula(formula.sex), data = df, subset=(sex==1), ties = "efron")
   model_sex2<- coxph(as.formula(formula.sex), data = df, subset=(sex==2), ties = "efron")
   model_edu0<- coxph(as.formula(formula.edu), data = df, subset=(edu2_base==0), ties = "efron")
   model_edu1<- coxph(as.formula(formula.edu), data = df, subset=(edu2_base==1), ties = "efron")
   model_region0<- coxph(as.formula(formula.region), data = df, subset=(region3_base==0), ties = "efron")
   model_region1<- coxph(as.formula(formula.region), data = df, subset=(region3_base==1), ties = "efron")
   model_region2<- coxph(as.formula(formula.region), data = df, subset=(region3_base==2), ties = "efron")
   
   # Bind results from all strata for the current exposure 
   
   exposure_names<-expo_list_name[i]
   
   df_result_expo <- rbind(df_result_expo,
                cbind(outcome_name[j], exposure_names,group="age",strata="<60",HRfun(model_age0)),
                cbind(outcome_name[j], exposure_names,group="age",strata=">=60",HRfun(model_age1)),
                cbind(outcome_name[j], exposure_names,group="sex",strata="male",HRfun(model_sex1)),
                cbind(outcome_name[j], exposure_names,group="sex",strata="female",HRfun(model_sex2)),
                cbind(outcome_name[j], exposure_names,group="edu",strata="<high school",HRfun(model_edu0)),
                cbind(outcome_name[j], exposure_names,group="edu",strata=">=high school",HRfun(model_edu1)),
                cbind(outcome_name[j], exposure_names,group="region",strata="north",HRfun(model_region0)),
                cbind(outcome_name[j], exposure_names,group="region",strata="central",HRfun(model_region1)),
                cbind(outcome_name[j], exposure_names,group="region",strata="south",HRfun(model_region2))
                )
   }
   # Bind results from current outcome
   df_result_all<-rbind(df_result_all, df_result_expo)

}  

   
###### 2. Loop for interaction for three outcomes separately using each dataset
   df_result_all<-data.frame()# Initialize results data frame for all outcomes
   outcome_name<-c('type 2 diabetes', 'dyslipidaemia', 'hypertension')
   
   for (j in 1:length(listn)){
     df <- read.csv(listn[j])
     
     df$smoke3<-as.factor(df$smoke3)
     df$drink3<-as.factor(df$drink3)
     df$vegefruit3<-as.factor(df$vegefruit3)
     df$region3<-as.factor(df$region3)
     
     ### Calculate IQR and scale exposure variable by IQR
     #HWE
     df$hmean90_2d_sum_p1y_piqr<-df$hmean90_2d_sum_p1y/IQR(df$hmean90_2d_sum_p1y)
     df$hmean95_2d_sum_p1y_piqr<-df$hmean95_2d_sum_p1y/IQR(df$hmean95_2d_sum_p1y)
     df$hmean90_3d_sum_p1y_piqr<-df$hmean90_3d_sum_p1y/IQR(df$hmean90_3d_sum_p1y)
     df$hmean95_3d_sum_p1y_piqr<-df$hmean95_3d_sum_p1y/IQR(df$hmean95_3d_sum_p1y)
     #HWD
     df$hmean90_2d_sumdays_p1y_piqr<-df$hmean90_2d_sumdays_p1y/IQR(df$hmean90_2d_sumdays_p1y)
     df$hmean95_2d_sumdays_p1y_piqr<-df$hmean95_2d_sumdays_p1y/IQR(df$hmean95_2d_sumdays_p1y)
     df$hmean90_3d_sumdays_p1y_piqr<-df$hmean90_3d_sumdays_p1y/IQR(df$hmean90_3d_sumdays_p1y)
     df$hmean95_3d_sumdays_p1y_piqr<-df$hmean95_3d_sumdays_p1y/IQR(df$hmean95_3d_sumdays_p1y)
     
     ### Create Baseline Stratification Variables 
     # Selects baseline (first) observation for each person and extracts stratification variables.
     base <- df %>% arrange(pid, n) %>% group_by(pid) %>% slice(1) %>%
       select(pid, age, edu2, region3)
     colnames(base) <- c("pid","age_base","edu2_base","region3_base")
     base$age2_base<-ifelse(base$age_base<60, 0, 1)
     base$region3_base<-as.factor(base$region3_base)
     
     ### Merge Baseline Strata into Main Data 
     # Joins baseline stratification variables back to the full dataset by 'pid'.
     df <- merge(df, base, by="pid")
     
     ###select exposure variable for loop
     expo_list<-df[, c('hmean90_2d_sum_p1y_piqr', 'hmean95_2d_sum_p1y_piqr', 'hmean90_3d_sum_p1y_piqr', 'hmean95_3d_sum_p1y_piqr',
                       'hmean90_2d_sumdays_p1y_piqr',  'hmean95_2d_sumdays_p1y_piqr',  'hmean90_3d_sumdays_p1y_piqr',  'hmean95_3d_sumdays_p1y_piqr')]
     expo_list_name<- c('hmean90_2d_sum_p1y_piqr', 'hmean95_2d_sum_p1y_piqr', 'hmean90_3d_sum_p1y_piqr', 'hmean95_3d_sum_p1y_piqr',
                        'hmean90_2d_sumdays_p1y_piqr',  'hmean95_2d_sumdays_p1y_piqr',  'hmean90_3d_sumdays_p1y_piqr',  'hmean95_3d_sumdays_p1y_piqr')
     
     df_result_expo<-data.frame() ## initial results data frame for each exposure variable
     
     ### loop for stratified analysis for each exposure matrix
     for (i in 1:ncol(expo_list)){
       
       df$expo_hmean<-expo_list[, i]
       
       ## Define Covariates for Stratified Models 
       # Re-defines covariate strings, excluding the stratification variable from the 'cov' string itself.
       cov <- "+ yr+ smoke3 + drink3 + vegefruit3 + pa2 + bmi2 + income2 + pm25 + ndvi"
       
       # Define model formulas for each stratified analysis
       formula.age_int <- as.formula(paste0("Surv(tstart, tstop, event) ~ expo_hmean*age2_base",cov,"+ sex + edu2 + region3"))
       formula.age <- as.formula(paste0("Surv(tstart, tstop, event) ~ expo_hmean+age2_base",cov,"+ sex + edu2 + region3"))
       formula.sex_int <- as.formula(paste0("Surv(tstart, tstop, event) ~ expo_hmean*sex",cov,"+ age + edu2+ region3"))
       formula.sex <- as.formula(paste0("Surv(tstart, tstop, event) ~ expo_hmean+sex",cov,"+ age + edu2+ region3"))
       formula.edu_int <- as.formula(paste0("Surv(tstart, tstop, event) ~ expo_hmean*edu2_base",cov,"+ age + sex + region3"))
       formula.edu <- as.formula(paste0("Surv(tstart, tstop, event) ~ expo_hmean+edu2_base",cov,"+ age + sex + region3"))
       formula.region_int <- as.formula(paste0("Surv(tstart, tstop, event) ~ expo_hmean*region3_base",cov,"+ age + sex+ edu2"))
       formula.region <- as.formula(paste0("Surv(tstart, tstop, event) ~ expo_hmean+region3_base",cov,"+ age + sex+ edu2"))
       
       # Fit models for each stratum
       model_age_int<- coxph(as.formula(formula.age_int), data = df, ties = "efron")
       model_age<- coxph(as.formula(formula.age), data = df, ties = "efron")
       model_sex_int<- coxph(as.formula(formula.sex_int), data = df,  ties = "efron")
       model_sex<- coxph(as.formula(formula.sex), data = df, ties = "efron")
       model_edu_int<- coxph(as.formula(formula.edu_int), data = df, ties = "efron")
       model_edu<- coxph(as.formula(formula.edu), data = df,  ties = "efron")
       model_region_int<- coxph(as.formula(formula.region_int), data = df,  ties = "efron")
       model_region<- coxph(as.formula(formula.region), data = df, ties = "efron")
       
       
       # Bind results from all strata for the current exposure 
       
       exposure_names<-expo_list_name[i]
       
       df_result_expo <- rbind(df_result_expo,
                               cbind(outcome_name[j], exposure_names,group="age", intfun(model_age_int, model_age)),
                               
                               cbind(outcome_name[j], exposure_names,group="sex", intfun(model_sex_int, model_sex)),
                               
                               cbind(outcome_name[j], exposure_names,group="edu", intfun(model_edu_int, model_edu)),
                               
                               cbind(outcome_name[j], exposure_names,group="region", intfun(model_region_int, model_region))
       )
     }
     # Bind results from current outcome
     df_result_all<-rbind(df_result_all, df_result_expo)
     
   }  
   
  