# -----------------------------------------------------------------------------
# SCRIPT DESCRIPTION
# -----------------------------------------------------------------------------
# This script investigate the combined effects of heatwave and NDVI/PM2.5 on
# cardiometabolic disorders (including type 2 diabetes, dyslipidaemie, 
# hypertension) using cox regression model (coxph).
#
# The key functionalities are:
# 1.  Categorize heatwave exposure metrics [heatwave events (HWE) and heatwave 
#     days (HWD)], NDVI, PM2.5 into two groups separately based on their median 
#     values. The combination of Heatwave and NDVI/PM2.5 are created and 
#     classified into four groups, with the reference group consist of low 
#     heatwave with either high NDVI or low pm2.5. 
# 2.  Systematically fits a series of coxph for different combinations between 
#     heatwave exposure metrics [heatwave events (HWE) and heatwave days (HWD)]
#     and NDVI/PM2.5.
# 3.  tested the additive interaction and multiplicative interaction between 
#     heatwave exposure metrics and NDVI/PM2.5.
# 4.  Extracts, formats, and stores the coefficients and their 95% confidence
#     intervals (CIs) for the exposure variable from each model.
# -----------------------------------------------------------------------------
# DATA STRUCTURES
# -----------------------------------------------------------------------------
# 
# OUTCOME VARIABLE:
#   - event
#
# PREDICTOR VARIABLES:
#   - combination of heatwave and ndvi/pm2.5: 
#     ndvi_HWE_gr2d90, ndvi_HWE_gr2d95, ndvi_HWE_gr3d90, ndvi_HWE_gr3d95,
#     ndvi_HWD_gr2d90, ndvi_HWD_gr2d95, ndvi_HWD_gr3d90, ndvi_HWD_gr3d95,
#     
#     pm25_HWE_gr2d90, pm25_HWE_gr2d95, pm25_HWE_gr3d90, pm25_HWE_gr3d95,
#     pm25_HWD_gr2d90, pm25_HWD_gr2d95, pm25_HWD_gr3d90, pm25_HWD_gr3d95;
#     
#   - Covariates: age, sex, edu2, yr, smoke3, drink3, vegefruit3, pa2, bmi2, 
#                 income2, region3, pm25, ndvi
#
# EFFECT MODIFIERS (Stratification Variables):
#   - ndvi2
#   - pm25_2
# -----------------------------------------------------------------------------



# --- SCRIPT START ---

###### Load Required Libraries and Data
library(summarytools)
library("gdata")
library(readxl)
library("dplyr")
library(survival)
library("survminer")
library(interactionR)

# Helper Function ----
comHRfun <- function(model) {
  hr<-round(summary(model)$conf.int, digits = 2) 
  hr<-ifelse(hr==1, '1.00', hr)
  pv<-round(summary(model)$coefficients[, 5], digits = 4)
  pv<-ifelse(pv<0.001, '<0.001', as.character(round(pv,3)))
  
  hr1<-paste0(hr[1,1], ' (' , hr[1,3], ', ', hr[1,4], ')'  )
  pv1<-paste0(pv[1])
  hr2<-paste0(hr[2,1], ' (' , hr[2,3], ', ', hr[2,4], ')'  )
  pv2<-paste0(pv[2])
  hr3<-paste0(hr[3,1], ' (' , hr[3,3], ', ', hr[3,4], ')'  )
  pv3<-paste0(pv[3])
  result<-cbind(hr1, pv1,  hr2, pv2,  hr3, pv3)
  return(result)
}

intHRfun_NDVI <- function(model) {
  table_object = interactionR(model, exposure_names = c("ndvi2", "expo_hmean1"), ci.type = "delta", ci.level = 0.95, em = F, recode = F)
  reri<-round(table_object$dframe[10, c(2:4) ],  digits = 2)
  reri_r<-paste0(reri[1,1], ' (' , reri[1,2], ', ', reri[1,3], ')'  )
  
  multiscale<-round(table_object$dframe[9, c(2:4) ],  digits = 2)
  multiscale_r<-paste0(multiscale[1,1], ' (' , multiscale[1,2], ', ', multiscale[1,3], ')'  )
  result<-cbind(reri_r, multiscale_r)
  return(result)
}

intHRfun_PM25 <- function(model) {
  table_object = interactionR(model, exposure_names = c("pm25_2", "expo_hmean1"), ci.type = "delta", ci.level = 0.95, em = F, recode = F)
  reri<-round(table_object$dframe[10, c(2:4) ],  digits = 2)
  reri_r<-paste0(reri[1,1], ' (' , reri[1,2], ', ', reri[1,3], ')'  )
  
  multiscale<-round(table_object$dframe[9, c(2:4) ],  digits = 2)
  multiscale_r<-paste0(multiscale[1,1], ' (' , multiscale[1,2], ', ', multiscale[1,3], ')'  )
  result<-cbind(reri_r, multiscale_r)
  return(result)
}

############ 1 combine effects of heatwave and NDVI
###### import data for analysis 
setwd('D:/hw')
listn <- list.files(getwd(), pattern=".csv", full.names=FALSE)

###### Loop for analysis for three outcomes separately using each dataset
df_result_all<-data.frame()# Initialize results data frame for all outcomes
outcome_name<-c('type 2 diabetes', 'dyslipidaemia', 'hypertension')

for (j in 1:length(listn)){
  df <- read.csv(listn[j])
  
  df$smoke3<-as.factor(df$smoke3)
  df$drink3<-as.factor(df$drink3)
  df$vegefruit3<-as.factor(df$vegefruit3)
  df$region3<-as.factor(df$region3)

  ###### categorize heatwave exposure into two groups based on the median value
  ##hwe
  df$hwe_exp_2d90[df$hmean90_2d_sum_p1y<quantile(df$hmean90_2d_sum_p1y,1/2)]<-0
  df$hwe_exp_2d90[df$hmean90_2d_sum_p1y>=quantile(df$hmean90_2d_sum_p1y,1/2)]<-1
  
  df$hwe_exp_2d95[df$hmean95_2d_sum_p1y<quantile(df$hmean95_2d_sum_p1y,1/2)]<-0
  df$hwe_exp_2d95[df$hmean95_2d_sum_p1y>=quantile(df$hmean95_2d_sum_p1y,1/2)]<-1
  
  df$hwe_exp_3d90[df$hmean90_3d_sum_p1y<quantile(df$hmean90_3d_sum_p1y,1/2)]<-0
  df$hwe_exp_3d90[df$hmean90_3d_sum_p1y>=quantile(df$hmean90_3d_sum_p1y,1/2)]<-1
  
  df$hwe_exp_3d95[df$hmean95_3d_sum_p1y<quantile(df$hmean95_3d_sum_p1y,1/2)]<-0
  df$hwe_exp_3d95[df$hmean95_3d_sum_p1y>=quantile(df$hmean95_3d_sum_p1y,1/2)]<-1
  
  ##hwd
  df$hwd_exp_2d90[df$hmean90_2d_sumdays_p1y<quantile(df$hmean90_2d_sumdays_p1y,1/2)]<-0
  df$hwd_exp_2d90[df$hmean90_2d_sumdays_p1y>=quantile(df$hmean90_2d_sumdays_p1y,1/2)]<-1
  
  df$hwd_exp_2d95[df$hmean95_2d_sumdays_p1y<quantile(df$hmean95_2d_sumdays_p1y,1/2)]<-0
  df$hwd_exp_2d95[df$hmean95_2d_sumdays_p1y>=quantile(df$hmean95_2d_sumdays_p1y,1/2)]<-1
  
  df$hwd_exp_3d90[df$hmean90_3d_sumdays_p1y<quantile(df$hmean90_3d_sumdays_p1y,1/2)]<-0
  df$hwd_exp_3d90[df$hmean90_3d_sumdays_p1y>=quantile(df$hmean90_3d_sumdays_p1y,1/2)]<-1
  
  df$hwd_exp_3d95[df$hmean95_3d_sumdays_p1y<quantile(df$hmean95_3d_sumdays_p1y,1/2)]<-0
  df$hwd_exp_3d95[df$hmean95_3d_sumdays_p1y>=quantile(df$hmean95_3d_sumdays_p1y,1/2)]<-1
  
  ###### categorize ndvi into two groups based on the median value
  base <- df %>% arrange(pid, n) %>% group_by(pid) %>% slice(1) %>%
    select(pid, ndvi)
  colnames(base) <- c("pid","ndvi_base")
  base$ndvi2<-ifelse(base$ndvi_base>=round(quantile(base$ndvi_base, 1/2), 2), 0, 1)
  
  ### Merge Baseline Strata into Main Data 
  # Joins baseline stratification variables back to the full dataset by 'pid'.
  df <- merge(df, base, by="pid")
  
  #df$ndvi2[df$ndvi>=round(quantile(df$ndvi, 1/2), 2)]<-0
  #df$ndvi2[df$ndvi<round(quantile(df$ndvi, 1/2), 2)]<-1
  
  ###### combination of heatwave exposure and ndvi    
  df$ndvi_HWE_gr2d90[df$ndvi2==0 & df$hwe_exp_2d90==0]<-0
  df$ndvi_HWE_gr2d90[df$ndvi2==1 & df$hwe_exp_2d90==0]<-1
  df$ndvi_HWE_gr2d90[df$ndvi2==0 & df$hwe_exp_2d90==1]<-2
  df$ndvi_HWE_gr2d90[df$ndvi2==1 & df$hwe_exp_2d90==1]<-3
  
  df$ndvi_HWE_gr2d95[df$ndvi2==0 & df$hwe_exp_2d95==0]<-0
  df$ndvi_HWE_gr2d95[df$ndvi2==1 & df$hwe_exp_2d95==0]<-1
  df$ndvi_HWE_gr2d95[df$ndvi2==0 & df$hwe_exp_2d95==1]<-2
  df$ndvi_HWE_gr2d95[df$ndvi2==1 & df$hwe_exp_2d95==1]<-3
  
  df$ndvi_HWE_gr3d90[df$ndvi2==0 & df$hwe_exp_3d90==0]<-0
  df$ndvi_HWE_gr3d90[df$ndvi2==1 & df$hwe_exp_3d90==0]<-1
  df$ndvi_HWE_gr3d90[df$ndvi2==0 & df$hwe_exp_3d90==1]<-2
  df$ndvi_HWE_gr3d90[df$ndvi2==1 & df$hwe_exp_3d90==1]<-3
  
  df$ndvi_HWE_gr3d95[df$ndvi2==0 & df$hwe_exp_3d95==0]<-0
  df$ndvi_HWE_gr3d95[df$ndvi2==1 & df$hwe_exp_3d95==0]<-1
  df$ndvi_HWE_gr3d95[df$ndvi2==0 & df$hwe_exp_3d95==1]<-2
  df$ndvi_HWE_gr3d95[df$ndvi2==1 & df$hwe_exp_3d95==1]<-3
  
  # HW days   
  df$ndvi_HWD_gr2d90[df$ndvi2==0 & df$hwd_exp_2d90==0]<-0
  df$ndvi_HWD_gr2d90[df$ndvi2==1 & df$hwd_exp_2d90==0]<-1
  df$ndvi_HWD_gr2d90[df$ndvi2==0 & df$hwd_exp_2d90==1]<-2
  df$ndvi_HWD_gr2d90[df$ndvi2==1 & df$hwd_exp_2d90==1]<-3
  
  df$ndvi_HWD_gr2d95[df$ndvi2==0 & df$hwd_exp_2d95==0]<-0
  df$ndvi_HWD_gr2d95[df$ndvi2==1 & df$hwd_exp_2d95==0]<-1
  df$ndvi_HWD_gr2d95[df$ndvi2==0 & df$hwd_exp_2d95==1]<-2
  df$ndvi_HWD_gr2d95[df$ndvi2==1 & df$hwd_exp_2d95==1]<-3
  
  df$ndvi_HWD_gr3d90[df$ndvi2==0 & df$hwd_exp_3d90==0]<-0
  df$ndvi_HWD_gr3d90[df$ndvi2==1 & df$hwd_exp_3d90==0]<-1
  df$ndvi_HWD_gr3d90[df$ndvi2==0 & df$hwd_exp_3d90==1]<-2
  df$ndvi_HWD_gr3d90[df$ndvi2==1 & df$hwd_exp_3d90==1]<-3
  
  df$ndvi_HWD_gr3d95[df$ndvi2==0 & df$hwd_exp_3d95==0]<-0
  df$ndvi_HWD_gr3d95[df$ndvi2==1 & df$hwd_exp_3d95==0]<-1
  df$ndvi_HWD_gr3d95[df$ndvi2==0 & df$hwd_exp_3d95==1]<-2
  df$ndvi_HWD_gr3d95[df$ndvi2==1 & df$hwd_exp_3d95==1]<-3
  
  ###### select each combination exposure for loop
  expo_list<-df[, c('ndvi_HWE_gr2d90', 'ndvi_HWE_gr2d95', 'ndvi_HWE_gr3d90', 'ndvi_HWE_gr3d95',
                    'ndvi_HWD_gr2d90', 'ndvi_HWD_gr2d95', 'ndvi_HWD_gr3d90', 'ndvi_HWD_gr3d95')]
  
  ###### select each categorical heatwave exposure for loop
  expo_list1<-df[, c('hwe_exp_2d90', 'hwe_exp_2d95', 'hwe_exp_3d90', 'hwe_exp_3d95',
                     'hwd_exp_2d90', 'hwd_exp_2d95', 'hwd_exp_3d90', 'hwd_exp_3d95')]
  
  df_result_expo<-data.frame() ## initial results data frame for each exposure variable
  
  ### loop for model for each exposure matrix
  for (i in 1:ncol(expo_list)){
    
    # specify the exposure for analysis
    df$expo_hmean<-expo_list[, i]
    df$expo_hmean1<-expo_list1[, i]
    
    # Define Covariates 
    cov <- "+ age + sex + edu2 + yr+ smoke3 + drink3 + vegefruit3 + pa2 + bmi2 + income2 + region3+ pm25 "
    
    # Define model formulas 
    formula_com<- as.formula(paste0("Surv(tstart, tstop, event) ~ as.factor(expo_hmean)",cov)) #for the combination exposure
    formula_int<- as.formula(paste0("Surv(tstart, tstop, event) ~ expo_hmean1*ndvi2",cov))#for the interaction
    
    # Fit models 
    model_com<- coxph(as.formula(formula_com), data = df, ties = "efron")#for the combination exposure
    model_int<- coxph(as.formula(formula_int), data = df, ties = "efron")#for the interaction
    
    # Bind results for the current exposure 
    exposure_names<-names(expo_list[i])
    outcome_names<-outcome_name[j]
    
    df_result_expo <- rbind(df_result_expo,
                            cbind(outcome_names, exposure_names, comHRfun(model_com), intHRfun_NDVI(model_int))
                           )
    
  }

  df_result_all<-rbind(df_result_all, df_result_expo)
  
}


############ 2 combine effects of heatwave and pm25
###### import data for analysis 
setwd('D:/hw')
listn <- list.files(getwd(), pattern=".csv", full.names=FALSE)

###### Loop for analysis for three outcomes separately using each dataset
df_result_all<-data.frame()# Initialize results data frame for all outcomes
outcome_name<-c('type 2 diabetes', 'dyslipidaemia', 'hypertension')

for (j in 1:length(listn)){
  df <- read.csv(listn[j])
  
  df$smoke3<-as.factor(df$smoke3)
  df$drink3<-as.factor(df$drink3)
  df$vegefruit3<-as.factor(df$vegefruit3)
  df$region3<-as.factor(df$region3)
  
  ###### categorize heatwave exposure into two groups based on the median value
  ##hwe
  df$hwe_exp_2d90[df$hmean90_2d_sum_p1y<quantile(df$hmean90_2d_sum_p1y,1/2)]<-0
  df$hwe_exp_2d90[df$hmean90_2d_sum_p1y>=quantile(df$hmean90_2d_sum_p1y,1/2)]<-1
  
  df$hwe_exp_2d95[df$hmean95_2d_sum_p1y<quantile(df$hmean95_2d_sum_p1y,1/2)]<-0
  df$hwe_exp_2d95[df$hmean95_2d_sum_p1y>=quantile(df$hmean95_2d_sum_p1y,1/2)]<-1
  
  df$hwe_exp_3d90[df$hmean90_3d_sum_p1y<quantile(df$hmean90_3d_sum_p1y,1/2)]<-0
  df$hwe_exp_3d90[df$hmean90_3d_sum_p1y>=quantile(df$hmean90_3d_sum_p1y,1/2)]<-1
  
  df$hwe_exp_3d95[df$hmean95_3d_sum_p1y<quantile(df$hmean95_3d_sum_p1y,1/2)]<-0
  df$hwe_exp_3d95[df$hmean95_3d_sum_p1y>=quantile(df$hmean95_3d_sum_p1y,1/2)]<-1
  
  ##hwd
  df$hwd_exp_2d90[df$hmean90_2d_sumdays_p1y<quantile(df$hmean90_2d_sumdays_p1y,1/2)]<-0
  df$hwd_exp_2d90[df$hmean90_2d_sumdays_p1y>=quantile(df$hmean90_2d_sumdays_p1y,1/2)]<-1
  
  df$hwd_exp_2d95[df$hmean95_2d_sumdays_p1y<quantile(df$hmean95_2d_sumdays_p1y,1/2)]<-0
  df$hwd_exp_2d95[df$hmean95_2d_sumdays_p1y>=quantile(df$hmean95_2d_sumdays_p1y,1/2)]<-1
  
  df$hwd_exp_3d90[df$hmean90_3d_sumdays_p1y<quantile(df$hmean90_3d_sumdays_p1y,1/2)]<-0
  df$hwd_exp_3d90[df$hmean90_3d_sumdays_p1y>=quantile(df$hmean90_3d_sumdays_p1y,1/2)]<-1
  
  df$hwd_exp_3d95[df$hmean95_3d_sumdays_p1y<quantile(df$hmean95_3d_sumdays_p1y,1/2)]<-0
  df$hwd_exp_3d95[df$hmean95_3d_sumdays_p1y>=quantile(df$hmean95_3d_sumdays_p1y,1/2)]<-1
  
  ###### categorize ndvi into two groups based on the median value
  base <- df %>% arrange(pid, n) %>% group_by(pid) %>% slice(1) %>%
    select(pid, pm25)
  colnames(base) <- c("pid","pm25_base")
  base$pm25_2<-ifelse(base$pm25_base<round(quantile(base$pm25_base, 1/2),2), 0, 1)
  
  ### Merge Baseline Strata into Main Data 
  # Joins baseline stratification variables back to the full dataset by 'pid'.
  df <- merge(df, base, by="pid")
  
  ###### categorize PM25 into two groups based on the median value
  #df$pm25_2[df$pm25<round(quantile(df$pm25, 1/2), 2)]<-0
  #df$pm25_2[df$pm25>=round(quantile(df$pm25, 1/2), 2)]<-1
  #df$pm25_2[df$pm25<23.78]<-0
  #df$pm25_2[df$pm25>=23.78]<-1
  
  ###### combination of heatwave exposure and pm25    
  df$pm25_HWE_gr2d90[df$pm25_2==0 & df$hwe_exp_2d90==0]<-0
  df$pm25_HWE_gr2d90[df$pm25_2==1 & df$hwe_exp_2d90==0]<-1
  df$pm25_HWE_gr2d90[df$pm25_2==0 & df$hwe_exp_2d90==1]<-2
  df$pm25_HWE_gr2d90[df$pm25_2==1 & df$hwe_exp_2d90==1]<-3
  
  df$pm25_HWE_gr2d95[df$pm25_2==0 & df$hwe_exp_2d95==0]<-0
  df$pm25_HWE_gr2d95[df$pm25_2==1 & df$hwe_exp_2d95==0]<-1
  df$pm25_HWE_gr2d95[df$pm25_2==0 & df$hwe_exp_2d95==1]<-2
  df$pm25_HWE_gr2d95[df$pm25_2==1 & df$hwe_exp_2d95==1]<-3
  
  df$pm25_HWE_gr3d90[df$pm25_2==0 & df$hwe_exp_3d90==0]<-0
  df$pm25_HWE_gr3d90[df$pm25_2==1 & df$hwe_exp_3d90==0]<-1
  df$pm25_HWE_gr3d90[df$pm25_2==0 & df$hwe_exp_3d90==1]<-2
  df$pm25_HWE_gr3d90[df$pm25_2==1 & df$hwe_exp_3d90==1]<-3
  
  df$pm25_HWE_gr3d95[df$pm25_2==0 & df$hwe_exp_3d95==0]<-0
  df$pm25_HWE_gr3d95[df$pm25_2==1 & df$hwe_exp_3d95==0]<-1
  df$pm25_HWE_gr3d95[df$pm25_2==0 & df$hwe_exp_3d95==1]<-2
  df$pm25_HWE_gr3d95[df$pm25_2==1 & df$hwe_exp_3d95==1]<-3
  
  # HW days   
  df$pm25_HWD_gr2d90[df$pm25_2==0 & df$hwd_exp_2d90==0]<-0
  df$pm25_HWD_gr2d90[df$pm25_2==1 & df$hwd_exp_2d90==0]<-1
  df$pm25_HWD_gr2d90[df$pm25_2==0 & df$hwd_exp_2d90==1]<-2
  df$pm25_HWD_gr2d90[df$pm25_2==1 & df$hwd_exp_2d90==1]<-3
  
  df$pm25_HWD_gr2d95[df$pm25_2==0 & df$hwd_exp_2d95==0]<-0
  df$pm25_HWD_gr2d95[df$pm25_2==1 & df$hwd_exp_2d95==0]<-1
  df$pm25_HWD_gr2d95[df$pm25_2==0 & df$hwd_exp_2d95==1]<-2
  df$pm25_HWD_gr2d95[df$pm25_2==1 & df$hwd_exp_2d95==1]<-3
  
  df$pm25_HWD_gr3d90[df$pm25_2==0 & df$hwd_exp_3d90==0]<-0
  df$pm25_HWD_gr3d90[df$pm25_2==1 & df$hwd_exp_3d90==0]<-1
  df$pm25_HWD_gr3d90[df$pm25_2==0 & df$hwd_exp_3d90==1]<-2
  df$pm25_HWD_gr3d90[df$pm25_2==1 & df$hwd_exp_3d90==1]<-3
  
  df$pm25_HWD_gr3d95[df$pm25_2==0 & df$hwd_exp_3d95==0]<-0
  df$pm25_HWD_gr3d95[df$pm25_2==1 & df$hwd_exp_3d95==0]<-1
  df$pm25_HWD_gr3d95[df$pm25_2==0 & df$hwd_exp_3d95==1]<-2
  df$pm25_HWD_gr3d95[df$pm25_2==1 & df$hwd_exp_3d95==1]<-3
  
  df_result_expo<-data.frame() ## initial results data frame for each combination
  
  ###### select each combination exposure for loop
  expo_list<-df[, c('pm25_HWE_gr2d90', 'pm25_HWE_gr2d95', 'pm25_HWE_gr3d90', 'pm25_HWE_gr3d95',
                    'pm25_HWD_gr2d90', 'pm25_HWD_gr2d95', 'pm25_HWD_gr3d90', 'pm25_HWD_gr3d95')]
  
  ###### select each categorical heatwave exposure for loop
  expo_list1<-df[, c('hwe_exp_2d90', 'hwe_exp_2d95', 'hwe_exp_3d90', 'hwe_exp_3d95',
                     'hwd_exp_2d90', 'hwd_exp_2d95', 'hwd_exp_3d90', 'hwd_exp_3d95')]
  
  df_result_expo<-data.frame() ## initial results data frame for each exposure variable
  
  ### loop for model for each exposure matrix
  for (i in 1:ncol(expo_list)){
    
    # specify the combination exposure for analysis
    df$expo_hmean<-expo_list[, i]
    df$expo_hmean1<-expo_list1[, i]
    
    # Define Covariates 
    cov <- "+ age + sex + edu2 + yr+ smoke3 + drink3 + vegefruit3 + pa2 + bmi2 + income2 + region3+ ndvi "
    
    # Define model formulas 
    formula_com<- as.formula(paste0("Surv(tstart, tstop, event) ~ as.factor(expo_hmean)",cov)) #for the combination exposure
    formula_int<- as.formula(paste0("Surv(tstart, tstop, event) ~ expo_hmean1*pm25_2",cov))#for the interaction
    
    # Fit models 
    model_com<- coxph(as.formula(formula_com), data = df, ties = "efron")#for the combination exposure
    model_int<- coxph(as.formula(formula_int), data = df, ties = "efron")#for the interaction
    
    # Bind results for the current exposure 
    exposure_names<-names(expo_list[i])
    
    df_result_expo <- rbind(df_result_expo,
                            cbind(outcome_name[j], exposure_names, comHRfun(model_com), intHRfun_PM25(model_int))
    )
    
  }
  
  df_result_all<-rbind(df_result_all, df_result_expo)
  
}
