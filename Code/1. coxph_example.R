# -----------------------------------------------------------------------------
# SCRIPT DESCRIPTION
# -----------------------------------------------------------------------------
# This script investigate the association between heatwave exposure
# and cardiometabolic disorders (including type 2 diabetes,
# dyslipidaemie, hypertension) using cox regression model (coxph).
#
# The key functionalities are:
# 1.  Systematically fits a series of coxph for different heatwave exposure 
#     metrics, including both heatwave events (HWE) and heatwave days (HWD).
# 2.  Fit a full adjusted model for each exposure metric.
# 3.  Scales the  exposure variable by its Interquartile Range (IQR) to
#     standardize the interpretation of the resulting coefficients.
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
  
  pv<-summary(model)$coefficients[, 5]
  pv<-ifelse(pv<0.001, '<0.001', as.character(round(pv,3)))
  result<-cbind(hr_result, pv[1])
  return(result)
}


###### import data for analysis 
setwd('D:/hw')
listn<-list.files(getwd(),  pattern=".csv", full.names="FALSE")

###### Loop for analysis for three outcomes separately using each dataset
df_result_all<-data.frame()# Initialize results data frame for all outcomes
outcome_name<-c('type 2 diabetes', 'dyslipidaemia', 'hypertension')

for (j in 1:length(listn)){
  df <- read.csv(listn[j])

  df$smoke3<-as.factor(df$smoke3)
  df$drink3<-as.factor(df$drink3)
  df$vegefruit3<-as.factor(df$vegefruit3)
  df$region3<-as.factor(df$region3)
  
  ###### Calculate IQR and scale exposure variable by IQR
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
  
  ###### select exposure variable for loop
  expo_list<-df[, c('hmean90_2d_sum_p1y_piqr', 'hmean95_2d_sum_p1y_piqr', 'hmean90_3d_sum_p1y_piqr', 'hmean95_3d_sum_p1y_piqr',
                    'hmean90_2d_sumdays_p1y_piqr',  'hmean95_2d_sumdays_p1y_piqr',  'hmean90_3d_sumdays_p1y_piqr',  'hmean95_3d_sumdays_p1y_piqr')]
  expo_list_name<- c('hmean90_2d_sum_p1y_piqr', 'hmean95_2d_sum_p1y_piqr', 'hmean90_3d_sum_p1y_piqr', 'hmean95_3d_sum_p1y_piqr',
                     'hmean90_2d_sumdays_p1y_piqr',  'hmean95_2d_sumdays_p1y_piqr',  'hmean90_3d_sumdays_p1y_piqr',  'hmean95_3d_sumdays_p1y_piqr')
  
  df_result_expo<-data.frame() ## initial results data frame for each exposure variable
  
  ### loop for model for each exposure metric
  for (i in 1:ncol(expo_list)){
    
    df$expo_hmean<-expo_list[, i]
    
    ## Define Covariates 
    cov <- "+ age + sex + edu2 + yr+ smoke3 + drink3 + vegefruit3 + pa2 + bmi2 + income2 + region3+ pm25 + ndvi"
    
    # Define model formulas 
    formula<- as.formula(paste0("Surv(tstart, tstop, event) ~ expo_hmean",cov))
    
    # Fit models 
    model<- coxph(as.formula(formula), data = df, ties = "efron")
    
    # Bind results for the current exposure 
    
    exposure_names<-expo_list_name[i]
    
    df_result_expo <- rbind(df_result_expo,
                            cbind(outcome_name[j], exposure_names, HRfun(model))
                           )
  }
  # Bind results from current outcome
  df_result_all<-rbind(df_result_all, df_result_expo)
  
}  
