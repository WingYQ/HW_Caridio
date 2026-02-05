# -----------------------------------------------------------------------------
# SCRIPT DESCRIPTION
# -----------------------------------------------------------------------------
# This script conducts a detailed analysis of the exposure-response (E-R)
# relationship between various heatwave exposures and  cardiometabolic disorders 
# (including type 2 diabetes, dyslipidaemie, hypertension).
# It is a multi-stage process designed to select the best model form and then
# visualize the final relationships.
#
# The key workflow is as follows:
# 1.  knot Selection: For exposures, it fits multiple models with varying knot 
#     location and uses the Akaike Information Criterion (AIC) to select
#     the optimal knot location that best balances model fit and complexity.
# 2.  Linearity Assessment: For each heatwave exposure metric, it fits both a
#     linear model and a non-linear model (using a natural cubic spline). It
#     then uses a Likelihood Ratio Test (LRT) to determine if the non-linear
#     fit is statistically superior.
# 3.  Concentration-Response Curve Generation: Using the optimal model form
#     (linear or non-linear with the best knot), it fits the final models and
#     generates C-R plots for heatwave events (HWE) and heatwave days (HWD).
#
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

###### Load Required Libraries and Data
library("gdata")
library("dplyr")
library(survival)
library("survminer")
library(splines) ## for natural cubic spline
library(dlnm)
library(ggplot2)
library(rms)

############ 1. Loop to Calculate AIC for Different knots
# Iterates through chosen exposures and different knots
# to find the optimal knots based on AIC.

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
  
  ###### defind knot location
  knotsn<-data.frame(knots1=c(0.1,0.5,0.9), knots2=c(0.2, 0.5, 0.8), knots3=c(0.25, 0.5, 0.75), knots4=c(0.10, 0.75, 0.90))
  knotsn_name<-c('0.1, 0.5, 0.9', '0.2, 0.5, 0.8', '0.25, 0.5, 0.75', '0.10, 0.75, 0.90')
  
  ### loop for model for each exposure matrix
  for (i in 1:ncol(expo_list)){
    
    df$expo_hmean<-expo_list[, i]
    
    df_result_knot<-data.frame()## initial results data frame for each knot
    
    for (k in 1:ncol(knotsn)){
      
      ## Define Covariates 
      cov <- "+ age + sex + edu2 + yr+ smoke3 + drink3 + vegefruit3 + pa2 + bmi2 + income2 + region3+ pm25 + ndvi"
      
      # Define model formulas 
      formula<- as.formula(paste0("Surv(tstart, tstop, event) ~ ns(expo_hmean, knots=quantile(expo_hmean, probs = knotsn[, k]))",cov))
      
      # Fit models 
      model<- coxph(as.formula(formula), data = df, ties = "efron")
      
      # bind result for all knot for 
      exposure_names<-expo_list_name[i]
      
      df_result_knot<- rbind(df_result_knot, cbind(outcome=outcome_name[j], 
                                                   exposure=exposure_names, 
                                                   knot=knotsn_name[k], 
                                                   best_knot_number=k, 
                                                   aic=AIC(model)
      )
      )
    }
    # Bind results for the current exposure
    df_result_expo <- rbind(df_result_expo, df_result_knot)
    
  }
  # Bind results from current outcome
  df_result_all<-rbind(df_result_all, df_result_expo)
  
}  


###### select best knot accroding to minimum AIC
#create exposure and outcome id so that the results were presented in the given order
df_exposure_id<-data.frame(exposure_id=c(1:8),
                  exposure=c('hmean90_2d_sum_p1y_piqr', 'hmean95_2d_sum_p1y_piqr', 'hmean90_3d_sum_p1y_piqr', 'hmean95_3d_sum_p1y_piqr',
                             'hmean90_2d_sumdays_p1y_piqr',  'hmean95_2d_sumdays_p1y_piqr',  'hmean90_3d_sumdays_p1y_piqr',  'hmean95_3d_sumdays_p1y_piqr')
                  )

df_result_all1<-merge(df_result_all, df_exposure_id, by='exposure', all.x = TRUE)

df_outcome_id<-data.frame(outcome_id=c(1:3), outcome=c('type 2 diabetes', 'dyslipidaemia', 'hypertension'))

df_result_all2<-merge(df_result_all1, df_outcome_id,  by='outcome', all.x = TRUE)

#Groups results by outcome and exposure, then filters for the minimum AIC to find the optimal knot.
best_knot<- df_result_all2 %>%
  group_by(outcome_id, exposure_id) %>%
  filter(aic == min(aic)) %>%
  select(outcome, exposure, best_knot_number) %>%
  ungroup() %>%
  arrange(outcome_id, exposure_id)
best_knot

############ 2 linearity assessment
###### import data for analysis 
setwd('D:/hw')
listn <- list.files(getwd(), pattern=".csv", full.names=FALSE)

###### extract the best knot for each exposure
knotsn<-as.numeric(best_knot$best_knot_number)
knotsn_df<-data.frame(knots1=c(0.1,0.5,0.9), knots2=c(0.2, 0.5, 0.8), knots3=c(0.25, 0.5, 0.75), knots4=c(0.10, 0.75, 0.90))

###### Loop for analysis for three outcomes separately using each dataset
df_result_all<-data.frame()# Initialize results data frame for all outcomes
outcome_name<-c('type 2 diabetes', 'dyslipidaemia', 'hypertension')

for (j in 1:length(listn)){
  df <- read.csv(listn[j])
  
  df$smoke3<-as.factor(df$smoke3)
  df$drink3<-as.factor(df$drink3)
  df$vegefruit3<-as.factor(df$vegefruit3)
  df$region3<-as.factor(df$region3)
  
  ###### select exposure variable for loop
  expo_list<-df[, c('hmean90_2d_sum_p1y', 'hmean95_2d_sum_p1y', 'hmean90_3d_sum_p1y', 'hmean95_3d_sum_p1y',
                    'hmean90_2d_sumdays_p1y',  'hmean95_2d_sumdays_p1y',  'hmean90_3d_sumdays_p1y',  'hmean95_3d_sumdays_p1y')]
  expo_list_name<- c('hmean90_2d_sum_p1y', 'hmean95_2d_sum_p1y', 'hmean90_3d_sum_p1y', 'hmean95_3d_sum_p1y',
                     'hmean90_2d_sumdays_p1y',  'hmean95_2d_sumdays_p1y',  'hmean90_3d_sumdays_p1y',  'hmean95_3d_sumdays_p1y')
  
  
  df_result_expo<-data.frame() ## initial results data frame for each exposure variable
  
  ### loop for model for each exposure metric
  for (i in 1:ncol(expo_list)){
    
    df$expo_hmean<-expo_list[, i]
    
    # Define covariates 
    cov <- "+ age + sex + edu2 + yr+ smoke3 + drink3 + vegefruit3 + pa2 + bmi2 + income2 + region3+ pm25 + ndvi"
      
    # Define model formulas 
    formula_knot<- as.formula(paste0("Surv(tstart, tstop, event) ~ ns(expo_hmean, knots=quantile(expo_hmean, probs = knotsn_df[knotsn[(j-1)*8+1]]))",cov))
    formula<- as.formula(paste0("Surv(tstart, tstop, event) ~ expo_hmean",cov))
                                                                                    
    # Fit models 
    model_knot<- coxph(as.formula(formula_knot), data = df, ties = "efron")
    model<- coxph(as.formula(formula), data = df, ties = "efron")
      
    # extract p-value
    p<-anova(model, model_knot)[2, 4]
    p<-ifelse(p<0.001, '<0.001', round(p,3))
      
    # bind result for the current exposure
    exposure_names<-expo_list_name[i]
      
    df_result_expo<- rbind(df_result_expo, cbind(outcome=outcome_name[j], 
                                                 exposure=exposure_names, 
                                                 p_nolinear=p)
                           )

  }
  # Bind results from current outcome
  df_result_all<-rbind(df_result_all, df_result_expo)
  
}  


############ 3 dose response curve
###### import data for analysis 
setwd('D:/hw')
listn<-list.files(getwd(),  pattern=".csv", full.names="FALSE")

knotsn<-as.numeric(best_knot$best_knot_number)
p_c<-ifelse(substr(df_result_all$p_nolinear, 1, 4)=='0.00' | df_result_all$p_nolinear=='<0.001', 
            df_result_all$p_nolinear,
            round(as.numeric(df_result_all$p_nolinear),2))

p_c<-ifelse(p_c=='<0.001', 'p nonlinear <0.001', paste0('p nonlinear = ', p_c))

knotsn_df<-data.frame(knots1=c(0.1,0.5,0.9), knots2=c(0.2, 0.5, 0.8), knots3=c(0.25, 0.5, 0.75), knots4=c(0.10, 0.75, 0.90))

###### 3.1 plotting HWE
# Setup Plotting Environment 
attach(mtcars)
par(mfrow=c(3,4), mar=c(5, 4, 4, 5))

#loop for plotting for each outcome
for (j in 1:length(listn)){
  df <- read.csv(listn[j])
  
  df$smoke3<-as.factor(df$smoke3)
  df$drink3<-as.factor(df$drink3)
  df$vegefruit3<-as.factor(df$vegefruit3)
  df$region3<-as.factor(df$region3)
  
  ###### select exposure variable for loop
  expo_list<-df[, c('hmean90_2d_sum_p1y', 'hmean95_2d_sum_p1y', 'hmean90_3d_sum_p1y','hmean95_3d_sum_p1y',
                    'hmean90_2d_sumdays_p1y','hmean95_2d_sumdays_p1y', 'hmean90_3d_sumdays_p1y','hmean95_3d_sumdays_p1y')]
  expo_list_name<- c('HWE_90th_2d', 'HWE_95th_2d', 'HWE_90th_3d', 'HWE_95th_3d',
                     'HWD_90th_2d', 'HWD_95th_2d', 'HWD_90th_3d', 'HWD_95th_3d')
  
  #loop for plotting for exposure
  for (i in 1:4){
    
    df$expo_hmean<-expo_list[, i]
    
    outnames<-c('A. Type 2 diabetes','B. Dyslipidaemia','C. Hypertension')
    
    # Create basis for non-linear exposure, using
    pp2<-onebasis(df$expo_hmean, fun="ns", knots=quantile(df$expo_hmean, probs = knotsn_df[knotsn[(j-1)*8+i]]))
    # Define Covariates 
    cov <- "+ age + sex + edu2 + yr+ smoke3 + drink3 + vegefruit3 + pa2 + bmi2 + income2 + region3+ pm25 + ndvi"
    # Define model formulas 
    formula<- as.formula(paste0("Surv(tstart, tstop, event) ~ pp2", cov))
    # Fit models 
    model<- coxph(as.formula(formula), data = df, ties = "efron")
    # Predict coefficients and CIs, centered at the median exposure
    pred <- crosspred(pp2, model, model.link="log", cen = median(df$expo_hmean) )
    # Plot the concentration-response curve
    plot(pred, xlab=expo_list_name[i], ylab="Hazard Ratio (95%CI)"  )
    # Overlay a histogram of the exposure variable
    par(new = TRUE)
    hist(df$expo_hmean, col = rgb(0.9, 0.9, 0.9, 0.3), border = NA, axes = FALSE, xlab = "", ylab = "", main = "", breaks = 30)
    axis(side = 4) # Add a right-hand y-axis for the histogram frequency
    mtext("Frequency", side = 4, line = 2, cex = 0.7) # Label the right-hand y-axis
    if(expo_list_name[i]=='HWE_90th_2d'){ 
      mtext(outnames[j], side = 3, line = 2, adj = 0)
    }
    mtext(p_c[(j-1)*8+i], side = 3, line = 0, cex = 0.8, padj = 0)
  }
    
}
  

###### 3.2 plotting HWD
# Setup Plotting Environment 
attach(mtcars)
par(mfrow=c(3,4), mar=c(5, 4, 4, 5))

#loop for plotting for each outcome
for (j in 1:length(listn)){
  df <- read.csv(listn[j])
  
  ###### select exposure variable for loop
  expo_list<-df[, c('hmean90_2d_sum_p1y', 'hmean95_2d_sum_p1y', 'hmean90_3d_sum_p1y','hmean95_3d_sum_p1y',
                    'hmean90_2d_sumdays_p1y','hmean95_2d_sumdays_p1y', 'hmean90_3d_sumdays_p1y','hmean95_3d_sumdays_p1y')]
  expo_list_name<- c('HWE_90th_2d', 'HWE_95th_2d', 'HWE_90th_3d', 'HWE_95th_3d',
                     'HWD_90th_2d', 'HWD_95th_2d', 'HWD_90th_3d', 'HWD_95th_3d')
  
  #loop for plotting for exposure
  for (i in 5:8){
    
    df$expo_hmean<-expo_list[, i]
    
    outnames<-c('A. Type 2 diabetes','B. Dyslipidaemia','C. Hypertension')
    
    # Create basis for non-linear exposure, using
    pp2<-onebasis(df$expo_hmean, fun="ns", knots=quantile(df$expo_hmean, probs = knotsn_df[knotsn[(j-1)*8+i]]))
    # Define Covariates 
    cov <- "+ age + sex + edu2 + yr+ smoke3 + drink3 + vegefruit3 + pa2 + bmi2 + income2 + region3+ pm25 + ndvi"
    # Define model formulas 
    formula<- as.formula(paste0("Surv(tstart, tstop, event) ~ pp2", cov))
    # Fit models 
    model<- coxph(as.formula(formula), data = df, ties = "efron")
    # Predict coefficients and CIs, centered at the median exposure
    pred <- crosspred(pp2, model, model.link="log", cen = median(df$expo_hmean) )
    # Plot the concentration-response curve
    plot(pred, xlab=expo_list_name[i], ylab="Hazard Ratio (95%CI)"  )
    # Overlay a histogram of the exposure variable
    par(new = TRUE)
    hist(df$expo_hmean, col = rgb(0.9, 0.9, 0.9, 0.3), border = NA, axes = FALSE, xlab = "", ylab = "", main = "", breaks = 30)
    axis(side = 4) # Add a right-hand y-axis for the histogram frequency
    mtext("Frequency", side = 4, line = 2, cex = 0.7) # Label the right-hand y-axis
    if(expo_list_name[i]=='HWE_90th_2d'){ 
      mtext(outnames[j], side = 3, line = 2, adj = 0)
      }
    mtext(p_c[(j-1)*8+i], side = 3, line = 0, cex = 0.8, padj = 0)
  }
  
}
