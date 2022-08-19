workflow
================
14 August, 2022

# overview

This is the main document for my master thesis. I wrote “eval=FALSE” so
that it does not overwrite existing data.

# Current Study

The key aim of this paper is to assess if there is a meaningful
difference between using multilevel models and regression-based models
when analyzing dyadic data. A further aim is to investigate and seek a
better way to address categorical variable in the dyadic analysis
framework. To this end, a simulation study can answer these questions.
Convergence and coverage rates, bias, and Type I error rates will be
used to evaluate model performance.

## Simulation design

To address these questions, I generated data in R based on that
hypothetical research scenario. In that hypothetical scenario,
researchers are interested in the effects of education and gender, and
their interaction on health outcomes (DV). Furthermore, in this
scenario, researchers collected sibling data because sibling data allow
researchers to address between- and within-familial effects in this
relationship. In this scenario, MLM or regression-based approach can be
used.

I will be generating dyadic data to address research questions, using
discord (Garrison, Trattner, et al), \[in the final manuscripts, R
packages info and R info will be added here\]. With this generated
dataset, we investigate and compare the multilevel model and the
sibling-comparison model that Garrison and Rodgers (2016) proposed. The
four factors are manipulated for this study: 1) sample size, 2)
intraclass correlation, 3) the main effect of the categorical variable,
and 4) moderation effect of categorical variable.

### Sample Size

Sample size affects estimations. In dyadic analysis, it has been
suggested that at least 50 dyads are needed (when there are no
singletons) to get a reliable and valid estimation in multilevel
modeling (Du & Wang, 2016).

I choose the sample size as 30 dyads (total 60 people), 120 dyads (total
240 people), and 510 dyads (total 1020 people) to address how small,
medium and large sample sizes affect the estimations. Further, Du & Wang
(2016) included 500 dyads (1000 people) as the largest sample size in
their simulation study. We chose this size because these dyad numbers
are multiple of 6; we have models that have categorical variable with
two levels, and three levels (it is discussed in the following section),
and we tried to make a balanced sample. the effect of unbalanced sample
design is not a focus of the interest in this study. Thus, the dyads
numbers total should be multiple of 6.

### Intraclass Correlation (ICC).

As discussed in the introduction, ICC influences the estimates.
Generally, ICC denotes the how much variance can be explained by dyads.
We identify the ICC from the null model, as discussed in introduction .

It is worth noting here is that as group size gets smaller, ICC tends to
be higher (Hox, 2010) in MLM. In Du & Wang (2016), 0.1, 0.2, 0.3, 0.5,
and 0.7 of ICC were used to simulate dyadic analysis using MLM. In this
work, Du & Wang (2016) found that when ICC is equal to or less than 0.2
(ICC≤0.2), more convergence issues may occurs in dyadic analysis using
MLM. Thus, we included 0.2 as a minimum level of ICC because of this
finding.

Further, 0.4 and 0.8 of ICC can be found in the literature as well. For
example, 0.43 of ICC was observed among opposite-gender sibling pairs in
weight (Raskdyad2ind et al., 2018). Furthermore, ICC sometimes can be
over 0.8; 0.8 of ICC is observed in romantic dyads (McIsaac et al.,
2008).

### Main effect and Interaction effect of Categorical Variable.

The main effect of the categorical variable is also included in our
study. In our hypothetical scenario, researchers try to predict health
outcome (Y) with education (X), Gender (S), and interaction between
income and gender (XS) at dyad2individual levels.

$ Y =  + 1X + 2S + β3XS + e$

At the dyad2individual level, I set the standardized beta coefficient of
a categorical variable (Gender) to be 0.1, 0.3, and 0.5, as conventional
wisdom suggests as small, medium, and large effect size (Cohen, 1992) .
Similarly, I also set the beta weight of interaction term (gender \*
education) to be 0.1, 0.3, and 0.5.

### Other Settings

I set that the other continuous IV (e.g., education) to have a 0.3
standardized beta weight (moderate association), because our research
question in this current study is focused on categorical variable.
Furthermore, to pursue parsimoniously, I assumed that there are no
singletons nor missing data in our simulation study because this issue
is not our focus of interest. The further limitations derived from these
settings is discussed in the discussion section.

### condition summary

## Condition

1)  sample size: cluster size 30, 120, and 510.

2)  ICC (from the null model)

0.2, 0.4, 0.8

3)  y = beta0 + beta1*edu + beta2*gender+ beta3*edu*gender

beta1= 0.3 beta2= 0.1, 0.3 or 0.5 ( and 0 when beta3 is 0 ) beta3= 0,
0.1, 0.3, or 0.5

# Code

## libraries

``` r
library(tidyverse)
library(MASS)
library(jtools) #to get nicer output
#library(lme4)
library(lmerTest)
#lme4's lmer() is masked by lmerTest's lmer() function, and
#lmerTest's lmer function is more conservative when determining model convergence. 
library(car)

library(multilevel)
library(nlme)
library(simstudy)
library(MuMIn) # To calculate R^2 for MLM models 
#library(r2mlm) # To calculate R^2 for MLM models -want to discuss 

set.seed(2022)
```

## Equations

### ICC

![ρ=(σ\_(u_0)^2)/(σ\_(u_0)^2+σ_e^2 )](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%CF%81%3D%28%CF%83_%28u_0%29%5E2%29%2F%28%CF%83_%28u_0%29%5E2%2B%CF%83_e%5E2%20%29 "ρ=(σ_(u_0)^2)/(σ_(u_0)^2+σ_e^2 )")
In the null model. If it is not null model, I call it pseudo ICC.

### when Generate the data, we follow the following regression equation:

y = beta1*(edu btw + edu within) + beta2 (gender btw + gender within ) +
beta3 (gender with * edu with ) + error_btw + error within

## functions

``` r
# Data generation function (dyad structure)
source('scripts/gen_dyadic.R')
source('scripts/gen_rmvn.R')

# function to convert dyad structure into dyad2individual structure
source('scripts/dyad_to_long.R')

# function to convert plain dyad structure data to the dyad structure data that is suitable for kinship-discordant model

source('scripts/reconstruct.R')

# the function for sibling comparison model
source('scripts/model_sib.R')

# the function for null model
source('scripts/model_null.R')

# the function for mlm model
source('scripts/model_mlm.R')

# the function for ICC weight calculator
source('scripts/ICC_calculator.R')



# the function for Catching  warning message

source('scripts/myTryCatch.R')

# the function for convergence (need more suspection)

source('scripts/hasConverged.R')
```

# condition

1)  sample size: cluster size 30, 120, and 510. (3 conditions in sample
    size)

2)  ICC (from the null model)

0.2, 0.4, 0.8 (3 conditions in ICC)

3)  y = beta0 + beta1*edu + beta2*gender+ beta3*edu*gender

beta1= 0.3 ( 1 condition in beta1) beta2= 0, 0.1, 0.3 or 0.5 (3
conditions in beta2) beta3= 0, 0.1, 0.3, or 0.5 (4 conditions in beta3)

but if if beta3 is 0 beta 2 is all 0.

#### Conditions

``` r
### Conditions

# cluster number 
cn <- c(30,120,510)
# ICC
target_ICC <- c(0.2,0.4,0.8)

#beta1
beta1 <- 0.3

# beta2 
beta2 <- c(0.1,0.3,0.5)
beta2_2 <- c(0, 0.1, 0.3, 0.5) # only when the beta3 is  0 

# beta3

beta3 <- c(0,0.1,0.3,0.5)


# 
# ncond <- length(cn)*length(target_ICC)*length(beta1)*length(beta2)*length(beta3)

conditions <- expand.grid(cn,target_ICC,beta1,beta2,beta3)
names(conditions) <- c('cn','target_ICC','beta1','beta2','beta3')

# when the beta3 is 0, beta2 have additional conditions (0)

conditions <- rbind(conditions, expand.grid(cn=c(30,120,510), target_ICC=c(0.2,0.4,0.8),
                              beta1=0.3, beta2=0, beta3=0))

ncond <- nrow(conditions)

#conditions$nrep <- 1 #number of simulations per condition

conditions$wtn <- 100 # variance weight



conditions$condition <- 1:nrow(conditions)


conditions$btw <- ICC_calculator(  
   beta1 = conditions$beta1,
   beta2 = conditions$beta2,
   beta3 = conditions$beta3,
   ICC = conditions$target_ICC,
   within = conditions$wtn)

 # within ICC weights 
# arbitrary choice. it will be multiplied by 0.01 ~0.09 so, and err variance is 1, so make it comparable. 

# conditions_summary<-conditions

## save the condition summary table 
# 
write.csv(conditions,"conditions/Conditions_summary.csv")
```

# data generation & get results

``` r
# nuber of replication for each condition 
nrep <- 1000

# measuring a time

start_time <- Sys.time()

# two loops one for reps

results <- data.frame(nrep = 1:nrep)

#setseed

set.seed(2022)

for (j in 1:ncond) {
print(j)
  
# make a result data frame 

# one row per one individual data set's results. 

for (i in 1:nrep) {
  

  # data generation
  
   sim_data <- gen_dyadic(
    cn = conditions$cn[j],
    beta1 = conditions$beta1[j],
    beta2 = conditions$beta2[j],
    beta3 = conditions$beta3[j],
    con_icc_wtn = conditions$wtn[j],
    con_icc_btw = conditions$btw[j],
    cat_icc_wtn = conditions$wtn[j],
    cat_icc_btw = conditions$btw[j])
  
# fit the model 
   
sib_data <- reconstruct(sim_data)# reconstruct dyadic data to dyadic data, but suitable for kinship discordant model 
  

   ## Sib0
   ## most simple model in kinship-discordant model. 
   ##lm(y_diff  ~ y_mean + con_mean + con_diff + gender_composition_two_eff, data=dat3)


 results$warning_sib0[i] <- tryc(myTryCatch(
   sib0 <- model_sib(sib_data,
             gender_composition=2,
             effect_code=TRUE,
             interaction=FALSE,

   ## sib 1?
## gender composition is thre, and no interaction model 
## sib_model<- lm(y_diff  ~ y_mean + con_mean + con_diff + ev1 + ev2, data=dat3) // ev = effect code variable 
# 
results$warning_sib1[i] <- tryc(myTryCatch(
   sib1 <- model_sib(sib_data,
             gender_composition=3,
             effect_code=TRUE,
             interaction=FALSE,
             return_data=FALSE)))



   ## sib 2?
## lm(y_diff  ~ y_mean + con_mean + con_diff + gender_composition_two_effsim
#                     , data=dat3)



results$warning_sib2[i] <- tryc(myTryCatch(
                                sib2 <- model_sib(sib_data,
              gender_composition = 2,
              effect_code = TRUE,
              interaction = TRUE,
              return_data = FALSE)))
 
## what is sib 3?
## lm(y_diff  ~ y_mean + con_mean + con_diff + ev1 + ev2 +
 #                      ev1*con_diff + ev2*con_diff, data=dat3)

results$warning_sib3[i]<-tryc(myTryCatch(sib3<-model_sib(sib_data,
              gender_composition=3,
              effect_code=TRUE,
              interaction=TRUE,
              return_data=FALSE)))


# change the data into individual structure so that it can fit into multilevel models
ind_sim_data<-dyad2ind(sim_data)

# null model (Multilevel model)
results$warning_MLM0[i]<-tryc(
  myTryCatch(
    MLM0<-model_null(ind_sim_data))
  )

# MLM1: NO level 2 model without interaction
#lmer(y~1 + catwtn + con +  (1|pid), data=data)

results$warning_MLM1[i]<-tryc(myTryCatch(
                          MLM1<-no2model(ind_sim_data,
                          interaction=FALSE))
 
# MLM2: No level 2 model with interaction 
#lmer(y~1 + catwtn + con + catwtn*con + (1|pid), data=data)


results$warning_MLM2[i]<-tryc(
                          myTryCatch(
                                MLM2<-no2model(ind_sim_data,
                    interaction=TRUE)))
# MLM2<-no2model(sim_data,
#                interaction=TRUE)

# 
# MLM3<- yes2model(sim_data,
#                 gender_composition=2,
#                 interaction=FALSE)
# 
results$warning_MLM3[i]<-tryc(myTryCatch(
                MLM3<-yes2model(ind_sim_data,
                gender_composition=2,
                interaction=FALSE)))

# MLM4<-yes2model(sim_data,
#           gender_composition=2,
#           interaction=TRUE)

# See the model cheatsheet 


# 
 results$warning_MLM4[i]<-tryc(myTryCatch(
MLM4<-yes2model(ind_sim_data,
          gender_composition=2,
          interaction=TRUE)
 ))


# 
# MLM5<- yes2model(sim_data,
#                 gender_composition=3,
#                 interaction=FALSE)
# 
# 
# See the model cheatsheet
 
results$warning_MLM5[i]<-tryc(myTryCatch(
MLM5<- yes2model(data=ind_sim_data,
                gender_composition=3,
                interaction=FALSE)))
# 

# MLM6<- yes2model(sim_data,
#                 gender_composition=3,
#                 interaction=TRUE)
# See the model cheatsheet
# 
results$warning_MLM6[i]<-tryc(myTryCatch(
  MLM6 <-yes2model(ind_sim_data,
                gender_composition=3,
                interaction=TRUE)))

# hasConverged :this function determines whether the model has converged or not. 

results$hasConverged_MLM0[i]<-hasConverged(MLM0)
results$hasConverged_MLM1[i]<-hasConverged(MLM1)
results$hasConverged_MLM2[i]<-hasConverged(MLM2)
results$hasConverged_MLM3[i]<-hasConverged(MLM3)
results$hasConverged_MLM4[i]<-hasConverged(MLM4)
results$hasConverged_MLM5[i]<-hasConverged(MLM5)
results$hasConverged_MLM6[i]<-hasConverged(MLM6)


# descriptive statistics of DV

results$ind_y_mean[i]<-mean(ind_sim_data$y,na.rm=TRUE) 
results$ind_y_var[i]<-var(ind_sim_data$y,na.rm=TRUE)
results$dyad_y_1_mean[i]<-mean(sib_data$y_1,na.rm=TRUE)
results$dyad_y_2_mean[i]<-mean(sib_data$y_2,na.rm=TRUE)
results$dyad_con_1_mean[i]<-mean(sib_data$con_1,na.rm=TRUE)
results$dyad_con_2_mean[i]<-mean(sib_data$con_2,na.rm=TRUE)
results$ind_con_mean[i]<-mean(ind_sim_data$con,na.rm=TRUE)
results$ind_con_var[i]<-var(ind_sim_data$y,na.rm=TRUE)

results$n_three_level_gender_mimi[i]<-sum(sib_data$gender_composition_three=="mimi")
results$n_three_level_gender_mixed[i]<-sum(sib_data$gender_composition_three=="mixed")
results$n_three_level_gender_oneone[i]<-sum(sib_data$gender_composition_three=="oneone")


## ICC from the null model

results$ICC_MLM0[i]<-unlist(jtools::summ(MLM0)[3])[3]

## pseudo-ICC from non-null multilevel model 


results$pseudo_ICC_MLM1[i]<-unlist(jtools::summ(MLM1)[3])[3]
results$pseudo_ICC_MLM2[i]<-unlist(jtools::summ(MLM2)[3])[3]
results$pseudo_ICC_MLM3[i]<-unlist(jtools::summ(MLM3)[3])[3]
results$pseudo_ICC_MLM4[i]<-unlist(jtools::summ(MLM4)[3])[3]
results$pseudo_ICC_MLM5[i]<-unlist(jtools::summ(MLM5)[3])[3]
results$pseudo_ICC_MLM6[i]<-unlist(jtools::summ(MLM6)[3])[3]
# deviance?

results$Dev_sib0[i]<-deviance(sib0)
results$Dev_sib1[i]<-deviance(sib1)
results$Dev_sib2[i]<-deviance(sib2)
results$Dev_sib3[i]<-deviance(sib3)


#REML MLM cannot use deviance
#deviance() is deprecated for REML fits; use REMLcrit for the REML
results$REMLcriterion_MLM0[i] <- tryc(summary(MLM0)[[3]]$cmp[7])
results$REMLcriterion_MLM1[i] <- tryc(summary(MLM1)[[3]]$cmp[7])
results$REMLcriterion_MLM2[i] <- tryc(summary(MLM2)[[3]]$cmp[7])
results$REMLcriterion_MLM3[i] <- tryc(summary(MLM3)[[3]]$cmp[7])
results$REMLcriterion_MLM4[i] <- tryc(summary(MLM4)[[3]]$cmp[7])
results$REMLcriterion_MLM5[i] <- tryc(summary(MLM5)[[3]]$cmp[7])
results$REMLcriterion_MLM6[i]<- tryc(summary(MLM6)[[3]]$cmp[7])

# AIC


results$AIC_sib0[i]<-AIC(sib0)
results$AIC_sib1[i]<-AIC(sib1)
results$AIC_sib2[i]<-AIC(sib2)
results$AIC_sib3[i]<-AIC(sib3)


#BIC 
results$BIC_sib0[i]<-BIC(sib0)
results$BIC_sib1[i]<-BIC(sib1)
results$BIC_sib2[i]<-BIC(sib2)
results$BIC_sib3[i]<-BIC(sib3)


# residual standard error 

#Residual standard error: ResSE on ResSE_df degrees of freedom
results$ResSE_sib0[i]<- tryc(summary(sib0)$sigma)
results$ResSE_sib1[i]<- tryc(summary(sib1)$sigma)
results$ResSE_sib2[i]<- tryc(summary(sib2)$sigma)
results$ResSE_sib3[i]<- tryc(summary(sib3)$sigma)
## residual standard error df 
results$ResSE_df_sib0[i]<- tryc(summary(sib0)[["df"]][[2]])
results$ResSE_df_sib1[i]<- tryc(summary(sib1)[["df"]][[2]])
results$ResSE_df_sib2[i]<- tryc(summary(sib2)[["df"]][[2]])
results$ResSE_df_sib3[i]<- tryc(summary(sib3)[["df"]][[2]])

# F statistics

## F-statistic: F_sib on F_numdf and F_denDF

results$F_sib0[i]<- tryc(summary(sib0)[["fstatistic"]][["value"]])
results$F_sib1[i]<- tryc(summary(sib1)[["fstatistic"]][["value"]])
results$F_sib2[i]<- tryc(summary(sib2)[["fstatistic"]][["value"]])
results$F_sib3[i]<- tryc(summary(sib3)[["fstatistic"]][["value"]])

results$F_numdf_sib0[i]<- tryc(summary(sib0)[["fstatistic"]][["numdf"]])
results$F_numdf_sib1[i]<- tryc(summary(sib1)[["fstatistic"]][["numdf"]])
results$F_numdf_sib2[i]<- tryc(summary(sib2)[["fstatistic"]][["numdf"]])
results$F_numdf_sib3[i]<- tryc(summary(sib3)[["fstatistic"]][["numdf"]])

results$F_dendf_sib0[i]<- tryc(summary(sib0)[["fstatistic"]][["dendf"]])
results$F_dendf_sib1[i]<- tryc(summary(sib1)[["fstatistic"]][["dendf"]])
results$F_dendf_sib2[i]<- tryc(summary(sib2)[["fstatistic"]][["dendf"]])
results$F_dendf_sib3[i]<- tryc(summary(sib3)[["fstatistic"]][["dendf"]])

# # p-value
# 
 results$f_p_sib0[i]<- pf(summary(sib0)$fstatistic[1L], summary(sib0)$fstatistic[2L], summary(sib0)$fstatistic[3L], lower.tail = FALSE)
# 
 results$f_p_sib1[i]<- pf(summary(sib1)$fstatistic[1L], summary(sib1)$fstatistic[2L], summary(sib1)$fstatistic[3L], lower.tail = FALSE)
# 
 results$f_p_sib2[i]<- pf(summary(sib2)$fstatistic[1L], summary(sib2)$fstatistic[2L], summary(sib2)$fstatistic[3L], lower.tail = FALSE)
# 
 results$f_p_sib3[i]<- pf(summary(sib3)$fstatistic[1L], summary(sib3)$fstatistic[2L], summary(sib3)$fstatistic[3L], lower.tail = FALSE)

## Type2 Wald Chisquare test - not recommended in dyadic study 
 # so commentised and did not include. 
 
 
# null model does not have Iv so... no wald test.
# results$wald_catwtn_chisq_MLM1[i]<-car::Anova(MLM1)['catwtn','Chisq']
# results$wald_con_chisq_MLM1[i]<-car::Anova(MLM1)['con','Chisq']
# results$wald_catwtn_df_MLM1[i]<-car::Anova(MLM1)['catwtn','Df']
# results$wald_con_df_MLM1[i]<-car::Anova(MLM1)['con','Df']
# results$wald_catwtn_p_MLM1[i]<-car::Anova(MLM1)['catwtn','Pr(>Chisq)']
# results$wald_con_p_MLM1[i]<-car::Anova(MLM1)['con','Pr(>Chisq)']
# 
# results$wald_catwtn_chisq_MLM2[i]<-car::Anova(MLM2)['catwtn','Chisq']
# results$wald_con_chisq_MLM2[i]<-car::Anova(MLM2)['con','Chisq']
# results$wald_catwtnxcon_chisq_MLM2[i]<-car::Anova(MLM2)['catwtn:con','Chisq']
# results$wald_catwtn_df_MLM2[i]<-car::Anova(MLM2)['catwtn','Df']
# results$wald_con_df_MLM2[i]<-car::Anova(MLM2)['con','Df']
# results$wald_catwtnxcon_df_MLM2[i]<-car::Anova(MLM2)['catwtn:con','Df']
# results$wald_catwtn_p_MLM2[i]<-car::Anova(MLM2)['catwtn','Pr(>Chisq)']
# results$wald_con_p_MLM2[i]<-car::Anova(MLM2)['con','Pr(>Chisq)']
# results$wald_catwtnxcon_p_MLM2[i]<-car::Anova(MLM2)['catwtn:con','Pr(>Chisq)']
# 
# results$wald_catbtw_chisq_MLM3[i]<-car::Anova(MLM3)['catbtw','Chisq']
# results$wald_con_chisq_MLM3[i]<-car::Anova(MLM3)['con','Chisq']
# results$wald_catbtw_df_MLM3[i]<-car::Anova(MLM3)['catbtw','Df']
# results$wald_con_df_MLM3[i]<-car::Anova(MLM3)['con','Df']
# results$wald_catbtw_p_MLM3[i]<-car::Anova(MLM3)['catbtw','Pr(>Chisq)']
# results$wald_con_p_MLM3[i]<-car::Anova(MLM3)['con','Pr(>Chisq)']
# 
# 
# results$wald_catbtw_chisq_MLM4[i]<-car::Anova(MLM4)['catbtw','Chisq']
# results$wald_con_chisq_MLM4[i]<-car::Anova(MLM4)['con','Chisq']
# results$wald_conxcatbtw_chisq_MLM4[i]<-car::Anova(MLM4)['con:catbtw','Chisq']
# results$wald_catbtw_df_MLM4[i]<-car::Anova(MLM4)['catbtw','Df']
# results$wald_con_df_MLM4[i]<-car::Anova(MLM4)['con','Df']
# results$wald_conxcatbtw_df_MLM4[i]<-car::Anova(MLM4)['con:catbtw','Df']
# results$wald_catbtw_p_MLM4[i]<-car::Anova(MLM4)['catbtw','Pr(>Chisq)']
# results$wald_con_p_MLM4[i]<-car::Anova(MLM4)['con','Pr(>Chisq)']
# results$wald_conxcatbtw_p_MLM4[i]<-car::Anova(MLM4)['con:catbtw','Pr(>Chisq)']
# 
# 
# results$wald_con_chisq_MLM5[i]<-car::Anova(MLM5)['con','Chisq']
# results$wald_ev1_chisq_MLM5[i]<-car::Anova(MLM5)['ev1','Chisq']
# results$wald_ev2_chisq_MLM5[i]<-car::Anova(MLM5)['ev2','Chisq']
# results$wald_con_df_MLM5[i]<-car::Anova(MLM5)['con','Df']
# results$wald_ev1_df_MLM5[i]<-car::Anova(MLM5)['ev1','Df']
# results$wald_ev2_df_MLM5[i]<-car::Anova(MLM5)['ev2','Df']
# results$wald_con_p_MLM5[i]<-car::Anova(MLM5)['con','Pr(>Chisq)']
# results$wald_ev1_p_MLM5[i]<-car::Anova(MLM5)['ev1','Pr(>Chisq)']
# results$wald_ev2_p_MLM5[i]<-car::Anova(MLM5)['ev2','Pr(>Chisq)']
# 
# 
# results$wald_con_chisq_MLM6[i]<-car::Anova(MLM6)['con','Chisq']
# results$wald_ev1_chisq_MLM6[i]<-car::Anova(MLM6)['ev1','Chisq']
# results$wald_ev2_chisq_MLM6[i]<-car::Anova(MLM6)['ev2','Chisq']
# results$wald_conxev1_chisq_MLM6[i]<-car::Anova(MLM6)['con:ev1','Chisq']
# results$wald_conxev2_chisq_MLM6[i]<-car::Anova(MLM6)['con:ev2','Chisq']
# results$wald_con_df_MLM6[i]<-car::Anova(MLM6)['con','Df']
# results$wald_ev1_df_MLM6[i]<-car::Anova(MLM6)['ev1','Df']
# results$wald_ev2_df_MLM6[i]<-car::Anova(MLM6)['ev2','Df']
# results$wald_conxev1_df_MLM6[i]<-car::Anova(MLM6)['con:ev1','Df']
# results$wald_conxev2_df_MLM6[i]<-car::Anova(MLM6)['con:ev2','Df']
# results$wald_con_p_MLM6[i]<-car::Anova(MLM6)['con','Pr(>Chisq)']
# results$wald_ev1_p_MLM6[i]<-car::Anova(MLM6)['ev1','Pr(>Chisq)']
# results$wald_ev2_p_MLM6[i]<-car::Anova(MLM6)['ev2','Pr(>Chisq)']
# results$wald_conxev1_p_MLM6[i]<-car::Anova(MLM6)['con:ev1','Pr(>Chisq)']
# results$wald_conxev2_p_MLM6[i]<-car::Anova(MLM6)['con:ev2','Pr(>Chisq)']

# residual standard deviation
results$ResSD_MLM0[i]<- summary(MLM0)$sigma
results$ResSD_MLM1[i]<- summary(MLM1)$sigma
results$ResSD_MLM2[i]<- summary(MLM2)$sigma
results$ResSD_MLM3[i]<- summary(MLM3)$sigma
results$ResSD_MLM4[i]<- summary(MLM4)$sigma
results$ResSD_MLM5[i]<- summary(MLM5)$sigma
results$ResSD_MLM6[i]<- summary(MLM6)$sigma

## MLM's random effect. 

results$ran_pid_intercept_MLM0[i]<-attr(VarCorr(MLM0)$pid,"stddev")
results$ran_pid_intercept_MLM1[i]<-attr(VarCorr(MLM1)$pid,"stddev")
results$ran_pid_intercept_MLM2[i]<-attr(VarCorr(MLM2)$pid,"stddev")
results$ran_pid_intercept_MLM3[i]<-attr(VarCorr(MLM3)$pid,"stddev")
results$ran_pid_intercept_MLM4[i]<-attr(VarCorr(MLM4)$pid,"stddev")
results$ran_pid_intercept_MLM5[i]<-attr(VarCorr(MLM5)$pid,"stddev")
results$ran_pid_intercept_MLM6[i]<-attr(VarCorr(MLM6)$pid,"stddev")

### MLM singularity 

results$singularity_MLM0[i]<-isSingular(MLM0)
results$singularity_MLM1[i]<-isSingular(MLM1)
results$singularity_MLM2[i]<-isSingular(MLM2)
results$singularity_MLM3[i]<-isSingular(MLM3)
results$singularity_MLM4[i]<-isSingular(MLM4)
results$singularity_MLM5[i]<-isSingular(MLM5)
results$singularity_MLM6[i]<-isSingular(MLM6)

# R squred

## Rsquare for kinship discordant model 

results$rsquared_sib0[i]<-summary(sib0)$r.squared
results$adj_rsquared_sib0[i]<-summary(sib0)$adj.r.squared
results$rsquared_sib1[i]<-summary(sib1)$r.squared
results$adj_rsquared_sib1[i]<-summary(sib1)$adj.r.squared
results$rsquared_sib2[i]<-summary(sib2)$r.squared
results$adj_rsquared_sib2[i]<-summary(sib2)$adj.r.squared
results$rsquared_sib3[i]<-summary(sib3)$r.squared
results$adj_rsquared_sib3[i]<-summary(sib3)$adj.r.squared

## r square for MLM 

results$marg_rsquared_MLM0[i]<-MuMIn::r.squaredGLMM(MLM0)[,'R2m']
results$cond_rsquared_MLM0[i]<-MuMIn::r.squaredGLMM(MLM0)[,'R2c']


results$marg_rsquared_MLM1[i]<-MuMIn::r.squaredGLMM(MLM1)[,'R2m']
results$cond_rsquared_MLM1[i]<-MuMIn::r.squaredGLMM(MLM1)[,'R2c']


results$marg_rsquared_MLM2[i]<-MuMIn::r.squaredGLMM(MLM2)[,'R2m']
results$cond_rsquared_MLM2[i]<-MuMIn::r.squaredGLMM(MLM2)[,'R2c']


results$marg_rsquared_MLM3[i]<-MuMIn::r.squaredGLMM(MLM3)[,'R2m']
results$cond_rsquared_MLM3[i]<-MuMIn::r.squaredGLMM(MLM3)[,'R2c']


results$marg_rsquared_MLM4[i]<-MuMIn::r.squaredGLMM(MLM4)[,'R2m']
results$cond_rsquared_MLM4[i]<-MuMIn::r.squaredGLMM(MLM4)[,'R2c']


results$marg_rsquared_MLM5[i]<-MuMIn::r.squaredGLMM(MLM5)[,'R2m']
results$cond_rsquared_MLM5[i]<-MuMIn::r.squaredGLMM(MLM5)[,'R2c']


results$marg_rsquared_MLM6[i]<-MuMIn::r.squaredGLMM(MLM6)[,'R2m']
results$cond_rsquared_MLM6[i]<-MuMIn::r.squaredGLMM(MLM6)[,'R2c']


# RMSE (NOT scaled RMSE)
#however, see: https://stackoverflow.com/questions/43123462/how-to-obtain-rmse-out-of-lm-result
# Now I think we do not need this becuase we know the true parameters 
results$RMSE_sib0[i]<-sqrt(mean(residuals(sib0)^2))
results$RMSE_sib1[i]<-sqrt(mean(residuals(sib1)^2))
results$RMSE_sib2[i]<-sqrt(mean(residuals(sib2)^2))
results$RMSE_sib3[i]<-sqrt(mean(residuals(sib3)^2))

results$RMSE_MLM0[i]<-sqrt(mean(residuals(MLM0)^2))
results$RMSE_MLM1[i]<-sqrt(mean(residuals(MLM1)^2))
results$RMSE_MLM2[i]<-sqrt(mean(residuals(MLM2)^2))
results$RMSE_MLM3[i]<-sqrt(mean(residuals(MLM3)^2))
results$RMSE_MLM4[i]<-sqrt(mean(residuals(MLM4)^2))
results$RMSE_MLM5[i]<-sqrt(mean(residuals(MLM5)^2))
results$RMSE_MLM6[i]<-sqrt(mean(residuals(MLM6)^2))

# coefficients 

## sib0

results$coef_intercept_sib0[i]<-tryc(coef(summary(sib0))["(Intercept)","Estimate"])
results$coef_intercept_se_sib0[i]<-tryc(coef(summary(sib0))["(Intercept)","Std. Error"])
results$coef_intercept_t_sib0[i]<-tryc(coef(summary(sib0))["(Intercept)","t value"])
results$coef_intercept_p_sib0[i]<-tryc(coef(summary(sib0))["(Intercept)","Pr(>|t|)"])
results$coef_y_mean_sib0[i]<-tryc(coef(summary(sib0))["y_mean","Estimate"])
results$coef_y_mean_se_sib0[i]<-tryc(coef(summary(sib0))["y_mean","Std. Error"])
results$coef_y_mean_t_sib0[i]<-tryc(coef(summary(sib0))["y_mean","t value"])
results$coef_y_mean_p_sib0[i]<-tryc(coef(summary(sib0))["y_mean","Pr(>|t|)"])
results$coef_con_mean_sib0[i]<-tryc(coef(summary(sib0))["con_mean","Estimate"])
results$coef_con_mean_se_sib0[i]<-tryc(coef(summary(sib0))["con_mean","Std. Error"])
results$coef_con_mean_t_sib0[i]<-tryc(coef(summary(sib0))["con_mean","t value"])
results$coef_con_mean_p_sib0[i]<-tryc(coef(summary(sib0))["con_mean","Pr(>|t|)"])
results$coef_con_diff_sib0[i]<-tryc(coef(summary(sib0))["con_diff","Estimate"])
results$coef_con_diff_se_sib0[i]<-tryc(coef(summary(sib0))["con_diff","Std. Error"])
results$coef_con_diff_t_sib0[i]<-tryc(coef(summary(sib0))["con_diff","t value"])
results$coef_con_diff_p_sib0[i]<-tryc(coef(summary(sib0))["con_diff","Pr(>|t|)"])
results$coef_gender_composition_two_eff_sib0[i]<-tryc(coef(summary(sib0))["gender_composition_two_eff","Estimate"])
results$coef_gender_composition_two_eff_se_sib0[i]<-tryc(coef(summary(sib0))["gender_composition_two_eff","Std. Error"])
results$coef_gender_composition_two_eff_t_sib0[i]<-tryc(coef(summary(sib0))["gender_composition_two_eff","t value"])
results$coef_gender_composition_two_eff_p_sib0[i]<-tryc(coef(summary(sib0))["gender_composition_two_eff","Pr(>|t|)"])

## sib1

results$coef_intercept_sib1[i]<-tryc(coef(summary(sib1))["(Intercept)","Estimate"])
results$coef_intercept_se_sib1[i]<-tryc(coef(summary(sib1))["(Intercept)","Std. Error"])
results$coef_intercept_t_sib1[i]<-tryc(coef(summary(sib1))["(Intercept)","t value"])
results$coef_intercept_p_sib1[i]<-tryc(coef(summary(sib1))["(Intercept)","Pr(>|t|)"])

results$coef_y_mean_sib1[i]<-tryc(coef(summary(sib1))["y_mean","Estimate"])
results$coef_y_mean_se_sib1[i]<-tryc(coef(summary(sib1))["y_mean","Std. Error"])
results$coef_y_mean_t_sib1[i]<-tryc(coef(summary(sib1))["y_mean","t value"])
results$coef_y_mean_p_sib1[i]<-tryc(coef(summary(sib1))["y_mean","Pr(>|t|)"])

results$coef_con_mean_sib1[i]<-tryc(coef(summary(sib1))["con_mean","Estimate"])
results$coef_con_mean_se_sib1[i]<-tryc(coef(summary(sib1))["con_mean","Std. Error"])
results$coef_con_mean_t_sib1[i]<-tryc(coef(summary(sib1))["con_mean","t value"])
results$coef_con_mean_p_sib1[i]<-tryc(coef(summary(sib1))["con_mean","Pr(>|t|)"])
results$coef_con_diff_sib1[i]<-tryc(coef(summary(sib1))["con_diff","Estimate"])
results$coef_con_diff_se_sib1[i]<-tryc(coef(summary(sib1))["con_diff","Std. Error"])
results$coef_con_diff_t_sib1[i]<-tryc(coef(summary(sib1))["con_diff","t value"])
results$coef_con_diff_p_sib1[i]<-tryc(coef(summary(sib1))["con_diff","Pr(>|t|)"])

results$coef_ev1_sib1[i]<-tryc(coef(summary(sib1))["ev1","Estimate"])
results$coef_ev1_se_sib1[i]<-tryc(coef(summary(sib1))["ev1","Std. Error"])
results$coef_ev1_t_sib1[i]<-tryc(coef(summary(sib1))["ev1","t value"])
results$coef_ev1_p_sib1[i]<-tryc(coef(summary(sib1))["ev1","Pr(>|t|)"])
results$coef_ev2_sib1[i]<-tryc(coef(summary(sib1))["ev2","Estimate"])
results$coef_ev2_se_sib1[i]<-tryc(coef(summary(sib1))["ev2","Std. Error"])
results$coef_ev2_t_sib1[i]<-tryc(coef(summary(sib1))["ev2","t value"])
results$coef_ev2_p_sib1[i]<-tryc(coef(summary(sib1))["ev2","Pr(>|t|)"])

## sib2 coefficient


results$coef_intercept_sib2[i]<-tryc(coef(summary(sib2))["(Intercept)","Estimate"])
results$coef_intercept_se_sib2[i]<-tryc(coef(summary(sib2))["(Intercept)","Std. Error"])
results$coef_intercept_t_sib2[i]<-tryc(coef(summary(sib2))["(Intercept)","t value"])
results$coef_intercept_p_sib2[i]<-tryc(coef(summary(sib2))["(Intercept)","Pr(>|t|)"])
results$coef_y_mean_sib2[i]<-tryc(coef(summary(sib2))["y_mean","Estimate"])
results$coef_y_mean_se_sib2[i]<-tryc(coef(summary(sib2))["y_mean","Std. Error"])
results$coef_y_mean_t_sib2[i]<-tryc(coef(summary(sib2))["y_mean","t value"])
results$coef_y_mean_p_sib2[i]<-tryc(coef(summary(sib2))["y_mean","Pr(>|t|)"])
results$coef_con_mean_sib2[i]<-tryc(coef(summary(sib2))["con_mean","Estimate"])
results$coef_con_mean_se_sib2[i]<-tryc(coef(summary(sib2))["con_mean","Std. Error"])
results$coef_con_mean_t_sib2[i]<-tryc(coef(summary(sib2))["con_mean","t value"])
results$coef_con_mean_p_sib2[i]<-tryc(coef(summary(sib2))["con_mean","Pr(>|t|)"])
results$coef_con_diff_sib2[i]<-tryc(coef(summary(sib2))["con_diff","Estimate"])
results$coef_con_diff_se_sib2[i]<-tryc(coef(summary(sib2))["con_diff","Std. Error"])
results$coef_con_diff_t_sib2[i]<-tryc(coef(summary(sib2))["con_diff","t value"])
results$coef_con_diff_p_sib2[i]<-tryc(coef(summary(sib2))["con_diff","Pr(>|t|)"])
results$coef_gender_composition_two_eff_sib2[i]<-tryc(coef(summary(sib2))["gender_composition_two_eff","Estimate"])
results$coef_gender_composition_two_eff_se_sib2[i]<-tryc(coef(summary(sib2))["gender_composition_two_eff","Std. Error"])
results$coef_gender_composition_two_eff_t_sib2[i]<-tryc(coef(summary(sib2))["gender_composition_two_eff","t value"])
results$coef_gender_composition_two_eff_p_sib2[i]<-tryc(coef(summary(sib2))["gender_composition_two_eff","Pr(>|t|)"])

results$coef_con_diffxgender_composition_two_eff_sib2[i]<-tryc(coef(summary(sib2))["con_diff:gender_composition_two_eff","Estimate"])
results$coef_con_diffxgender_composition_two_eff_se_sib2[i]<-tryc(coef(summary(sib2))["con_diff:gender_composition_two_eff","Std. Error"])
results$coef_con_diffxgender_composition_two_eff_t_sib2[i]<-tryc(coef(summary(sib2))["con_diff:gender_composition_two_eff","t value"])
results$coef_con_diffxgender_composition_two_eff_p_sib2[i]<-tryc(coef(summary(sib2))["con_diff:gender_composition_two_eff","Pr(>|t|)"])

## sib3 coefficient


results$coef_intercept_sib3[i]<-tryc(coef(summary(sib3))["(Intercept)","Estimate"])
results$coef_intercept_se_sib3[i]<-tryc(coef(summary(sib3))["(Intercept)","Std. Error"])
results$coef_intercept_t_sib3[i]<-tryc(coef(summary(sib3))["(Intercept)","t value"])
results$coef_intercept_p_sib3[i]<-tryc(coef(summary(sib3))["(Intercept)","Pr(>|t|)"])
results$coef_y_mean_sib3[i]<-tryc(coef(summary(sib3))["y_mean","Estimate"])
results$coef_y_mean_se_sib3[i]<-tryc(coef(summary(sib3))["y_mean","Std. Error"])
results$coef_y_mean_t_sib3[i]<-tryc(coef(summary(sib3))["y_mean","t value"])
results$coef_y_mean_p_sib3[i]<-tryc(coef(summary(sib3))["y_mean","Pr(>|t|)"])
results$coef_con_mean_sib3[i]<-tryc(coef(summary(sib3))["con_mean","Estimate"])
results$coef_con_mean_se_sib3[i]<-tryc(coef(summary(sib3))["con_mean","Std. Error"])
results$coef_con_mean_t_sib3[i]<-tryc(coef(summary(sib3))["con_mean","t value"])
results$coef_con_mean_p_sib3[i]<-tryc(coef(summary(sib3))["con_mean","Pr(>|t|)"])
results$coef_con_diff_sib3[i]<-tryc(coef(summary(sib3))["con_diff","Estimate"])
results$coef_con_diff_se_sib3[i]<-tryc(coef(summary(sib3))["con_diff","Std. Error"])
results$coef_con_diff_t_sib3[i]<-tryc(coef(summary(sib3))["con_diff","t value"])
results$coef_con_diff_p_sib3[i]<-tryc(coef(summary(sib3))["con_diff","Pr(>|t|)"])

results$coef_ev1_sib3[i]<-tryc(coef(summary(sib3))["ev1","Estimate"])
results$coef_ev1_se_sib3[i]<-tryc(coef(summary(sib3))["ev1","Std. Error"])
results$coef_ev1_t_sib3[i]<-tryc(coef(summary(sib3))["ev1","t value"])
results$coef_ev1_p_sib3[i]<-tryc(coef(summary(sib3))["ev1","Pr(>|t|)"])
results$coef_ev2_sib3[i]<-tryc(coef(summary(sib3))["ev2","Estimate"])
results$coef_ev2_se_sib3[i]<-tryc(coef(summary(sib3))["ev2","Std. Error"])
results$coef_ev2_t_sib3[i]<-tryc(coef(summary(sib3))["ev2","t value"])
results$coef_ev2_p_sib3[i]<-tryc(coef(summary(sib3))["ev2","Pr(>|t|)"])


results$coef_con_diffxev1_sib3[i]<-tryc(coef(summary(sib3))["con_diff:ev1","Estimate"])
results$coef_con_diffxev1_sib3se[i]<-tryc(coef(summary(sib3))["con_diff:ev1","Std. Error"])
results$coef_con_diffxev1_t_sib3[i]<-tryc(coef(summary(sib3))["con_diff:ev1","t value"])
results$coef_con_diffxev1_p_sib3[i]<-tryc(coef(summary(sib3))["con_diff:ev1","Pr(>|t|)"])
results$coef_con_diffxev2_sib3[i]<-tryc(coef(summary(sib3))["con_diff:ev2","Estimate"])
results$coef_con_diffxev2_se_sib3[i]<-tryc(coef(summary(sib3))["con_diff:ev2","Std. Error"])
results$coef_con_diffxev2_t_sib3[i]<-tryc(coef(summary(sib3))["con_diff:ev2","t value"])
results$coef_con_diffxev2_p_sib3[i]<-tryc(coef(summary(sib3))["con_diff:ev2","Pr(>|t|)"])



## MLM fixed effect

### MLM0

results$coef_fixed_Intercept_MLM0[i] <- tryc(summary(MLM0)[['coefficients']]['(Intercept)','Estimate'])
results$coef_fixed_Intercept_se_MLM0[i] <- tryc(summary(MLM0)[['coefficients']]['(Intercept)','Std. Error'])
results$coef_fixed_Intercept_df_MLM0[i] <- tryc(summary(MLM0)[['coefficients']]['(Intercept)','df'])
results$coef_fixed_Intercept_t_MLM0[i] <- tryc(summary(MLM0)[['coefficients']]['(Intercept)','t value'])
results$coef_fixed_Intercept_p_MLM0[i] <- tryc(summary(MLM0)[['coefficients']]['(Intercept)','Pr(>|t|)'])


## MLM 1


results$coef_fixed_Intercept_MLM1[i] <- tryc(summary(MLM1)[['coefficients']]['(Intercept)','Estimate'])
results$coef_fixed_Intercept_se_MLM1[i] <- tryc(summary(MLM1)[['coefficients']]['(Intercept)','Std. Error'])
results$coef_fixed_Intercept_df_MLM1[i] <- tryc(summary(MLM1)[['coefficients']]['(Intercept)','df'])
results$coef_fixed_Intercept_t_MLM1[i] <- tryc(summary(MLM1)[['coefficients']]['(Intercept)','t value'])
results$coef_fixed_Intercept_p_MLM1[i] <- tryc(summary(MLM1)[['coefficients']]['(Intercept)','Pr(>|t|)'])

results$coef_fixed_con_MLM1[i] <- tryc(summary(MLM1)[['coefficients']]['con','Estimate'])
results$coef_fixed_con_se_MLM1[i] <- tryc(summary(MLM1)[['coefficients']]['con','Std. Error'])
results$coef_fixed_con_df_MLM1[i] <- tryc(summary(MLM1)[['coefficients']]['con','df'])
results$coef_fixed_con_t_MLM1[i] <- tryc(summary(MLM1)[['coefficients']]['con','t value'])
results$coef_fixed_con_p_MLM1[i] <- tryc(summary(MLM1)[['coefficients']]['con','Pr(>|t|)'])


results$coef_fixed_catwtn_MLM1[i] <- tryc(summary(MLM1)[['coefficients']]['catwtn','Estimate'])
results$coef_fixed_catwtn_se_MLM1[i] <- tryc(summary(MLM1)[['coefficients']]['catwtn','Std. Error'])
results$coef_fixed_catwtn_df_MLM1[i] <- tryc(summary(MLM1)[['coefficients']]['catwtn','df'])
results$coef_fixed_catwtn_t_MLM1[i] <- tryc(summary(MLM1)[['coefficients']]['catwtn','t value'])
results$coef_fixed_catwtn_p_MLM1[i] <- tryc(summary(MLM1)[['coefficients']]['catwtn','Pr(>|t|)'])

### MLM2



results$coef_fixed_Intercept_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['(Intercept)','Estimate'])
results$coef_fixed_Intercept_se_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['(Intercept)','Std. Error'])
results$coef_fixed_Intercept_df_MLM2 <- tryc(summary(MLM2)[['coefficients']]['(Intercept)','df'])
results$coef_fixed_Intercept_t_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['(Intercept)','t value'])
results$coef_fixed_Intercept_p_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['(Intercept)','Pr(>|t|)'])

results$coef_fixed_con_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['con','Estimate'])
results$coef_fixed_con_se_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['con','Std. Error'])
results$coef_fixed_con_df_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['con','df'])
results$coef_fixed_con_t_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['con','t value'])
results$coef_fixed_con_p_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['con','Pr(>|t|)'])


results$coef_fixed_catwtn_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['catwtn','Estimate'])
results$coef_fixed_catwtn_se_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['catwtn','Std. Error'])
results$coef_fixed_catwtn_df_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['catwtn','df'])
results$coef_fixed_catwtn_t_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['catwtn','t value'])
results$coef_fixed_catwtn_p_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['catwtn','Pr(>|t|)'])



results$coef_fixed_catwtnxcon_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['catwtn:con','Estimate'])
results$coef_fixed_catwtnxcon_se_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['catwtn:con','Std. Error'])
results$coef_fixed_catwtnxcon_df_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['catwtn:con','df'])
results$coef_fixed_catwtnxcon_t_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['catwtn:con','t value'])
results$coef_fixed_catwtnxcon_p_MLM2[i] <- tryc(summary(MLM2)[['coefficients']]['catwtn:con','Pr(>|t|)'])



### MLM3


results$coef_fixed_Intercept_MLM3[i] <- tryc(summary(MLM3)[['coefficients']]['(Intercept)','Estimate'])
results$coef_fixed_Intercept_se_MLM3[i] <- tryc(summary(MLM3)[['coefficients']]['(Intercept)','Std. Error'])
results$coef_fixed_Intercept_df_MLM3[i] <- tryc(summary(MLM3)[['coefficients']]['(Intercept)','df'])
results$coef_fixed_Intercept_t_MLM3[i] <- tryc(summary(MLM3)[['coefficients']]['(Intercept)','t value'])
results$coef_fixed_Intercept_p_MLM3[i] <- tryc(summary(MLM3)[['coefficients']]['(Intercept)','Pr(>|t|)'])

results$coef_fixed_con_MLM3[i] <- tryc(summary(MLM3)[['coefficients']]['con','Estimate'])
results$coef_fixed_con_se_MLM3[i] <- tryc(summary(MLM3)[['coefficients']]['con','Std. Error'])
results$coef_fixed_con_df_MLM3[i] <- tryc(summary(MLM3)[['coefficients']]['con','df'])
results$coef_fixed_con_t_MLM3[i] <- tryc(summary(MLM3)[['coefficients']]['con','t value'])
results$coef_fixed_con_p_MLM3[i] <- tryc(summary(MLM3)[['coefficients']]['con','Pr(>|t|)'])


results$coef_fixed_catbtw_MLM3[i] <- tryc(summary(MLM3)[['coefficients']]['catbtw','Estimate'])
results$coef_fixed_catbtw_se_MLM3[i] <- tryc(summary(MLM3)[['coefficients']]['catbtw','Std. Error'])
results$coef_fixed_catbtw_df_MLM3[i] <- tryc(summary(MLM3)[['coefficients']]['catbtw','df'])
results$coef_fixed_catbtw_t_MLM3[i] <- tryc(summary(MLM3)[['coefficients']]['catbtw','t value'])
results$coef_fixed_catbtw_p_MLM3[i] <- tryc(summary(MLM3)[['coefficients']]['catbtw','Pr(>|t|)'])

### MLM4


results$coef_fixed_Intercept_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['(Intercept)','Estimate'])
results$coef_fixed_Intercept_se_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['(Intercept)','Std. Error'])
results$coef_fixed_Intercept_df_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['(Intercept)','df'])
results$coef_fixed_Intercept_t_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['(Intercept)','t value'])
results$coef_fixed_Intercept_p_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['(Intercept)','Pr(>|t|)'])

results$coef_fixed_con_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['con','Estimate'])
results$coef_fixed_con_se_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['con','Std. Error'])
results$coef_fixed_con_df_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['con','df'])
results$coef_fixed_con_t_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['con','t value'])
results$coef_fixed_con_p_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['con','Pr(>|t|)'])


results$coef_fixed_catbtw_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['catbtw','Estimate'])
results$coef_fixed_catbtw_se_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['catbtw','Std. Error'])
results$coef_fixed_catbtw_df_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['catbtw','df'])
results$coef_fixed_catbtw_t_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['catbtw','t value'])
results$coef_fixed_catbtw_p_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['catbtw','Pr(>|t|)'])


results$coef_fixed_conxcatbtw_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['con:catbtw','Estimate'])
results$coef_fixed_conxcatbtw_se_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['con:catbtw','Std. Error'])
results$coef_fixed_conxcatbtw_df_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['con:catbtw','df'])
results$coef_fixed_conxcatbtw_t_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['con:catbtw','t value'])
results$coef_fixed_conxcatbtw_p_MLM4[i] <- tryc(summary(MLM4)[['coefficients']]['con:catbtw','Pr(>|t|)'])

### MLM5

results$coef_fixed_Intercept_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['(Intercept)','Estimate'])
results$coef_fixed_Intercept_se_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['(Intercept)','Std. Error'])
results$coef_fixed_Intercept_df_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['(Intercept)','df'])
results$coef_fixed_Intercept_t_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['(Intercept)','t value'])
results$coef_fixed_Intercept_p_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['(Intercept)','Pr(>|t|)'])

results$coef_fixed_con_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['con','Estimate'])
results$coef_fixed_con_se_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['con','Std. Error'])
results$coef_fixed_con_df_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['con','df'])
results$coef_fixed_con_t_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['con','t value'])
results$coef_fixed_con_p_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['con','Pr(>|t|)'])


results$coef_fixed_ev1_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['ev1','Estimate'])
results$coef_fixed_ev1_se_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['ev1','Std. Error'])
results$coef_fixed_ev1_df_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['ev1','df'])
results$coef_fixed_ev1_t_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['ev1','t value'])
results$coef_fixed_ev1_p_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['ev1','Pr(>|t|)'])


results$coef_fixed_ev2_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['ev2','Estimate'])
results$coef_fixed_ev2_se_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['ev2','Std. Error'])
results$coef_fixed_ev2_df_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['ev2','df'])
results$coef_fixed_ev2_t_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['ev2','t value'])
results$coef_fixed_ev2_p_MLM5[i] <- tryc(summary(MLM5)[['coefficients']]['ev2','Pr(>|t|)'])

## MLM6 
results$coef_fixed_Intercept_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['(Intercept)','Estimate'])
results$coef_fixed_Intercept_se_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['(Intercept)','Std. Error'])
results$coef_fixed_Intercept_df_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['(Intercept)','df'])
results$coef_fixed_Intercept_t_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['(Intercept)','t value'])
results$coef_fixed_Intercept_p_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['(Intercept)','Pr(>|t|)'])

results$coef_fixed_con_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['con','Estimate'])
results$coef_fixed_con_se_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['con','Std. Error'])
results$coef_fixed_con_df_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['con','df'])
results$coef_fixed_con_t_MLM6[i] <- 
  tryc(summary(MLM6)[['coefficients']]['con','t value'])
results$coef_fixed_con_p_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['con','Pr(>|t|)'])


results$coef_fixed_ev1_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['ev1','Estimate'])
results$coef_fixed_ev1_se_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['ev1','Std. Error'])
results$coef_fixed_ev1_df_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['ev1','df'])
results$coef_fixed_ev1_t_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['ev1','t value'])
results$coef_fixed_ev1_p_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['ev1','Pr(>|t|)'])


results$coef_fixed_ev2_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['ev2','Estimate'])
results$coef_fixed_ev2_se_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['ev2','Std. Error'])
results$coef_fixed_ev2_df_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['ev2','df'])
results$coef_fixed_ev2_t_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['ev2','t value'])
results$coef_fixed_ev2_p_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['ev2','Pr(>|t|)'])


results$coef_fixed_conxev1_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['con:ev1','Estimate'])
results$coef_fixed_conxev1_se_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['con:ev1','Std. Error'])
results$coef_fixed_conxev1_df_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['con:ev1','df'])
results$coef_fixed_conxev1_t_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['con:ev1','t value'])
results$coef_fixed_conxev1_p_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['con:ev1','Pr(>|t|)'])


results$coef_fixed_conxev2_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['con:ev2','Estimate'])
results$coef_fixed_conxev2_se_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['con:ev2','Std. Error'])
results$coef_fixed_conxev2_df_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['con:ev2','df'])
results$coef_fixed_conxev2_t_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['con:ev2','t value'])
results$coef_fixed_conxev2_p_MLM6[i] <- tryc(summary(MLM6)[['coefficients']]['con:ev2','Pr(>|t|)'])

}

### save the result

write.csv(results,paste0("model_results/results_",conditions$condition[j],".csv"))

}

end_time<-Sys.time()

end_time-start_time
```

# post-generation process

### Convergence

``` r
## empty holder to keep data 

##  warning_model (indicates whether this model throws warning)
## hasConverged_model  (indiciates whetehr this model has converged or not
## these two information was cathed by by grep() functions. 

### following three lines are for making a empty holder with column names.

results_1 <- read.csv("model_results/results_1.csv")
results_1 <- results_1[,grep("hasConverged",names(results_1))]
convergence <- data.frame(results_1[0,])

# get convergence information. 

for(i in 1:ncond){
# read the data 
results <- read.csv(paste0("model_results/results_",i,".csv"))

#mean() You can use mean() to get the proportion of TRUE of a logical vector.


convergence[i,]<-results %>%
  summarise(across(contains("hasConverged"), ~mean(.x==1)))


convergence$condition[i] <- i 



}
#moved outside the loop so it goes faster
write.csv(convergence,"convergence/convergence.csv")

# binary variable that indicate the converge or not. 
# and when making summary statistics or report something, exclude
# non-convergence results. 
# in a loop. data collection process 
```

### warning message

as said earlier, keeping warning message might be worthy (just in case)

``` r
## empty holder to keep data 

### following three lines are for making a empty holder with column names. 
results_1 <- read.csv("model_results/results_1.csv")
results_1 <- results_1[,grep("warning",names(results_1))]
warning <- data.frame(results_1[0,])

# keep the data if the warning was generated. 

for(i in 1:ncond){
# read the data 
results <- read.csv(paste0("model_results/results_",i,".csv"))

#mean() You can use mean() to get the proportion of TRUE of a logical vector.
 warning[i,]<-results %>%
   summarise(across(contains("warning"), ~ mean(.x!='Nothing')))



warning$condition[i] <- i 



}
#moved outside the loop so it goes faster
write.csv(warning,"convergence/warning.csv")

# binary variable that indicate the converge or not. 
# and when making summary statistics or report something, exclude
# non-convergence results. 
# in a loop. data collection process 
```

## type 1 error rate

``` r
## empty holderto keep data 


# to make empty holder with col names 

results_1<-read.csv("model_results/results_1.csv")
results_1<-results_1[,grep("_p_",names(results_1))]
non_sig_rate<-data.frame(results_1[0,])

nMLM<-6 # how many MLM models we have? 0 to 6

## get information 

for(i in 1:ncond){
  
# read the data 
results <- read.csv(paste0("model_results/results_",i,".csv"))

# exclude non-convergence model  

for(k in 0:nMLM){  
  results[
    ## rows to replace
    results[[paste0("hasConverged_MLM", k)]]==0,
    ## column to replace
    grep(paste0("MLM",k),names(results))
  ] <- NA
}
# These two lines generates error 
# results <- results %>%
#   filter(has.)

# mean() You can use mean() to get the proportion of TRUE of a logical vector.

non_sig_rate[i,]<-  results %>%
  summarise(across(contains("_p_"), ~ mean(.x > 0.05, na.rm=TRUE)))
 # type1error$conditions[i] <- i 
non_sig_rate$condition[i] <- i

}

write.csv(non_sig_rate,"type1error/non_sig_rate.csv")
```

# get p-value of parameters

(excluding results from non-convergence model )

``` r
## empty holder to keep data 

# identical(dput(colnames(ho)),dput(colnames(type1error)))

# to make empty holder with col names 

results_1<-read.csv("model_results/results_1.csv")
results_1<-results_1[,grep("_p_",names(results_1))]
sig_rate<-data.frame(results_1[0,])

nMLM<-6 # how many MLM models we have? 0 to 6

for(i in 1:ncond){
  
# read the data 
results <- read.csv(paste0("model_results/results_",i,".csv"))

# exclude non-convergence model  

for(k in 0:nMLM){  
  results[
    ## rows to replace
    results[[paste0("hasConverged_MLM", k)]]==0,
    ## column to replace
    grep(paste0("MLM",k),names(results))
  ] <- NA
}
# These two lines generates error 
# results <- results %>%
#   filter(has.)

# mean() You can use mean() to get the proportion of TRUE of a logical vector.

sig_rate[i,]<-  results %>%
  summarise(across(contains("_p_"), ~ mean(.x <= 0.05, na.rm=TRUE)))
 # type1error$conditions[i] <- i 
sig_rate$condition[i] <- i

}

write.csv(sig_rate,"type1error/sig_rate.csv")
```

sig rate table with 0.01 alpha level

``` r
## empty holder

# identical(dput(colnames(ho)),dput(colnames(type1error)))

# to make empty holder with col names 

results_1<-read.csv("model_results/results_1.csv")
results_1<-results_1[,grep("_p_",names(results_1))]
sig_rate3<-data.frame(results_1[0,])

nMLM<-6 # how many MLM models we have? 0 to 6

for(i in 1:ncond){
  
# read the data 
results <- read.csv(paste0("model_results/results_",i,".csv"))

# exclude non-convergence model  

for(k in 0:nMLM){  
  results[
    ## rows to replace
    results[[paste0("hasConverged_MLM", k)]]==0,
    ## column to replace
    grep(paste0("MLM",k),names(results))
  ] <- NA
}
# These two lines generates error 
# results <- results %>%
#   filter(has.)

# mean() You can use mean() to get the proportion of TRUE of a logical vector.

sig_rate3[i,]<-  results %>%
  summarise(across(contains("_p_"), ~ mean(.x <= 0.01, na.rm=TRUE)))
 # type1error$conditions[i] <- i 
sig_rate3$condition[i] <- i

}

write.csv(sig_rate3,"type1error/sig_rate_0.01alpha.csv")
```

## descriptive statistics

### Quantiles

``` r
## empty holder to keep data
quantiles<-as.data.frame(matrix(nrow = ncond))

for(i in 1:ncond){
# read the data 
results <- read.csv(paste0("model_results/results_",i,".csv"))


# exclude non-convergence model  
## why 6? becuase we have 6 model in MLM
for(k in 0:nMLM){  
  results[
    ## rows to replace
    results[[paste0("hasConverged_MLM", k)]]==0,
    ## column to replace
    grep(paste0("MLM",k),names(results))
  ] <- NA
}

quantiles$conditions[i] <- i

start <- which( colnames(results)=="hasConverged_MLM6")+1

# becuase 1st~20th column is about convergence/error. 
for ( k in start:ncol(results)){

quantiles[i,paste0(names(results)[k],"_",c("01","05","25","50","75", "95","99"))]<- quantile(results[k],c(.01,.05,.25,.5,.75,.95,.99),na.rm=TRUE)

}

}
# fixed qualilites to not burnout
write.csv(quantiles,"descriptive_results/quantiles.csv")

# 
# 25%/75% included
```

### Mean and standard deviation

``` r
##  empty holder

mean_var<-as.data.frame(matrix(nrow = ncond))

for(i in 1:ncond){
# read the data 
results <- read.csv(paste0("model_results/results_",i,".csv"))


start <- which( colnames(results)=="hasConverged_MLM6")+1

# exclude non-convergence model  
## why 6? becuase we have 6 model in MLM
for(k in 0:nMLM){  
  results[
    ## rows to replace
    results[[paste0("hasConverged_MLM", k)]]==0,
    ## column to replace
    grep(paste0("MLM",k),names(results))
  ] <- NA
}

mean_var$condition[i] <- i

## what is K and why does it start with 20
# k is for columns in results dataset. 
# it starts with 20 because I think 1st~20th columns does not need mean value or 
# variance value because it is just an collection of warning message, or 
# indication of singular fit. 

# 20 self-updating way



for ( k in start:ncol(results)){

mean_var[i,paste0(names(results)[k],"mean")]<- mean(results[,k],na.rm=TRUE)

mean_var[i,paste0(names(results)[k],"var")]<- var(results[,k],na.rm=TRUE)

}


}
# outside the loop otherwise you're just overwriting everything and making the results much slower
write.csv(mean_var,"descriptive_results/mean_var.csv")
```

the variance table to check

``` r
mean_var <- read.csv("descriptive_results/mean_var.csv")
var<- mean_var[,grep("var$",names(mean_var))]
write.csv(var,"descriptive_results/var.csv")
```

## some useful, selected tables.

includes mean of parameters, and other information such as AIC and
criterion.

``` r
mean_var <- read.csv("descriptive_results/mean_var.csv")
mean<- mean_var[,grep("mean$|condition",names(mean_var))]
mean_all<-mean[,grep("AIC|BIC|deviance|criterion|rsquared|ICC_MLM0|condition",names(mean))]
write.csv(mean_all,"descriptive_results/mean_all.csv")

 ## will split into MLM and sib becuase the table is too big. 

#mean_all<-read.csv("descriptive_results/mean_all.csv")
```

#### Pseudo ICC

``` r
mean_var <- read.csv("descriptive_results/mean_var.csv")
mean<- mean_var[,grep("mean$|condition",names(mean_var))]
psuedo_ICC<-mean[,grep("ICC|condition",names(mean))]
write.csv(psuedo_ICC,"appendix/psuedo_ICC.csv")
```

#### singularity rate

``` r
# make an empty tray

results_1<-read.csv("model_results/results_1.csv")
 results_1<-results_1[,grep("singularity",names(results_1))]
 singularity<-data.frame(results_1[0,])

nMLM<-6 # how many MLM models we have? 0 to 6 

for(i in 1:ncond){
# read the data 
results <- read.csv(paste0("model_results/results_",i,".csv"))

# exclude non-convergence model  
## why 6? becuase we have 6 model in MLM
for(k in 0:nMLM){  
  results[
    ## rows to replace
    results[[paste0("hasConverged_MLM", k)]]==0,
    ## column to replace
    grep(paste0("MLM",k),names(results))
  ] <- NA
}



# mean() You can use mean() to get the proportion of TRUE of a logical vector.
 # type1error[i,]<-results %>%
 #  summarise(across(contains("_p_"), ~ mean(.x > 0.05)))
# 
singularity[i,]<-  results %>%
  summarise(across(contains("singularity"), ~ mean(.x==TRUE, na.rm=TRUE)))
 singularity$condition[i]<-i
}
# one csv file in type 1 error rate 
# fixed
write.csv(singularity,"appendix/singularity.csv")
```

siguarity rate that including non-convergence models.

``` r
# make an empty tray to keep data 
# below three lines are to make a empty tray with col names. 

results_1<-read.csv("model_results/results_1.csv")
 results_1<-results_1[,grep("singularity",names(results_1))]
 singularity2<-data.frame(results_1[0,])

nMLM<-6 # how many MLM models we have? 0 to 6 

for(i in 1:ncond){
# read the data 
results <- read.csv(paste0("model_results/results_",i,".csv"))




# mean() You can use mean() to get the proportion of TRUE of a logical vector.
 # type1error[i,]<-results %>%
 #  summarise(across(contains("_p_"), ~ mean(.x > 0.05)))
# 
singularity2[i,]<-  results %>%
  summarise(across(contains("singularity"), ~ mean(.x==TRUE, na.rm=TRUE)))
 # type1error$conditions[i] <- i 

}
# one csv file in type 1 error rate 
# fixed
write.csv(singularity2,"appendix/singularity_including_nonconvergence_model.csv")
```

#### gender composition varaible’s descriptive statistics

``` r
#descriptive gender compositiion 

mean_var <- read.csv("descriptive_results/mean_var.csv")
mean<- mean_var[,grep("mean$|condition",names(mean_var))]
gender_composition<-mean[,grep("gender",names(mean))]
write.csv(gender_composition,"appendix/gender_composition.csv")
```

##### Additional tables

##### RMES, ResSE, REsSD.

Not relavent to our research question.

``` r
# Make sure to note that this RMSE is not standardized(scaled)
# Now I think we do not need RMSE

 
mean_var <- read.csv("descriptive_results/mean_var.csv")
mean<- mean_var[,grep("mean$|condition",names(mean_var))]
SD<-mean[,grep("RMSE|ResSE|ResSD|condition",names(mean))]
write.csv(SD,"appendix/SD.csv")

#make a nicer table 

SD_sib<-mean[,grep("sib|condition",names(SD))]
SD_MLM<-mean[,grep("MLM|condition",names(SD))]

write.csv(SD,"appendix/SD_sib.csv")
write.csv(SD,"appendix/SD_MLM.csv")
```

check the reference group is equal across all the model -chekced. check
the digit.
