Singularity_issue
================
Yoo Ri Hwang
2022-07-01

# overview

The results of this research project showed that some of multilevel
models have singularity issues. In this document, I tried to look
closely how singularity rate are related to other conditions, such as CN
(sample size), beta2 (gender effects), beta3 (interaction effect), and
ICC.

## pacakge

``` r
library(tidyverse)
```

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.1 --

    ## v ggplot2 3.3.6     v purrr   0.3.4
    ## v tibble  3.1.4     v dplyr   1.0.8
    ## v tidyr   1.2.0     v stringr 1.4.0
    ## v readr   2.1.2     v forcats 0.5.1

    ## Warning: package 'ggplot2' was built under R version 4.1.3

    ## Warning: package 'tidyr' was built under R version 4.1.3

    ## Warning: package 'readr' was built under R version 4.1.3

    ## Warning: package 'dplyr' was built under R version 4.1.3

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(apaTables)
```

## Make informative table

### condition information

``` r
conditions<-read.csv("conditions/Conditions_summary.csv")
conditions<-conditions %>% 
  rename(
    CN = cn,
    ICC = target_ICC
    )

# this is for merging 
conditions_con<-conditions %>%
  dplyr::select("CN","beta1","beta2","beta3","ICC","condition")
```

load data

``` r
singularity<-read.csv("appendix/singularity.csv")
psuedo_ICC<-read.csv("appendix/psuedo_ICC.csv")
sing_psuedo<-read.csv("appendix/singularity_with_psuedoICC.csv")
mlm56<-sing_psuedo[,grep("CN|beta|condition|ICC|MLM5|MLM6",names(sing_psuedo))]
# 
# singularity<-merge(singularity,conditions_con, key="condition")
# pseudo_ICC<-merge(psuedo_ICC,conditions_con, key="condition")
# 

# write.csv(singularity,"appendix/singularity.csv")
# write.csv(psuedo_ICC,"appendix/psuedo_ICC.csv")
# 
# singularity_with_psuedoICC<-merge(singularity, pseudo_ICC, merge="condition")
# write.csv(sing_psuedo,"appendix/singularity_with_psuedoICC.csv")
```

## Regression

### Regression as a whole

put singularity rate as a DV, and other conditions as IVs.

``` r
# data reshaping 

sing_psuedo$X<-NULL
singularity$X<-NULL

singularity_long<-singularity %>%
  pivot_longer(
    cols=c('singularity_MLM0','singularity_MLM1','singularity_MLM2','singularity_MLM3','singularity_MLM4','singularity_MLM5','singularity_MLM6'),
    names_to = "model_type",
    values_to = "singularity_rate"
  )

## regression

### reminder: beta1 only have 0.3 in every condition 

sin<-lm(singularity_rate ~ factor(model_type)+ICC+beta2+beta3+CN,data=singularity_long)
summary(sin)
```

    ## 
    ## Call:
    ## lm(formula = singularity_rate ~ factor(model_type) + ICC + beta2 + 
    ##     beta3 + CN, data = singularity_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.51598 -0.11927 -0.00875  0.14260  0.44302 
    ## 
    ## Coefficients:
    ##                                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        -1.803e-01  2.708e-02  -6.658 5.11e-11 ***
    ## factor(model_type)singularity_MLM1 -3.043e-03  2.631e-02  -0.116 0.907971    
    ## factor(model_type)singularity_MLM2 -1.108e-02  2.631e-02  -0.421 0.673902    
    ## factor(model_type)singularity_MLM3  1.941e-01  2.631e-02   7.378 4.00e-13 ***
    ## factor(model_type)singularity_MLM4  1.956e-01  2.631e-02   7.434 2.68e-13 ***
    ## factor(model_type)singularity_MLM5  6.469e-01  2.631e-02  24.583  < 2e-16 ***
    ## factor(model_type)singularity_MLM6  5.911e-01  2.631e-02  22.463  < 2e-16 ***
    ## ICC                                -1.155e-02  2.819e-02  -0.410 0.682128    
    ## beta2                               7.291e-01  4.034e-02  18.075  < 2e-16 ***
    ## beta3                               1.451e-01  3.661e-02   3.964 8.04e-05 ***
    ## CN                                 -1.241e-04  3.376e-05  -3.677 0.000252 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2013 on 808 degrees of freedom
    ## Multiple R-squared:  0.6837, Adjusted R-squared:  0.6798 
    ## F-statistic: 174.7 on 10 and 808 DF,  p-value: < 2.2e-16

``` r
sin_factor<-lm(singularity_rate ~ factor(model_type)+factor(ICC)+factor(beta2)+factor(beta3)+factor(CN),data=singularity_long)
summary(sin_factor)
```

    ## 
    ## Call:
    ## lm(formula = singularity_rate ~ factor(model_type) + factor(ICC) + 
    ##     factor(beta2) + factor(beta3) + factor(CN), data = singularity_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.47831 -0.12429 -0.01349  0.16939  0.42185 
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        -0.185085   0.032933  -5.620 2.63e-08 ***
    ## factor(model_type)singularity_MLM1 -0.003043   0.025694  -0.118  0.90576    
    ## factor(model_type)singularity_MLM2 -0.011077   0.025694  -0.431  0.66651    
    ## factor(model_type)singularity_MLM3  0.194133   0.025694   7.556 1.14e-13 ***
    ## factor(model_type)singularity_MLM4  0.195624   0.025694   7.614 7.50e-14 ***
    ## factor(model_type)singularity_MLM5  0.646889   0.025694  25.177  < 2e-16 ***
    ## factor(model_type)singularity_MLM6  0.591085   0.025694  23.005  < 2e-16 ***
    ## factor(ICC)0.4                     -0.011791   0.016821  -0.701  0.48352    
    ## factor(ICC)0.8                     -0.008827   0.016821  -0.525  0.59990    
    ## factor(beta2)0.1                    0.061942   0.030324   2.043  0.04141 *  
    ## factor(beta2)0.3                    0.303533   0.030324  10.010  < 2e-16 ***
    ## factor(beta2)0.5                    0.348387   0.030324  11.489  < 2e-16 ***
    ## factor(beta3)0.1                    0.006190   0.020216   0.306  0.75955    
    ## factor(beta3)0.3                    0.042294   0.020216   2.092  0.03674 *  
    ## factor(beta3)0.5                    0.065337   0.020216   3.232  0.00128 ** 
    ## factor(CN)120                      -0.036433   0.016821  -2.166  0.03061 *  
    ## factor(CN)510                      -0.068287   0.016821  -4.060 5.39e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1965 on 802 degrees of freedom
    ## Multiple R-squared:  0.7007, Adjusted R-squared:  0.6947 
    ## F-statistic: 117.3 on 16 and 802 DF,  p-value: < 2.2e-16

``` r
#apa.reg.table(sin_factor, filename="tables/singularity_regression.doc")

## remove ICC


sin<-lm(singularity_rate ~ factor(model_type)+beta2+beta3+CN,data=singularity_long)
summary(sin)
```

    ## 
    ## Call:
    ## lm(formula = singularity_rate ~ factor(model_type) + beta2 + 
    ##     beta3 + CN, data = singularity_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.51521 -0.11970 -0.00763  0.14406  0.43917 
    ## 
    ## Coefficients:
    ##                                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        -1.857e-01  2.366e-02  -7.850 1.32e-14 ***
    ## factor(model_type)singularity_MLM1 -3.043e-03  2.630e-02  -0.116  0.90792    
    ## factor(model_type)singularity_MLM2 -1.108e-02  2.630e-02  -0.421  0.67374    
    ## factor(model_type)singularity_MLM3  1.941e-01  2.630e-02   7.381 3.89e-13 ***
    ## factor(model_type)singularity_MLM4  1.956e-01  2.630e-02   7.438 2.61e-13 ***
    ## factor(model_type)singularity_MLM5  6.469e-01  2.630e-02  24.596  < 2e-16 ***
    ## factor(model_type)singularity_MLM6  5.911e-01  2.630e-02  22.474  < 2e-16 ***
    ## beta2                               7.291e-01  4.032e-02  18.084  < 2e-16 ***
    ## beta3                               1.451e-01  3.660e-02   3.966 7.97e-05 ***
    ## CN                                 -1.241e-04  3.374e-05  -3.679  0.00025 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2012 on 809 degrees of freedom
    ## Multiple R-squared:  0.6837, Adjusted R-squared:  0.6801 
    ## F-statistic: 194.3 on 9 and 809 DF,  p-value: < 2.2e-16

``` r
sin_factor<-lm(singularity_rate ~ factor(model_type)+factor(beta2)+factor(beta3)+factor(CN),data=singularity_long)
summary(sin_factor)
```

    ## 
    ## Call:
    ## lm(formula = singularity_rate ~ factor(model_type) + factor(beta2) + 
    ##     factor(beta3) + factor(CN), data = singularity_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.47644 -0.12314 -0.01191  0.16778  0.42872 
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                        -0.191958   0.031440  -6.106 1.59e-09 ***
    ## factor(model_type)singularity_MLM1 -0.003043   0.025671  -0.119  0.90568    
    ## factor(model_type)singularity_MLM2 -0.011077   0.025671  -0.432  0.66622    
    ## factor(model_type)singularity_MLM3  0.194133   0.025671   7.562 1.08e-13 ***
    ## factor(model_type)singularity_MLM4  0.195624   0.025671   7.621 7.11e-14 ***
    ## factor(model_type)singularity_MLM5  0.646889   0.025671  25.200  < 2e-16 ***
    ## factor(model_type)singularity_MLM6  0.591085   0.025671  23.026  < 2e-16 ***
    ## factor(beta2)0.1                    0.061942   0.030296   2.045  0.04123 *  
    ## factor(beta2)0.3                    0.303533   0.030296  10.019  < 2e-16 ***
    ## factor(beta2)0.5                    0.348387   0.030296  11.499  < 2e-16 ***
    ## factor(beta3)0.1                    0.006190   0.020198   0.306  0.75933    
    ## factor(beta3)0.3                    0.042294   0.020198   2.094  0.03657 *  
    ## factor(beta3)0.5                    0.065337   0.020198   3.235  0.00127 ** 
    ## factor(CN)120                      -0.036433   0.016805  -2.168  0.03045 *  
    ## factor(CN)510                      -0.068287   0.016805  -4.063 5.31e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1963 on 804 degrees of freedom
    ## Multiple R-squared:  0.7005, Adjusted R-squared:  0.6953 
    ## F-statistic: 134.3 on 14 and 804 DF,  p-value: < 2.2e-16

### regression with predictors such as (with/out interaction, the type of gender composition)

1)  is singularity rate related to whether the model contains
    interaction term(s) or not?
2)  is singularity rate related to how gender variable was coded?

``` r
singularity_long <- singularity_long %>%
  mutate(
    interaction=case_when(
      model_type=="singularity_MLM0" ~ "NO",
      model_type=="singularity_MLM1" ~ "NO",
      model_type=="singularity_MLM2" ~ "YES",
      model_type=="singularity_MLM3" ~ "NO",
      model_type=="singularity_MLM4" ~ "YES",
      model_type=="singularity_MLM5" ~ "NO",
      model_type=="singularity_MLM6" ~ "YES"
    ),
    gender=case_when(
      model_type=="singularity_MLM0" ~ "no_gender_variable",
      model_type=="singularity_MLM1" ~ "individual", # wheter individual is male/female
      model_type=="singularity_MLM2" ~ "individual",
      model_type=="singularity_MLM3" ~ "same_mixed",
      model_type=="singularity_MLM4" ~ "same_mixed",
      model_type=="singularity_MLM5" ~ "mm_ff_mixed",
      model_type=="singularity_MLM6" ~ "mm_ff_mixed"
    )
  )

# set reference group 
singularity_long$gender <- relevel(factor(singularity_long$gender), ref="no_gender_variable")

#regression (numeric)
sin<-lm(singularity_rate ~ factor(interaction)+factor(gender)+ICC+beta2+beta3+CN,data=singularity_long)
summary(sin)
```

    ## 
    ## Call:
    ## lm(formula = singularity_rate ~ factor(interaction) + factor(gender) + 
    ##     ICC + beta2 + beta3 + CN, data = singularity_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.51256 -0.11940 -0.00943  0.14920  0.42551 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)               -1.803e-01  2.709e-02  -6.655 5.20e-11 ***
    ## factor(interaction)YES    -2.078e-02  1.520e-02  -1.367 0.171900    
    ## factor(gender)individual   3.331e-03  2.403e-02   0.139 0.889785    
    ## factor(gender)mm_ff_mixed  6.294e-01  2.403e-02  26.189  < 2e-16 ***
    ## factor(gender)same_mixed   2.053e-01  2.403e-02   8.542  < 2e-16 ***
    ## ICC                       -1.155e-02  2.821e-02  -0.410 0.682262    
    ## beta2                      7.291e-01  4.036e-02  18.067  < 2e-16 ***
    ## beta3                      1.451e-01  3.663e-02   3.962 8.09e-05 ***
    ## CN                        -1.241e-04  3.377e-05  -3.675 0.000253 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2014 on 810 degrees of freedom
    ## Multiple R-squared:  0.6827, Adjusted R-squared:  0.6795 
    ## F-statistic: 217.8 on 8 and 810 DF,  p-value: < 2.2e-16

``` r
#regression (Factor)

sin_factor<-lm(singularity_rate ~ factor(interaction)+factor(gender)+factor(ICC)+factor(beta2)+factor(beta3)+factor(CN),data=singularity_long)
summary(sin_factor)
```

    ## 
    ## Call:
    ## lm(formula = singularity_rate ~ factor(interaction) + factor(gender) + 
    ##     factor(ICC) + factor(beta2) + factor(beta3) + factor(CN), 
    ##     data = singularity_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.48250 -0.12319 -0.01312  0.17433  0.41548 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)               -0.185085   0.032951  -5.617 2.68e-08 ***
    ## factor(interaction)YES    -0.020782   0.014842  -1.400  0.16184    
    ## factor(gender)individual   0.003331   0.023468   0.142  0.88715    
    ## factor(gender)mm_ff_mixed  0.629378   0.023468  26.819  < 2e-16 ***
    ## factor(gender)same_mixed   0.205270   0.023468   8.747  < 2e-16 ***
    ## factor(ICC)0.4            -0.011791   0.016830  -0.701  0.48375    
    ## factor(ICC)0.8            -0.008827   0.016830  -0.524  0.60009    
    ## factor(beta2)0.1           0.061942   0.030340   2.042  0.04152 *  
    ## factor(beta2)0.3           0.303533   0.030340  10.004  < 2e-16 ***
    ## factor(beta2)0.5           0.348387   0.030340  11.483  < 2e-16 ***
    ## factor(beta3)0.1           0.006190   0.020227   0.306  0.75967    
    ## factor(beta3)0.3           0.042294   0.020227   2.091  0.03684 *  
    ## factor(beta3)0.5           0.065337   0.020227   3.230  0.00129 ** 
    ## factor(CN)120             -0.036433   0.016830  -2.165  0.03069 *  
    ## factor(CN)510             -0.068287   0.016830  -4.058 5.44e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1966 on 804 degrees of freedom
    ## Multiple R-squared:  0.6996, Adjusted R-squared:  0.6944 
    ## F-statistic: 133.8 on 14 and 804 DF,  p-value: < 2.2e-16

``` r
#apa.reg.table(sin_factor, filename="tables/singularity_regression.doc")
```

remove ICC

``` r
#regression (numeric)
sin<-lm(singularity_rate ~ factor(interaction)+factor(gender)+ICC+beta2+beta3+CN,data=singularity_long)
summary(sin)
```

    ## 
    ## Call:
    ## lm(formula = singularity_rate ~ factor(interaction) + factor(gender) + 
    ##     ICC + beta2 + beta3 + CN, data = singularity_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.51256 -0.11940 -0.00943  0.14920  0.42551 
    ## 
    ## Coefficients:
    ##                             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)               -1.803e-01  2.709e-02  -6.655 5.20e-11 ***
    ## factor(interaction)YES    -2.078e-02  1.520e-02  -1.367 0.171900    
    ## factor(gender)individual   3.331e-03  2.403e-02   0.139 0.889785    
    ## factor(gender)mm_ff_mixed  6.294e-01  2.403e-02  26.189  < 2e-16 ***
    ## factor(gender)same_mixed   2.053e-01  2.403e-02   8.542  < 2e-16 ***
    ## ICC                       -1.155e-02  2.821e-02  -0.410 0.682262    
    ## beta2                      7.291e-01  4.036e-02  18.067  < 2e-16 ***
    ## beta3                      1.451e-01  3.663e-02   3.962 8.09e-05 ***
    ## CN                        -1.241e-04  3.377e-05  -3.675 0.000253 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2014 on 810 degrees of freedom
    ## Multiple R-squared:  0.6827, Adjusted R-squared:  0.6795 
    ## F-statistic: 217.8 on 8 and 810 DF,  p-value: < 2.2e-16

``` r
#regression (Factor)

sin_factor<-lm(singularity_rate ~ factor(interaction)+factor(gender)+factor(ICC)+factor(beta2)+factor(beta3)+factor(CN),data=singularity_long)
summary(sin_factor)
```

    ## 
    ## Call:
    ## lm(formula = singularity_rate ~ factor(interaction) + factor(gender) + 
    ##     factor(ICC) + factor(beta2) + factor(beta3) + factor(CN), 
    ##     data = singularity_long)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.48250 -0.12319 -0.01312  0.17433  0.41548 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)               -0.185085   0.032951  -5.617 2.68e-08 ***
    ## factor(interaction)YES    -0.020782   0.014842  -1.400  0.16184    
    ## factor(gender)individual   0.003331   0.023468   0.142  0.88715    
    ## factor(gender)mm_ff_mixed  0.629378   0.023468  26.819  < 2e-16 ***
    ## factor(gender)same_mixed   0.205270   0.023468   8.747  < 2e-16 ***
    ## factor(ICC)0.4            -0.011791   0.016830  -0.701  0.48375    
    ## factor(ICC)0.8            -0.008827   0.016830  -0.524  0.60009    
    ## factor(beta2)0.1           0.061942   0.030340   2.042  0.04152 *  
    ## factor(beta2)0.3           0.303533   0.030340  10.004  < 2e-16 ***
    ## factor(beta2)0.5           0.348387   0.030340  11.483  < 2e-16 ***
    ## factor(beta3)0.1           0.006190   0.020227   0.306  0.75967    
    ## factor(beta3)0.3           0.042294   0.020227   2.091  0.03684 *  
    ## factor(beta3)0.5           0.065337   0.020227   3.230  0.00129 ** 
    ## factor(CN)120             -0.036433   0.016830  -2.165  0.03069 *  
    ## factor(CN)510             -0.068287   0.016830  -4.058 5.44e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1966 on 804 degrees of freedom
    ## Multiple R-squared:  0.6996, Adjusted R-squared:  0.6944 
    ## F-statistic: 133.8 on 14 and 804 DF,  p-value: < 2.2e-16

### regression per model (closer look)MLM3 MLM4, MLM5, MM6

Because MLM3, MLM4, MLM5, MLM6 showed high singularity rate, we closely
look each model’s signualrity rate to get a clear picture.

#### regression on singularity rate of MLM3

ICC were not sig, and removed from the regression equation

``` r
# 
# 
# m3<-lm(singularity_MLM4~CN+beta2+beta3+ICC,data=singularity)
# m3_factor<-lm(singularity_MLM4~factor(CN)+factor(beta2)+factor(beta3)+factor(ICC),data=singularity)
# 
# summary(m3)
# summary(m3_factor)



m3<-lm(singularity_MLM4~CN+beta2+beta3,data=singularity)
m3_factor<-lm(singularity_MLM4~factor(CN)+factor(beta2)+factor(beta3),data=singularity)

summary(m3)
```

    ## 
    ## Call:
    ## lm(formula = singularity_MLM4 ~ CN + beta2 + beta3, data = singularity)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.11511 -0.04486 -0.01087  0.04858  0.11109 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  6.426e-02  1.208e-02   5.319 5.34e-07 ***
    ## CN          -3.161e-04  2.513e-05 -12.578  < 2e-16 ***
    ## beta2        4.838e-01  3.003e-02  16.109  < 2e-16 ***
    ## beta3        4.140e-01  2.726e-02  15.189  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.05663 on 113 degrees of freedom
    ## Multiple R-squared:  0.8657, Adjusted R-squared:  0.8621 
    ## F-statistic: 242.8 on 3 and 113 DF,  p-value: < 2.2e-16

``` r
summary(m3_factor)
```

    ## 
    ## Call:
    ## lm(formula = singularity_MLM4 ~ factor(CN) + factor(beta2) + 
    ##     factor(beta3), data = singularity)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.11605 -0.04091  0.00639  0.03654  0.09166 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       0.08160    0.01902   4.291 3.89e-05 ***
    ## factor(CN)120    -0.07377    0.01203  -6.134 1.44e-08 ***
    ## factor(CN)510    -0.16736    0.01203 -13.916  < 2e-16 ***
    ## factor(beta2)0.1  0.05392    0.02168   2.487   0.0144 *  
    ## factor(beta2)0.3  0.15378    0.02168   7.093 1.44e-10 ***
    ## factor(beta2)0.5  0.24986    0.02168  11.524  < 2e-16 ***
    ## factor(beta3)0.1  0.01189    0.01445   0.823   0.4126    
    ## factor(beta3)0.3  0.10141    0.01445   7.016 2.10e-10 ***
    ## factor(beta3)0.5  0.20159    0.01445  13.947  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.05311 on 108 degrees of freedom
    ## Multiple R-squared:  0.8871, Adjusted R-squared:  0.8787 
    ## F-statistic: 106.1 on 8 and 108 DF,  p-value: < 2.2e-16

#### regression on singularity rate of MLM4

ICC were not sig, and removed from the regression equation

``` r
# m4<-lm(singularity_MLM4~CN+beta2+beta3+ICC,data=singularity)
# m4_factor<-lm(singularity_MLM4~factor(CN)+factor(beta2)+factor(beta3)+factor(ICC),data=singularity)
# 
# summary(m4)
# summary(m4_factor)

m4<-lm(singularity_MLM4~CN+beta2+beta3,data=singularity)
m4_factor<-lm(singularity_MLM4~factor(CN)+factor(beta2)+factor(beta3),data=singularity)

summary(m4)
```

    ## 
    ## Call:
    ## lm(formula = singularity_MLM4 ~ CN + beta2 + beta3, data = singularity)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.11511 -0.04486 -0.01087  0.04858  0.11109 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  6.426e-02  1.208e-02   5.319 5.34e-07 ***
    ## CN          -3.161e-04  2.513e-05 -12.578  < 2e-16 ***
    ## beta2        4.838e-01  3.003e-02  16.109  < 2e-16 ***
    ## beta3        4.140e-01  2.726e-02  15.189  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.05663 on 113 degrees of freedom
    ## Multiple R-squared:  0.8657, Adjusted R-squared:  0.8621 
    ## F-statistic: 242.8 on 3 and 113 DF,  p-value: < 2.2e-16

``` r
summary(m4_factor)
```

    ## 
    ## Call:
    ## lm(formula = singularity_MLM4 ~ factor(CN) + factor(beta2) + 
    ##     factor(beta3), data = singularity)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.11605 -0.04091  0.00639  0.03654  0.09166 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       0.08160    0.01902   4.291 3.89e-05 ***
    ## factor(CN)120    -0.07377    0.01203  -6.134 1.44e-08 ***
    ## factor(CN)510    -0.16736    0.01203 -13.916  < 2e-16 ***
    ## factor(beta2)0.1  0.05392    0.02168   2.487   0.0144 *  
    ## factor(beta2)0.3  0.15378    0.02168   7.093 1.44e-10 ***
    ## factor(beta2)0.5  0.24986    0.02168  11.524  < 2e-16 ***
    ## factor(beta3)0.1  0.01189    0.01445   0.823   0.4126    
    ## factor(beta3)0.3  0.10141    0.01445   7.016 2.10e-10 ***
    ## factor(beta3)0.5  0.20159    0.01445  13.947  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.05311 on 108 degrees of freedom
    ## Multiple R-squared:  0.8871, Adjusted R-squared:  0.8787 
    ## F-statistic: 106.1 on 8 and 108 DF,  p-value: < 2.2e-16

#### regression on singularity rate of MLM5

ICC were non sig, removed from reg equation

``` r
m5<-lm(singularity_MLM5~CN+beta2+beta3,data=singularity)
m5_factor<-lm(singularity_MLM5~factor(CN)+factor(beta2)+factor(beta3),data=singularity)

summary(m5)
```

    ## 
    ## Call:
    ## lm(formula = singularity_MLM5 ~ CN + beta2 + beta3, data = singularity)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.27564 -0.13238 -0.07274  0.15238  0.33739 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  7.809e-02  4.105e-02   1.902   0.0597 .  
    ## CN          -4.459e-05  8.539e-05  -0.522   0.6026    
    ## beta2        2.024e+00  1.020e-01  19.838   <2e-16 ***
    ## beta3        1.787e-01  9.262e-02   1.929   0.0562 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1924 on 113 degrees of freedom
    ## Multiple R-squared:  0.7866, Adjusted R-squared:  0.7809 
    ## F-statistic: 138.8 on 3 and 113 DF,  p-value: < 2.2e-16

``` r
summary(m5_factor)
```

    ## 
    ## Call:
    ## lm(formula = singularity_MLM5 ~ factor(CN) + factor(beta2) + 
    ##     factor(beta3), data = singularity)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.27796 -0.03844  0.01630  0.05296  0.21123 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       0.01265    0.03589   0.352   0.7252    
    ## factor(CN)120    -0.01110    0.02270  -0.489   0.6257    
    ## factor(CN)510    -0.02385    0.02270  -1.051   0.2958    
    ## factor(beta2)0.1  0.17379    0.04092   4.247  4.6e-05 ***
    ## factor(beta2)0.3  0.92098    0.04092  22.508  < 2e-16 ***
    ## factor(beta2)0.5  0.95823    0.04092  23.419  < 2e-16 ***
    ## factor(beta3)0.1  0.01281    0.02728   0.470   0.6395    
    ## factor(beta3)0.3  0.06593    0.02728   2.417   0.0173 *  
    ## factor(beta3)0.5  0.06533    0.02728   2.395   0.0183 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1002 on 108 degrees of freedom
    ## Multiple R-squared:  0.9447, Adjusted R-squared:  0.9406 
    ## F-statistic: 230.4 on 8 and 108 DF,  p-value: < 2.2e-16

#### regression on singularity rate of MLM6

ICC were non sig, removed from reg equation

``` r
### regression on singularity rate of MLM6
m6<-lm(singularity_MLM6~CN+beta2+beta3,data=singularity)
m6_factor<-lm(singularity_MLM6~factor(CN)+factor(beta2)+factor(beta3),data=singularity)

summary(m6)
```

    ## 
    ## Call:
    ## lm(formula = singularity_MLM6 ~ CN + beta2 + beta3, data = singularity)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.24640 -0.14217 -0.06807  0.15861  0.37002 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  2.824e-02  4.145e-02   0.681   0.4970    
    ## CN          -3.455e-05  8.621e-05  -0.401   0.6894    
    ## beta2        2.263e+00  1.030e-01  21.968   <2e-16 ***
    ## beta3       -1.791e-01  9.351e-02  -1.916   0.0579 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1943 on 113 degrees of freedom
    ## Multiple R-squared:  0.8108, Adjusted R-squared:  0.8057 
    ## F-statistic: 161.4 on 3 and 113 DF,  p-value: < 2.2e-16

``` r
summary(m6_factor)
```

    ## 
    ## Call:
    ## lm(formula = singularity_MLM6 ~ factor(CN) + factor(beta2) + 
    ##     factor(beta3), data = singularity)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.3402 -0.0329 -0.0036  0.0739  0.3435 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)       0.015385   0.037064   0.415 0.678907    
    ## factor(CN)120    -0.018821   0.023442  -0.803 0.423813    
    ## factor(CN)510    -0.022000   0.023442  -0.939 0.350079    
    ## factor(beta2)0.1  0.117093   0.042260   2.771 0.006586 ** 
    ## factor(beta2)0.3  0.922759   0.042260  21.835  < 2e-16 ***
    ## factor(beta2)0.5  1.015259   0.042260  24.024  < 2e-16 ***
    ## factor(beta3)0.1  0.002259   0.028173   0.080 0.936233    
    ## factor(beta3)0.3 -0.008222   0.028173  -0.292 0.770965    
    ## factor(beta3)0.5 -0.097963   0.028173  -3.477 0.000731 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1035 on 108 degrees of freedom
    ## Multiple R-squared:  0.9487, Adjusted R-squared:  0.9448 
    ## F-statistic: 249.4 on 8 and 108 DF,  p-value: < 2.2e-16
