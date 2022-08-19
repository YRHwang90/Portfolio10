table
================
Yoo Ri Hwang
5/19/2022

# Overview

This document is to generate tables so that the author can
see/understand the results more easily. Tables will be saved as csv file
in certain local path. The code that saves the tables are all
commentized so that it does not overwrite.

## package

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

\#condition table

``` r
#number of conditions

ncond <-117

#condition table

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

# the Main table

## Convergence rate

``` r
convergence<-read.csv("convergence/convergence.csv")


convergence<-convergence[,grep("Converged|condition", names(convergence))]
convergence<-convergence %>% 
  rename(
    MLM0 = hasConverged_MLM0,
    MLM1 = hasConverged_MLM1,
    MLM2 = hasConverged_MLM2,
    MLM3 = hasConverged_MLM3,
    MLM4 = hasConverged_MLM4,
    MLM5 = hasConverged_MLM5,
    MLM6 = hasConverged_MLM6
    )


convergence_table<-merge(conditions_con, convergence, key="condition")
```

need to trim the table after I see the results.

##### Convergence raw number

``` r
# 
# identical(dput(colnames(ho)),dput(colnames(convergence)))
# to make a empty holder with column names 
### following three lines are for making a empty holder with column names.

results_1 <- read.csv("model_results/results_1.csv")
results_1 <- results_1[,grep("hasConverged",names(results_1))]
convergence_raw <- data.frame(results_1[0,])

nMLM<-6 # the number of MLM models ~ from 0 to which number?

for(i in 1:ncond){
# read the data 
results <- read.csv(paste0("model_results/results_",i,".csv"))

#mean() You can use mean() to get the proportion of TRUE of a logical vector.


for( k in 0:nMLM){

  convergence_raw[i,paste0("hasConverged_MLM",k)]<- sum(results[,paste0("hasConverged_MLM",k)])

}

convergence_raw$condition[i] <- i 



}
#moved outside the loop so it goes faster
#write.csv(convergence_raw,"convergence/convergence_raw.csv")

  
con_raw<-convergence_raw
```

# type one err rate

``` r
# 
# 
 sig_rate<-read.csv("type1error/sig_rate.csv")
# sig_rate2<-read.csv("type1error/sig_rate_0.1alpha.csv")
# # merge with conditions_con
# 
 sig_rate<-merge(sig_rate,conditions_con,key="condition")
# sig_rate2<-merge(sig_rate2,conditions_con,key="condition")
# write.csv(sig_rate,"type1error/sig_rate_with_conditions.csv")
# write.csv(sig_rate2,"type1error/sig_rate_0.1alpha.csv")

# the condtion where the beta 2 is 0, the gender's main effect should be non-sig 

beta2<-sig_rate %>%
  filter(beta2==0)
#these guys should be non-sig. 
beta2<-beta2[,grep("_gender|ev|catbtw|ICC|beta|condition|CN",names(beta2))]

# the condition wehre bata 3 is 0, the interaction shoould be non-sig
beta3<-sig_rate %>%
  filter(beta3==0)
#these guys should be non-sig. 
beta3<-beta3[,grep("fxg|fxe|nxe|nxc|ICC|beta|condition|CN",names(beta3))]
# 
# write.csv(beta2,"type1error/beta2.csv")
# write.csv(beta3,"type1error/beta3.csv")
```

``` r
# just for check by model 
sig_rate2<-read.csv("type1error/sig_rate_without_singularity.csv")
sig_rate2<-merge(sig_rate2, conditions_con, key="condition")

sib0<-sig_rate[,grep("ICC|beta|CN|sib0|condition",names(sig_rate))]
 
sib1<-sig_rate[,grep("ICC|beta|CN|sib1|condition",names(sig_rate))]

sib2<-sig_rate[,grep("ICC|beta|CN|sib2|condition",names(sig_rate))]

sib3<-sig_rate[,grep("ICC|beta|CN|sib3|condition",names(sig_rate))]

MLM0<-sig_rate[,grep("ICC|beta|CN|MLM0|condition",names(sig_rate))]

MLM1<-sig_rate[,grep("ICC|beta|CN|MLM1|condition",names(sig_rate))]
MLM1_2<-sig_rate2[,grep("ICC|beta|CN|MLM1|condition",names(sig_rate2))]

MLM2<-sig_rate[,grep("ICC|beta|CN|MLM2|condition",names(sig_rate))]

MLM3<-sig_rate[,grep("ICC|beta|CN|MLM3|condition",names(sig_rate))]

MLM4<-sig_rate[,grep("ICC|beta|CN|MLM4|condition",names(sig_rate))]

MLM5<-sig_rate[,grep("ICC|beta|CN|MLM5|condition",names(sig_rate))]

MLM6<-sig_rate[,grep("ICC|beta|CN|MLM6|condition",names(sig_rate))]
MLM6_2<-sig_rate2[,grep("ICC|beta|CN|MLM6|condition",names(sig_rate2))]
# 
# write.csv(sib0,"type1error/sib0_sig_rate.csv")
# 
# write.csv(sib1,"type1error/sib1_sig_rate.csv")
# 
# write.csv(sib2,"type1error/sib2_sig_rate.csv")
# 
# write.csv(sib3,"type1error/sib3_sig_rate.csv")
# 
# 
# write.csv(MLM0,"type1error/MLM0_sig_rate.csv")
# 
# write.csv(MLM1,"type1error/MLM1_sig_rate.csv")
# 
# 
# write.csv(MLM2,"type1error/MLM2_sig_rate.csv")
# 
# 
# write.csv(MLM3,"type1error/MLM3_sig_rate.csv")
# 
# write.csv(MLM4,"type1error/MLM4_sig_rate.csv")
# 
# 
# write.csv(MLM5,"type1error/MLM5_sig_rate.csv")
# 
# write.csv(MLM6,"type1error/MLM6_sig_rate.csv")
```

# Appendix

# Appendix F : singularity rate

``` r
singular<-read.csv("appendix/singularity.csv")
```

singularity raw table

``` r
# load the data 

results_1<-read.csv("model_results/results_1.csv")
 results_1<-results_1[,grep("singularity",names(results_1))]
 singularity_raw<-data.frame(results_1[0,])
 
 # read the data 

for(i in 1:ncond){
# read the data 
results <- read.csv(paste0("model_results/results_",i,".csv"))

#mean() You can use mean() to get the proportion of TRUE of a logical vector.
for(k in 0:nMLM){  
  results[
    ## rows to replace
    results[[paste0("hasConverged_MLM", k)]]==0,
    ## column to replace
    grep(paste0("MLM",k),names(results))
  ] <- NA
}

# get the numbers 

for( n in 0:nMLM){

  singularity_raw[i,paste0("singularity_MLM",n)]<- sum(results[,paste0("singularity_MLM",n)], na.rm=T)

}

singularity_raw$condition[i] <- i 



}
#moved outside the loop so it goes faster
# write.csv(singularity_raw,"appendix/singularity_raw.csv")
```

#### For Check - please ignore this part.

``` r
## empty holder

# identical(dput(colnames(ho)),dput(colnames(type1error)))

# to make empty holder with col names 

results_1<-read.csv("model_results/results_1.csv")
results_1<-results_1[,grep("_p_",names(results_1))]
sig_rate<-data.frame(results_1[0,])

nMLM<-6 # how many MLM models we have? 0 to 6


ncond<-117

for(i in 1:ncond){
  
# read the data 
results <- read.csv(paste0("model_results/results_",i,".csv"))

# exclude non-convergence model  and singulairty 

nMLM<-6 # how many MLM models we have? 0 to 6

for(k in 0:nMLM){  
  results[
    ## rows to replace
    results[[paste0("hasConverged_MLM", k)]]==0 |results[[paste0("singularity_MLM", k)]]==TRUE ,
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

# 
# write.csv(sig_rate,"type1error/sig_rate_without_singularity.csv")
```
