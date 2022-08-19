closer_look_MLM5&MLM6
================
Yoo Ri Hwang
2022-07-01

## Overview

MLM5: gender composition (ff mm mixed) without interaction MLM6 gender
composition (FF mm mixed) with interaction

Under certain conditions, it shows 100% of singularity and 0 pseudo ICC
in MLM5 and MLM6. In this document, we attempt to get closer look into
how generated data are fitted into MLM5 and MLM6.

So, this is NOT main documents.

## load package

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
source('scripts/gen_dyadic.R')
source('scripts/gen_rmvn.R')
source('scripts/dyad_to_long.R')
source('scripts/reconstruct.R')
source('scripts/model_sib.R')
source('scripts/model_null.R')
source('scripts/model_mlm.R')
source('scripts/ICC_calculator.R')
source('scripts/myTryCatch.R')
```

    ## 
    ## Attaching package: 'NLP'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     annotate

``` r
source('scripts/hasConverged.R')
```

### Table

``` r
#overall singularity rate and psuedoICC across all multilevel model 
table<-read.csv("appendix/singularity_with_psuedoICC.csv")

# pick mlm5 and mlm 6

mlm56<-table[,grep("CN|beta|condition|ICC|MLM5|MLM6",names(table))]
```

## generate raw dataset, and fit into MLM5 and MLM6

### data generated per each condition (one dataset per one condition)

117 conditions

``` r
set.seed(5)
# condition information


conditions<-read.csv("conditions/Conditions_summary.csv")
ncond<-nrow(conditions)

# empty holder to put generated dataset


dyad_data<-vector(mode = "list", length = ncond)

### generate the data 

for (j in 1:ncond) {

  
    sim_data <- gen_dyadic(
    cn = conditions$cn[j],
    beta1 = conditions$beta1[j],
    beta2 = conditions$beta2[j],
    beta3 = conditions$beta3[j],
    con_icc_wtn = conditions$wtn[j],
    con_icc_btw = conditions$btw[j],
    cat_icc_wtn = conditions$wtn[j],
    cat_icc_btw = conditions$btw[j])
    
    dyad_data[[j]]<-sim_data
}
```

    ## Loading required package: MASS

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r
## holder  for generated data so that I can keep the data
ind_data<-vector(mode = "list", length = ncond)


# dyad_data has dyadic structure. let's convert it into individual structure so that I can run multilevel model

for(i in 1:ncond){
ind_data[[i]]<-dyad2ind(dyad_data[[i]])
}
```

### fit that into MLM5 and MLM6

``` r
# empty holder for generated data 

MLM5<-vector("list",length=ncond)
MLM6<-vector("list",length=ncond)
sum_MLM5<-vector("list",length=ncond)
sum_MLM6<-vector("list",length=ncond)

## fit the data into the models (MLM5, MLM6)

for(i in 1:ncond){

MLM5[[i]]<- yes2model(data=ind_data[[i]],
                gender_composition=3,
                interaction=FALSE)

MLM6[[i]]<- yes2model(data=ind_data[[i]],
                gender_composition=3,
                interaction=TRUE)
print(paste0(i," MLM5 cn:",conditions$cn[i]," ticc: ",conditions$target_ICC[i]," b2:",conditions$beta2[i]," b3:",conditions$beta3[i]))
print(sum_MLM5[[i]]<-summary(MLM5[[i]]))

print(paste0(i," MLM6 cn:",conditions$cn[i]," ticc: ",conditions$target_ICC[i]," b2:",conditions$beta2[i]," b3:",conditions$beta3[i]))
print(sum_MLM6[[i]]<-summary(MLM6[[i]]))
}
```

    ## Loading required package: lme4

    ## Warning: package 'lme4' was built under R version 4.1.3

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

    ## [1] "1 MLM5 cn:30 ticc: 0.2 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 219.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.74255 -0.48680  0.07618  0.53904  1.83007 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.9725   0.9861  
    ##  Residual             1.3724   1.1715  
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.16378    0.25364   0.646
    ## con          0.34330    0.01713  20.037
    ## ev1          1.45668    0.38077   3.826
    ## ev2         -0.95563    0.37837  -2.526
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.049              
    ## ev1  0.148 -0.112       
    ## ev2  0.155  0.008 -0.652
    ## [1] "1 MLM6 cn:30 ticc: 0.2 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 230.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.71032 -0.47868  0.08018  0.55049  1.81236 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 1.013    1.006   
    ##  Residual             1.417    1.190   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.175639   0.260855   0.673
    ## con          0.339479   0.020815  16.309
    ## ev1          1.470439   0.390444   3.766
    ## ev2         -0.967766   0.387667  -2.496
    ## con:ev1     -0.010703   0.034915  -0.307
    ## con:ev2      0.003517   0.027750   0.127
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.033                            
    ## ev1      0.161 -0.149                     
    ## ev2      0.140  0.061 -0.651              
    ## con:ev1 -0.133  0.485 -0.115  0.066       
    ## con:ev2  0.068 -0.167  0.083  0.012 -0.711
    ## [1] "2 MLM5 cn:120 ticc: 0.2 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 897
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.19711 -0.53956 -0.03114  0.59985  2.77792 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.7182   0.8475  
    ##  Residual             1.7503   1.3230  
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.131461   0.120309   1.093
    ## con          0.290459   0.009505  30.558
    ## ev1          1.301227   0.178842   7.276
    ## ev2         -0.948643   0.177785  -5.336
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.001              
    ## ev1  0.124 -0.109       
    ## ev2  0.124  0.004 -0.623
    ## [1] "2 MLM6 cn:120 ticc: 0.2 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 910
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.19409 -0.53474 -0.03501  0.58123  2.77522 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.7267   0.8524  
    ##  Residual             1.7572   1.3256  
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.119707   0.121854   0.982
    ## con          0.292018   0.010187  28.666
    ## ev1          1.283674   0.180786   7.101
    ## ev2         -0.937240   0.179190  -5.230
    ## con:ev1      0.012666   0.015465   0.819
    ## con:ev2     -0.008243   0.014878  -0.554
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.043                            
    ## ev1      0.136 -0.123                     
    ## ev2      0.110  0.032 -0.625              
    ## con:ev1 -0.120  0.201 -0.119  0.079       
    ## con:ev2  0.033  0.091  0.082 -0.017 -0.648
    ## [1] "3 MLM5 cn:510 ticc: 0.2 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 3890.3
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.71301 -0.61992 -0.00703  0.62405  2.76918 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.5087   0.7132  
    ##  Residual             2.1548   1.4679  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.129912   0.058633   2.216
    ## con          0.298671   0.004745  62.943
    ## ev1          1.110888   0.086783  12.801
    ## ev2         -0.657631   0.086779  -7.578
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.042              
    ## ev1  0.131 -0.012       
    ## ev2  0.131  0.007 -0.633
    ## [1] "3 MLM6 cn:510 ticc: 0.2 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 3905.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.71192 -0.61843 -0.00977  0.61830  2.76959 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.5074   0.7123  
    ##  Residual             2.1583   1.4691  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.130896   0.058652   2.232
    ## con          0.299026   0.005080  58.864
    ## ev1          1.107990   0.086845  12.758
    ## ev2         -0.653763   0.086912  -7.522
    ## con:ev1     -0.006855   0.007407  -0.925
    ## con:ev2      0.007188   0.007728   0.930
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.044                            
    ## ev1      0.130 -0.010                     
    ## ev2      0.132  0.020 -0.633              
    ## con:ev1 -0.011  0.086  0.034 -0.027       
    ## con:ev2  0.019  0.206 -0.026  0.052 -0.649
    ## [1] "4 MLM5 cn:30 ticc: 0.4 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 228.3
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.04397 -0.47004  0.00548  0.54713  2.52868 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.7998   0.8943  
    ##  Residual             1.7803   1.3343  
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) -0.10569    0.25616  -0.413
    ## con          0.31664    0.01639  19.324
    ## ev1          1.33592    0.38222   3.495
    ## ev2         -0.07375    0.38398  -0.192
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.060              
    ## ev1  0.156  0.040       
    ## ev2  0.147 -0.103 -0.655
    ## [1] "4 MLM6 cn:30 ticc: 0.4 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 238.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.06445 -0.49419 -0.00443  0.50971  2.46634 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.597    0.7727  
    ##  Residual             1.967    1.4024  
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) -0.08282    0.24945  -0.332
    ## con          0.31642    0.01784  17.738
    ## ev1          1.36223    0.37252   3.657
    ## ev2         -0.06428    0.37199  -0.173
    ## con:ev1      0.02660    0.02626   1.013
    ## con:ev2     -0.02150    0.02689  -0.800
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.039                            
    ## ev1      0.154  0.078                     
    ## ev2      0.150 -0.117 -0.653              
    ## con:ev1  0.079  0.113  0.089  0.009       
    ## con:ev2 -0.116  0.180  0.008 -0.043 -0.648
    ## [1] "5 MLM5 cn:120 ticc: 0.4 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 908.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.45350 -0.66102 -0.05297  0.70022  2.56298 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.3959   0.6292  
    ##  Residual             2.1105   1.4527  
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.202088   0.116373   1.737
    ## con          0.300916   0.007841  38.379
    ## ev1          1.493722   0.172043   8.682
    ## ev2         -0.341097   0.172142  -1.981
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.088              
    ## ev1  0.137  0.036       
    ## ev2  0.130 -0.049 -0.637
    ## [1] "5 MLM6 cn:120 ticc: 0.4 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 922.3
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.55532 -0.64758 -0.04676  0.68118  2.55505 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.381    0.6172  
    ##  Residual             2.132    1.4602  
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.202051   0.116323   1.737
    ## con          0.303266   0.008223  36.879
    ## ev1          1.505582   0.172876   8.709
    ## ev2         -0.338712   0.172086  -1.968
    ## con:ev1      0.003271   0.011418   0.287
    ## con:ev2      0.007953   0.012760   0.623
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.079                            
    ## ev1      0.140  0.040                     
    ## ev2      0.127 -0.038 -0.636              
    ## con:ev1  0.043 -0.052  0.105 -0.040       
    ## con:ev2 -0.036  0.263 -0.036  0.044 -0.620
    ## [1] "6 MLM5 cn:510 ticc: 0.4 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 3779.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.9610 -0.6099  0.0090  0.6470  2.7640 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.4675   0.6837  
    ##  Residual             1.9231   1.3868  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.202997   0.055735   3.642
    ## con          0.307337   0.003761  81.711
    ## ev1          1.538318   0.082651  18.612
    ## ev2         -0.458121   0.082631  -5.544
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.007              
    ## ev1  0.133  0.022       
    ## ev2  0.134 -0.005 -0.635
    ## [1] "6 MLM6 cn:510 ticc: 0.4 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 3797.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.96745 -0.60771  0.00924  0.64973  2.76118 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.4692   0.685   
    ##  Residual             1.9258   1.388   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.202633   0.055810   3.631
    ## con          0.307361   0.003830  80.243
    ## ev1          1.537905   0.082746  18.586
    ## ev2         -0.458295   0.082735  -5.539
    ## con:ev1     -0.002050   0.005623  -0.365
    ## con:ev2      0.002524   0.005535   0.456
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.004                            
    ## ev1      0.134  0.023                     
    ## ev2      0.133 -0.008 -0.635              
    ## con:ev1  0.023  0.105  0.012 -0.008       
    ## con:ev2 -0.008  0.061 -0.009 -0.009 -0.586
    ## [1] "7 MLM5 cn:30 ticc: 0.8 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 222.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.24497 -0.58235  0.01505  0.70729  1.40470 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 1.068    1.034   
    ##  Residual             1.408    1.187   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)   0.5914     0.2982   1.983
    ## con           0.3012     0.0108  27.899
    ## ev1           2.4199     0.4551   5.318
    ## ev2           0.8266     0.4554   1.815
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.060              
    ## ev1  0.215 -0.039       
    ## ev2  0.221  0.053 -0.715
    ## [1] "7 MLM6 cn:30 ticc: 0.8 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 232.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.3405 -0.5173 -0.0701  0.5886  1.3979 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.9284   0.9636  
    ##  Residual             1.4262   1.1942  
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.56059    0.28744   1.950
    ## con          0.30416    0.01180  25.778
    ## ev1          2.45045    0.43824   5.592
    ## ev2          0.75348    0.44000   1.712
    ## con:ev1      0.02197    0.01942   1.131
    ## con:ev2     -0.02987    0.01558  -1.917
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.047                            
    ## ev1      0.213 -0.031                     
    ## ev2      0.224  0.046 -0.715              
    ## con:ev1 -0.029  0.430  0.018 -0.028       
    ## con:ev2  0.054 -0.194 -0.035  0.082 -0.668
    ## [1] "8 MLM5 cn:120 ticc: 0.8 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 865.3
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.26482 -0.62375 -0.02651  0.63256  2.27937 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.4431   0.6656  
    ##  Residual             1.6546   1.2863  
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 0.606155   0.106779   5.677
    ## con         0.305083   0.004472  68.224
    ## ev1         2.266149   0.156888  14.444
    ## ev2         0.483452   0.156846   3.082
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.072              
    ## ev1  0.110  0.039       
    ## ev2  0.115 -0.032 -0.616
    ## [1] "8 MLM6 cn:120 ticc: 0.8 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 878.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.61975 -0.57966 -0.06735  0.59689  2.21162 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.4257   0.6525  
    ##  Residual             1.6584   1.2878  
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.594640   0.106351   5.591
    ## con          0.308007   0.004878  63.136
    ## ev1          2.278400   0.156160  14.590
    ## ev2          0.450416   0.157004   2.869
    ## con:ev1     -0.005114   0.007088  -0.721
    ## con:ev2      0.012722   0.007533   1.689
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.081                            
    ## ev1      0.106  0.044                     
    ## ev2      0.122 -0.067 -0.616              
    ## con:ev1  0.044  0.076 -0.046  0.054       
    ## con:ev2 -0.064  0.249  0.051 -0.115 -0.667
    ## [1] "9 MLM5 cn:510 ticc: 0.8 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 3859
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.61307 -0.59088 -0.00431  0.59102  2.75924 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.6663   0.8163  
    ##  Residual             1.9521   1.3972  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.63610    0.05922  10.741
    ## con          0.29935    0.00237 126.306
    ## ev1          2.27994    0.08779  25.972
    ## ev2          0.45943    0.08773   5.237
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.010              
    ## ev1  0.124 -0.082       
    ## ev2  0.123  0.075 -0.628
    ## [1] "9 MLM6 cn:510 ticc: 0.8 b2:0.1 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 3878.3
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.60714 -0.59151 -0.00254  0.58814  2.75555 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.6715   0.8195  
    ##  Residual             1.9531   1.3975  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error t value
    ## (Intercept)  0.6357576  0.0596719  10.654
    ## con          0.2993739  0.0024843 120.506
    ## ev1          2.2794480  0.0883205  25.809
    ## ev2          0.4595774  0.0881476   5.214
    ## con:ev1      0.0002277  0.0036575   0.062
    ## con:ev2     -0.0001059  0.0036919  -0.029
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.011                            
    ## ev1      0.129 -0.101                     
    ## ev2      0.123  0.093 -0.628              
    ## con:ev1 -0.101  0.114 -0.078  0.010       
    ## con:ev2  0.093  0.140  0.010  0.053 -0.629

    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')

    ## [1] "10 MLM5 cn:30 ticc: 0.2 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 294
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.26671 -0.60513 -0.00345  0.68267  2.21200 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             7.849    2.802   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.55504    0.42130   1.317
    ## con          0.29662    0.03139   9.449
    ## ev1          3.99878    0.64124   6.236
    ## ev2         -1.78457    0.62293  -2.865
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.211              
    ## ev1  0.230  0.239       
    ## ev2  0.177 -0.034 -0.675
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "10 MLM6 cn:30 ticc: 0.2 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 302
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.21230 -0.60437 -0.00242  0.66700  2.17077 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             8.039    2.835   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.65101    0.45791   1.422
    ## con          0.29925    0.03596   8.321
    ## ev1          4.25914    0.72947   5.839
    ## ev2         -1.95453    0.66388  -2.944
    ## con:ev1      0.04348    0.05434   0.800
    ## con:ev2     -0.03682    0.05429  -0.678
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.316                            
    ## ev1      0.338  0.290                     
    ## ev2      0.070 -0.068 -0.710              
    ## con:ev1  0.306  0.188  0.455 -0.310       
    ## con:ev2 -0.066  0.185 -0.283  0.244 -0.685
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "11 MLM5 cn:120 ticc: 0.2 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1145.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.37626 -0.68065  0.02927  0.71504  2.07746 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             6.733    2.595   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.67759    0.17844   3.797
    ## con          0.27360    0.01528  17.908
    ## ev1          3.60954    0.26639  13.550
    ## ev2         -1.88374    0.26564  -7.091
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.019              
    ## ev1  0.146  0.079       
    ## ev2  0.144 -0.023 -0.646
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "11 MLM6 cn:120 ticc: 0.2 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1152.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.4367 -0.6536  0.0415  0.7648  2.2885 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             6.676    2.584   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.70554    0.17856   3.951
    ## con          0.28894    0.01725  16.753
    ## ev1          3.65194    0.26672  13.692
    ## ev2         -1.91663    0.26508  -7.230
    ## con:ev1      0.02887    0.02544   1.135
    ## con:ev2      0.01135    0.02664   0.426
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.042                            
    ## ev1      0.155  0.093                     
    ## ev2      0.137 -0.045 -0.647              
    ## con:ev1  0.094  0.119  0.101 -0.052       
    ## con:ev2 -0.043  0.250 -0.050  0.007 -0.686
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "12 MLM5 cn:510 ticc: 0.2 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4688.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.3912 -0.6600 -0.0063  0.6465  2.7421 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             5.733    2.394   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.350928   0.076910   4.563
    ## con          0.307000   0.006681  45.949
    ## ev1          3.922631   0.112657  34.819
    ## ev2         -2.140078   0.112658 -18.996
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.010              
    ## ev1  0.099  0.006       
    ## ev2  0.100 -0.008 -0.602
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "12 MLM6 cn:510 ticc: 0.2 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4703.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.38205 -0.65924 -0.00522  0.66743  2.73975 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             5.739    2.396   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.350500   0.076958   4.554
    ## con          0.306152   0.006909  44.312
    ## ev1          3.923170   0.112724  34.803
    ## ev2         -2.140759   0.112739 -18.989
    ## con:ev1     -0.009344   0.010214  -0.915
    ## con:ev2      0.005378   0.010109   0.532
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.011                            
    ## ev1      0.099  0.006                     
    ## ev2      0.100 -0.010 -0.602              
    ## con:ev1  0.006  0.125 -0.005  0.007       
    ## con:ev2 -0.010  0.096  0.007 -0.017 -0.613
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "13 MLM5 cn:30 ticc: 0.4 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 290.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5281 -0.7679 -0.0835  0.7075  1.7089 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             7.376    2.716   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.78757    0.37948   2.075
    ## con          0.31659    0.02986  10.601
    ## ev1          4.71379    0.56931   8.280
    ## ev2         -1.75128    0.56422  -3.104
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.096              
    ## ev1  0.165  0.134       
    ## ev2  0.153 -0.012 -0.651
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "13 MLM6 cn:30 ticc: 0.4 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 299
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.53636 -0.70866  0.00299  0.68738  1.84004 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             7.555    2.749   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.730198   0.390707   1.869
    ## con          0.306842   0.033308   9.212
    ## ev1          4.638165   0.591414   7.843
    ## ev2         -1.709169   0.578015  -2.957
    ## con:ev1     -0.025540   0.048189  -0.530
    ## con:ev2     -0.005558   0.051728  -0.107
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.141                            
    ## ev1      0.192  0.139                     
    ## ev2      0.128 -0.014 -0.661              
    ## con:ev1  0.145  0.064  0.225 -0.154       
    ## con:ev2 -0.014  0.265 -0.140  0.114 -0.670
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "14 MLM5 cn:120 ticc: 0.4 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1136.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.19542 -0.70312  0.00236  0.69139  2.65897 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             6.475    2.545   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.88364    0.16772   5.268
    ## con          0.32214    0.01334  24.149
    ## ev1          4.37419    0.24549  17.818
    ## ev2         -1.16447    0.24476  -4.758
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.043              
    ## ev1  0.093  0.086       
    ## ev2  0.088 -0.038 -0.593
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "14 MLM6 cn:120 ticc: 0.4 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1148.5
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.23669 -0.71316  0.01357  0.68108  2.62123 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             6.521    2.554   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.881188   0.169182   5.209
    ## con          0.322001   0.013535  23.789
    ## ev1          4.371398   0.248200  17.612
    ## ev2         -1.160749   0.246152  -4.716
    ## con:ev1     -0.002324   0.020165  -0.115
    ## con:ev2     -0.007876   0.018548  -0.425
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.058                            
    ## ev1      0.104  0.103                     
    ## ev2      0.080 -0.047 -0.595              
    ## con:ev1  0.101  0.147  0.122 -0.064       
    ## con:ev2 -0.050 -0.089 -0.069  0.023 -0.537
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "15 MLM5 cn:510 ticc: 0.4 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4829
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.73307 -0.66945 -0.01087  0.68671  2.58989 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.00     0.000   
    ##  Residual             6.58     2.565   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.756005   0.083658   9.037
    ## con          0.305442   0.006228  49.041
    ## ev1          4.502353   0.123539  36.445
    ## ev2         -1.362452   0.123549 -11.028
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.023              
    ## ev1  0.122 -0.038       
    ## ev2  0.120  0.040 -0.624
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "15 MLM6 cn:510 ticc: 0.4 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4844.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.72500 -0.65656 -0.00951  0.68662  2.58738 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             6.592    2.568   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error t value
    ## (Intercept)  7.565e-01  8.385e-02   9.022
    ## con          3.052e-01  6.416e-03  47.566
    ## ev1          4.504e+00  1.239e-01  36.359
    ## ev2         -1.363e+00  1.237e-01 -11.018
    ## con:ev1     -1.510e-03  9.331e-03  -0.162
    ## con:ev2      5.174e-05  9.498e-03   0.005
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.022                            
    ## ev1      0.123 -0.047                     
    ## ev2      0.120  0.044 -0.624              
    ## con:ev1 -0.047  0.079 -0.053  0.012       
    ## con:ev2  0.044  0.129  0.012  0.010 -0.607
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "16 MLM5 cn:30 ticc: 0.8 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 300.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.54473 -0.74834 -0.02558  0.73627  1.98106 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             8.692    2.948   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.96883    0.46627   4.223
    ## con          0.30805    0.01733  17.779
    ## ev1          7.73999    0.71681  10.798
    ## ev2          0.42615    0.71829   0.593
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.022              
    ## ev1  0.219  0.115       
    ## ev2  0.219  0.131 -0.688
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "16 MLM6 cn:30 ticc: 0.8 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 309.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6019 -0.7434  0.1024  0.6905  1.8959 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             8.694    2.949   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  2.03484    0.48722   4.176
    ## con          0.31119    0.02246  13.853
    ## ev1          7.91437    0.73485  10.770
    ## ev2          0.21203    0.74874   0.283
    ## con:ev1      0.04179    0.03087   1.354
    ## con:ev2     -0.02764    0.03768  -0.733
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.198                            
    ## ev1      0.182  0.014                     
    ## ev2      0.235  0.200 -0.706              
    ## con:ev1  0.016 -0.081  0.201 -0.250       
    ## con:ev2  0.184  0.485 -0.208  0.273 -0.739
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "17 MLM5 cn:120 ticc: 0.8 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1139.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.19871 -0.68617  0.04627  0.72736  2.35183 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             6.543    2.558   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 1.886123   0.174137   10.83
    ## con         0.295017   0.007652   38.55
    ## ev1         7.121442   0.258545   27.54
    ## ev2         1.088569   0.258583    4.21
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.033              
    ## ev1  0.136  0.056       
    ## ev2  0.133 -0.058 -0.638
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "17 MLM6 cn:120 ticc: 0.8 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1153.3
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.21778 -0.67190  0.05261  0.70682  2.28916 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             6.575    2.564   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  1.876545   0.175352  10.702
    ## con          0.291401   0.008611  33.839
    ## ev1          7.094745   0.260913  27.192
    ## ev2          1.101141   0.259574   4.242
    ## con:ev1     -0.010196   0.013794  -0.739
    ## con:ev2      0.001147   0.012048   0.095
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.055                            
    ## ev1      0.144  0.099                     
    ## ev2      0.129 -0.076 -0.638              
    ## con:ev1  0.092  0.353  0.108 -0.038       
    ## con:ev2 -0.080 -0.030 -0.043 -0.001 -0.680
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "18 MLM5 cn:510 ticc: 0.8 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4848.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.4934 -0.7283  0.0079  0.6529  2.8181 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             6.699    2.588   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 2.128050   0.085341  24.936
    ## con         0.300770   0.003597  83.625
    ## ev1         6.953324   0.126522  54.958
    ## ev2         1.206103   0.126677   9.521
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.019              
    ## ev1  0.133 -0.017       
    ## ev2  0.134  0.052 -0.635
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "18 MLM6 cn:510 ticc: 0.8 b2:0.3 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4865.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.49686 -0.71846  0.00062  0.65463  2.81337 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.00     0.00    
    ##  Residual             6.71     2.59    
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  2.129576   0.085598  24.879
    ## con          0.300673   0.003776  79.621
    ## ev1          6.952244   0.126745  54.852
    ## ev2          1.209443   0.127132   9.513
    ## con:ev1     -0.003231   0.005436  -0.594
    ## con:ev2      0.002210   0.005732   0.386
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.035                            
    ## ev1      0.130 -0.028                     
    ## ev2      0.139  0.066 -0.636              
    ## con:ev1 -0.029  0.050  0.013 -0.043       
    ## con:ev2  0.064  0.200 -0.041  0.074 -0.630
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "19 MLM5 cn:30 ticc: 0.2 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 343.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.76660 -0.66433 -0.00968  0.78779  1.87996 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             19.23    4.386   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.85569    0.61494   1.392
    ## con          0.30386    0.06312   4.814
    ## ev1          6.08327    0.91124   6.676
    ## ev2         -3.28239    0.91640  -3.582
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.127              
    ## ev1  0.151  0.020       
    ## ev2  0.166 -0.108 -0.653
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "19 MLM6 cn:30 ticc: 0.2 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 349.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.76602 -0.64782 -0.01741  0.77496  1.84022 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             19.83    4.454   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.87670    0.62957   1.393
    ## con          0.28040    0.07692   3.645
    ## ev1          6.15259    0.93742   6.563
    ## ev2         -3.33132    0.94679  -3.519
    ## con:ev1     -0.07231    0.13172  -0.549
    ## con:ev2      0.03343    0.09730   0.344
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.142                            
    ## ev1      0.146 -0.054                     
    ## ev2      0.174 -0.042 -0.660              
    ## con:ev1 -0.047  0.544 -0.143  0.111       
    ## con:ev2 -0.050 -0.316  0.152 -0.182 -0.707
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "20 MLM5 cn:120 ticc: 0.2 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1320.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.07720 -0.61488  0.02443  0.66562  2.15707 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             14.15    3.761   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)   0.7867     0.2535   3.103
    ## con           0.2739     0.0215  12.742
    ## ev1           6.5407     0.3747  17.454
    ## ev2          -3.2266     0.3747  -8.610
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.015              
    ## ev1  0.124  0.025       
    ## ev2  0.125 -0.027 -0.626
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "20 MLM6 cn:120 ticc: 0.2 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1329.5
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.10771 -0.64950 -0.02533  0.65133  2.27972 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             14.18    3.765   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.78366    0.25391   3.086
    ## con          0.27901    0.02197  12.703
    ## ev1          6.54823    0.37513  17.456
    ## ev2         -3.24177    0.37536  -8.636
    ## con:ev1      0.01224    0.03251   0.377
    ## con:ev2      0.02325    0.03159   0.736
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.014                            
    ## ev1      0.124  0.028                     
    ## ev2      0.125 -0.031 -0.626              
    ## con:ev1  0.027  0.129  0.006  0.009       
    ## con:ev2 -0.032  0.047  0.009 -0.035 -0.591
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "21 MLM5 cn:510 ticc: 0.2 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5680
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.32132 -0.66866 -0.00784  0.70680  2.44891 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.0     
    ##  Residual             15.21    3.9     
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.71960    0.12882   5.586
    ## con          0.31268    0.01105  28.288
    ## ev1          6.62937    0.19067  34.769
    ## ev2         -3.37503    0.19087 -17.682
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.062              
    ## ev1  0.132 -0.024       
    ## ev2  0.136  0.052 -0.635
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "21 MLM6 cn:510 ticc: 0.2 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5692.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.32486 -0.66599 -0.00162  0.71135  2.46699 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             15.23    3.903   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.716611   0.129288   5.543
    ## con          0.310064   0.011798  26.281
    ## ev1          6.629075   0.191224  34.666
    ## ev2         -3.379397   0.192381 -17.566
    ## con:ev1     -0.007677   0.017433  -0.440
    ## con:ev2     -0.001904   0.017702  -0.108
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.077                            
    ## ev1      0.127 -0.029                     
    ## ev2      0.144  0.073 -0.637              
    ## con:ev1 -0.029  0.124  0.051 -0.065       
    ## con:ev2  0.072  0.167 -0.064  0.118 -0.647
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "22 MLM5 cn:30 ticc: 0.4 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 333.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.75990 -0.46205 -0.07017  0.63119  1.91282 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             15.59    3.948   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.78281    0.52942   1.479
    ## con          0.28987    0.03804   7.620
    ## ev1          7.80121    0.78029   9.998
    ## ev2         -2.66340    0.79641  -3.344
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.084              
    ## ev1  0.104 -0.098       
    ## ev2  0.129  0.222 -0.619
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "22 MLM6 cn:30 ticc: 0.4 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 341
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.76524 -0.45341 -0.09982  0.60057  1.90434 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             16.13    4.017   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.827118   0.558422   1.481
    ## con          0.293361   0.041257   7.111
    ## ev1          7.757048   0.805887   9.625
    ## ev2         -2.589311   0.845904  -3.061
    ## con:ev1     -0.003665   0.061964  -0.059
    ## con:ev2      0.016893   0.060876   0.278
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.122                            
    ## ev1      0.057 -0.138                     
    ## ev2      0.195  0.253 -0.630              
    ## con:ev1 -0.133  0.170  0.020 -0.134       
    ## con:ev2  0.260  0.120 -0.143  0.281 -0.647
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "23 MLM5 cn:120 ticc: 0.4 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1319.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.99625 -0.62041 -0.00435  0.66266  2.28507 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             14.08    3.752   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.34053    0.25360   5.286
    ## con          0.32718    0.01969  16.617
    ## ev1          7.79859    0.37430  20.835
    ## ev2         -2.46722    0.37548  -6.571
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.075              
    ## ev1  0.120 -0.056       
    ## ev2  0.131  0.097 -0.628
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "23 MLM6 cn:120 ticc: 0.4 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1329.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.99901 -0.58165  0.00276  0.70307  2.32268 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             14.13    3.759   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.32326    0.25578   5.173
    ## con          0.32222    0.02022  15.933
    ## ev1          7.81786    0.37600  20.792
    ## ev2         -2.50823    0.38097  -6.584
    ## con:ev1     -0.01359    0.02950  -0.461
    ## con:ev2     -0.01574    0.02973  -0.529
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.086                            
    ## ev1      0.109 -0.065                     
    ## ev2      0.147  0.116 -0.630              
    ## con:ev1 -0.066  0.087  0.036 -0.075       
    ## con:ev2  0.117  0.110 -0.075  0.158 -0.601
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "24 MLM5 cn:510 ticc: 0.4 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5620.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.36892 -0.69581 -0.02486  0.69678  2.20929 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             14.34    3.786   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  1.342881   0.124332   10.80
    ## con          0.303084   0.009007   33.65
    ## ev1          7.850876   0.183897   42.69
    ## ev2         -2.413345   0.183895  -13.12
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.036              
    ## ev1  0.129  0.007       
    ## ev2  0.128 -0.005 -0.630
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "24 MLM6 cn:510 ticc: 0.4 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5634.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.3643 -0.6992 -0.0254  0.6966  2.2037 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.00    
    ##  Residual             14.36    3.79    
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error t value
    ## (Intercept)  1.3427740  0.1244652  10.788
    ## con          0.3027959  0.0093303  32.453
    ## ev1          7.8508872  0.1842486  42.610
    ## ev2         -2.4139737  0.1841814 -13.107
    ## con:ev1      0.0001391  0.0136622   0.010
    ## con:ev2     -0.0015434  0.0138006  -0.112
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.038                            
    ## ev1      0.129  0.012                     
    ## ev2      0.128 -0.001 -0.630              
    ## con:ev1  0.012  0.098  0.043 -0.025       
    ## con:ev2 -0.002  0.127 -0.024  0.034 -0.615
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "25 MLM5 cn:30 ticc: 0.8 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 305.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.86095 -0.51332 -0.04137  0.59670  2.09485 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.00    
    ##  Residual             9.364    3.06    
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  3.47552    0.40004   8.688
    ## con          0.29342    0.01829  16.045
    ## ev1         11.70110    0.55905  20.930
    ## ev2          2.24158    0.55894   4.010
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.157              
    ## ev1  0.006 -0.036       
    ## ev2 -0.005  0.030 -0.501
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "25 MLM6 cn:30 ticc: 0.8 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 316
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.84108 -0.49977 -0.06406  0.63083  2.06840 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             9.696    3.114   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  3.464393   0.408977   8.471
    ## con          0.297775   0.024064  12.374
    ## ev1         11.706463   0.579730  20.193
    ## ev2          2.214407   0.580089   3.817
    ## con:ev1     -0.003085   0.030484  -0.101
    ## con:ev2      0.010730   0.042007   0.255
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.178                            
    ## ev1      0.007  0.012                     
    ## ev2      0.008 -0.090 -0.508              
    ## con:ev1  0.013 -0.312 -0.189  0.161       
    ## con:ev2 -0.073  0.600  0.117 -0.195 -0.747
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "26 MLM5 cn:120 ticc: 0.8 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1292.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.63011 -0.51539  0.00611  0.48842  2.43015 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             12.49    3.534   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  3.36669    0.23252  14.479
    ## con          0.29816    0.01052  28.336
    ## ev1         11.74262    0.33752  34.791
    ## ev2          2.07615    0.33919   6.121
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.099              
    ## ev1  0.085 -0.089       
    ## ev2  0.062  0.133 -0.583
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "26 MLM6 cn:120 ticc: 0.8 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1304.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.56536 -0.50829  0.01118  0.50536  2.45515 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.00    
    ##  Residual             12.53    3.54    
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  3.379599   0.235216  14.368
    ## con          0.300363   0.010716  28.030
    ## ev1         11.692042   0.343733  34.015
    ## ev2          2.090073   0.339960   6.148
    ## con:ev1      0.005297   0.015489   0.342
    ## con:ev2      0.010626   0.015712   0.676
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.088                            
    ## ev1      0.093 -0.110                     
    ## ev2      0.062  0.137 -0.580              
    ## con:ev1 -0.111  0.062 -0.160  0.025       
    ## con:ev2  0.135  0.102  0.024  0.010 -0.585
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "27 MLM5 cn:510 ticc: 0.8 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5701.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.20219 -0.82460 -0.00605  0.82261  2.34813 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             15.51    3.939   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  3.394181   0.133247  25.473
    ## con          0.296928   0.005574  53.267
    ## ev1         11.651918   0.199225  58.486
    ## ev2          1.899696   0.199364   9.529
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.007              
    ## ev1  0.156 -0.025       
    ## ev2  0.157  0.045 -0.658
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "27 MLM6 cn:510 ticc: 0.8 b2:0.5 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5716.5
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.21445 -0.83380 -0.00614  0.81705  2.35846 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             15.54    3.942   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept) 3.394e+00  1.336e-01  25.414
    ## con         2.990e-01  6.277e-03  47.641
    ## ev1         1.165e+01  1.995e-01  58.371
    ## ev2         1.904e+00  1.998e-01   9.526
    ## con:ev1     5.180e-03  9.616e-03   0.539
    ## con:ev2     2.225e-04  9.317e-03   0.024
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.015                            
    ## ev1      0.155 -0.040                     
    ## ev2      0.159  0.057 -0.658              
    ## con:ev1 -0.039  0.227 -0.013 -0.018       
    ## con:ev2  0.058  0.137 -0.018  0.052 -0.682
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')

    ## [1] "28 MLM5 cn:30 ticc: 0.2 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 251.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.0102 -0.5359 -0.1246  0.6882  2.1909 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.0472   0.2173  
    ##  Residual             3.6200   1.9026  
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.08132    0.28925   0.281
    ## con          0.25463    0.02451  10.388
    ## ev1          1.10290    0.42954   2.568
    ## ev2         -0.65446    0.43095  -1.519
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.202              
    ## ev1  0.200 -0.077       
    ## ev2  0.162  0.111 -0.690
    ## [1] "28 MLM6 cn:30 ticc: 0.2 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 254.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.0342 -0.6187 -0.1775  0.6078  2.3772 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev. 
    ##  pid      (Intercept) 1.883e-15 4.339e-08
    ##  Residual             3.485e+00 1.867e+00
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) -0.09233    0.29194  -0.316
    ## con          0.28342    0.03922   7.226
    ## ev1          0.69201    0.45750   1.513
    ## ev2         -0.45596    0.42719  -1.067
    ## con:ev1      0.12814    0.05919   2.165
    ## con:ev2     -0.08642    0.06292  -1.373
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.264                            
    ## ev1      0.291 -0.240                     
    ## ev2      0.097  0.118 -0.696              
    ## con:ev1 -0.249  0.184 -0.382  0.211       
    ## con:ev2  0.107  0.357  0.186 -0.151 -0.764
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "29 MLM5 cn:120 ticc: 0.2 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 998.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.85856 -0.57060  0.02599  0.56641  2.12224 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.9173   0.9578  
    ##  Residual             2.8234   1.6803  
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)   0.1477     0.1443   1.024
    ## con           0.3015     0.0110  27.420
    ## ev1           1.3461     0.2123   6.342
    ## ev2          -0.7846     0.2123  -3.696
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.030              
    ## ev1  0.113  0.004       
    ## ev2  0.113  0.011 -0.615
    ## [1] "29 MLM6 cn:120 ticc: 0.2 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 967.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.2222 -0.5674  0.0349  0.5425  2.2697 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.8077   0.8987  
    ##  Residual             2.3258   1.5250  
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.14345    0.13272   1.081
    ## con          0.30374    0.01009  30.096
    ## ev1          1.39106    0.19536   7.120
    ## ev2         -0.83498    0.19540  -4.273
    ## con:ev1      0.08404    0.01469   5.720
    ## con:ev2     -0.09255    0.01429  -6.477
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.031                            
    ## ev1      0.113  0.006                     
    ## ev2      0.114  0.010 -0.616              
    ## con:ev1  0.006  0.082  0.033 -0.026       
    ## con:ev2  0.011  0.003 -0.026  0.037 -0.545
    ## [1] "30 MLM5 cn:510 ticc: 0.2 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4167.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.41195 -0.67680  0.00943  0.61762  2.66866 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.5968   0.7726  
    ##  Residual             2.8891   1.6997  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.051244   0.066947   0.765
    ## con          0.297622   0.005242  56.772
    ## ev1          1.335328   0.099484  13.423
    ## ev2         -0.687785   0.099399  -6.919
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.023              
    ## ev1  0.139  0.042       
    ## ev2  0.139  0.004 -0.639
    ## [1] "30 MLM6 cn:510 ticc: 0.2 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4056.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.54198 -0.61063 -0.00323  0.60438  2.77076 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.6745   0.8213  
    ##  Residual             2.4352   1.5605  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.070415   0.064556   1.091
    ## con          0.297584   0.005148  57.809
    ## ev1          1.403508   0.095968  14.625
    ## ev2         -0.737191   0.095799  -7.695
    ## con:ev1      0.080719   0.007617  10.597
    ## con:ev2     -0.080073   0.007611 -10.521
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.036                            
    ## ev1      0.141  0.045                     
    ## ev2      0.136 -0.001 -0.640              
    ## con:ev1  0.045  0.128  0.063 -0.046       
    ## con:ev2 -0.001  0.126 -0.046  0.033 -0.629

    ## boundary (singular) fit: see help('isSingular')

    ## [1] "31 MLM5 cn:30 ticc: 0.4 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 260.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.2465 -0.6504  0.1242  0.8447  1.9299 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             4.293    2.072   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.52700    0.30460   1.730
    ## con          0.31564    0.02017  15.650
    ## ev1          1.67808    0.46426   3.615
    ## ev2         -0.87304    0.46646  -1.872
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.016              
    ## ev1  0.189 -0.128       
    ## ev2  0.184  0.160 -0.694
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "31 MLM6 cn:30 ticc: 0.4 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 264.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.3858 -0.4714  0.1054  0.6730  1.6763 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.5476   0.740   
    ##  Residual             3.5547   1.885   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.33463    0.32868   1.018
    ## con          0.30303    0.02735  11.081
    ## ev1          1.69566    0.49929   3.396
    ## ev2         -1.06447    0.49763  -2.139
    ## con:ev1      0.06737    0.04399   1.531
    ## con:ev2     -0.09847    0.04028  -2.445
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.035                            
    ## ev1      0.202 -0.254                     
    ## ev2      0.193  0.245 -0.696              
    ## con:ev1 -0.240  0.365 -0.187  0.021       
    ## con:ev2  0.252  0.115  0.022  0.135 -0.742

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "32 MLM5 cn:120 ticc: 0.4 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 971.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.25369 -0.74984  0.01541  0.60558  2.80569 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             3.216    1.793   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 0.308671   0.122566   2.518
    ## con         0.291359   0.009242  31.525
    ## ev1         1.358859   0.181749   7.477
    ## ev2         0.055942   0.181155   0.309
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.095              
    ## ev1  0.142  0.092       
    ## ev2  0.138  0.044 -0.629
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "32 MLM6 cn:120 ticc: 0.4 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 968.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.34335 -0.71436  0.00814  0.60832  2.87009 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             3.019    1.738   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.323212   0.120552   2.681
    ## con          0.294034   0.009438  31.153
    ## ev1          1.505779   0.180113   8.360
    ## ev2         -0.067251   0.178113  -0.378
    ## con:ev1      0.049934   0.014473   3.450
    ## con:ev2     -0.053697   0.013458  -3.990
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.143                            
    ## ev1      0.155  0.122                     
    ## ev2      0.124  0.020 -0.641              
    ## con:ev1  0.119  0.229  0.204 -0.156       
    ## con:ev2  0.021  0.023 -0.166  0.150 -0.633
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "33 MLM5 cn:510 ticc: 0.4 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4173.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.79838 -0.59781 -0.00664  0.63357  2.80388 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.5039   0.7099  
    ##  Residual             2.9870   1.7283  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.20436    0.06656   3.070
    ## con          0.29630    0.00446  66.437
    ## ev1          1.39749    0.09900  14.116
    ## ev2         -0.34539    0.09908  -3.486
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.027              
    ## ev1  0.144  0.024       
    ## ev2  0.142 -0.046 -0.645
    ## [1] "33 MLM6 cn:510 ticc: 0.4 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4064.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.89840 -0.54852  0.01167  0.58491  2.90901 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.6029   0.7765  
    ##  Residual             2.5145   1.5857  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.242201   0.064315   3.766
    ## con          0.298411   0.004547  65.630
    ## ev1          1.427942   0.095709  14.920
    ## ev2         -0.340423   0.095634  -3.560
    ## con:ev1      0.071294   0.006899  10.335
    ## con:ev2     -0.071864   0.006668 -10.777
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.024                            
    ## ev1      0.144  0.042                     
    ## ev2      0.142 -0.049 -0.644              
    ## con:ev1  0.041  0.199  0.048 -0.005       
    ## con:ev2 -0.050  0.103 -0.006 -0.012 -0.653

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "34 MLM5 cn:30 ticc: 0.8 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 254.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.18938 -0.69705  0.04551  0.78288  2.23681 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             3.766    1.941   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.54853    0.25698   2.135
    ## con          0.29006    0.01165  24.908
    ## ev1          3.07975    0.36609   8.413
    ## ev2          0.30060    0.36572   0.822
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.178              
    ## ev1  0.071  0.050       
    ## ev2  0.066  0.022 -0.563
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "34 MLM6 cn:30 ticc: 0.8 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 265
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.94264 -0.62738  0.01463  0.72195  2.23896 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev. 
    ##  pid      (Intercept) 1.656e-15 4.070e-08
    ##  Residual             3.736e+00 1.933e+00
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.56035    0.25794   2.172
    ## con          0.29103    0.01211  24.037
    ## ev1          3.21368    0.37593   8.549
    ## ev2          0.17770    0.37308   0.476
    ## con:ev1      0.02549    0.01820   1.400
    ## con:ev2     -0.02498    0.01756  -1.423
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.205                            
    ## ev1      0.085  0.083                     
    ## ev2      0.064  0.024 -0.577              
    ## con:ev1  0.080  0.172  0.242 -0.166       
    ## con:ev2  0.024  0.071 -0.170  0.212 -0.625
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "35 MLM5 cn:120 ticc: 0.8 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 986.3
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.34016 -0.74680  0.02796  0.71741  2.80366 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.1625   0.4031  
    ##  Residual             3.2580   1.8050  
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 0.553500   0.134601   4.112
    ## con         0.308585   0.005588  55.228
    ## ev1         2.507227   0.191280  13.108
    ## ev2         0.196364   0.191033   1.028
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.291              
    ## ev1  0.144 -0.051       
    ## ev2  0.130 -0.002 -0.635
    ## [1] "35 MLM6 cn:120 ticc: 0.8 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 997.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.37533 -0.73473  0.01371  0.71079  2.76812 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.2904   0.5389  
    ##  Residual             3.0996   1.7606  
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.496452   0.138803   3.577
    ## con          0.313936   0.006122  51.282
    ## ev1          2.324488   0.210472  11.044
    ## ev2          0.289542   0.202973   1.427
    ## con:ev1      0.022080   0.009873   2.236
    ## con:ev2     -0.010463   0.008198  -1.276
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.334                            
    ## ev1      0.197 -0.187                     
    ## ev2      0.095  0.060 -0.648              
    ## con:ev1 -0.175  0.373 -0.389  0.225       
    ## con:ev2  0.065 -0.154  0.261 -0.297 -0.646
    ## [1] "36 MLM5 cn:510 ticc: 0.8 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4133.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.4123 -0.6504  0.0030  0.6564  3.0791 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.1795   0.4237  
    ##  Residual             3.1402   1.7721  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 0.678891   0.063293  10.726
    ## con         0.296954   0.002416 122.903
    ## ev1         2.463313   0.094361  26.105
    ## ev2         0.297060   0.094464   3.145
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.077              
    ## ev1  0.158  0.049       
    ## ev2  0.148 -0.067 -0.656
    ## [1] "36 MLM6 cn:510 ticc: 0.8 b2:0.1 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4126
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3777 -0.6712 -0.0132  0.6689  3.0841 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.2604   0.5103  
    ##  Residual             2.9860   1.7280  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.706509   0.063582  11.112
    ## con          0.297241   0.002548 116.668
    ## ev1          2.506712   0.095281  26.309
    ## ev2          0.280172   0.094641   2.960
    ## con:ev1      0.017798   0.003814   4.667
    ## con:ev2     -0.017912   0.003748  -4.779
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.078                            
    ## ev1      0.164  0.079                     
    ## ev2      0.145 -0.071 -0.655              
    ## con:ev1  0.079  0.161  0.122 -0.040       
    ## con:ev2 -0.072  0.111 -0.040  0.022 -0.638

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "37 MLM5 cn:30 ticc: 0.2 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 287.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.59257 -0.57432 -0.05316  0.65063  2.08305 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             6.967    2.639   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.63497    0.35396   1.794
    ## con          0.28593    0.03347   8.543
    ## ev1          4.19964    0.52576   7.988
    ## ev2         -2.55871    0.53755  -4.760
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.084              
    ## ev1  0.125 -0.158       
    ## ev2  0.087  0.259 -0.628
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "37 MLM6 cn:30 ticc: 0.2 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 292.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7393 -0.4081 -0.0305  0.6650  1.9004 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             6.748    2.598   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.43743    0.36531   1.197
    ## con          0.26704    0.03560   7.501
    ## ev1          4.29690    0.53901   7.972
    ## ev2         -2.76169    0.54031  -5.111
    ## con:ev1      0.04870    0.05221   0.933
    ## con:ev2     -0.10332    0.05417  -1.907
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.030                            
    ## ev1      0.120 -0.235                     
    ## ev2      0.127  0.305 -0.625              
    ## con:ev1 -0.237  0.103 -0.188 -0.034       
    ## con:ev2  0.296  0.207 -0.033  0.173 -0.657
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "38 MLM5 cn:120 ticc: 0.2 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1175.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.54446 -0.72912 -0.01418  0.70430  2.44858 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             7.651    2.766   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.38124    0.19562   1.949
    ## con          0.26272    0.01668  15.752
    ## ev1          4.15127    0.29187  14.223
    ## ev2         -2.07543    0.29189  -7.110
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.094              
    ## ev1  0.165 -0.021       
    ## ev2  0.161  0.024 -0.664
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "38 MLM6 cn:120 ticc: 0.2 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1179.3
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.48573 -0.62984  0.01135  0.62294  2.61223 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             7.478    2.735   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.35263    0.19368   1.821
    ## con          0.27519    0.01739  15.827
    ## ev1          4.05257    0.29088  13.932
    ## ev2         -2.02463    0.28950  -6.994
    ## con:ev1      0.07376    0.02708   2.724
    ## con:ev2     -0.04460    0.02419  -1.844
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.104                            
    ## ev1      0.170 -0.056                     
    ## ev2      0.157  0.031 -0.664              
    ## con:ev1 -0.054  0.273 -0.125  0.061       
    ## con:ev2  0.033 -0.046  0.068 -0.078 -0.628
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "39 MLM5 cn:510 ticc: 0.2 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4968.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5278 -0.6898  0.0374  0.6921  3.4313 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             7.549    2.748   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.42031    0.09125   4.606
    ## con          0.31246    0.00770  40.580
    ## ev1          3.85167    0.13567  28.389
    ## ev2         -2.08046    0.13566 -15.335
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.022              
    ## ev1  0.141  0.032       
    ## ev2  0.140 -0.030 -0.642
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "39 MLM6 cn:510 ticc: 0.2 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4928.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5898 -0.6716  0.0392  0.6467  3.5248 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             7.168    2.677   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.447596   0.088998   5.029
    ## con          0.305256   0.008006  38.128
    ## ev1          3.882219   0.132375  29.327
    ## ev2         -2.083631   0.132227 -15.758
    ## con:ev1      0.077410   0.011451   6.760
    ## con:ev2     -0.084165   0.012344  -6.818
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.019                            
    ## ev1      0.143  0.040                     
    ## ev2      0.140 -0.035 -0.642              
    ## con:ev1  0.041  0.032  0.045 -0.012       
    ## con:ev2 -0.033  0.245 -0.011 -0.006 -0.645
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "40 MLM5 cn:30 ticc: 0.4 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 300.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.10982 -0.53423  0.02168  0.72035  1.98032 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             8.771    2.962   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.67061    0.42325   1.584
    ## con          0.36917    0.02815  13.116
    ## ev1          4.13597    0.61805   6.692
    ## ev2         -1.32981    0.61821  -2.151
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.230              
    ## ev1  0.128  0.095       
    ## ev2  0.172 -0.098 -0.658
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "40 MLM6 cn:30 ticc: 0.4 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 301.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.26575 -0.46262  0.07294  0.59056  2.07498 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             7.775    2.788   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.98386    0.41237   2.386
    ## con          0.32441    0.03036  10.685
    ## ev1          3.79103    0.59692   6.351
    ## ev2         -0.53733    0.63823  -0.842
    ## con:ev1      0.06806    0.03929   1.732
    ## con:ev2     -0.15234    0.05079  -2.999
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.314                            
    ## ev1      0.066  0.170                     
    ## ev2      0.255 -0.277 -0.665              
    ## con:ev1  0.190 -0.251 -0.204  0.248       
    ## con:ev2 -0.256  0.477  0.205 -0.408 -0.681
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "41 MLM5 cn:120 ticc: 0.4 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1197.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.29474 -0.75133 -0.03916  0.76476  2.18226 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.0     
    ##  Residual             8.409    2.9     
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.72660    0.20988   3.462
    ## con          0.28475    0.01507  18.897
    ## ev1          4.44994    0.32269  13.790
    ## ev2         -1.26869    0.32138  -3.948
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.001              
    ## ev1  0.177  0.196       
    ## ev2  0.178 -0.175 -0.691
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "41 MLM6 cn:120 ticc: 0.4 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1198
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.34425 -0.71776  0.00503  0.74121  2.21052 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             8.077    2.842   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.92641    0.21406   4.328
    ## con          0.29103    0.01619  17.974
    ## ev1          4.59650    0.32614  14.094
    ## ev2         -1.21330    0.32169  -3.772
    ## con:ev1      0.07624    0.02534   3.008
    ## con:ev2     -0.07445    0.02312  -3.221
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.049                            
    ## ev1      0.211  0.270                     
    ## ev2      0.172 -0.228 -0.690              
    ## con:ev1  0.263  0.288  0.214 -0.050       
    ## con:ev2 -0.240  0.027 -0.054 -0.114 -0.666
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "42 MLM5 cn:510 ticc: 0.4 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4932.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6494 -0.7229  0.0288  0.6732  3.4005 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             7.283    2.699   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.792066   0.089180   8.882
    ## con          0.292084   0.006394  45.679
    ## ev1          4.636584   0.131909  35.150
    ## ev2         -1.432135   0.131938 -10.855
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.069              
    ## ev1  0.133 -0.007       
    ## ev2  0.135  0.022 -0.635
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "42 MLM6 cn:510 ticc: 0.4 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4905.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6594 -0.6995  0.0209  0.6430  3.5247 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             7.003    2.646   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.785821   0.087536   8.977
    ## con          0.295981   0.006724  44.018
    ## ev1          4.698278   0.129685  36.229
    ## ev2         -1.495345   0.129958 -11.506
    ## con:ev1      0.064136   0.010093   6.355
    ## con:ev2     -0.052891   0.009966  -5.307
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.079                            
    ## ev1      0.132  0.000                     
    ## ev2      0.137  0.035 -0.636              
    ## con:ev1  0.000  0.169  0.071 -0.058       
    ## con:ev2  0.035  0.133 -0.059  0.095 -0.652
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "43 MLM5 cn:30 ticc: 0.8 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 289
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.64643 -0.56762  0.04066  0.61336  1.83881 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             7.003    2.646   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  2.16195    0.34932   6.189
    ## con          0.32144    0.01816  17.700
    ## ev1          7.60458    0.49871  15.249
    ## ev2          1.41652    0.50423   2.809
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.160              
    ## ev1  0.066  0.023       
    ## ev2  0.085  0.149 -0.555
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "43 MLM6 cn:30 ticc: 0.8 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 299.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.68057 -0.55968  0.00278  0.64142  1.79253 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             7.167    2.677   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  2.14844    0.37262   5.766
    ## con          0.31886    0.02178  14.641
    ## ev1          7.69444    0.52525  14.649
    ## ev2          1.29652    0.56406   2.299
    ## con:ev1      0.02312    0.02806   0.824
    ## con:ev2     -0.01588    0.03695  -0.430
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.297                            
    ## ev1     -0.009 -0.109                     
    ## ev2      0.193  0.328 -0.599              
    ## con:ev1 -0.119 -0.264  0.243 -0.328       
    ## con:ev2  0.293  0.518 -0.268  0.425 -0.704
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "44 MLM5 cn:120 ticc: 0.8 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1122.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.9194 -0.5723 -0.0155  0.7261  2.6639 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             6.078    2.465   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 2.175159   0.161829  13.441
    ## con         0.304494   0.007277  41.843
    ## ev1         6.962031   0.236340  29.458
    ## ev2         1.294039   0.236896   5.462
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.073              
    ## ev1  0.085 -0.124       
    ## ev2  0.065  0.142 -0.586
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "44 MLM6 cn:120 ticc: 0.8 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1133.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.92495 -0.62054  0.00722  0.67978  2.63603 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.000   
    ##  Residual             6.036    2.457   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  2.126377   0.163670  12.992
    ## con          0.300121   0.007752  38.716
    ## ev1          6.958786   0.239424  29.065
    ## ev2          1.260177   0.237323   5.310
    ## con:ev1      0.013567   0.010773   1.259
    ## con:ev2     -0.023260   0.012186  -1.909
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.034                            
    ## ev1      0.096 -0.147                     
    ## ev2      0.071  0.167 -0.586              
    ## con:ev1 -0.155 -0.049 -0.140  0.003       
    ## con:ev2  0.154  0.300  0.003  0.076 -0.642
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "45 MLM5 cn:510 ticc: 0.8 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4960.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.82484 -0.73364  0.01818  0.65807  2.84642 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.00     0.000   
    ##  Residual             7.48     2.735   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 2.036106   0.090673  22.456
    ## con         0.295629   0.003738  79.088
    ## ev1         7.130551   0.134566  52.989
    ## ev2         1.317812   0.134539   9.795
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.042              
    ## ev1  0.137  0.020       
    ## ev2  0.138  0.002 -0.640
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "45 MLM6 cn:510 ticc: 0.8 b2:0.3 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4964.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8338 -0.6938  0.0240  0.6692  2.8486 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.000    0.00    
    ##  Residual             7.399    2.72    
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  2.039958   0.090197  22.617
    ## con          0.293568   0.003899  75.293
    ## ev1          7.122475   0.133884  53.199
    ## ev2          1.335890   0.133944   9.973
    ## con:ev1      0.015512   0.005629   2.756
    ## con:ev2     -0.020942   0.005903  -3.548
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.039                            
    ## ev1      0.137  0.016                     
    ## ev2      0.138 -0.010 -0.639              
    ## con:ev1  0.017  0.058 -0.025  0.013       
    ## con:ev2 -0.009  0.193  0.013 -0.041 -0.630
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "46 MLM5 cn:30 ticc: 0.2 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 328.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.2393 -0.6231 -0.1245  0.7496  2.1907 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.00    
    ##  Residual             14.52    3.81    
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)   0.7505     0.5313   1.413
    ## con           0.2952     0.0405   7.289
    ## ev1           6.4903     0.7920   8.195
    ## ev2          -3.2781     0.7917  -4.140
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.072              
    ## ev1  0.157 -0.036       
    ## ev2  0.156 -0.024 -0.654
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "46 MLM6 cn:30 ticc: 0.2 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 334.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.23879 -0.54476 -0.03637  0.55028  2.17171 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.0     0.000   
    ##  Residual             14.6     3.821   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.71462    0.53599   1.333
    ## con          0.31548    0.04734   6.665
    ## ev1          6.33910    0.80355   7.889
    ## ev2         -3.16337    0.79910  -3.959
    ## con:ev1      0.09338    0.07870   1.187
    ## con:ev2     -0.07572    0.06332  -1.196
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.116                            
    ## ev1      0.165 -0.094                     
    ## ev2      0.149  0.011 -0.658              
    ## con:ev1 -0.085  0.460 -0.149  0.097       
    ## con:ev2  0.012 -0.157  0.120 -0.108 -0.699
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "47 MLM5 cn:120 ticc: 0.2 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1358.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.96973 -0.71940  0.02792  0.68999  2.33081 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             16.64    4.079   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.62147    0.28149   2.208
    ## con          0.28236    0.02512  11.242
    ## ev1          6.82192    0.41754  16.338
    ## ev2         -3.31929    0.41750  -7.950
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.085              
    ## ev1  0.146 -0.014       
    ## ev2  0.144 -0.001 -0.646
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "47 MLM6 cn:120 ticc: 0.2 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1361.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.00494 -0.69123  0.00184  0.63185  2.38722 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             16.27    4.034   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.61938    0.27863   2.223
    ## con          0.27616    0.02659  10.386
    ## ev1          6.71383    0.41491  16.181
    ## ev2         -3.22314    0.41485  -7.770
    ## con:ev1      0.09900    0.03760   2.633
    ## con:ev2     -0.09098    0.04141  -2.197
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.094                            
    ## ev1      0.146 -0.013                     
    ## ev2      0.146 -0.025 -0.647              
    ## con:ev1 -0.013 -0.001 -0.099  0.071       
    ## con:ev2 -0.024  0.273  0.065 -0.097 -0.647
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "48 MLM5 cn:510 ticc: 0.2 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5700.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.62619 -0.73663 -0.01608  0.74624  2.30321 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.00    
    ##  Residual             15.52    3.94    
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.78468    0.13060   6.008
    ## con          0.27900    0.01097  25.444
    ## ev1          6.50982    0.19386  33.580
    ## ev2         -3.27419    0.19384 -16.891
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.038              
    ## ev1  0.139  0.024       
    ## ev2  0.139  0.019 -0.639
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "48 MLM6 cn:510 ticc: 0.2 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5688.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.62843 -0.69100 -0.00479  0.66657  2.36438 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             15.18    3.896   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.79163    0.12932   6.121
    ## con          0.28286    0.01142  24.774
    ## ev1          6.57654    0.19215  34.226
    ## ev2         -3.33511    0.19206 -17.365
    ## con:ev1      0.08055    0.01715   4.696
    ## con:ev2     -0.07267    0.01669  -4.354
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.053                            
    ## ev1      0.140  0.029                     
    ## ev2      0.139  0.017 -0.641              
    ## con:ev1  0.029  0.171  0.067 -0.055       
    ## con:ev2  0.017  0.094 -0.056  0.061 -0.635
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "49 MLM5 cn:30 ticc: 0.4 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 315.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.93134 -0.50184  0.03005  0.43228  2.57490 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             11.32    3.365   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.23607    0.43451   2.845
    ## con          0.31041    0.03524   8.808
    ## ev1          7.27845    0.61877  11.763
    ## ev2         -2.24442    0.66097  -3.396
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.025              
    ## ev1 -0.003  0.120       
    ## ev2  0.009 -0.369 -0.506
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "49 MLM6 cn:30 ticc: 0.4 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 322.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.94092 -0.36434 -0.02757  0.36822  2.49134 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             11.43    3.381   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.36426    0.48389   2.819
    ## con          0.29927    0.03708   8.071
    ## ev1          7.23447    0.65404  11.061
    ## ev2         -1.96440    0.72705  -2.702
    ## con:ev1      0.05807    0.04834   1.201
    ## con:ev2     -0.04525    0.05796  -0.781
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.136                            
    ## ev1     -0.128  0.172                     
    ## ev2      0.171 -0.442 -0.534              
    ## con:ev1  0.178 -0.230 -0.023  0.290       
    ## con:ev2 -0.425  0.284  0.269 -0.399 -0.561
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "50 MLM5 cn:120 ticc: 0.4 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1327.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.22995 -0.68312 -0.00264  0.76545  2.34264 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             14.57    3.818   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.49785    0.25510   5.872
    ## con          0.32295    0.01913  16.885
    ## ev1          7.54865    0.37576  20.089
    ## ev2         -2.08346    0.37625  -5.537
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.012              
    ## ev1  0.114  0.039       
    ## ev2  0.114  0.064 -0.611
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "50 MLM6 cn:120 ticc: 0.4 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1332.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.36296 -0.68659  0.00038  0.76237  2.30313 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.0     0.000   
    ##  Residual             14.3     3.781   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.45020    0.25419   5.705
    ## con          0.30937    0.02081  14.867
    ## ev1          7.62758    0.37351  20.421
    ## ev2         -2.17097    0.37450  -5.797
    ## con:ev1      0.04711    0.03059   1.540
    ## con:ev2     -0.08150    0.03185  -2.559
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.056                            
    ## ev1      0.108  0.022                     
    ## ev2      0.116  0.068 -0.614              
    ## con:ev1  0.022  0.109  0.067 -0.087       
    ## con:ev2  0.065  0.224 -0.084  0.094 -0.668
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "51 MLM5 cn:510 ticc: 0.4 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5766.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.56043 -0.72171 -0.00326  0.72730  2.37994 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             16.56    4.069   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  1.475538   0.136527   10.81
    ## con          0.308663   0.009586   32.20
    ## ev1          7.709606   0.203595   37.87
    ## ev2         -2.190614   0.203752  -10.75
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.002              
    ## ev1  0.150  0.005       
    ## ev2  0.150  0.040 -0.650
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "51 MLM6 cn:510 ticc: 0.4 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5766.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.66531 -0.69649  0.00578  0.71754  2.36595 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             16.38    4.047   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.45559    0.13597  10.706
    ## con          0.30167    0.01073  28.116
    ## ev1          7.73497    0.20260  38.179
    ## ev2         -2.22867    0.20288 -10.985
    ## con:ev1      0.04512    0.01589   2.839
    ## con:ev2     -0.05989    0.01647  -3.635
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.024                            
    ## ev1      0.148 -0.007                     
    ## ev2      0.152  0.046 -0.651              
    ## con:ev1 -0.007  0.131  0.017 -0.037       
    ## con:ev2  0.045  0.233 -0.036  0.051 -0.682
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "52 MLM5 cn:30 ticc: 0.8 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 333.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.08063 -0.65997  0.06103  0.54349  2.17009 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             15.59    3.948   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  2.99481    0.52100   5.748
    ## con          0.30346    0.02881  10.533
    ## ev1         11.79466    0.74394  15.854
    ## ev2          2.09462    0.75624   2.770
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.158              
    ## ev1  0.064 -0.013       
    ## ev2  0.090 -0.180 -0.554
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "52 MLM6 cn:30 ticc: 0.8 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 338.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.10270 -0.54574 -0.01482  0.56164  2.31534 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             14.92    3.863   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  2.89014    0.52177   5.539
    ## con          0.30122    0.02821  10.678
    ## ev1         11.65792    0.74592  15.629
    ## ev2          2.28340    0.77039   2.964
    ## con:ev1      0.07754    0.03950   1.963
    ## con:ev2     -0.00884    0.03971  -0.223
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.148                            
    ## ev1      0.031 -0.013                     
    ## ev2      0.122 -0.173 -0.580              
    ## con:ev1 -0.014 -0.028 -0.157  0.204       
    ## con:ev2 -0.181 -0.013  0.209 -0.265 -0.478
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "53 MLM5 cn:120 ticc: 0.8 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1350.5
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.23062 -0.67895  0.03672  0.72312  2.45529 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             15.97    3.996   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  3.48161    0.27288  12.759
    ## con          0.29255    0.01154  25.350
    ## ev1         11.82853    0.40412  29.270
    ## ev2          1.78427    0.40418   4.415
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.085              
    ## ev1  0.140 -0.065       
    ## ev2  0.128  0.067 -0.638
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "53 MLM6 cn:120 ticc: 0.8 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1361.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.24541 -0.65666  0.02705  0.71625  2.46181 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             15.94    3.993   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  3.44057    0.27395  12.559
    ## con          0.29532    0.01179  25.047
    ## ev1         11.74524    0.40850  28.752
    ## ev2          1.81960    0.40451   4.498
    ## con:ev1      0.02631    0.01781   1.477
    ## con:ev2     -0.02023    0.01647  -1.228
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.097                            
    ## ev1      0.150 -0.093                     
    ## ev2      0.122  0.075 -0.638              
    ## con:ev1 -0.092  0.187 -0.148  0.056       
    ## con:ev2  0.079 -0.034  0.060 -0.041 -0.585
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "54 MLM5 cn:510 ticc: 0.8 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5687.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.63006 -0.66928 -0.02739  0.75602  2.21290 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             15.31    3.912   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  3.459905   0.128267  26.974
    ## con          0.295009   0.005452  54.106
    ## ev1         11.868938   0.189431  62.656
    ## ev2          1.861466   0.189456   9.825
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.050              
    ## ev1  0.126  0.003       
    ## ev2  0.127  0.017 -0.628
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "54 MLM6 cn:510 ticc: 0.8 b2:0.5 b3:0.1"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5698.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.62672 -0.66629 -0.01218  0.72781  2.26990 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             15.26    3.906   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  3.462681   0.128142  27.022
    ## con          0.297131   0.005681  52.300
    ## ev1         11.893820   0.189468  62.775
    ## ev2          1.842295   0.189592   9.717
    ## con:ev1      0.019382   0.008494   2.282
    ## con:ev2     -0.011435   0.008289  -1.380
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.058                            
    ## ev1      0.126  0.009                     
    ## ev2      0.128  0.021 -0.629              
    ## con:ev1  0.009  0.157  0.058 -0.046       
    ## con:ev2  0.021  0.088 -0.047  0.067 -0.625
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "55 MLM5 cn:30 ticc: 0.2 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 303.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.1620 -0.6468  0.1089  0.4790  2.3843 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.9818   0.9909  
    ##  Residual             8.4578   2.9082  
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.54772    0.47462   3.261
    ## con          0.25880    0.03886   6.659
    ## ev1          2.90277    0.80756   3.594
    ## ev2         -0.36525    0.78416  -0.466
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.023              
    ## ev1  0.157 -0.459       
    ## ev2  0.182  0.404 -0.744
    ## [1] "55 MLM6 cn:30 ticc: 0.2 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 300
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.3650 -0.5096  0.1097  0.5501  2.5358 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 1.422    1.193   
    ##  Residual             6.683    2.585   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.30955    0.58952   0.525
    ## con          0.28096    0.04445   6.321
    ## ev1          1.99723    0.95251   2.097
    ## ev2         -0.74278    0.86039  -0.863
    ## con:ev1      0.20922    0.07357   2.844
    ## con:ev2     -0.21393    0.06148  -3.480
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.197                            
    ## ev1      0.378 -0.630                     
    ## ev2      0.089  0.533 -0.738              
    ## con:ev1 -0.615  0.447 -0.528  0.181       
    ## con:ev2  0.562 -0.063  0.196  0.190 -0.722
    ## [1] "56 MLM5 cn:120 ticc: 0.2 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1273.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.4164 -0.6525  0.1190  0.7417  2.3181 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.03603 0.1898  
    ##  Residual             11.53689 3.3966  
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.14986    0.22593   0.663
    ## con          0.29336    0.01863  15.749
    ## ev1          1.55398    0.33279   4.670
    ## ev2         -0.78502    0.33223  -2.363
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.016              
    ## ev1  0.100  0.099       
    ## ev2  0.103 -0.081 -0.607
    ## [1] "56 MLM6 cn:120 ticc: 0.2 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1184.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.87364 -0.49136  0.04604  0.62285  2.86283 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.4382   0.662   
    ##  Residual             7.1877   2.681   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.40021    0.18963   2.110
    ## con          0.27304    0.01579  17.288
    ## ev1          1.62894    0.27802   5.859
    ## ev2         -0.61442    0.27851  -2.206
    ## con:ev1      0.22648    0.02221  10.198
    ## con:ev2     -0.23769    0.02407  -9.876
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.027                            
    ## ev1      0.102  0.106                     
    ## ev2      0.107 -0.104 -0.607              
    ## con:ev1  0.111 -0.016  0.049  0.012       
    ## con:ev2 -0.100  0.211  0.011 -0.092 -0.606
    ## [1] "57 MLM5 cn:510 ticc: 0.2 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5366.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.2228 -0.6894 -0.0053  0.6502  3.6050 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.444   0.6664  
    ##  Residual             10.741   3.2774  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.301940   0.111673   2.704
    ## con          0.289114   0.009275  31.171
    ## ev1          1.664973   0.165174  10.080
    ## ev2         -0.765259   0.165133  -4.634
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.012              
    ## ev1  0.126 -0.025       
    ## ev2  0.126  0.012 -0.628
    ## [1] "57 MLM6 cn:510 ticc: 0.2 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5009.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8684 -0.6051  0.0138  0.5368  3.7786 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 1.330    1.153   
    ##  Residual             6.553    2.560   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.252406   0.099438   2.538
    ## con          0.285672   0.008009  35.670
    ## ev1          1.582382   0.147093  10.758
    ## ev2         -0.732330   0.147005  -4.982
    ## con:ev1      0.222409   0.011677  19.046
    ## con:ev2     -0.227444   0.011766 -19.330
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.014                            
    ## ev1      0.127 -0.027                     
    ## ev2      0.126  0.014 -0.628              
    ## con:ev1 -0.028  0.086 -0.032  0.016       
    ## con:ev2  0.014  0.108  0.015 -0.004 -0.600
    ## [1] "58 MLM5 cn:30 ticc: 0.4 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 329.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.15547 -0.55689 -0.00978  0.76642  1.95736 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  1.146   1.071   
    ##  Residual             13.562   3.683   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.83475    0.54038   1.545
    ## con          0.29555    0.03182   9.287
    ## ev1          1.79899    0.78974   2.278
    ## ev2         -0.06594    0.79111  -0.083
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.175              
    ## ev1  0.088  0.128       
    ## ev2  0.086  0.141 -0.586
    ## [1] "58 MLM6 cn:30 ticc: 0.4 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 319
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.49406 -0.30781  0.01376  0.50485  2.30737 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 1.329    1.153   
    ##  Residual             9.364    3.060   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.81141    0.48449   1.675
    ## con          0.26561    0.02820   9.419
    ## ev1          1.85168    0.69658   2.658
    ## ev2         -0.16161    0.69678  -0.232
    ## con:ev1      0.16213    0.03824   4.239
    ## con:ev2     -0.19003    0.04327  -4.392
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.133                            
    ## ev1      0.047  0.100                     
    ## ev2      0.048  0.122 -0.549              
    ## con:ev1  0.106 -0.118 -0.063 -0.096       
    ## con:ev2  0.115  0.231 -0.085 -0.041 -0.573
    ## [1] "59 MLM5 cn:120 ticc: 0.4 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1271.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.7643 -0.6154  0.0085  0.6085  2.5161 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.4769  0.6906  
    ##  Residual             11.0095  3.3181  
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.68229    0.23575   2.894
    ## con          0.30362    0.01424  21.329
    ## ev1          2.32966    0.35009   6.655
    ## ev2         -0.56040    0.34999  -1.601
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.059              
    ## ev1  0.139 -0.076       
    ## ev2  0.130  0.072 -0.638
    ## [1] "59 MLM6 cn:120 ticc: 0.4 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1238.5
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -3.14868 -0.52941 -0.00364  0.52655  2.79682 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 2.517    1.587   
    ##  Residual             7.400    2.720   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.47179    0.24198   1.950
    ## con          0.31239    0.01594  19.602
    ## ev1          2.01448    0.36098   5.581
    ## ev2         -0.46714    0.35757  -1.306
    ## con:ev1      0.17576    0.02414   7.280
    ## con:ev2     -0.15407    0.02430  -6.339
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.065                            
    ## ev1      0.151 -0.119                     
    ## ev2      0.124  0.092 -0.639              
    ## con:ev1 -0.118  0.195 -0.136  0.048       
    ## con:ev2  0.089  0.214  0.047  0.002 -0.702
    ## [1] "60 MLM5 cn:510 ticc: 0.4 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5397
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8859 -0.6668 -0.0110  0.6838  3.2434 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.1715  0.4141  
    ##  Residual             11.3352  3.3668  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.196381   0.112491   1.746
    ## con          0.317635   0.007185  44.210
    ## ev1          1.708097   0.166884  10.235
    ## ev2         -0.403233   0.166988  -2.415
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.041              
    ## ev1  0.128  0.068       
    ## ev2  0.134 -0.077 -0.635
    ## [1] "60 MLM6 cn:510 ticc: 0.4 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5225.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6804 -0.6122 -0.0005  0.6052  3.5206 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 1.568    1.252   
    ##  Residual             8.159    2.856   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.356592   0.111166   3.208
    ## con          0.316717   0.007178  44.126
    ## ev1          1.677180   0.164333  10.206
    ## ev2         -0.209285   0.165251  -1.266
    ## con:ev1      0.138134   0.010676  12.939
    ## con:ev2     -0.145494   0.010617 -13.704
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.043                            
    ## ev1      0.125  0.079                     
    ## ev2      0.141 -0.096 -0.635              
    ## con:ev1  0.078  0.143  0.014  0.031       
    ## con:ev2 -0.096  0.127  0.032 -0.104 -0.636
    ## [1] "61 MLM5 cn:30 ticc: 0.8 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 314.5
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.34225 -0.73920  0.09108  0.66419  2.13900 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 1.865    1.366   
    ##  Residual             9.404    3.067   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.39116    0.54509   2.552
    ## con          0.29379    0.01778  16.524
    ## ev1          2.65081    0.81431   3.255
    ## ev2          0.81085    0.83699   0.969
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.212              
    ## ev1  0.214 -0.148       
    ## ev2  0.120  0.272 -0.695
    ## [1] "61 MLM6 cn:30 ticc: 0.8 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 320
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.31181 -0.68543  0.03075  0.66359  2.18737 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 2.624    1.620   
    ##  Residual             8.240    2.871   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.87574    0.59582   1.470
    ## con          0.30789    0.02177  14.141
    ## ev1          1.92194    0.93637   2.053
    ## ev2          1.00893    0.86383   1.168
    ## con:ev1      0.07965    0.03618   2.201
    ## con:ev2     -0.06433    0.02981  -2.158
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.273                            
    ## ev1      0.299 -0.370                     
    ## ev2      0.070  0.328 -0.690              
    ## con:ev1 -0.350  0.458 -0.432  0.158       
    ## con:ev2  0.347 -0.092  0.176 -0.036 -0.719
    ## [1] "62 MLM5 cn:120 ticc: 0.8 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1268.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.20494 -0.71867  0.06745  0.62268  2.63717 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.603   0.7765  
    ##  Residual             10.692   3.2698  
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.09397    0.25549   4.282
    ## con          0.30861    0.00703  43.900
    ## ev1          2.86601    0.38326   7.478
    ## ev2          0.80471    0.38346   2.099
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.125              
    ## ev1  0.188  0.005       
    ## ev2  0.183 -0.033 -0.687
    ## [1] "62 MLM6 cn:120 ticc: 0.8 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1276.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.44002 -0.66264  0.04874  0.59239  2.41882 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  1.113   1.055   
    ##  Residual             10.015   3.165   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  1.108925   0.259362   4.276
    ## con          0.307439   0.007937  38.735
    ## ev1          2.980313   0.392461   7.594
    ## ev2          0.706182   0.391041   1.806
    ## con:ev1      0.028124   0.011789   2.386
    ## con:ev2     -0.029164   0.012114  -2.407
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.125                            
    ## ev1      0.191  0.034                     
    ## ev2      0.181 -0.009 -0.685              
    ## con:ev1  0.035  0.139  0.134 -0.073       
    ## con:ev2 -0.008  0.216 -0.071  0.103 -0.677
    ## [1] "63 MLM5 cn:510 ticc: 0.8 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5395.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3885 -0.6425 -0.0039  0.6501  3.0825 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.3009  0.5485  
    ##  Residual             11.1816  3.3439  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 0.971128   0.116570   8.331
    ## con         0.298622   0.003768  79.242
    ## ev1         2.837908   0.174233  16.288
    ## ev2         1.015605   0.174262   5.828
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.044              
    ## ev1  0.158 -0.018       
    ## ev2  0.160  0.026 -0.659
    ## [1] "63 MLM6 cn:510 ticc: 0.8 b2:0.1 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5373.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3552 -0.6112 -0.0138  0.6369  2.7782 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.8735  0.9346  
    ##  Residual             10.2245  3.1976  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.943110   0.117625   8.018
    ## con          0.296606   0.004189  70.811
    ## ev1          2.884157   0.175768  16.409
    ## ev2          0.942134   0.176197   5.347
    ## con:ev1      0.036566   0.006135   5.960
    ## con:ev2     -0.037802   0.006451  -5.860
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.055                            
    ## ev1      0.156 -0.022                     
    ## ev2      0.163  0.047 -0.660              
    ## con:ev1 -0.023  0.099  0.035 -0.042       
    ## con:ev2  0.045  0.242 -0.040  0.078 -0.673

    ## boundary (singular) fit: see help('isSingular')

    ## [1] "64 MLM5 cn:30 ticc: 0.2 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 328.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.24509 -0.52083 -0.09198  0.65882  2.74659 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             14.49    3.806   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.71262    0.49967   1.426
    ## con          0.27727    0.04038   6.867
    ## ev1          5.39199    0.71903   7.499
    ## ev2         -2.30464    0.71749  -3.212
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.123              
    ## ev1  0.071 -0.074       
    ## ev2  0.067 -0.036 -0.561
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "64 MLM6 cn:30 ticc: 0.2 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 319.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5625 -0.4832 -0.1631  0.5960  3.1524 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.1154  0.3396  
    ##  Residual             10.8241  3.2900  
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.61527    0.44024   1.398
    ## con          0.29699    0.03565   8.331
    ## ev1          4.82002    0.64135   7.515
    ## ev2         -1.82974    0.63562  -2.879
    ## con:ev1      0.21537    0.05338   4.035
    ## con:ev2     -0.19499    0.04951  -3.939
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.140                            
    ## ev1      0.084 -0.101                     
    ## ev2      0.059 -0.015 -0.574              
    ## con:ev1 -0.098  0.162 -0.195  0.143       
    ## con:ev2 -0.016 -0.051  0.153 -0.151 -0.563

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "65 MLM5 cn:120 ticc: 0.2 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1281.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.73443 -0.61363 -0.01555  0.69211  2.55504 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             12.01    3.465   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.81626    0.23014   3.547
    ## con          0.32522    0.02032  16.003
    ## ev1          4.18679    0.33688  12.428
    ## ev2         -1.63536    0.33691  -4.854
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.059              
    ## ev1  0.100 -0.023       
    ## ev2  0.100 -0.027 -0.603
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "65 MLM6 cn:120 ticc: 0.2 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1210.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1459 -0.5625 -0.0043  0.5625  3.0552 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.00     0.000   
    ##  Residual             8.51     2.917   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.80484    0.19392   4.150
    ## con          0.30376    0.01745  17.411
    ## ev1          4.25454    0.28391  14.985
    ## ev2         -1.67583    0.28391  -5.903
    ## con:ev1      0.22340    0.02416   9.245
    ## con:ev2     -0.21910    0.02635  -8.314
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.052                            
    ## ev1      0.098 -0.020                     
    ## ev2      0.098 -0.019 -0.600              
    ## con:ev1 -0.021 -0.059  0.037  0.002       
    ## con:ev2 -0.018  0.187  0.002  0.035 -0.574
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "66 MLM5 cn:510 ticc: 0.2 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5716.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.0150 -0.6692 -0.0020  0.6572  3.7320 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             15.77    3.971   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.52291    0.13273   3.940
    ## con          0.28733    0.01081  26.590
    ## ev1          3.89341    0.19745  19.718
    ## ev2         -1.59631    0.19736  -8.088
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.054              
    ## ev1  0.147  0.036       
    ## ev2  0.144 -0.018 -0.647
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "66 MLM6 cn:510 ticc: 0.2 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5457.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3826 -0.5903 -0.0297  0.5414  4.4722 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             12.08    3.475   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.60782    0.11630   5.226
    ## con          0.28867    0.01016  28.403
    ## ev1          4.14383    0.17349  23.885
    ## ev2         -1.75501    0.17292 -10.149
    ## con:ev1      0.25279    0.01498  16.873
    ## con:ev2     -0.23108    0.01535 -15.053
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.062                            
    ## ev1      0.151  0.049                     
    ## ev2      0.142 -0.018 -0.647              
    ## con:ev1  0.049  0.117  0.090 -0.050       
    ## con:ev2 -0.017  0.186 -0.048  0.044 -0.653
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "67 MLM5 cn:30 ticc: 0.4 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 329.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.34324 -0.76518  0.07525  0.82588  1.89819 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             14.61    3.822   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.35702    0.53187   2.551
    ## con          0.31499    0.03864   8.153
    ## ev1          5.15784    0.80350   6.419
    ## ev2         -1.00136    0.80551  -1.243
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.033              
    ## ev1  0.147 -0.154       
    ## ev2  0.158  0.169 -0.664
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "67 MLM6 cn:30 ticc: 0.4 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 328.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.22601 -0.71140 -0.04372  0.69374  2.26649 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             12.99    3.604   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.07495    0.52118   2.063
    ## con          0.33183    0.04297   7.723
    ## ev1          4.92131    0.77082   6.385
    ## ev2         -1.03895    0.79291  -1.310
    ## con:ev1      0.17253    0.06020   2.866
    ## con:ev2     -0.09724    0.06983  -1.393
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.108                            
    ## ev1      0.127 -0.220                     
    ## ev2      0.207  0.282 -0.668              
    ## con:ev1 -0.232 -0.026 -0.053 -0.100       
    ## con:ev2  0.264  0.394 -0.088  0.261 -0.704
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "68 MLM5 cn:120 ticc: 0.4 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1297.3
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.59269 -0.69250 -0.01295  0.60650  2.45385 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.0     0.000   
    ##  Residual             12.8     3.578   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.50311    0.23794   6.317
    ## con          0.33279    0.01761  18.896
    ## ev1          4.96204    0.35421  14.009
    ## ev2         -0.72906    0.36409  -2.002
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.078              
    ## ev1  0.085 -0.190       
    ## ev2  0.120  0.296 -0.623
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "68 MLM6 cn:120 ticc: 0.4 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1283.5
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.42503 -0.59596  0.00763  0.57038  2.55078 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.0     0.000   
    ##  Residual             11.6     3.406   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.06389    0.24312   4.376
    ## con          0.31530    0.01747  18.043
    ## ev1          5.15719    0.34619  14.897
    ## ev2         -1.35375    0.37357  -3.624
    ## con:ev1      0.10547    0.02424   4.352
    ## con:ev2     -0.13073    0.02702  -4.838
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.157                            
    ## ev1      0.019 -0.240                     
    ## ev2      0.235  0.365 -0.635              
    ## con:ev1 -0.246 -0.055 -0.014 -0.176       
    ## con:ev2  0.362  0.253 -0.171  0.368 -0.613
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "69 MLM5 cn:510 ticc: 0.4 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5767.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1014 -0.6684 -0.0012  0.7039  3.4907 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             16.57    4.071   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  1.107655   0.138522   7.996
    ## con          0.304648   0.008788  34.668
    ## ev1          4.920896   0.207374  23.730
    ## ev2         -0.969857   0.207368  -4.677
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.005              
    ## ev1  0.161  0.008       
    ## ev2  0.161 -0.001 -0.661
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "69 MLM6 cn:510 ticc: 0.4 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5636.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.2213 -0.6119  0.0225  0.6318  3.6315 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.0     0.000   
    ##  Residual             14.4     3.794   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  1.117361   0.129108   8.654
    ## con          0.302125   0.008653  34.917
    ## ev1          4.925679   0.193276  25.485
    ## ev2         -0.962198   0.193276  -4.978
    ## con:ev1      0.136855   0.012900  10.609
    ## con:ev2     -0.151227   0.012788 -11.826
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.003                            
    ## ev1      0.161  0.008                     
    ## ev2      0.161 -0.004 -0.661              
    ## con:ev1  0.008  0.149  0.003 -0.001       
    ## con:ev2 -0.004  0.125 -0.001 -0.005 -0.639
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "70 MLM5 cn:30 ticc: 0.8 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 337.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.16085 -0.66608 -0.02668  0.66363  2.08141 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             16.78    4.096   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)   2.0004     0.6634   3.015
    ## con           0.3128     0.0205  15.256
    ## ev1           7.9440     1.0056   7.900
    ## ev2           0.9598     0.9957   0.964
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.217              
    ## ev1  0.248  0.180       
    ## ev2  0.187 -0.114 -0.719
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "70 MLM6 cn:30 ticc: 0.8 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 347.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.20116 -0.68364 -0.06908  0.65626  1.98524 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             17.25    4.154   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  1.886083   0.698382   2.701
    ## con          0.306788   0.023367  13.129
    ## ev1          7.683289   1.107676   6.936
    ## ev2          1.066771   1.031752   1.034
    ## con:ev1     -0.017682   0.034670  -0.510
    ## con:ev2      0.001073   0.035798   0.030
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.253                            
    ## ev1      0.325  0.241                     
    ## ev2      0.124 -0.126 -0.725              
    ## con:ev1  0.258  0.136  0.378 -0.205       
    ## con:ev2 -0.122  0.226 -0.185  0.141 -0.681
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "71 MLM5 cn:120 ticc: 0.8 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1331.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -3.12693 -0.60760 -0.00201  0.72291  2.57595 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             14.71    3.835   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  2.49800    0.26558   9.406
    ## con          0.29869    0.00988  30.234
    ## ev1          7.77300    0.37723  20.605
    ## ev2          2.01688    0.37720   5.347
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.263              
    ## ev1  0.105 -0.016       
    ## ev2  0.111  0.008 -0.615
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "71 MLM6 cn:120 ticc: 0.8 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1325.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.2465 -0.6142  0.0666  0.6060  2.6768 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             13.64    3.693   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  2.504528   0.256031   9.782
    ## con          0.302139   0.009753  30.978
    ## ev1          8.175775   0.376174  21.734
    ## ev2          1.625187   0.376759   4.314
    ## con:ev1      0.059800   0.014588   4.099
    ## con:ev2     -0.055735   0.013932  -4.001
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.266                            
    ## ev1      0.108  0.026                     
    ## ev2      0.112  0.017 -0.612              
    ## con:ev1  0.026  0.159  0.260 -0.150       
    ## con:ev2  0.018  0.028 -0.157  0.266 -0.598
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "72 MLM5 cn:510 ticc: 0.8 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5699.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -3.00562 -0.72227  0.00301  0.68144  2.84216 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             15.48    3.935   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 2.375388   0.131715   18.03
    ## con         0.305054   0.004865   62.70
    ## ev1         7.620252   0.196224   38.83
    ## ev2         2.179060   0.196544   11.09
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.029              
    ## ev1  0.147 -0.022       
    ## ev2  0.149  0.061 -0.649
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "72 MLM6 cn:510 ticc: 0.8 b2:0.3 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5680.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.61465 -0.71586  0.02471  0.68972  2.82310 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             14.98    3.871   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  2.326109   0.129911  17.905
    ## con          0.304971   0.005002  60.968
    ## ev1          7.665762   0.193255  39.667
    ## ev2          2.089859   0.194008  10.772
    ## con:ev1      0.035050   0.007540   4.648
    ## con:ev2     -0.042300   0.007237  -5.845
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.038                            
    ## ev1      0.143 -0.029                     
    ## ev2      0.154  0.067 -0.650              
    ## con:ev1 -0.029  0.181  0.014 -0.042       
    ## con:ev2  0.069  0.064 -0.044  0.081 -0.626
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "73 MLM5 cn:30 ticc: 0.2 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 344.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.31892 -0.68968 -0.01921  0.64564  2.54642 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             18.99    4.357   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.14109    0.57444   1.986
    ## con          0.28472    0.04462   6.381
    ## ev1          6.33265    0.82181   7.706
    ## ev2         -3.08980    0.82095  -3.764
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.152              
    ## ev1  0.069 -0.046       
    ## ev2  0.063 -0.003 -0.564
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "73 MLM6 cn:30 ticc: 0.2 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 341.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.54719 -0.45438 -0.02299  0.43249  2.59941 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             16.46    4.058   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.04359    0.53939   1.935
    ## con          0.30741    0.04566   6.732
    ## ev1          5.80295    0.78800   7.364
    ## ev2         -2.63578    0.77712  -3.392
    ## con:ev1      0.20081    0.07263   2.765
    ## con:ev2     -0.19887    0.06373  -3.120
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.191                            
    ## ev1      0.092 -0.118                     
    ## ev2      0.053  0.020 -0.575              
    ## con:ev1 -0.109  0.333 -0.239  0.146       
    ## con:ev2  0.020 -0.037  0.164 -0.176 -0.666
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "74 MLM5 cn:120 ticc: 0.2 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1396.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.44228 -0.78005  0.01205  0.65020  3.05592 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             19.56    4.422   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.99817    0.29321   3.404
    ## con          0.29215    0.02675  10.923
    ## ev1          7.48738    0.43557  17.190
    ## ev2         -3.29529    0.43007  -7.662
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.005              
    ## ev1  0.099 -0.162       
    ## ev2  0.102  0.034 -0.601
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "74 MLM6 cn:120 ticc: 0.2 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1365.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5142 -0.6027  0.0412  0.5831  3.4938 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             16.53    4.065   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.72863    0.27378   2.661
    ## con          0.28194    0.02611  10.797
    ## ev1          7.15691    0.40393  17.718
    ## ev2         -3.16526    0.39846  -7.944
    ## con:ev1      0.24268    0.03731   6.505
    ## con:ev2     -0.22068    0.04019  -5.491
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.019                            
    ## ev1      0.120 -0.161                     
    ## ev2      0.081  0.067 -0.603              
    ## con:ev1 -0.166  0.029 -0.131  0.076       
    ## con:ev2  0.063  0.239  0.070  0.027 -0.641
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "75 MLM5 cn:510 ticc: 0.2 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6176.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.68118 -0.69046  0.01133  0.71850  2.77524 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.0     0.00    
    ##  Residual             24.8     4.98    
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.88855    0.16663   5.332
    ## con          0.30192    0.01365  22.123
    ## ev1          6.74225    0.24837  27.146
    ## ev2         -3.24546    0.24828 -13.072
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.006              
    ## ev1  0.148  0.027       
    ## ev2  0.148  0.004 -0.648
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "75 MLM6 cn:510 ticc: 0.2 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6030.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.91372 -0.50979 -0.00929  0.58752  2.97883 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             21.26    4.611   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.91659    0.15439   5.937
    ## con          0.29405    0.01320  22.270
    ## ev1          6.83818    0.23010  29.718
    ## ev2         -3.31114    0.22997 -14.398
    ## con:ev1      0.22392    0.01929  11.611
    ## con:ev2     -0.23623    0.01974 -11.967
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.015                            
    ## ev1      0.148  0.027                     
    ## ev2      0.147  0.001 -0.649              
    ## con:ev1  0.028  0.091  0.032 -0.026       
    ## con:ev2  0.001  0.157 -0.026  0.014 -0.627
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "76 MLM5 cn:30 ticc: 0.4 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 354.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.4528 -0.6592 -0.1117  0.6638  1.9681 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             23.13    4.809   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.95921    0.63228   1.517
    ## con          0.26016    0.05862   4.438
    ## ev1          8.26890    0.92570   8.933
    ## ev2         -1.75571    0.90660  -1.937
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.134              
    ## ev1  0.034 -0.205       
    ## ev2  0.058 -0.035 -0.546
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "76 MLM6 cn:30 ticc: 0.4 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 350.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.57820 -0.50070 -0.08381  0.35713  2.08080 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             19.75    4.444   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.67303    0.60007   1.122
    ## con          0.26723    0.05549   4.816
    ## ev1          8.06882    0.85861   9.398
    ## ev2         -1.66503    0.85367  -1.950
    ## con:ev1      0.26412    0.08126   3.251
    ## con:ev2     -0.22318    0.08115  -2.750
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.084                            
    ## ev1      0.033 -0.191                     
    ## ev2      0.017  0.008 -0.526              
    ## con:ev1 -0.187  0.099 -0.051  0.086       
    ## con:ev2  0.008  0.095  0.085  0.086 -0.599
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "77 MLM5 cn:120 ticc: 0.4 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1445.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -3.14825 -0.63208  0.04048  0.64332  2.44210 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.0     
    ##  Residual             24.01    4.9     
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.56694    0.33224   4.716
    ## con          0.35297    0.02455  14.378
    ## ev1          8.17870    0.48807  16.757
    ## ev2         -2.52633    0.48812  -5.176
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.110              
    ## ev1  0.122  0.018       
    ## ev2  0.126 -0.023 -0.626
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "77 MLM6 cn:120 ticc: 0.4 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1428.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.4255 -0.5152  0.0243  0.5976  2.4127 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             21.56    4.643   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.62790    0.31526   5.164
    ## con          0.34458    0.02463  13.992
    ## ev1          7.94814    0.46503  17.091
    ## ev2         -2.20145    0.46697  -4.714
    ## con:ev1      0.15804    0.03659   4.319
    ## con:ev2     -0.19060    0.03657  -5.212
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.118                            
    ## ev1      0.119  0.007                     
    ## ev2      0.131 -0.044 -0.627              
    ## con:ev1  0.007  0.140 -0.102  0.078       
    ## con:ev2 -0.044  0.138  0.078 -0.136 -0.640
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "78 MLM5 cn:510 ticc: 0.4 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6047.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3625 -0.7006  0.0457  0.7110  3.3792 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             21.83    4.672   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.48171    0.15153   9.779
    ## con          0.30306    0.01122  27.012
    ## ev1          7.96733    0.22328  35.683
    ## ev2         -2.06141    0.22291  -9.248
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.035              
    ## ev1  0.115 -0.059       
    ## ev2  0.113  0.007 -0.615
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "78 MLM6 cn:510 ticc: 0.4 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5956.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.4618 -0.6398  0.0288  0.6322  3.6129 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             19.74    4.443   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.41020    0.14447   9.761
    ## con          0.30097    0.01113  27.050
    ## ev1          7.78060    0.21316  36.501
    ## ev2         -1.93950    0.21232  -9.135
    ## con:ev1      0.15924    0.01619   9.837
    ## con:ev2     -0.14902    0.01667  -8.941
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.049                            
    ## ev1      0.120 -0.063                     
    ## ev2      0.109  0.012 -0.616              
    ## con:ev1 -0.064  0.080 -0.089  0.058       
    ## con:ev2  0.011  0.163  0.056 -0.037 -0.624
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "79 MLM5 cn:30 ticc: 0.8 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 328
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.88532 -0.46383  0.09979  0.62530  2.34790 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             13.88    3.725   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  4.51001    0.48992   9.206
    ## con          0.28230    0.01803  15.656
    ## ev1         13.51255    0.69032  19.574
    ## ev2          2.22632    0.67033   3.321
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.120              
    ## ev1 -0.108 -0.250       
    ## ev2 -0.089 -0.075 -0.379
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "79 MLM6 cn:30 ticc: 0.8 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 335.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.79221 -0.62136  0.09679  0.73786  2.40742 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             13.56    3.682   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  4.35671    0.54541   7.988
    ## con          0.28148    0.02053  13.708
    ## ev1         13.41599    0.71663  18.721
    ## ev2          2.36561    0.70919   3.336
    ## con:ev1      0.04052    0.02628   1.542
    ## con:ev2     -0.02881    0.02544  -1.132
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.311                            
    ## ev1     -0.208 -0.352                     
    ## ev2     -0.238 -0.232 -0.244              
    ## con:ev1 -0.361 -0.282  0.095  0.272       
    ## con:ev2 -0.244 -0.375  0.278  0.199 -0.110
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "80 MLM5 cn:120 ticc: 0.8 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1421.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.28425 -0.75689  0.05136  0.70578  2.52869 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             21.58    4.645   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  3.81507    0.30804  12.385
    ## con          0.30504    0.01295  23.560
    ## ev1         12.53293    0.45213  27.720
    ## ev2          3.04198    0.45179   6.733
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.020              
    ## ev1  0.100 -0.054       
    ## ev2  0.102  0.037 -0.605
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "80 MLM6 cn:120 ticc: 0.8 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1425
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.17866 -0.73852 -0.03605  0.63326  2.61303 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             20.93    4.575   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  3.75833    0.30397  12.364
    ## con          0.30268    0.01316  22.993
    ## ev1         12.51467    0.44544  28.095
    ## ev2          3.02210    0.44591   6.777
    ## con:ev1      0.05586    0.01836   3.043
    ## con:ev2     -0.04055    0.02012  -2.015
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.026                            
    ## ev1      0.101 -0.057                     
    ## ev2      0.104  0.053 -0.604              
    ## con:ev1 -0.060 -0.040 -0.015 -0.010       
    ## con:ev2  0.050  0.220 -0.009  0.058 -0.601
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "81 MLM5 cn:510 ticc: 0.8 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6151.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.4559 -0.6616 -0.0129  0.6640  3.4084 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             24.15    4.915   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  3.926720   0.160288   24.50
    ## con          0.298066   0.006394   46.62
    ## ev1         12.473604   0.235841   52.89
    ## ev2          2.699575   0.235898   11.44
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.068              
    ## ev1  0.118  0.009       
    ## ev2  0.117  0.024 -0.620
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "81 MLM6 cn:510 ticc: 0.8 b2:0.5 b3:0.3"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6130.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.4999 -0.6606  0.0008  0.6698  3.4584 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             23.35    4.832   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  3.914957   0.157625  24.837
    ## con          0.298641   0.006556  45.556
    ## ev1         12.413273   0.232357  53.423
    ## ev2          2.755456   0.232215  11.866
    ## con:ev1      0.051923   0.009775   5.312
    ## con:ev2     -0.053864   0.009584  -5.620
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.062                            
    ## ev1      0.117 -0.005                     
    ## ev2      0.116  0.016 -0.619              
    ## con:ev1 -0.005  0.150 -0.060  0.021       
    ## con:ev2  0.016  0.094  0.022 -0.047 -0.624
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')

    ## [1] "82 MLM5 cn:30 ticc: 0.2 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 377.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.71675 -0.74675 -0.05911  0.55711  2.26895 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             34.39    5.864   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.64708    0.76655   0.844
    ## con          0.30610    0.05218   5.866
    ## ev1          1.49621    1.06677   1.403
    ## ev2         -1.64620    1.05711  -1.557
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.049              
    ## ev1 -0.072  0.164       
    ## ev2 -0.085 -0.096 -0.420
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "82 MLM6 cn:30 ticc: 0.2 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 326.5
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.69954 -0.45298 -0.03755  0.41965  2.02123 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 7.293    2.700   
    ##  Residual             7.375    2.716   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.22315    0.61699   1.982
    ## con          0.32110    0.04138   7.759
    ## ev1          2.22949    0.85427   2.610
    ## ev2         -1.69328    0.84526  -2.003
    ## con:ev1      0.30552    0.05190   5.887
    ## con:ev2     -0.44860    0.04916  -9.125
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.010                            
    ## ev1     -0.060  0.106                     
    ## ev2     -0.090 -0.032 -0.418              
    ## con:ev1  0.117 -0.341  0.096 -0.066       
    ## con:ev2 -0.036 -0.496 -0.069 -0.015  0.011
    ## [1] "83 MLM5 cn:120 ticc: 0.2 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1480
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.29945 -0.63032  0.03774  0.59499  3.16086 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  2.471   1.572   
    ##  Residual             25.448   5.045   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.10836    0.36550   0.296
    ## con          0.29839    0.02765  10.792
    ## ev1          0.84098    0.53614   1.569
    ## ev2         -0.45844    0.53584  -0.856
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.005              
    ## ev1  0.101  0.036       
    ## ev2  0.102 -0.013 -0.604
    ## [1] "83 MLM6 cn:120 ticc: 0.2 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1389.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.41024 -0.45789  0.02817  0.52847  2.90197 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  6.905   2.628   
    ##  Residual             12.748   3.570   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.23126    0.34188   0.676
    ## con          0.31210    0.02361  13.218
    ## ev1          0.96692    0.50127   1.929
    ## ev2         -0.48302    0.50107  -0.964
    ## con:ev1      0.38177    0.03458  11.039
    ## con:ev2     -0.32612    0.03389  -9.623
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.000                            
    ## ev1      0.102  0.034                     
    ## ev2      0.101 -0.015 -0.604              
    ## con:ev1  0.034  0.099  0.023 -0.013       
    ## con:ev2 -0.016  0.042 -0.013 -0.011 -0.573
    ## [1] "84 MLM5 cn:510 ticc: 0.2 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6218.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.5736 -0.6347  0.0061  0.6941  2.8349 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.8402  0.9166  
    ##  Residual             25.0111  5.0011  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.53283    0.17636   3.021
    ## con          0.29894    0.01245  24.011
    ## ev1          1.60221    0.26414   6.066
    ## ev2         -0.25536    0.26414  -0.967
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.023              
    ## ev1  0.163 -0.006       
    ## ev2  0.163  0.008 -0.663
    ## [1] "84 MLM6 cn:510 ticc: 0.2 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 5969.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.5982 -0.6026  0.0023  0.6056  3.1787 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  3.88    1.970   
    ##  Residual             16.52    4.064   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.50082    0.16820   2.978
    ## con          0.28957    0.01230  23.551
    ## ev1          1.69154    0.25195   6.714
    ## ev2         -0.36962    0.25203  -1.467
    ## con:ev1      0.30338    0.01805  16.807
    ## con:ev2     -0.28906    0.01878 -15.388
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.025                            
    ## ev1      0.163 -0.006                     
    ## ev2      0.164  0.015 -0.663              
    ## con:ev1 -0.006  0.106  0.019 -0.018       
    ## con:ev2  0.015  0.219 -0.017  0.032 -0.664

    ## boundary (singular) fit: see help('isSingular')

    ## [1] "85 MLM5 cn:30 ticc: 0.4 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 358.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -3.14194 -0.66880 -0.06335  0.66571  2.18937 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             24.25    4.924   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.46144    0.66621   2.194
    ## con          0.28715    0.03481   8.249
    ## ev1          1.24565    0.96907   1.285
    ## ev2          1.59288    0.97039   1.641
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.156              
    ## ev1  0.107 -0.032       
    ## ev2  0.121  0.061 -0.616
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "85 MLM6 cn:30 ticc: 0.4 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 347
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.88106 -0.63588 -0.00415  0.66005  1.78565 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 10.02    3.165   
    ##  Residual             10.72    3.274   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.40163    0.75043   1.868
    ## con          0.33509    0.03550   9.438
    ## ev1          2.05746    1.10014   1.870
    ## ev2          0.78301    1.11117   0.705
    ## con:ev1      0.30937    0.05187   5.964
    ## con:ev2     -0.23302    0.05150  -4.525
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.146                            
    ## ev1      0.102 -0.022                     
    ## ev2      0.130  0.067 -0.618              
    ## con:ev1 -0.022  0.092  0.121 -0.098       
    ## con:ev2  0.069  0.072 -0.100  0.182 -0.585
    ## [1] "86 MLM5 cn:120 ticc: 0.4 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1491.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.51228 -0.63505  0.05271  0.59157  2.61360 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  3.289   1.814   
    ##  Residual             26.021   5.101   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.39136    0.40887   0.957
    ## con          0.26134    0.02075  12.593
    ## ev1          1.31211    0.61240   2.143
    ## ev2          1.43091    0.62018   2.307
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.086              
    ## ev1  0.170 -0.018       
    ## ev2  0.183  0.159 -0.666
    ## [1] "86 MLM6 cn:120 ticc: 0.4 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1448.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.69918 -0.52470 -0.02768  0.56278  3.06461 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  9.143   3.024   
    ##  Residual             16.070   4.009   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.12095    0.42645   0.284
    ## con          0.31908    0.02305  13.841
    ## ev1          1.98778    0.63509   3.130
    ## ev2          0.50686    0.64830   0.782
    ## con:ev1      0.29237    0.03822   7.650
    ## con:ev2     -0.24334    0.03116  -7.809
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.113                            
    ## ev1      0.146 -0.002                     
    ## ev2      0.205  0.115 -0.675              
    ## con:ev1 -0.002  0.451  0.090 -0.113       
    ## con:ev2  0.129 -0.128 -0.141  0.195 -0.703

    ## boundary (singular) fit: see help('isSingular')

    ## [1] "87 MLM5 cn:510 ticc: 0.4 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6193.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.93790 -0.66997  0.02606  0.67884  2.95878 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev. 
    ##  pid      (Intercept) 1.978e-14 1.407e-07
    ##  Residual             2.519e+01 5.019e+00
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 0.737402   0.164439   4.484
    ## con         0.306409   0.008797  34.829
    ## ev1         2.110440   0.243048   8.683
    ## ev2         0.232279   0.243454   0.954
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.034              
    ## ev1  0.126 -0.017       
    ## ev2  0.128  0.060 -0.628
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "87 MLM6 cn:510 ticc: 0.4 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6104.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -3.01577 -0.61971  0.03715  0.62502  2.87411 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  3.026   1.74    
    ##  Residual             19.981   4.47    
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.630834   0.167649   3.763
    ## con          0.299903   0.009150  32.775
    ## ev1          2.239190   0.247450   9.049
    ## ev2          0.001027   0.248571   0.004
    ## con:ev1      0.126394   0.013290   9.510
    ## con:ev2     -0.140161   0.013755 -10.190
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.050                            
    ## ev1      0.121 -0.026                     
    ## ev2      0.134  0.073 -0.629              
    ## con:ev1 -0.027  0.075  0.029 -0.055       
    ## con:ev2  0.072  0.173 -0.054  0.094 -0.627
    ## [1] "88 MLM5 cn:30 ticc: 0.8 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 378.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.00339 -0.61880  0.00669  0.59124  1.83768 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  5.98    2.445   
    ##  Residual             28.85    5.372   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.43146    0.91885   1.558
    ## con          0.32469    0.01948  16.666
    ## ev1          4.17333    1.40444   2.972
    ## ev2          1.47724    1.32789   1.112
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.255              
    ## ev1  0.224  0.327       
    ## ev2  0.141 -0.033 -0.630
    ## [1] "88 MLM6 cn:30 ticc: 0.8 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 380.5
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.84848 -0.57880  0.00584  0.51678  2.06126 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  7.111   2.667   
    ##  Residual             24.298   4.929   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  2.12400    0.99188   2.141
    ## con          0.31275    0.02345  13.337
    ## ev1          6.48577    1.58358   4.096
    ## ev2         -0.26109    1.41206  -0.185
    ## con:ev1      0.09637    0.03208   3.004
    ## con:ev2     -0.09557    0.03915  -2.441
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.356                            
    ## ev1      0.344  0.236                     
    ## ev2      0.019  0.038 -0.693              
    ## con:ev1  0.275 -0.094  0.499 -0.404       
    ## con:ev2  0.033  0.471 -0.296  0.323 -0.726

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "89 MLM5 cn:120 ticc: 0.8 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1502.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6454 -0.6516  0.0047  0.7101  3.2825 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev. 
    ##  pid      (Intercept) 6.671e-14 2.583e-07
    ##  Residual             3.030e+01 5.504e+00
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.87258    0.38231   4.898
    ## con          0.31505    0.00832  37.865
    ## ev1          3.71218    0.56335   6.589
    ## ev2          1.47605    0.56830   2.597
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.142              
    ## ev1  0.144 -0.001       
    ## ev2  0.124  0.132 -0.641
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "89 MLM6 cn:120 ticc: 0.8 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1496.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6167 -0.6091 -0.0090  0.7215  3.5373 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             27.96    5.288   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  1.798358   0.369854   4.862
    ## con          0.316762   0.008561  37.002
    ## ev1          3.398384   0.548481   6.196
    ## ev2          1.458264   0.546002   2.671
    ## con:ev1      0.057067   0.012248   4.659
    ## con:ev2     -0.038638   0.013224  -2.922
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.094                            
    ## ev1      0.134 -0.045                     
    ## ev2      0.121  0.118 -0.630              
    ## con:ev1 -0.046  0.033 -0.120 -0.007       
    ## con:ev2  0.113  0.250 -0.006 -0.006 -0.649
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "90 MLM5 cn:510 ticc: 0.8 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6277.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -3.09881 -0.68409 -0.00723  0.70122  2.94470 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             27.31    5.226   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 1.066406   0.172049   6.198
    ## con         0.304632   0.004138  73.614
    ## ev1         3.316259   0.254630  13.024
    ## ev2         1.441600   0.254676   5.661
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.043              
    ## ev1  0.131  0.002       
    ## ev2  0.132 -0.019 -0.633
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "90 MLM6 cn:510 ticc: 0.8 b2:0.1 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6253
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1681 -0.6691 -0.0085  0.6963  2.9035 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             26.27    5.126   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  1.083313   0.168838   6.416
    ## con          0.305041   0.004319  70.621
    ## ev1          3.239438   0.250027  12.956
    ## ev2          1.537208   0.250271   6.142
    ## con:ev1      0.036637   0.006520   5.620
    ## con:ev2     -0.038739   0.006332  -6.118
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.050                            
    ## ev1      0.130 -0.002                     
    ## ev2      0.133 -0.024 -0.633              
    ## con:ev1 -0.002  0.184 -0.046  0.040       
    ## con:ev2 -0.025  0.101  0.041 -0.063 -0.645
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "91 MLM5 cn:30 ticc: 0.2 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 366
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.47935 -0.56858  0.06442  0.63518  1.99795 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  1.513   1.230   
    ##  Residual             26.826   5.179   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.75785    0.80402   2.186
    ## con          0.31476    0.05039   6.247
    ## ev1          6.28127    1.23249   5.096
    ## ev2         -1.97886    1.21450  -1.629
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.048              
    ## ev1  0.178 -0.172       
    ## ev2  0.190  0.025 -0.681
    ## [1] "91 MLM6 cn:30 ticc: 0.2 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 354.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.54035 -0.58213 -0.00731  0.46928  2.56192 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  1.08    1.039   
    ##  Residual             19.84    4.454   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.11809    0.70410   1.588
    ## con          0.33833    0.05246   6.449
    ## ev1          5.63394    1.06840   5.273
    ## ev2         -1.65624    1.05590  -1.569
    ## con:ev1      0.35134    0.07545   4.657
    ## con:ev2     -0.25597    0.08426  -3.038
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.001                            
    ## ev1      0.199 -0.181                     
    ## ev2      0.166  0.107 -0.682              
    ## con:ev1 -0.191  0.048 -0.125  0.053       
    ## con:ev2  0.099  0.361  0.047  0.067 -0.713

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "92 MLM5 cn:120 ticc: 0.2 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1521.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.2052 -0.7858  0.0840  0.6932  3.7041 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             33.07    5.751   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.82547    0.38576   2.140
    ## con          0.26277    0.02737   9.601
    ## ev1          3.35047    0.56571   5.923
    ## ev2         -1.15145    0.56556  -2.036
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.089              
    ## ev1  0.115  0.023       
    ## ev2  0.113  0.002 -0.615
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "92 MLM6 cn:120 ticc: 0.2 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1431.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.66596 -0.63852  0.05181  0.60685  2.74318 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             21.83    4.672   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.88982    0.31360   2.837
    ## con          0.26003    0.02231  11.655
    ## ev1          3.87890    0.46212   8.394
    ## ev2         -1.54296    0.46148  -3.343
    ## con:ev1      0.34432    0.03161  10.894
    ## con:ev2     -0.24881    0.03229  -7.706
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.091                            
    ## ev1      0.116  0.023                     
    ## ev2      0.112  0.007 -0.617              
    ## con:ev1  0.024  0.005  0.103 -0.065       
    ## con:ev2  0.007  0.065 -0.063  0.090 -0.537
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "93 MLM5 cn:510 ticc: 0.2 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6426.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.7843 -0.6482  0.0236  0.6954  3.5248 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             31.72    5.632   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.09894    0.18612   5.904
    ## con          0.32402    0.01368  23.694
    ## ev1          4.36663    0.27620  15.810
    ## ev2         -1.52974    0.27616  -5.539
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.018              
    ## ev1  0.136  0.022       
    ## ev2  0.136  0.014 -0.637
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "93 MLM6 cn:510 ticc: 0.2 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6146
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.1429 -0.5484  0.0063  0.5997  4.2840 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.0     0.000   
    ##  Residual             23.8     4.878   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.12080    0.16130   6.948
    ## con          0.31112    0.01254  24.806
    ## ev1          4.40519    0.23928  18.410
    ## ev2         -1.56836    0.23927  -6.555
    ## con:ev1      0.31381    0.01810  17.333
    ## con:ev2     -0.30417    0.01911 -15.913
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.008                            
    ## ev1      0.135  0.016                     
    ## ev2      0.135  0.009 -0.637              
    ## con:ev1  0.017  0.058  0.004 -0.014       
    ## con:ev2  0.008  0.212 -0.013 -0.001 -0.639
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "94 MLM5 cn:30 ticc: 0.4 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 358.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.6504 -0.6327 -0.2724  0.6553  2.1358 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.886   0.9413  
    ##  Residual             23.393   4.8366  
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.34675    0.65519   2.056
    ## con          0.29279    0.04283   6.836
    ## ev1          6.31711    0.95225   6.634
    ## ev2          0.40865    0.93456   0.437
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.152              
    ## ev1 -0.041 -0.274       
    ## ev2  0.030  0.199 -0.526
    ## [1] "94 MLM6 cn:30 ticc: 0.4 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 347.5
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.91071 -0.58564  0.04869  0.42985  2.83086 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.8945  0.9458  
    ##  Residual             17.1057  4.1359  
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.44767    0.60210   0.744
    ## con          0.26016    0.03964   6.563
    ## ev1          6.41015    0.83390   7.687
    ## ev2         -0.29321    0.88637  -0.331
    ## con:ev1      0.24655    0.05325   4.630
    ## con:ev2     -0.20776    0.06334  -3.280
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.218                            
    ## ev1     -0.059 -0.301                     
    ## ev2      0.114  0.315 -0.532              
    ## con:ev1 -0.311 -0.145  0.010 -0.133       
    ## con:ev2  0.290  0.346 -0.119  0.382 -0.633

    ## boundary (singular) fit: see help('isSingular')

    ## [1] "95 MLM5 cn:120 ticc: 0.4 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1535.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.6623 -0.6787 -0.0144  0.5971  2.8999 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             35.12    5.926   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.38344    0.39353   3.515
    ## con          0.28945    0.02461  11.760
    ## ev1          6.07267    0.57830  10.501
    ## ev2         -1.47937    0.57603  -2.568
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.056              
    ## ev1  0.106 -0.089       
    ## ev2  0.101  0.012 -0.603
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "95 MLM6 cn:120 ticc: 0.4 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1486.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.68834 -0.69528 -0.04301  0.65896  3.14111 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  1.242   1.115   
    ##  Residual             26.287   5.127   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.21131    0.35791   3.384
    ## con          0.28473    0.02208  12.895
    ## ev1          5.58398    0.52699  10.596
    ## ev2         -1.16459    0.52307  -2.226
    ## con:ev1      0.22545    0.03112   7.244
    ## con:ev2     -0.23213    0.03193  -7.271
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.057                            
    ## ev1      0.114 -0.085                     
    ## ev2      0.093  0.011 -0.606              
    ## con:ev1 -0.089 -0.010 -0.115  0.080       
    ## con:ev2  0.011  0.063  0.078 -0.046 -0.528

    ## boundary (singular) fit: see help('isSingular')

    ## [1] "96 MLM5 cn:510 ticc: 0.4 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6528.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3604 -0.6863 -0.0312  0.7065  3.0952 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             35.04    5.919   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.62896    0.19654   8.288
    ## con          0.31354    0.01092  28.701
    ## ev1          5.80754    0.29121  19.943
    ## ev2         -0.57534    0.29125  -1.975
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.069              
    ## ev1  0.137 -0.012       
    ## ev2  0.139  0.020 -0.640
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "96 MLM6 cn:510 ticc: 0.4 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6320.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.7628 -0.6535 -0.0389  0.6439  3.3624 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  1.199   1.095   
    ##  Residual             27.056   5.202   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.55060    0.18035   8.598
    ## con          0.30429    0.01059  28.745
    ## ev1          6.07960    0.26758  22.720
    ## ev2         -0.93641    0.26822  -3.491
    ## con:ev1      0.21651    0.01560  13.877
    ## con:ev2     -0.23625    0.01592 -14.838
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.077                            
    ## ev1      0.136 -0.008                     
    ## ev2      0.142  0.037 -0.640              
    ## con:ev1 -0.009  0.117  0.065 -0.055       
    ## con:ev2  0.037  0.174 -0.054  0.094 -0.647

    ## boundary (singular) fit: see help('isSingular')

    ## [1] "97 MLM5 cn:30 ticc: 0.8 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 376
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8244 -0.6047  0.1510  0.6740  1.7353 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             32.71    5.719   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  2.29293    0.80119   2.862
    ## con          0.28481    0.02606  10.927
    ## ev1          9.17368    1.23728   7.414
    ## ev2          2.90218    1.18843   2.442
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.120              
    ## ev1  0.114  0.279       
    ## ev2  0.150  0.025 -0.622
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "97 MLM6 cn:30 ticc: 0.8 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 383.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.74530 -0.55012  0.07452  0.68699  1.74426 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.9914  0.9957  
    ##  Residual             32.2750  5.6811  
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  2.53200    0.85924   2.947
    ## con          0.27969    0.03269   8.556
    ## ev1          9.38831    1.28396   7.312
    ## ev2          2.80555    1.25125   2.242
    ## con:ev1      0.05248    0.04582   1.145
    ## con:ev2     -0.05090    0.05361  -0.949
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.016                            
    ## ev1      0.156  0.246                     
    ## ev2      0.083 -0.101 -0.622              
    ## con:ev1  0.263 -0.025  0.160 -0.100       
    ## con:ev2 -0.090  0.421 -0.083 -0.075 -0.720

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "98 MLM5 cn:120 ticc: 0.8 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1512.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3896 -0.7268 -0.0171  0.6114  2.9059 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             31.65    5.626   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  3.52711    0.39548   8.919
    ## con          0.30071    0.01021  29.439
    ## ev1          9.45181    0.58820  16.069
    ## ev2          3.12007    0.59093   5.280
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.146              
    ## ev1  0.169 -0.114       
    ## ev2  0.129  0.149 -0.661
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "98 MLM6 cn:120 ticc: 0.8 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1518.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.4274 -0.6698 -0.0375  0.6030  2.7473 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             30.86    5.555   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  3.27023    0.40098   8.156
    ## con          0.31137    0.01097  28.386
    ## ev1          8.99577    0.61139  14.714
    ## ev2          3.30467    0.58760   5.624
    ## con:ev1      0.04761    0.01778   2.678
    ## con:ev2     -0.03506    0.01449  -2.420
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.205                            
    ## ev1      0.213 -0.222                     
    ## ev2      0.101  0.181 -0.659              
    ## con:ev1 -0.209  0.387 -0.304  0.117       
    ## con:ev2  0.201 -0.194  0.138 -0.075 -0.640
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "99 MLM5 cn:510 ticc: 0.8 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6471.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.0005 -0.6392 -0.0140  0.6777  3.2202 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.00    
    ##  Residual             33.06    5.75    
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 2.900088   0.185490   15.63
    ## con         0.309557   0.005564   55.64
    ## ev1         9.241176   0.272155   33.96
    ## ev2         3.227953   0.272020   11.87
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.041              
    ## ev1  0.106  0.032       
    ## ev2  0.105 -0.003 -0.607
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "99 MLM6 cn:510 ticc: 0.8 b2:0.3 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6433.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.1038 -0.6445 -0.0065  0.6719  3.3479 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             31.42    5.605   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  2.941896   0.181022  16.252
    ## con          0.312488   0.005663  55.185
    ## ev1          9.383069   0.266040  35.269
    ## ev2          3.127767   0.265512  11.780
    ## con:ev1      0.059928   0.008440   7.101
    ## con:ev2     -0.050982   0.008293  -6.148
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.052                            
    ## ev1      0.109  0.041                     
    ## ev2      0.103 -0.006 -0.608              
    ## con:ev1  0.040  0.149  0.075 -0.047       
    ## con:ev2 -0.006  0.099 -0.048  0.044 -0.626
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "100 MLM5 cn:30 ticc: 0.2 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 390.5
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.87898 -0.74606 -0.03419  0.56637  2.10887 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.00    
    ##  Residual             44.22    6.65    
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.17270    1.05773   0.163
    ## con          0.33380    0.07771   4.295
    ## ev1          6.16691    1.63155   3.780
    ## ev2         -2.12930    1.60816  -1.324
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.109              
    ## ev1  0.233  0.176       
    ## ev2  0.211 -0.050 -0.711
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "100 MLM6 cn:30 ticc: 0.2 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 387.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.92695 -0.59541 -0.07544  0.59327  2.26764 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             40.13    6.334   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)   0.6013     1.0393   0.579
    ## con           0.2449     0.1039   2.357
    ## ev1           7.2246     1.6250   4.446
    ## ev2          -2.7886     1.5562  -1.792
    ## con:ev1       0.3761     0.1383   2.719
    ## con:ev2      -0.4390     0.1796  -2.445
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.153                            
    ## ev1      0.285  0.158                     
    ## ev2      0.162 -0.044 -0.721              
    ## con:ev1  0.185 -0.172  0.265 -0.167       
    ## con:ev2 -0.038  0.570 -0.123  0.092 -0.764
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "101 MLM5 cn:120 ticc: 0.2 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1547.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.0353 -0.6329  0.0252  0.6496  2.4278 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             36.99    6.082   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.58910    0.42317   3.755
    ## con          0.31634    0.03376   9.371
    ## ev1          6.51279    0.62250  10.462
    ## ev2         -1.91758    0.62494  -3.068
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.153              
    ## ev1  0.143 -0.001       
    ## ev2  0.156  0.088 -0.644
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "101 MLM6 cn:120 ticc: 0.2 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1516.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.3503 -0.5262 -0.0338  0.5425  2.6437 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             31.54    5.616   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.45473    0.39679   3.666
    ## con          0.33108    0.03527   9.387
    ## ev1          7.29590    0.58714  12.426
    ## ev2         -2.83398    0.59803  -4.739
    ## con:ev1      0.32877    0.05413   6.074
    ## con:ev2     -0.31050    0.05237  -5.929
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.207                            
    ## ev1      0.128  0.014                     
    ## ev2      0.180  0.115 -0.655              
    ## con:ev1  0.014  0.231  0.191 -0.174       
    ## con:ev2  0.117  0.138 -0.183  0.262 -0.685
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "102 MLM5 cn:510 ticc: 0.2 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6613.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.5763 -0.7207  0.0299  0.7262  2.7124 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.0     0.000   
    ##  Residual             38.1     6.173   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.24729    0.20189   6.178
    ## con          0.31737    0.01717  18.488
    ## ev1          7.42722    0.29800  24.924
    ## ev2         -3.27793    0.29803 -10.999
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.045              
    ## ev1  0.124 -0.013       
    ## ev2  0.123  0.020 -0.626
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "102 MLM6 cn:510 ticc: 0.2 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6362.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.0622 -0.6206  0.0106  0.6119  3.0842 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             29.49    5.431   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.18287    0.17765   6.658
    ## con          0.31113    0.01559  19.954
    ## ev1          7.23219    0.26258  27.542
    ## ev2         -3.14764    0.26231 -12.000
    ## con:ev1      0.35271    0.02275  15.501
    ## con:ev2     -0.35846    0.02308 -15.528
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.044                            
    ## ev1      0.125 -0.021                     
    ## ev2      0.122  0.017 -0.625              
    ## con:ev1 -0.021  0.089 -0.055  0.023       
    ## con:ev2  0.017  0.130  0.023 -0.029 -0.612
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "103 MLM5 cn:30 ticc: 0.4 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 368.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.2351 -0.5085 -0.1004  0.5563  2.2803 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  9.28    3.046   
    ##  Residual             21.81    4.670   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  3.10276    0.83554   3.713
    ## con          0.33799    0.05238   6.453
    ## ev1          9.57689    1.15465   8.294
    ## ev2         -2.05389    1.14661  -1.791
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.120              
    ## ev1 -0.060 -0.158       
    ## ev2 -0.093  0.106 -0.421
    ## [1] "103 MLM6 cn:30 ticc: 0.4 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 352.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6027 -0.5659  0.1244  0.3909  2.6538 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  3.096   1.760   
    ##  Residual             17.156   4.142   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  2.65862    0.64276   4.136
    ## con          0.34388    0.04232   8.126
    ## ev1          8.86630    0.90043   9.847
    ## ev2         -1.72666    0.87587  -1.971
    ## con:ev1      0.20808    0.05742   3.624
    ## con:ev2     -0.28345    0.05660  -5.008
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.115                            
    ## ev1     -0.027 -0.131                     
    ## ev2     -0.105  0.100 -0.429              
    ## con:ev1 -0.136 -0.117 -0.218  0.088       
    ## con:ev2  0.102 -0.158  0.087 -0.051 -0.347

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "104 MLM5 cn:120 ticc: 0.4 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1596.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3195 -0.6385  0.0430  0.6260  2.2389 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             45.52    6.747   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  2.23439    0.50624   4.414
    ## con          0.29716    0.03152   9.427
    ## ev1          9.00471    0.76417  11.784
    ## ev2         -1.10423    0.76654  -1.441
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.093              
    ## ev1  0.197  0.010       
    ## ev2  0.203  0.079 -0.692
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "104 MLM6 cn:120 ticc: 0.4 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1591.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.2573 -0.5861 -0.0190  0.5825  2.3392 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             43.42    6.589   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  2.21001    0.50034   4.417
    ## con          0.31531    0.03539   8.911
    ## ev1          9.39704    0.75486  12.449
    ## ev2         -1.53747    0.76410  -2.012
    ## con:ev1      0.19493    0.05365   3.633
    ## con:ev2     -0.15187    0.05352  -2.838
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.150                            
    ## ev1      0.183  0.006                     
    ## ev2      0.218  0.106 -0.698              
    ## con:ev1  0.006  0.197  0.135 -0.138       
    ## con:ev2  0.107  0.190 -0.140  0.200 -0.692
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "105 MLM5 cn:510 ticc: 0.4 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6611.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1498 -0.6702  0.0174  0.7088  3.1710 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             38.01    6.165   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.56372    0.20212   7.737
    ## con          0.27965    0.01327  21.074
    ## ev1          8.43392    0.29871  28.235
    ## ev2         -1.40135    0.29859  -4.693
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.048              
    ## ev1  0.128  0.035       
    ## ev2  0.127  0.021 -0.627
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "105 MLM6 cn:510 ticc: 0.4 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6452.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.3463 -0.6053 -0.0036  0.6224  3.5274 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             32.18    5.673   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  1.61015    0.18646   8.635
    ## con          0.28794    0.01289  22.335
    ## ev1          8.77779    0.27602  31.801
    ## ev2         -1.70663    0.27567  -6.191
    ## con:ev1      0.24684    0.01924  12.831
    ## con:ev2     -0.22352    0.01903 -11.748
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.069                            
    ## ev1      0.129  0.041                     
    ## ev2      0.126  0.020 -0.629              
    ## con:ev1  0.041  0.152  0.090 -0.072       
    ## con:ev2  0.020  0.121 -0.073  0.076 -0.638
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "106 MLM5 cn:30 ticc: 0.8 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 381.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.8937 -0.6208  0.1143  0.5950  2.5698 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             35.99    5.999   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  3.25743    0.78173   4.167
    ## con          0.24847    0.02905   8.554
    ## ev1         13.33703    1.15375  11.560
    ## ev2          4.50135    1.14927   3.917
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.018              
    ## ev1  0.065  0.201       
    ## ev2  0.065  0.181 -0.508
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "106 MLM6 cn:30 ticc: 0.8 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 386
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.9317 -0.6339  0.0239  0.6069  2.6220 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.00    
    ##  Residual             34.57    5.88    
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  3.04532    0.82842   3.676
    ## con          0.21926    0.03487   6.287
    ## ev1         14.03326    1.17975  11.895
    ## ev2          3.60062    1.20826   2.980
    ## con:ev1      0.08642    0.04613   1.873
    ## con:ev2     -0.11614    0.05909  -1.966
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.212                            
    ## ev1      0.020  0.039                     
    ## ev2      0.087  0.268 -0.556              
    ## con:ev1  0.042 -0.189  0.254 -0.341       
    ## con:ev2  0.231  0.514 -0.273  0.330 -0.723
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "107 MLM5 cn:120 ticc: 0.8 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1501.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.4166 -0.7070 -0.1146  0.5746  2.8026 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             30.22    5.497   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  4.13750    0.36776  11.251
    ## con          0.29891    0.01105  27.056
    ## ev1         13.66997    0.54065  25.284
    ## ev2          3.63879    0.54526   6.673
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.051              
    ## ev1  0.114 -0.010       
    ## ev2  0.119 -0.130 -0.609
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "107 MLM6 cn:120 ticc: 0.8 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1497.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.49591 -0.74581 -0.08609  0.58447  2.64685 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.0     0.000   
    ##  Residual             28.3     5.319   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  4.28553    0.36055  11.886
    ## con          0.30055    0.01082  27.766
    ## ev1         13.40447    0.52729  25.421
    ## ev2          3.98778    0.53402   7.468
    ## con:ev1      0.05205    0.01600   3.253
    ## con:ev2     -0.06280    0.01534  -4.094
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.066                            
    ## ev1      0.095 -0.008                     
    ## ev2      0.131 -0.122 -0.615              
    ## con:ev1 -0.008  0.125 -0.066  0.118       
    ## con:ev2 -0.127  0.006  0.124 -0.149 -0.570
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## [1] "108 MLM5 cn:510 ticc: 0.8 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6667.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.5052 -0.7206  0.0174  0.7269  3.6558 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             40.11    6.333   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  4.593383   0.206508   22.24
    ## con          0.285108   0.006967   40.92
    ## ev1         13.819892   0.304808   45.34
    ## ev2          3.966806   0.305201   13.00
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.009              
    ## ev1  0.121 -0.005       
    ## ev2  0.121 -0.051 -0.622
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "108 MLM6 cn:510 ticc: 0.8 b2:0.5 b3:0.5"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 6624.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.5683 -0.6717 -0.0014  0.6903  3.8202 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept)  0.00    0.000   
    ##  Residual             37.96    6.161   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  4.634321   0.201319  23.020
    ## con          0.283725   0.006902  41.106
    ## ev1         13.744065   0.296825  46.304
    ## ev2          4.096611   0.297401  13.775
    ## con:ev1      0.071180   0.009919   7.176
    ## con:ev2     -0.066556   0.010192  -6.530
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.021                            
    ## ev1      0.118  0.002                     
    ## ev2      0.123 -0.054 -0.623              
    ## con:ev1  0.002  0.045 -0.019  0.046       
    ## con:ev2 -0.054  0.122  0.045 -0.056 -0.587
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')
    ## 
    ## [1] "109 MLM5 cn:30 ticc: 0.2 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 210.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.17390 -0.57005  0.00255  0.68141  1.70470 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.6352   0.797   
    ##  Residual             1.2819   1.132   
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.01943    0.22277   0.087
    ## con          0.27626    0.01817  15.205
    ## ev1         -0.03984    0.33195  -0.120
    ## ev2          0.26038    0.33509   0.777
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.071              
    ## ev1  0.156  0.021       
    ## ev2  0.143 -0.138 -0.652
    ## [1] "109 MLM6 cn:30 ticc: 0.2 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 219.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.20175 -0.46025 -0.04394  0.48762  1.66341 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.5649   0.7516  
    ##  Residual             1.3230   1.1502  
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)  0.05990    0.22223   0.270
    ## con          0.26495    0.02351  11.272
    ## ev1         -0.10767    0.32891  -0.327
    ## ev2          0.28599    0.33331   0.858
    ## con:ev1     -0.01048    0.02833  -0.370
    ## con:ev2     -0.02809    0.04160  -0.675
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.062                            
    ## ev1      0.129  0.082                     
    ## ev2      0.166 -0.211 -0.649              
    ## con:ev1  0.101 -0.454 -0.001  0.142       
    ## con:ev2 -0.179  0.640  0.098 -0.166 -0.743
    ## [1] "110 MLM5 cn:120 ticc: 0.2 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 849
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.42827 -0.53535 -0.01377  0.58027  1.87830 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 1.177    1.085   
    ##  Residual             1.086    1.042   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.010020   0.120356   0.083
    ## con          0.299681   0.008964  33.431
    ## ev1         -0.131959   0.173382  -0.761
    ## ev2          0.201696   0.173153   1.165
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.014              
    ## ev1  0.048 -0.053       
    ## ev2  0.049  0.015 -0.550
    ## [1] "110 MLM6 cn:120 ticc: 0.2 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 861.5
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.47913 -0.52617  0.03173  0.55318  1.98020 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 1.126    1.061   
    ##  Residual             1.114    1.056   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.001831   0.119247   0.015
    ## con          0.300092   0.009094  33.001
    ## ev1         -0.136452   0.171576  -0.795
    ## ev2          0.208721   0.171473   1.217
    ## con:ev1      0.014657   0.012769   1.148
    ## con:ev2     -0.002874   0.013423  -0.214
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.014                            
    ## ev1      0.049 -0.053                     
    ## ev2      0.047  0.020 -0.550              
    ## con:ev1 -0.054 -0.020 -0.024  0.017       
    ## con:ev2  0.019  0.121  0.016  0.026 -0.555
    ## [1] "111 MLM5 cn:510 ticc: 0.2 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 3438.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6717 -0.5622  0.0320  0.5878  2.7486 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.9321   0.9655  
    ##  Residual             0.9867   0.9933  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept) -0.010860   0.055178  -0.197
    ## con          0.296478   0.003770  78.638
    ## ev1          0.003503   0.081280   0.043
    ## ev2          0.032765   0.081262   0.403
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.067              
    ## ev1  0.122  0.021       
    ## ev2  0.120 -0.005 -0.623
    ## [1] "111 MLM6 cn:510 ticc: 0.2 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 3453.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7699 -0.5541  0.0430  0.5831  2.7385 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.9336   0.9662  
    ##  Residual             0.9862   0.9931  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept) -0.008633   0.055249  -0.156
    ## con          0.297065   0.003932  75.546
    ## ev1          0.013311   0.081657   0.163
    ## ev2          0.026135   0.081495   0.321
    ## con:ev1      0.007735   0.005867   1.318
    ## con:ev2     -0.005651   0.005748  -0.983
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.075                            
    ## ev1      0.125  0.035                     
    ## ev2      0.119 -0.003 -0.624              
    ## con:ev1  0.034  0.151  0.091 -0.055       
    ## con:ev2 -0.003  0.093 -0.056  0.067 -0.625
    ## [1] "112 MLM5 cn:30 ticc: 0.4 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 208.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.81078 -0.48574 -0.05435  0.49685  1.70434 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 1.477    1.2155  
    ##  Residual             0.830    0.9111  
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) -0.32754    0.28598  -1.145
    ## con          0.30968    0.01397  22.168
    ## ev1         -0.04631    0.43632  -0.106
    ## ev2          0.05595    0.43457   0.129
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.013              
    ## ev1  0.189  0.135       
    ## ev2  0.187 -0.102 -0.691
    ## [1] "112 MLM6 cn:30 ticc: 0.4 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 216.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.75987 -0.53950 -0.05361  0.55144  1.66577 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 1.5243   1.2346  
    ##  Residual             0.7768   0.8814  
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) -0.22390    0.29273  -0.765
    ## con          0.31564    0.01577  20.019
    ## ev1          0.08205    0.44528   0.184
    ## ev2          0.02404    0.44124   0.054
    ## con:ev1      0.04580    0.02373   1.930
    ## con:ev2     -0.03076    0.02402  -1.281
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.028                            
    ## ev1      0.206  0.183                     
    ## ev2      0.180 -0.155 -0.692              
    ## con:ev1  0.185  0.177  0.146 -0.032       
    ## con:ev2 -0.153  0.211 -0.031 -0.078 -0.692
    ## [1] "113 MLM5 cn:120 ticc: 0.4 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 836.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.43779 -0.56078 -0.00805  0.59801  2.46034 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.8988   0.948   
    ##  Residual             1.1350   1.065   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept) -0.071106   0.112137  -0.634
    ## con          0.290968   0.006927  42.008
    ## ev1          0.131566   0.162948   0.807
    ## ev2         -0.136717   0.162994  -0.839
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.023              
    ## ev1  0.077  0.027       
    ## ev2  0.076 -0.036 -0.579
    ## [1] "113 MLM6 cn:120 ticc: 0.4 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 850.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.42215 -0.58891 -0.01295  0.59535  2.41390 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.904    0.9508  
    ##  Residual             1.144    1.0697  
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept) -0.070689   0.112604  -0.628
    ## con          0.290589   0.006996  41.539
    ## ev1          0.128508   0.163682   0.785
    ## ev2         -0.135776   0.163556  -0.830
    ## con:ev1     -0.002213   0.010006  -0.221
    ## con:ev2     -0.002486   0.010144  -0.245
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.022                            
    ## ev1      0.078  0.031                     
    ## ev2      0.076 -0.037 -0.579              
    ## con:ev1  0.031  0.032  0.042 -0.006       
    ## con:ev2 -0.037  0.071 -0.006 -0.005 -0.554
    ## [1] "114 MLM5 cn:510 ticc: 0.4 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 3499.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.55010 -0.55942  0.01834  0.52808  3.09219 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 1.069    1.034   
    ##  Residual             1.009    1.004   
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.087217   0.058477   1.491
    ## con          0.299753   0.003464  86.542
    ## ev1         -0.010385   0.086731  -0.120
    ## ev2          0.054673   0.086822   0.630
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.008              
    ## ev1  0.133  0.029       
    ## ev2  0.134 -0.054 -0.636
    ## [1] "114 MLM6 cn:510 ticc: 0.4 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 3512
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.56553 -0.54399  0.02387  0.52313  3.11077 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 1.0863   1.0423  
    ##  Residual             0.9952   0.9976  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error t value
    ## (Intercept)  0.0877645  0.0587925   1.493
    ## con          0.3015595  0.0035834  84.154
    ## ev1         -0.0040783  0.0870975  -0.047
    ## ev2          0.0526205  0.0872707   0.603
    ## con:ev1      0.0101516  0.0052170   1.946
    ## con:ev2     -0.0007986  0.0053388  -0.150
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.017                            
    ## ev1      0.131  0.038                     
    ## ev2      0.137 -0.063 -0.636              
    ## con:ev1  0.039  0.082  0.010  0.025       
    ## con:ev2 -0.063  0.147  0.025 -0.058 -0.617
    ## [1] "115 MLM5 cn:30 ticc: 0.8 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 206.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.71427 -0.54497 -0.01491  0.60325  1.55277 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.9785   0.9892  
    ##  Residual             0.9425   0.9708  
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 0.144777   0.231654   0.625
    ## con         0.295991   0.008179  36.190
    ## ev1         0.443671   0.338185   1.312
    ## ev2         0.052535   0.336638   0.156
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.188              
    ## ev1  0.136 -0.138       
    ## ev2  0.092  0.100 -0.620
    ## [1] "115 MLM6 cn:30 ticc: 0.8 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 218
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.74152 -0.56745 -0.01422  0.62266  1.72990 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 1.1431   1.0692  
    ##  Residual             0.8542   0.9242  
    ## Number of obs: 60, groups:  pid, 30
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept)  0.134835   0.242244   0.557
    ## con          0.293539   0.008620  34.052
    ## ev1          0.456362   0.361760   1.262
    ## ev2          0.089733   0.352491   0.255
    ## con:ev1      0.002202   0.011397   0.193
    ## con:ev2     -0.020702   0.013289  -1.558
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.156                            
    ## ev1      0.154 -0.092                     
    ## ev2      0.081  0.067 -0.620              
    ## con:ev1 -0.104 -0.191 -0.228  0.102       
    ## con:ev2  0.064  0.244  0.085 -0.095 -0.552
    ## [1] "116 MLM5 cn:120 ticc: 0.8 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 826.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5155 -0.5620 -0.0211  0.5527  2.1044 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.936    0.9675  
    ##  Residual             1.047    1.0233  
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept) -0.038985   0.116395  -0.335
    ## con          0.293785   0.004594  63.945
    ## ev1          0.085439   0.175159   0.488
    ## ev2         -0.057473   0.173583  -0.331
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con  0.050              
    ## ev1  0.141  0.176       
    ## ev2  0.128 -0.115 -0.643
    ## [1] "116 MLM6 cn:120 ticc: 0.8 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 842.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5162 -0.5561 -0.0156  0.5279  2.0890 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 0.9584   0.979   
    ##  Residual             1.0477   1.024   
    ## Number of obs: 240, groups:  pid, 120
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error t value
    ## (Intercept) -0.0467562  0.1204562  -0.388
    ## con          0.2930377  0.0051794  56.578
    ## ev1          0.0695144  0.1815669   0.383
    ## ev2         -0.0478659  0.1770089  -0.270
    ## con:ev1     -0.0022275  0.0076252  -0.292
    ## con:ev2      0.0001539  0.0079811   0.019
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con      0.065                            
    ## ev1      0.180  0.222                     
    ## ev2      0.108 -0.169 -0.646              
    ## con:ev1  0.227  0.114  0.210 -0.070       
    ## con:ev2 -0.161  0.243 -0.065 -0.052 -0.680
    ## [1] "117 MLM5 cn:510 ticc: 0.8 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 3460.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.30379 -0.60140  0.02131  0.53661  2.85689 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 1.0467   1.0231  
    ##  Residual             0.9597   0.9796  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept) -0.021447   0.056929  -0.377
    ## con          0.297511   0.002068 143.848
    ## ev1         -0.048719   0.083889  -0.581
    ## ev2          0.002945   0.083952   0.035
    ## 
    ## Correlation of Fixed Effects:
    ##     (Intr) con    ev1   
    ## con -0.051              
    ## ev1  0.117  0.033       
    ## ev2  0.121 -0.051 -0.621
    ## [1] "117 MLM6 cn:510 ticc: 0.8 b2:0 b3:0"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con + (1 | pid)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 3475.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.32657 -0.59921  0.01796  0.53025  2.87722 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  pid      (Intercept) 1.0359   1.0178  
    ##  Residual             0.9619   0.9808  
    ## Number of obs: 1020, groups:  pid, 510
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error t value
    ## (Intercept) -0.020309   0.056858  -0.357
    ## con          0.298263   0.002116 140.949
    ## ev1         -0.050333   0.083703  -0.601
    ## ev2          0.003168   0.084068   0.038
    ## con:ev1      0.005598   0.003099   1.806
    ## con:ev2     -0.001145   0.003098  -0.370
    ## 
    ## Correlation of Fixed Effects:
    ##         (Intr) con    ev1    ev2    con:v1
    ## con     -0.056                            
    ## ev1      0.114  0.036                     
    ## ev2      0.126 -0.062 -0.622              
    ## con:ev1  0.036  0.099 -0.028  0.044       
    ## con:ev2 -0.062  0.098  0.044 -0.094 -0.601
