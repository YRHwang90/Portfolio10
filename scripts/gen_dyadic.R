# Data generation function (dyad structure)

# cat is for continuous variable.  (gender)
# con is for binary variable.  (education)
# data generation equation
#y = beta1*edu + beta2 *gender  + beta3 * interaction *within*whtin  + err bt + err within

## Name issue : con = continuous variable = V1
              # cat = categorical variable  = V2

gen_dyadic <- function(full=FALSE,
                       researcher=FALSE, # true releases only what the researcher would have
                       cn=120, # Clusternumber
                       con_mu=c(0,0,0,0), #mean of continuous variable
                       con_vcovmatrix=matrix(c(1,0,0,0,
                                            0,1,0,1,
                                            0,0,1,0,
                                            0,1,0,1),
                                          nrow = 4,
                                          byrow = TRUE,
                                          dimnames = list(c("con_wtn_1","con_btw_1","con_wtn_2","con_btw_2"),
                                                          c("con_wtn_1","con_btw_1","con_wtn_2","con_btw_2"))),
                      con_icc_wtn = .3,
                      con_icc_btw = .7,
                      cat_icc_wtn = .3,
                      cat_icc_btw = .7, # if con/cat things are different, it may cuase some problem in interaction
                      beta1 = .3, #beta weight for continuous variable on individual regression
                      beta2 = .3, #beta weight for the categorical variable on individual regression
                      beta3 = 0 # beta weight for interaction

){
  ## adjusting between and within variance in continuous variable
require(MASS)
require(tidyverse)



  ## Continuous variable
  dat <- MASS::mvrnorm(cn, # cluster  number
                       mu = con_mu, # Mu
                       Sigma = con_vcovmatrix) %>% # Covariance matrix
    as.data.frame() %>%
    mutate(
      pid = 1:cn, ## id column
      #con_1 is education score of sibling 1
      con_1 = sqrt(con_icc_wtn) * con_wtn_1 + sqrt(con_icc_btw) * con_btw_1,
      # con_2 is education score of sibling 2
      con_2 = sqrt(con_icc_wtn) * con_wtn_2 + sqrt(con_icc_btw) * con_btw_2,
      ## categorical within (individual gender, effect coded)
      catwtn_1 = sample(rep(c(-1, 1), cn / 2)),
      catwtn_2 = sample(rep(c(-1, 1), cn / 2)),
      ## categorical between (Gender composition) -- we are assumed that mixed sex sibing is beneficial to the health, compared to the
      ## same-sex siblign. this shold be addressed in limitation.
      catbtw_1 = case_when(
        catwtn_1 == 1 & catwtn_2 == 1 ~ 1,
        catwtn_1 == 1 & catwtn_2 == -1 ~ -1,
        catwtn_1 == -1 & catwtn_2 == 1 ~ -1,
        catwtn_1 == -1 & catwtn_2 == -1 ~ 1,
      ),
      catbtw_2 = catbtw_1,

      cat_1 = sqrt(cat_icc_btw) * catbtw_1 + sqrt(cat_icc_wtn) * catwtn_1,
      cat_2 = sqrt(cat_icc_btw) * catbtw_2 + sqrt(con_icc_wtn) * catwtn_2,

      er_1 = rnorm(cn, mean=0, sd=1), # error within for sibling 1
      er_2 = rnorm(cn, mean=0, sd=1), # error within for sibling 2
      erbt_1 = rnorm(cn, mean=0, sd=1), # error between
      erbt_2 = erbt_1,

# individual regression
# error variance
# interaction is within level. should be clarified in the manuscripts
      y_1 = beta1 * con_1 + beta2 * cat_1 + beta3 * sqrt(con_icc_wtn)* con_wtn_1* catwtn_1 + erbt_1 + er_1,
      y_2 = beta1 * con_2 + beta2 * cat_2 + beta3 * sqrt(con_icc_wtn)* con_wtn_2* catwtn_2 + erbt_2 + er_2,

    )%>%
    rename(conwtn_1 = con_wtn_1,
           conwtn_2 = con_wtn_2)
if(full){ # redundant code, does not have to pay attention this part.
  return(dat)
  } else if(researcher){
    dat %>%
      mutate(
      gender_1 = case_when(
        catwtn_1 == 1 ~ "male(1)",
        catwtn_1 == -1 ~ "female(-1)"),
      gender_2 = case_when(
        catwtn_2 == 1 ~ "male(1)",
        catwtn_2 == -1 ~ "female(-1)")) %>%
      subset(select=c(pid,y_1,y_2,
                      con_1,con_2,
                      gender_1,gender_2)) %>%
  rename(health_1 = y_1,
         health_2 = y_2,
         edu_1 = con_1,
         edu_2 = con_2) %>%
      return()
} else {
  dat %>%
    subset(select=c(pid,y_1,y_2,
                    con_1, con_2,
                    cat_1, cat_2,
                    catwtn_1, catwtn_2,
                    catbtw_1, catbtw_2,
                    conwtn_1, conwtn_2
                    )) %>%
  return()
  }
  }
