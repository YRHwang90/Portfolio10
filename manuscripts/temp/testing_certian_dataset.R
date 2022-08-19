results <- as.data.frame(conditions[108,])
j=108
sim_data <- gen_dyadic(
  cn= conditions$cn[j],
  beta1 = conditions$beta1[j],
  beta2 = conditions$beta2[j],
  beta3 = conditions$beta3[j],
  con_icc_wtn = conditions$wtn[j],
  con_icc_btw = conditions$btw[j],
  cat_icc_wtn = conditions$wtn[j],
  cat_icc_btw = conditions$btw[j])

sib_data<-reconstruct(sim_data)

tryc(unlist(myTryCatch(
  sib0<-model_sib(sib_data,
                  gender_composition=2,
                  effect_code=TRUE,
                  interaction=FALSE,
                  return_data=FALSE)
)$warning))


tryc(unlist(myTryCatch
            (sib1<-model_sib(sib_data,
                             gender_composition=3,
                             effect_code=TRUE,
                             interaction=FALSE,
                             return_data=FALSE)
            )$warning))

tryc(unlist(myTryCatch
                                     (sib2<-model_sib(sib_data,
                                                      gender_composition=2,
                                                      effect_code=TRUE,
                                                      interaction=TRUE,
                                                      return_data=FALSE)
                                     )$warning))


ind_sim_data<-dyad2ind(sim_data)
#
tryc(unlist(
  myTryCatch
  (MLM0<-model_null(ind_sim_data))
  $warning))

tryc(unlist(myTryCatch
                                     (MLM1<-no2model(sim_data,
                                                     interaction=FALSE))

tryc(unlist(myTryCatch(MLM2<-no2model(sim_data,
                                                               interaction=TRUE))$warning))
tryc(unlist(myTryCatch(MLM3<-yes2model(sim_data,
                                                                gender_composition=2,



tryc(unlist(myTryCatch(
  MLM4<-yes2model(sim_data,
                  gender_composition=2,
                  interaction=TRUE))$warning))

tryc(unlist(myTryCatch(
  MLM5<- yes2model(sim_data,
                   gender_composition=3,
                   interaction=FALSE))$warning))

tryc(unlist(myTryCatch(

  MLM6<- yes2model(sim_data,
                   gender_composition=3,
                   interaction=TRUE))$warning))


summary(sib4)$optinfo$conv$lme4$messages
