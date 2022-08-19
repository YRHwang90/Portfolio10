




sib_data<-reconstruct(sim_data)

sib0 <- model_sib(sib_data,
                  gender_composition=2,
                  effect_code=TRUE,
                  interaction=FALSE,
                  return_data=FALSE)


## what is sib 1?
## gender composition is thre, and no interaction model
## sib_model<- lm(y_diff  ~ y_mean + con_mean + con_diff + ev1 + ev2, data=dat3) // ev = effect code variable

sib1<-model_sib(sib_data,
                gender_composition=3,
                effect_code=TRUE,
                interaction=FALSE,
                return_data=FALSE)



## what is sib 2?
## lm(y_diff  ~ y_mean + con_mean + con_diff + gender_composition_two_effsim
#                     , data=dat3)
sib2<-model_sib(sib_data,
                gender_composition=2,
                effect_code=TRUE,
                interaction=TRUE,
                return_data=FALSE)


## what is sib 3?
## lm(y_diff  ~ y_mean + con_mean + con_diff + ev1 + ev2 +
#                      ev1*con_diff + ev2*con_diff, data=dat3)


sib3<-model_sib(sib_data,
                gender_composition=3,
                effect_code=TRUE,
                interaction=TRUE,
                return_data=FALSE)

ind_sim_data<-dyad2ind(sim_data)

MLM0<-model_null(ind_sim_data)


# MLM1: NO level 2 model without interaction
#lmer(y~1 + catwtn + con +  (1|pid), data=data)

# why does the mlm model require wide data?
# No, it requires dyadic structure except null model function

MLM1<-no2model(sim_data,
               interaction=FALSE)



# MLM2: No level 2 model with interaction
#lmer(y~1 + catwtn + con + catwtn*con + (1|pid), data=data)
MLM2<-no2model(sim_data,
               interaction=TRUE
)

MLM3<- yes2model(sim_data,
                 gender_composition=2,
                 interaction=FALSE)



MLM4<-yes2model(sim_data,
                gender_composition=2,
                interaction=TRUE)





MLM5<- yes2model(sim_data,
                 gender_composition=3,
                 interaction=FALSE)

MLM6<- yes2model(sim_data,
                 gender_composition=3,
                 interaction=TRUE)
