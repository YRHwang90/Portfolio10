# function for Discordant kinship model

model_sib<-function(data,
              gender_composition=2,
              effect_code=FALSE,
              interaction=FALSE,
              return_data=FALSE
){

  if(!interaction){
    if(gender_composition==2 & !effect_code){

      ## regression
      sib_model<- lm(y_diff  ~ y_mean + con_mean + con_diff + factor(gender_composition_two), data=data)

    # combine this two

}


    if(gender_composition==3 & !effect_code ){

      ## regression
      sib_model<- lm(y_diff  ~ y_mean + con_mean + con_diff + factor(gender_composition_three), data=data)

    # combine this two

}


    if(gender_composition==2 & effect_code ){

      ## regression
      sib_model<- lm(y_diff  ~ y_mean + con_mean + con_diff + gender_composition_two_eff, data=data)

}


    if(gender_composition==3 & effect_code ){

      ## regression
      sib_model<- lm(y_diff  ~ y_mean + con_mean + con_diff + ev1 + ev2, data=data)

}

  }

  if(interaction){

    if(gender_composition==2 & !effect_code){

      ## regression
      sib_model<- lm(y_diff  ~ y_mean + con_mean + con_diff + factor(gender_composition_two)
                     + factor(gender_composition_two)*con_diff, data=data)
    summary<-summary(sib_model)
    # combine this two
    out<-list(sib_model,summary)

}


    if(gender_composition==3 & !effect_code){

      ## regression
      sib_model<- lm(y_diff  ~ y_mean + con_mean + con_diff + factor(gender_composition_three)+
                       factor(gender_composition_three)*con_diff, data=data)
    summary<-summary(sib_model)
    # combine this two
    out<-list(sib_model,summary)

}


    if(gender_composition==2 & effect_code){

      ## regression
      sib_model<- lm(y_diff  ~ y_mean + con_mean + con_diff + gender_composition_two_eff + gender_composition_two_eff:con_diff
                     , data=data)
    summary<-summary(sib_model)
    # combine this two
    out<-list(sib_model,summary)
}


    if(gender_composition==3 & effect_code){

      ## regression
      sib_model<- lm(y_diff  ~ y_mean + con_mean + con_diff + ev1 + ev2 +
                       ev1*con_diff + ev2*con_diff, data=data)



}
  }

  if(return_data){
    return(data)
  }
  if(!return_data){
    # print it
    return(sib_model)
  }
}

