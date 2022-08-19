# The folowing function is to run MLMS
# activate dyad2ind() function beforehand

## First function:

##this model deals with the multilevel models.
## level 1: gender / education
## No level 2 variable

## The model is as follows:
# Yij= β0j + β1jXij+ β2jSij + β3jXijSij + eij
# Level 2 equation:
# β0j= γ00 + u0j
# β1j= γ10
# β2j= γ20
# β3j= γ30
# the fonts are breaking down everytime I corporate the changes in github.


# the level 1 interaction is,
#the whole education ( edu within + edu between)* gender within

# plug dyadic structure data

# when there is no level-2 variable.
no2model <- function(data,
                     interaction=TRUE

){
  require(lme4)


  if(interaction){
    result <- lmer(y~1 + catwtn + con + catwtn:con + (1|pid), data=data)}

  if(!interaction){
    result <- lmer(y~1 + catwtn + con +  (1|pid), data=data)
  }

  return(result)
}

## second function

##this model deals with the multilevel models.
## level 1:  education
## level 2: gender composition

## gender composition can be either have two levels or three levels.
## two level model(mixed sex vs same sex )
#
# Level 1:
#   Yij= β0j + β1jXij+ eij
#   Level 2:
#     β0j= γ00 + γ01Si+ u0j
#   β1i= γ10 + γ11Si
#
## Three level model (ff vs mm vs mixed sex)
#
# Level 1:
#   Yij= β0j + β1jXij+ eij
#   Level 2:
#     β0j= γ00 + γ01S1i+ γ02S2i+ u0j
#   β1i= γ10 + γ11S1i+ γ12S2i
# the fond is endlessly breaking.

# when there is level-2 varaible.
yes2model <- function(data,
                      gender_composition=2,
                      interaction=TRUE

){

  require(lme4)
  require(tidyverse)

    if(gender_composition==2){
      # convert dyadic structure into dyad2individual structure


      if(interaction){

      ## multilevel regression
      result <- lmer(y ~ 1 + con + catbtw + con:catbtw + (1|pid), data=data)
      }

      if(!interaction){

      result <- lmer(y ~ 1 + con + catbtw +  (1|pid), data=data)

      }

    }


    if(gender_composition==3){

      # Making effect coding variable
      # because there are three levels, we need two variable
      ## reference group = both -1


      if(interaction){
      ## multilevel regression
      result <- lmer(y ~ 1 + con + ev1 + ev2 + ev1:con + ev2:con
                             + (1|pid), data=data)
      }
      if(!interaction){
        result <- lmer(y ~ 1 + con + ev1 + ev2
                       + (1|pid), data=data)

      }

    }




  return(result)
}

