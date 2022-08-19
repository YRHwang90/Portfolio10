## the function to run null model

model_null <- function(data){
  lmer(y~1+(1|pid), data=data)
}
