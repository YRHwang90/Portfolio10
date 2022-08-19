# ICC calculator
# based on the current data generation equation
# put this in here to test a data generation


ICC_calculator <- function(
  beta1,
  beta2,
  beta3,
  ICC,

  # if it is ratio, it is determined by the numbers you put in the within argument
  within)
{

  a<-beta1
  b<-beta2
  c<-beta3
  d<-ICC
  y<-within

  x <- (1 - d * (a^2 * y + b^2 * y + c^2 * y + 2))/((d - 1) * (a^2 + b^2))

 return(x)
}

