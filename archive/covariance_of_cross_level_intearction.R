library(tidyverse)

sim_data<-function(cn=1000){

  require(tidyverse)

v1btw_1 <- rnorm(cn,0,1)
v1wtn_1 <- rnorm(cn,0,1)
v1wtn_2 <- rnorm(cn,0,1)

v1_1 <- v1btw_1 + v1wtn_1
v1_2 <- v1btw_1 + v1wtn_2

v2wtn_1 <- sample(rep(c(-1, 1),cn/2))
v2wtn_2 <- sample(rep(c(-1,1),cn/2))

data<-data.frame(v1_1,v1_2,v1btw_1,v1wtn_1,v1wtn_2,v2wtn_1,v2wtn_2)

data <- data %>%
  mutate(interaction_1= 99*v1btw_1*v2wtn_1,
         interaction_2= 99*v1btw_1*v2wtn_2,
         v1btw_2=v1btw_1,
         pid = 1:cn,
         y_1=v1_1+v2wtn_1+interaction_1,
         y_2=v1_2+v2wtn_2+ interaction_2)

return(data)
}

# make empty list

nrep=999
cov <- vector("list",length=nrep)

for(i in 1:nrep){


   sim<-sim_data(cn=1000)
   covp<-sim[,c('interaction_1','interaction_2')] %>%
     cov()
   cov[[i]]<-covp['interaction_1','interaction_2']}

df <- data.frame(matrix(unlist(cov), nrow=length(cov), byrow=TRUE))
colnames(df) <- c('cov')
df$cov2 <-as.numeric(df$cov)
# mean of these icc CHECK
meancov2 <- mean(df$cov2)


# Conclusion: as cross-class interaction coefficient goes up,
# the covariance between the members goes up althoung
# soley see the interaction term, the corr bettween interaction_1 and interaction_2 seems to be negligible.

