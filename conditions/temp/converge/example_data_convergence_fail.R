library(lme4)
read.table(textConnection("duration season  sites   effect
                          4d    mon s1  7305.91
                          4d    mon s2  856.297
                          4d    mon s3  649.93
                          4d    mon s1  10121.62
                          4d    mon s2  5137.85
                          4d    mon s3  3059.89
                          4d    mon s1  5384.3
                          4d    mon s2  5014.66
                          4d    mon s3  3378.15
                          4d    post    s1  6475.53
                          4d    post    s2  2923.15
                          4d    post    s3  554.05
                          4d    post    s1  7590.8
                          4d    post    s2  3888.01
                          4d    post    s3  600.07
                          4d    post    s1  6717.63
                          4d    post    s2  1542.93
                          4d    post    s3  1001.4
                          4d    pre s1  9290.84
                          4d    pre s2  2199.05
                          4d    pre s3  1149.99
                          4d    pre s1  5864.29
                          4d    pre s2  4847.92
                          4d    pre s3  4172.71
                          4d    pre s1  8419.88
                          4d    pre s2  685.18
                          4d    pre s3  4133.15
                          7d    mon s1  11129.86
                          7d    mon s2  1492.36
                          7d    mon s3  1375
                          7d    mon s1  10927.16
                          7d    mon s2  8131.14
                          7d    mon s3  9610.08
                          7d    mon s1  13732.55
                          7d    mon s2  13314.01
                          7d    mon s3  4075.65
                          7d    post    s1  11770.79
                          7d    post    s2  4254.88
                          7d    post    s3  753.2
                          7d    post    s1  11324.95
                          7d    post    s2  5133.76
                          7d    post    s3  2156.2
                          7d    post    s1  12103.76
                          7d    post    s2  3143.72
                          7d    post    s3  2603.23
                          7d    pre s1  13928.88
                          7d    pre s2  3208.28
                          7d    pre s3  8015.04
                          7d    pre s1  11851.47
                          7d    pre s2  6815.31
                          7d    pre s3  8478.77
                          7d    pre s1  13600.48
                          7d    pre s2  1219.46
                          7d    pre s3  6987.5
                          "),header=T)->dat1



#myTryCatch() function

myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  if(is.null(warn)){
    warn<-'null'
  }
  if(is.null(err)){
    err<-'null'
  }
  list(value=value, warning=warn, error=err)
}
##


hey<-function(effect=effect){lmer(effect ~ duration + (1+duration|sites) +(1+duration|season),
                                  data=dat1)

}

myTryCatch(m1<-lmer(effect ~ duration + (1+duration|sites) +(1+duration|season),
                           data=dat1))

myTryCatch(m1)

# eventhough the model has been failed to converge, it gives standard error estimation
summary(m1)[['coefficients']]['(Intercept)','Std. Error']
summary(m1)[['coefficients']]['duration7d','Std. Error']
summary(m1)[['coefficients']]

results2<-as.data.frame(matrix(data=NA, nrow=10, ncol=1))

# This seems to work fine.
hello<-unlist(myTryCatch(sib4<-lmer(effect ~ duration + (1+duration|sites) +(1+duration|season),
                                                   data=dat1))$warning)

unlist(myTryCatch(sib5<-sib4(data=dat1))$warning)


results2<-as.data.frame(matrix(data=NA, nrow=10, ncol=1))

results2[1,1]<-tryc(unlist(myTryCatch(sib5<-sib4(data=dat1)$warning)))


sib4<-function(data=dat1){
  lmer(effect ~ duration + (1+duration|sites) +(1+duration|season),
       data=dat1)

}
results2[1,1]<-isSingular(sib4)

results[2,1]<-unlist(myTryCatch(lmer(effect ~ (1|sites),
                                                   data=dat1))$error)


results[2,2]<-unlist(myTryCatch(lmer(effect ~ (1|sites),
                                     data=dat1))$warning)

hey<-myTryCatch(lmer(effect ~ (1|sites),
                     data=dat1)$warning)

lmer(effect ~ duration + (1+duration|sites) +(1+duration|season),
          data=dat1) %>%
  summary()

m1%>%
  summary()

error<-tryCatch(lmer(effect ~ duration + (1+duration|sites) +(1+duration|season),
              data=dat1),
         error=function(e) e$message)

lmer(effect ~ duration + (1+duration|sites) +(1+duration|season),
     data=dat1)
summary(sib5)


ca<-tryCatch(m_lme4 <- summary(lmer(effect ~ duration + (1+duration|sites) +(1+duration|season),
                                data=dat1)), error = function(e) e)

hasConverged <- function (mm) {

  if ( !inherits(mm, "merMod")) stop("Error: must pass a lmerMod object")

  retval <- NULL

  if(is.null(unlist(mm@optinfo$conv$lme4))) {
    retval = 1
  }
  else {
    if (isSingular(mm)) {
      retval = 0
    } else {
      retval = -1
    }
  }
  return(retval)
}
## Source of the function
## https://stackoverflow.com/questions/72090177/how-can-i-know-whether-the-model-is-converged-or-failed-to-converge-in-lme4-with
hasConverged(lmer(effect ~ duration + (1+duration|sites) +(1+duration|season),
                                     data=dat1))
#
# which returns 1 if the model converged normally
# ie not to a singular fit, 0 if it converges to a
# singular fit and -1 if it fails to converge.
# Another approach is to promote
# the warnings to errors as per the comment by @SamR:


ML<-lmer(effect ~ duration + (1+duration|sites) +(1+duration|season),
     data=dat1,REML=FALSE)

