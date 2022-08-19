type  1 error rate sig/nonsig rate.


results_1<-read.csv("model_results/results_1.csv")
results_1<-results_1[,grep("_p_",names(results_1))]
non_sig_rate<-data.frame(results_1[0,])

nMLM<-6 # how many MLM models we have? 0 to 6

for(i in 1:ncond){

  # read the data
  results <- read.csv(paste0("model_results/results_",i,".csv"))

  # exclude non-convergence model

  for(k in 0:nMLM){
    results[
      ## rows to replace
      results[[paste0("hasConverged_MLM", k)]]==0,
      ## column to replace
      grep(paste0("MLM",k),names(results))
    ] <- NA
  }
  # These two lines generates error
  # results <- results %>%
  #   filter(has.)

  # mean() You can use mean() to get the proportion of TRUE of a logical vector.

  non_sig_rate[i,]<-  results %>%
    summarise(across(contains("_p_"), ~ mean(.x > 0.05, na.rm=TRUE)))
  # type1error$conditions[i] <- i
  non_sig_rate[i] <- i

}

hey<-results %>%
  summarise(across(contains("_p_"), ~ mean(.x > 0.05, na.rm=TRUE)))

hey1<-non_sig_rate[3,]
rownames(hey1)<-NULL

identical(hey,hey1)
