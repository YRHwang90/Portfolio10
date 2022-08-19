
ac_1 <-c(0.1, 0.3, 0.03, 0.03)
ac_2 <-c(0.2, 0.4, 0.1, 0.008)
ac_3 <-c(0.8, 0.043, 0.7, 0.01)
ac_4 <-c(0.2, 0.73, 0.1, 0.1)
c_2<-c(1,2,5,23)
check_1<-c(0.01, 0.902,0,0.07)
check_2<-c(0.03, 0.042,0.002,0.00001)
check_3<-c(0.01,0,0.5,0.001)
check_4<-c(0.001, 0.042,0.02,0.2)
id<-1:4

df<-data.frame(id,ac_1, ac_2,ac_3,ac_4,c_2,check_1,check_2,check_3,check_4)
pivot_longer(df, -c(id,c_2)) %>%
  separate(name,into=c("type", "pos")) %>%
  pivot_wider(names_from=type, values_from = value) %>%
  mutate(ac=if_else(near(check,0.02), as.double(NA), ac)) %>%
  pivot_wider(names_from = pos, values_from = ac:check)

df
results_1 <- results_1[,grep("hasConverged",names(results_1))]
df
results


x <- c(1, 2, 3, -5, NA)

mean(x > 0,na.rm=TRUE)
