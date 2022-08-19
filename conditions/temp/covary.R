x1<-rnorm(999,0,1)
x2<-rnorm(999,0,1)
y <- x1+x2
iv1<-999999*x1
iv2<-999999*x2
yv<-iv1+iv2
x1x2<-x1*x2
iv1iv2<-iv1*iv2
data<-data.frame(x1,x2,y,yv,iv1,iv2,x1x2,iv1iv2)

cov(x1,x2) # nearly 0
cor(x1,x2) # nearly 0
cov(iv1,iv2) # very big,
cor(iv1,iv2) # nearly 0

summary(lm(y~x1+x2+x1*x2)) # interaction p=0.11

summary(lm(y~iv1+iv2+iv1:iv2)) #interaction significant.

summary(lm(yv~x1+x2+x1*x2)) # interaction non sig

summary(lm(yv~iv1+iv2+iv1:iv2)) #interaction non-sig

summary(lm(yv~iv1+x2+iv1:x2)) #interaction non-sig


cor(data)
cov(data)
