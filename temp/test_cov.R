n<-rnorm(9999,0,1)
m<-rnorm(9999,0,1)
cov(n,m)
cor(n,m)

bn<-999*n
bm<-999*m

cov(bn,bm) # cov becomes high after mulptiply by the same constant
cor(bn,bm)

sn<-0.00001*n
sm<-0.00001*m


cov(sn,sm) # cov becomes shrink dramatically after multiply by very small, less than 1
cor(sn,sm)

# correlation remains, becuase it is standardized, but covariance changes.

cov(nn,sm) # if the weight is different, cov goes down
         # if the weight is the same and is not under 1 and big enough, covary goe sup.

data<-data.frame(n,m,sn,sm,bn,bm)
data<-data %>%
  mutate (
     y = n+m,
     sy= sn+sm,
     by = bn+bm,
     wy = sn+bm,
     nm = n*m
  )
# covary varies
cov(n,m)
cov(sn,sm)
cov(sn,bm)
cov(bm,bn)

# cor never varies.
cor(n,m)
cor(sn,sm)
cor(sn,bm)
cor(bm,bn)
summary(lm(y~n + m + n:m, data)) ## No problem of interaction effect's type 1error rate.
summary(lm(sy~n + m + n:m, data))
summary(lm(wy~n + m  ,data))
summary(lm(wy~n + m  +n:m,data))# higher type 1 error rate of type 1 error rate.
summary(lm(wy~n ,data)) # non sig. becuase, small n
summary(lm(wy~m ,data)) # hightly sig, becuase small m.

cor(data)

# Now we begin to have imbalanced variance.

b<-rnorm(9999,0,999)
s<-rnorm(9999,0,1)
cov(b,s)
cor(b,s)

bb<-999*b
bs<-999*s

cov(bb,bs) # cov becomes high after mulptiply by the same constant
cor(bb,bs)

sb<-0.00001*b
ss<-0.00001*s

cov(sb,ss) # cov becomes small
cor(sb,ss)

cov(bb,ss)
cov(sb,bs)

data2<-data.frame(b,s,ss,sb,bb,bb)
data2<-data2 %>%
  mutate (
    y = s+b,
    sy= ss+sb,
    by = bs+bb,
    wy1 = ss+bb,
    wy2 = bs +sb,
    nm = n*m
  )

summary(lm(y~s+b,data2))
summary(lm(y~s+b+s:b,data2)) # got you baby
summary(lm(wy1~s+b+s:b,data2)) # criminal , more sig.
summary(lm(wy2~s+b+s:b,data2)) # not guilty, congrats. becase weights balanced you out.
summary(lm(y~ss+bb+ss:bb,data2))
summary(lm(y~bs+sb+bs:sb,data2))

# see the covariance of the interaction term and big one.

## okay we goes to the devil part
cat<-sample(rep(c(-1, 1), 500))
con<-rnorm(1000,0,1)

var(cat)
var(con) # SAME

dat_cat<-data.frame(cat,con)
dat_cat <- dat_cat %>%
  mutate(
    y = cat + con,
    y2 = 99999*cat + con,
    y3 = cat + 99999*con,
    y4= cat +2*con,
    catxcon=cat*con
  )

summary(lm(y~cat,dat_cat))
summary(lm(y~con,dat_cat))
summary(lm(y~cat+con,dat_cat)) # cat has more effect power, because..
summary(lm(y~cat+con+cat:con,dat_cat)) # not significant

summary(lm(y2~cat:con,dat_cat))

summary(lm(y2~cat,dat_cat))
summary(lm(y2~con,dat_cat))
summary(lm(y2~cat+con,dat_cat))
summary(lm(y2~cat+con+cat:con,dat_cat)) # not significant

summary(lm(y2~cat:con,dat_cat))

summary(lm(y3~cat,dat_cat))
summary(lm(y3~con,dat_cat))
summary(lm(y3~cat+con,dat_cat))
summary(lm(y3~cat+con+cat:con,dat_cat)) # not significant
summary(lm(y3~cat:con,dat_cat))

summary(lm(y4~cat,dat_cat))
summary(lm(y4~con,dat_cat))
summary(lm(y4~cat+con,dat_cat))
summary(lm(y4~cat+con+cat:con,dat_cat))
summary(lm(y4~cat:con,dat_cat))

cor(dat_cat)
cov(dat_cat)

##
#beta1,beta2 = 0.3, within = 100, between = 30

edu_w2<-rnorm(1000,0,1)
edu_w1<-rnorm(1000,0,1)
edu_b<-sample(rep(rnorm(500,0,1),each=2))
gender_w1<-sample(rep(c(-1, 1), 500))
gender_w2<-sample(rep(c(-1, 1), 500))
pid <- rep(1:50, each=2)
w <- 100
b <- 30

dat5<-data.frame(edu_w1,edu_w2,edu_b,gender_w1,gender_w2,w,b,pid)

dat5<-dat5 %>%
  mutate(
    gender_b= case_when(
        gender_w1==gender_w2 ~ 1,
        gender_w1!=gender_w2 ~ -1),
    edu_1=w*edu_w1 + b*edu_b,
    edu_2=w*edu_w2 + b*edu_b,
    gender_1=w*gender_w1 + b*gender_b,
    gender_2=w*gender_w1 + b*gender_b,
    y_1=0.3*edu_1 + 0.3*gender_1,
    y_2=0.3*edu_2 + 0.3*gender_2,
    edub_1=edu_b,
    edub_2=edu_b,
    genderb_1=gender_b,
    genderb_2=gender_b
  )

dat5<-dat5[,c(8,10:19)]
dat5$w <-NULL
dat5$b <- NULL
dat5$gender_b <- NULL

dat5 <- dat5 %>%
  dplyr::mutate(higher_one = case_when(
    y_1>y_2~"1",
    y_1<y_2~"2",
    y_1==y_2~"equal")
  )

## subset

# randomly assign tie siblings
p <- stats::rbinom(1,1,0.5)

data[data$higher_one=="equal",]<-ifelse(p==0,"1","2")

# splice the data

dat1<-dat5[dat5$higher_one=="1",]
dat2<-dat5[dat5$higher_one=="2",]

# change the column names
dat5 <- dat5 %>%
  dplyr::mutate(y_2_1=y_2,
                y_1_2=y_1,
                edu_1_2=edu_1,
                edu_2_1=edu_2,
                gender_1_2=gender_1,
                gender_2_1=gender_2
  ) %>%
  subset(
    select=-c(y_1,
              y_2,
              edu_1,
              edu_2,
              gender_1,
              gender_2

    )
  ) %>%
  dplyr::rename(y_1=y_2_1,
                y_2=y_1_2,
                edu_2=edu_1_2,
                edu_1=edu_2_1,
                gender_2=gender_1_2,
                gender_1=gender_2_1)

dat3<-rbind(dat1,dat2)

# calculate mean score and difference score

dat3<-dat3 %>%
  dplyr::mutate(
    edu_mean = (edu_1+edu_2)/2,
    edu_diff = edu_1-edu_2,
    y_mean = (y_1+y_2)/2,
    y_diff = y_1-y_2
  )

x<-rnorm(999,0,1)
y<-rnorm(999,0,1)
var(y)
var(x-y)
var(3*x)
dat3$higher_one<-NULL
var(dat3$edu_diff)
var(dat3$genderb_1)

dat3$edu_diffxgenderb_1 <-dat3$edu_diff * dat3$genderb_1

var(dat3$edu_diffxgenderb_1)



summary(lm(y_diff ~ edu_diff ,dat3)) # highly sig
summary(lm(y_diff ~ genderb_1,dat3)) # not significant because variance is low.
summary(lm(y_diff ~ edu_diff + genderb_1,dat3)) # still not sig.
summary(lm(y_diff ~ edu_diff + genderb_1 + edu_diffxgenderb_1,dat3)) # gender become less sig.
#


cor(dat3)
cov(dat3)
