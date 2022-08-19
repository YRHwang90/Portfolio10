# overview

# This documents is for getting and reporting statistics in my results section.
# (Power subsection)

# power of parameters in each models.
# 0.05 alpha level.




# power

## data management

sig_rate<-read.csv("type1error/sig_rate_with_conditions.csv")

sig_rate$beta2_factor<-as.factor(sig_rate$beta2)
sig_rate$beta3_factor<-as.factor(sig_rate$beta3)
sig_rate$ICC_factor<-as.factor(sig_rate$ICC)
sig_rate$CN_factor <- as.factor(sig_rate$CN)

sig_rate$X<-NULL
sig_rate$X.1<-NULL

sig_rate_long<-sig_rate%>%
  pivot_longer( cols=-c("beta1","condition","beta2","beta3","ICC","CN",
                        "beta2_factor","beta3_factor","ICC_factor","CN_factor"),
                names_to = c('coef', 'model'),
                names_sep = '(?<=_p_)')

sig_rate_long <- sig_rate_long %>%
  mutate(
    correspond = case_when(
      coef=="coef_con_diff_p_" ~ "edu_within",
      coef=="coef_con_diffxev1_p_" ~ "interaction",
      coef=="coef_con_diffxev2_p_" ~ "interaction",
      coef=="coef_con_diffxgender_composition_two_eff_p_" ~ "interaction",
      coef=="coef_con_mean_p_" ~ "edu_between",
      coef=="coef_ev1_p_" ~ "gender_between",
      coef=="coef_ev2_p_" ~ "gender_between",
      coef=="coef_fixed_catbtw_p_" ~ "gender_between",
      coef=="coef_fixed_catwtn_p_" ~ "gender_within",
      coef=="coef_fixed_catwtnxcon_p_"  ~ "interaction",
      coef=="coef_fixed_con_p_" ~ "edu_within",
      coef=="coef_fixed_conxcatbtw_p_" ~ "interaction",
      coef=="coef_fixed_conxev1_p_" ~ "interaction",
      coef=="coef_fixed_conxev2_p_" ~ "interaction",
      coef=="coef_fixed_ev1_p_" ~ "gender_between",
      coef=="coef_fixed_ev2_p_" ~ "gender_between",
      coef=="coef_fixed_Intercept_p_" ~ "intercept",
      coef=="coef_gender_composition_two_eff_p_"~"gender_between",
      coef=="coef_intercept_p_" ~ "intercept",
      coef=="coef_y_mean_p_" ~ "health_between",
      coef=="f_p_" ~ "F_test_dyadic"
    )
  )


sig_rate_long<-sig_rate_long %>%
  mutate(bi_model=case_when(
    model=='sib0' ~ 'REG',
    model=='sib1' ~ 'REG',
    model=='sib2' ~ 'REG',
    model=='sib3' ~ 'REG',
    model=='MLM0' ~ 'MLM',
    model=='MLM1' ~ 'MLM',
    model=='MLM2' ~ 'MLM',
    model=='MLM3' ~ 'MLM',
    model=='MLM4' ~ 'MLM',
    model=='MLM5' ~ 'MLM',
    model=='MLM6' ~ 'MLM',
  ))


sig_rate_long_edu_between<- sig_rate_long %>%
  subset(correspond=="edu_between")
sig_rate_long_edu_within<- sig_rate_long %>%
  subset(correspond=="edu_within")
sig_rate_long_F_test_dyadic<- sig_rate_long %>%
  subset(correspond=="F_test_dyadic")
sig_rate_long_gender_between<- sig_rate_long %>%
  subset(correspond=="gender_between")
sig_rate_long_gender_within<- sig_rate_long %>%
  subset(correspond=="gender_within")
sig_rate_long_health_between<- sig_rate_long %>%
  subset(correspond=="health_between")
sig_rate_long_interaction<- sig_rate_long %>%
  subset(correspond=="interaction")
sig_rate_long_intercept<- sig_rate_long %>%
  subset(correspond=="intercept")
# education's effect
sig_rate_long_edu_between<- sig_rate_long %>%
  subset(correspond=="edu_between")

sig_rate_long_edu_within<- sig_rate_long %>%
  subset(correspond=="edu_within")

## get number (education within)
edu_within_cn<-sig_rate_long_edu_within %>%
  group_by(CN_factor,bi_model) %>%
  summarize(Mean=mean(value))

print(edu_within_cn, digits=2)

sig_rate_long_edu_within %>%
  group_by(model) %>%
  summarize(mean = mean(value))

print ( sig_rate_long_edu_within %>%
  group_by(model) %>%
  summarize(mean = mean(value)), digits=2)

# not important
sig_rate_long_F_test_dyadic<- sig_rate_long %>%
  subset(correspond=="F_test_dyadic")

## power of gender
sig_rate_long_gender_between<- sig_rate_long %>%
  subset(correspond=="gender_between" & beta2!=0)
sig_rate_long_gender_within<- sig_rate_long %>%
  subset(correspond=="gender_within" & beta2!=0)

num_gen_bt<-sig_rate_long_gender_between%>%
  group_by(CN_factor,bi_model) %>%
  summarize(Mean=mean(value))

print(num_gen_bt,digits=2)

a<-sig_rate_long_gender_between %>%
  group_by(model) %>%
  summarize(mean = mean(value))

num_gen_w<-sig_rate_long_gender_within%>%
  group_by(CN_factor,bi_model) %>%
  summarize(Mean=mean(value))

print(num_gen_w,digits=2)

sig_rate_long_gender_within %>%
  group_by(model) %>%
  summarize(mean = mean(value))









## DV
sig_rate_long_health_between<- sig_rate_long %>%
  subset(correspond=="health_between")

##power of interaction
sig_rate_long_interaction<- sig_rate_long %>%
  subset(correspond=="interaction" & beta3!=0)

num_int<-sig_rate_long_interaction %>%
  group_by(CN_factor,bi_model) %>%
  summarize(Mean=mean(value))
print(num_int,digits=3)
print(sig_rate_long_interaction%>%
  group_by(model) %>%
  summarize(mean = mean(value)), digits=2)

a2<-sig_rate_long_interaction%>%
  group_by(model) %>%
  summarize(mean = mean(value))

b<-sig_rate_long_interaction %>%
  filter(model=="MLM4")

#power of gender composition variable (With two categories)

sig_rate_long_gender_between<- sig_rate_long %>%
  subset(correspond=="gender_between" & beta2!=0)


sig_rate_long_gender_betwee<-sig_rate_long_gender_between %>%
  subset(model!= 'MLM3' & model!= 'MLM4' & model!='sib0' & model!='sib2')
print(
gednder_betwen_table <- sig_rate_long_gender_betwee%>%
  group_by(coef,model) %>%
  summarize(Mean=mean(value)),digits=2)

