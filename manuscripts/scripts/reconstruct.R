## This function is to reconstruct dyadic data structure into dyadic data structure
## it makes indistinguishable dyads into distinguishable dyads
## for more details about discordnat kinship model, see Garrison & Rodgers (2016)
# to dyadic structure but more suitable for
# Mason's discordant kinship model

reconstruct<-function(data = data


){
  require(tidyverse)
  # make guidance column

  data <- data %>%
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

  dat1<-data[data$higher_one=="1",]
  dat2<-data[data$higher_one=="2",]

  # change the column names
  dat2 <- dat2 %>%
    dplyr::mutate(y_2_1=y_2,
                  y_1_2=y_1,
                  con_1_2=con_1,
                  con_2_1=con_2,
                  cat_1_2=cat_1,
                  cat_2_1=cat_2,
                  catwtn_1_2=catwtn_1,
                  catwtn_2_1=catwtn_2
    ) %>%
    subset(
      select=-c(y_1,
                y_2,
                con_1,
                con_2,
                cat_1,
                cat_2,
                catwtn_1,
                catwtn_2
      )
    ) %>%
    dplyr::rename(y_1=y_2_1,
                  y_2=y_1_2,
                  con_2=con_1_2,
                  con_1=con_2_1,
                  cat_2=cat_1_2,
                  cat_1=cat_2_1,
                  catwtn_1=catwtn_2_1,
                  catwtn_2=catwtn_1_2)

  dat3<-rbind(dat1,dat2)

  # calculate mean score and difference score

  dat3<-dat3 %>%
    dplyr::mutate(
      con_mean = (con_1+con_2)/2,
      con_diff = con_1-con_2,
      y_mean = (y_1+y_2)/2,
      y_diff = y_1-y_2
    )

  dat3$higher_one <- NULL

  dat3<-dat3 %>%
    mutate(
      gender_composition_two=case_when(
        catwtn_1==catwtn_2 ~ "same-sex",
        catwtn_1!=catwtn_2 ~ "mixed-sex"),

      gender_composition_three=case_when(
        catwtn_1==1 & catwtn_2==1 ~ "oneone",
        catwtn_1==-1 & catwtn_2==-1 ~ "mimi",
        catwtn_1!=catwtn_2 ~ "mixed"),

      gender_composition_two_eff=case_when(
        catwtn_1==catwtn_2 ~ 1,
        catwtn_1!=catwtn_2 ~ -1),
      # three levels : need two coding variable.
      ev1 = case_when(catwtn_1==catwtn_2 & catwtn_1==1 ~ 1,
                      catwtn_1==catwtn_2 & catwtn_1==-1 ~ 0,
                      catwtn_1!=catwtn_2 ~ -1), # reference
      ev2 = case_when(catwtn_1==catwtn_2 & catwtn_1==1 ~ 0,
                      catwtn_1==catwtn_2 & catwtn_1==-1 ~ 1,
                      catwtn_1!=catwtn_2 ~ -1) # reference
    )



return(dat3)





}
