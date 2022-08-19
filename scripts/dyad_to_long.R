## function to convert dyad structure into individual structure

dyad2ind<-function(data,
                   pid=pid){ #pid= pair id

  data <- data %>%
    mutate(
      ev1_1 = case_when(catwtn_1==catwtn_2 & catwtn_1==1 ~ 1,
                        catwtn_1==catwtn_2 & catwtn_1==-1 ~ 0, #reference group -1
                        catwtn_1!=catwtn_2 ~ -1),
      ev2_1 = case_when(catwtn_1==catwtn_2 & catwtn_1==1 ~ 0,
                        catwtn_1==catwtn_2 & catwtn_1==-1 ~ 1, # reference group -1
                        catwtn_1!=catwtn_2 ~ -1),
      ev1_2 = ev1_1,
      ev2_2 = ev2_1
    )


  dat_ind <-data %>%
    gather(key,value,-pid) %>%
    separate(key,c("key","indID"),sep="_") %>%
    spread(key,value)


  return(dat_ind)}
