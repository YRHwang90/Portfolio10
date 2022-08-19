poor<-gen_dyadic()
re_poor<-reconstruct(poor)
ind_poor<-dyad2ind(poor)

poor<-poor[,grep("cat",names(poor))]

re_poor<-re_poor[,grep("cat|eff|ev",names(re_poor))]

ind_poor<-ind_poor[,grep("cat|ev|pid",names(ind_poor))]
