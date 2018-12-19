
########RUN Random walk with restart 
library(readr)
library(dendextend)
Normalize <- function(x) { 
  y<-sweep(x, 2, apply(x, 2, sum), "/") 
  y[,which(apply(x,2,sum)==0)]<-0
  return(y)
}
condition<-function(x,y){
  c<-length(which((x-y)<=1e-10))
  return(c)
}
combine_no<-function(GG,GR,GD,DD,DG,DR,RR,RG,RD){ 
  #Normalize for the multiple network
  Normalize <- function(x) { 
    y<-sweep(x, 2, apply(x, 2, sum), "/") 
    y[,which(apply(x,2,sum)==0)]<-0
    return(y)
  }
  first<-rbind(DD,GD,RD) %>%   Normalize
  second<-rbind(DG,GG,RG) %>%   Normalize
  third<-rbind(DR,GR,RR) %>%   Normalize
  fs<-cbind(first,second)
  fst<-cbind(fs,third)
  return(fst)}
first_max<-function(x,n){
  max<-apply(x,2, function(x){names(x[order(x,decreasing=TRUE)[1:n]])})
  return(max)
}
##########################################################
library(readr)
gg<-read_rds("1080gene-coexpression_0weight.rds")
rr<-read_rds("2827_2827_rr.rds")
rg<-read_rds("2827_1080_rg.rds")
x<-read_rds("1080gne_705disease_01matrix.rds")
gd<-read_rds("T10%_1080gne_705disease_01matrix.rds")
max(gd) #[1] 1
dd<-read_rds("705disease-0weight-network.rds")
###################
GG<- Normalize(gg)
GD<- Normalize(gd)
DD<- Normalize (dd)
DG<-t(gd) %>% Normalize
RR<-Normalize(rr)
RG<-Normalize(rg)
GR<- t(rg) %>% Normalize
RD<-matrix(0,nrow(rr),ncol(dd))
colnames(RD)<-colnames(dd)
rownames(RD)<-rownames(rr)
DR<-t(RD)
#GG,GR,GD,DD,DG,DR,RR,RG,RD
w<-combine_no(GG,GR,GD,DD,DG,DR,RR,RG,RD) 
#this is new_normiazsed obsreved result
n<-nrow(w)
rwr<-matrix(0,n,n)
p0<-diag(n)
for (i in 1:n){
  f0<-p0[,i]
  rwr[,i]<-RWR(f0,w,0.75)
}
dim(rwr) #[1] 4550 4550
colnames(rwr)<-colnames(w)
rownames(rwr)<-rownames(w)
dr<-rwr[rownames(rr),colnames(dd)] 
dim(dr) #[1] 2765  705
write_rds(dr,"2827_705_dr.rds")

#############################Run Hypergeometric distribution test
library(readr)
#conver gene-rna 0/1 to list
gr<-read_rds("filter2765site_1080gene_01matrix.rds") %>% t
en<-function(x){rownames(as.data.frame(which(x==1)))}
rg<-apply(gr,2,en)
dr<-read_rds("1080gne_705disease_01matrix.rds") %>% t
dg<-apply(dr,1,en)
####total match gene number
length(unique(c(unlist(rg),unlist(dg))))
#[1] 1080
test<-function(rg,dg){
  m<-matrix(0,length(rg),length(dg))
  b<-vector("numeric",length(rg))
  c=d=c=t=b
  for( I in 1:length(dg)){
    for( i in 1:length(rg)){
      b[i]<-length(rg[[i]])
      d[i]<-length(dg[[I]])
      c[i]<-56-b[i]
      m[i,I]<-phyper(length(na.omit(match(rg[[i]],dg[[I]]))),b[i],c[i],d[i])}
  }
  return(m)}
tt<-test(rg,dg) #dim(tt) [1] 488 392
colnames(tt)<-names(dg)
rownames(tt)<-names(rg)
write_rds(tt,"hyhypergeometric_pvalue_result.rds")
