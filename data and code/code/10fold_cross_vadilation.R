gd<-read_rds("1080gne_705disease_01matrix.rds")
en<-function(x){rownames(as.data.frame(which(x==1)))}
gdf<-gd[,which(apply(gd,2,sum) >1)] #chose the disease containing more than one gene
nn<-apply(gd,2,sum)
n<-mean(nn[nn>0]) %>% ceiling() #259
gdl<-apply(gdf,2,en)
Random_choose<-function(x){
  if(length(x)>1) {
    y<-sample(unlist(x), ceiling(length(unlist(x))*0.1))
  }else{
    y<-x
  }
  return(y)
} # true random
convert<-function(x,test){
  for(i in 1:length(test))
    x[,names(test[i])][as.vector(test[[i]])]=0
  return(x)}
#######GD list
TEST<-list()
for (i in 1:n) {
  TEST[[i]]<-sapply(gdl,Random_choose)
}
############### GD matrix
dglist<-list()
for (i in 1:n) {
  dglist[[i]]<-convert(x=gd,test=TEST[[i]])
}
############## matrix rwr 
rg<-read_rds("filter2765site_1080gene_01matrix.rds") %>% t
gene_rna<-apply(rg,1,en)
ma<-function(x,y){M<-match(x,names(y))
s<-y[M]        
return(s)}
############## test list
l_10<-list()
for (i in 1:n){
l_10[[i]]<-sapply(TEST[[i]],ma,y=gene_rna)}
############################


##########################
library(readr)
gg<-read_rds("1080gene-coexpression_0weight.rds")
rr<-read_rds("2765rna_coexpression_0weight.rds")
rg<-read_rds("filter2765site_1080gene_01matrix.rds")
dd<-read_rds("705disease-0weight-network.rds")
#################################################
library(dendextend)
GG<- Normalize(gg)
DD<- Normalize (dd)
RR<-Normalize(rr)
RG<-Normalize(rg)
GR<- t(rg) %>% Normalize
RD<-matrix(0,nrow(rr),ncol(dd))
colnames(RD)<-colnames(dd)
rownames(RD)<-rownames(rr)
DR<-t(RD)
########
DG<-list()
for(i in 1:n){
  DG[[i]]<-t(dglist[[i]]) %>% Normalize}
GD<-list()
for (i in 1:n) {
GD[[i]]<- Normalize(dglist[[i]])}
##########
#GG,GR,GD,DD,DG,DR,RR,RG,RD
w<-list()
for (i in 1:n){
w[[i]]<-combine_no(GG,GR,GD[[i]],DD,DG[[i]],DR,RR,RG,RD)}
#this is new_normiazsed obsreved result
drname<-function(rwr,w,rr,dd){
  colnames(rwr)<-colnames(w)
  rownames(rwr)<-rownames(w)
  dr<-rwr[rownames(rr),colnames(dd)]
  return(dr)
}
rwr_function<-function(w){
nr<-nrow(w)
rwr<-matrix(0,nr,nr)
p0<-diag(nr)
for (i in 1:nr){
  f0<-p0[,i]
  rwr[,i]<-RWR(f0,w,0.75)
}
return(rwr)}
### 8 RWR results
disease_site<-list()
for (i in 1:length(w)){
disease_site[[i]]<-rwr_function(w[[i]])
}
##
###names!
dr_list<-list()
for (i in 1:length(disease_site)){
  dr_list[[i]]<-drname(disease_site[[i]],w[[i]],rr,dd) 
}

############ result
library(pROC)
dmatch_rwr<-function(result,test){
  v<-list()
  for (i in 1: ncol(result)){
    v[[i]]<-ifelse(is.na(match(result[,i],unlist(test[[i]])))=="TRUE",FALSE,TRUE)
  }
  return(v)}
exteact_rwr<-function(drfe,drf){ #length(result)=ncol(ttf)
  a<-list()
  for (i in 1:ncol(drfe)){
    a[[i]]<-as.vector(drf[drfe[,i],i])}
  return(a)
}
#test:l_10
#rwr:dr_list
drf=drfe=dmr=er=y=AUC=list()
auc<-vector("numeric",length(l_10))
for (i in 1:length(l_10)){
  drf[[i]]<-dr_list[[i]][,names(l_10[[i]])]
  drfe[[i]]<-extract_max(drf[[i]],nrow(drf[[i]]))
  dmr[[i]]<-dmatch_rwr(drfe[[i]],l_10[[i]])  %>% unlist()
  er[[i]]<-exteact_rwr(drfe[[i]],drf[[i]]) %>% unlist()
  y[[i]]<- prediction(er[[i]],as.numeric(dmr[[i]]))
  AUC[[i]]<- performance(y[[i]], measure = "auc")
  auc[i] <- AUC[[i]]@y.values[[1]]}
auc 