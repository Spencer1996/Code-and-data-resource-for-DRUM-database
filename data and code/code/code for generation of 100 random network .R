library(readr)
library(dendextend)
###############fucntions
###I change it !!!
Normalize <- function(x) { 
  y<-sweep(x, 2, apply(x, 2, sum), "/") 
  y[,which(apply(x,2,sum)==0)]<-0
  return(y)
}
#
# I change it!!!!
exchange<-function(newg,gg){newg[which(newg!=0)]<-gg[which(gg!=0)]
return(newg)}
#
combine<-function(GG,GR,GD,DD,DG,DR,RR,RG,RD){ 
  first<-cbind(DD,DG,DR)
  dim(first)
  second<-cbind(GD,GG,GR)
  third<-cbind(RD,RG,RR)
  fs<-rbind(first,second)
  fst<-rbind(fs,third)
  return(fst)}

#
rwr<-function(m,newd,newr){
  RWR<-function(f0,w,r){
    n<-ncol(w)
    fs0<-(1-r)*w%*%f0+(r*f0)
    fs0<-(1-r)*w%*%fs0+(r*f0)
    fs1<-fs0
    fs1<-(1-r)*w%*%fs1+(r*f0)
    while(condition(fs1,fs0)<n){
      fs0<-(1-r)*w%*%fs0+(r*f0)
      fs1<-(1-r)*w%*%fs1+(r*f0)
    }
    return(fs1)
  }
  n<-nrow(m)
  u<-nrow(newd)
  sub_rwr<-matrix(0,n,u)
  p0<-diag(u)%>%rbind(matrix(0,n-u,u)) #[1] 936 392
  for (I in 1:u){
    f0<-p0[,I]
    sub_rwr[,I]<-RWR(f0,m,0.75)
  }
  colnames(sub_rwr)<-colnames(newd)
  rownames(sub_rwr)<-rownames(m)
  sub<-sub_rwr[rownames(newr),]
  return(sub)
}
#
condition<-function(x,y){
  c<-length(which((x-y)<=1e-10))
  return(c)
}
##
RWR<-function(f0,w,r){
  n<-ncol(w)
  fs0<-(1-r)*w%*%f0+(r*f0)
  fs0<-(1-r)*w%*%fs0+(r*f0)
  fs1<-fs0
  fs1<-(1-r)*w%*%fs1+(r*f0)
  while(condition(fs1,fs0)<n){
    fs0<-(1-r)*w%*%fs0+(r*f0)
    fs1<-(1-r)*w%*%fs1+(r*f0)
  }
  return(fs1)
}
##
first10max<-function(x){ 
  max<-apply(x,2, function(x){names(x[order(x,decreasing=TRUE)[1:10]])})
  return(max)
}
##
match_list<-function(x,y){
  z<-x
  for(I in 1: length(y)){
    for (i in 1:length(x[[I]])) {
      z[[I]][i]<-length(na.omit(match(x[[I]][i],y[[I]])))}
  }
  return(z)
} 
#
mean1<-function(x){z<-x/sum(x)
return(z)}
#
###########input data
gg<-read_rds("1080gene-coexpression_0weight.rds")
gg1<-ifelse(gg>0,1,0)
rr<-read_rds("2765rna_coexpression_0weight.rds")
rr1<-ifelse(rr>0,1,0)
rg<-read_rds("filter2765site_1080gene_01matrix.rds")
gd<-read_rds("1080gne_705disease_01matrix.rds")
max(gd) #[1] 1
gdd<-ifelse(gd>0,1,0)
max(gdd)
dd<-read_rds("705disease-0weight-network.rds")
dd1<-ifelse(dd>0,1,0)
rd<-matrix(0,nrow(rr),ncol(dd))
colnames(rd)<-colnames(dd)
rownames(rd)<-rownames(rr)
ob<-read_rds("first10_RWR_matrix.rds")
#function:
extraction<-function(x){ 
  y<-apply(x,2, function(x){names(x[which(x>0)])})
  return(y)
}
#########################################################################################################################################
apply_exn<-function(x,n){
  exn<-function(x,n){
    if (length(which(x>0))<=n) {y<-names(x[which(x>0)])
    }else{y<-names(x[order(x,decreasing=TRUE)[1:n]])}
    return(y)}
  r<-list()
  for(i in 1:ncol(x)){
    r[[i]]<-exn(x[,i],n)
  } 
  names(r)<-colnames(x)
  return(r)
}
########### for lsit
match_list<-function(x,y){
  z<-x
  for(I in 1: length(y)){
    for (i in 1:length(x[[I]])) {
      z[[I]][i]<-length(na.omit(match(x[[I]][i],y[[I]])))}
  }
  return(z)
} 
library(dendextend)
trans<-function(x){t<-sapply(x,FUN=as.numeric)%>%  
  sapply(FUN=as.matrix) %>% sapply(FUN=t)
return(t)}
###################################################################################################################################################
extract_max<-function(x,n){  ###for matrix
  max<-apply(x,2, function(x){names(x[order(x,decreasing=TRUE)[1:n]])})
  return(max)
}

###
##functin match observe result to the random data
ml<-function(x,y){
  z<-matrix(0,nrow(x),ncol(x)) 
  for(I in 1: ncol(x)){
    for (i in 1:nrow(x)) {
      z[i,I]<-length(na.omit(mapply(match,x[i,I],y[,I])))}
  }
  colnames(z)<-colnames(y)
  return(z)
}
###
# random
for(i in 1:100){
  newg<-BiRewire::birewire.rewire.bipartite(gg1,verbose=FALSE) %>% exchange(gg=gg)  %>% mean1 
  newr<-BiRewire::birewire.rewire.bipartite(rr1,verbose=FALSE)  %>% exchange(gg=rr)  %>% mean1
  newd<-BiRewire::birewire.rewire.bipartite(dd1,verbose=FALSE) %>% exchange(gg=dd) %>% mean1
  newgr<-BiRewire::birewire.rewire.bipartite(rg,verbose=FALSE) %>% t %>% mean1
  newrg<-BiRewire::birewire.rewire.bipartite(rg,verbose=FALSE) %>%  mean1
  newdg<-BiRewire::birewire.rewire.bipartite(gdd,verbose=FALSE)  %>% t%>% mean1
  newgd<-BiRewire::birewire.rewire.bipartite(gdd,verbose=FALSE) %>% mean1
  DR<-t(rd)
  RD<-rd
  #GG,GR,GD,DD,DG,DR,RR,RG,RD
  #
  m=combine(newg,newgr,newgd,newd,newdg,DR,newr,newrg,RD) %>% Normalize %>% rwr(newd=newd,newr=newr) %>% extract_max(n=10) %>% ml(y=ob) 
  assign(paste("r",i,sep=""),m)
}

R<-matrix(0,nrow(ob),ncol(ob))
for (i in 1:100){
  R<-R+get(paste("r",i,sep=""))
}
p<-R/100
colnames(p)<-colnames(ob)
write_rds(p,"less0.01_first10RWR_matirx.rds")
##################################plot
library(dendextend)
library(readr)
l<-read_rds("obseved_data.rds")
p<-read_rds("less0.01_first10RWR_matirx.rds")
ob<-read_rds("first10_RWR_matrix.rds")#########top10
New_ROC_PR<-function(cut_off,ll,top){ #top is value of extraction value
  l<-ll[names(cut_off)]
  match_list<-function(r1,l){
    v<-vector("numeric",length(r1))
    for (i in 1: length(r1)){
      v[i]<-length(na.omit(match(r1[[i]],unlist(l[[i]]))))
    }
    return(v)}
  nr<-length(unique(unlist(l))) ############after filterring 0 cloumn of diease for p-value
  g<-match_list(r1=cut_off,l=l) 
  P<-sum(g)/sum((sapply(cut_off,length)))
  f<-rep(10,length(l))
  #(sapply(l,function(x){length(unlist(x))}))
  R<-sum(g)/sum(f*top/nr) #########modifications
  tpr<-R
  fpr<-sum(abs((sapply(cut_off,length))-g))/sum(abs(f-nr))#########modifications
  roc_pr<-matrix(c(R,fpr,P,tpr),2,2)
  #matrix(c("R","FPR","P","TPR"),2,2)
  #[,1]  [,2] 
  #[1,] "R"   "P"  
  #[2,] "FPR" "TPR"
  return(roc_pr)}
cutoff_rwr<-function(result,ob,cut_off){
  rrr<-result<=cut_off #rrr is logcal matrix
  data.list<-list() 
  for (i in 1:ncol(rrr)){
    a<-ob[,i][rrr[,i]]
    assign(paste("d",i,sep=""),a)
    data.list[[i]]= get(paste("d",i,sep="")) 
  } 
  names(data.list)<-colnames(rrr)
  data.list[which(sapply(data.list,length)!=0)] ##remove 0 list
}
########use 1000 p-value
cut<-seq(from=0.001,to=1,by=0.001)
roc_rwr<-matrix(0,length(cut),2)
mm<-function(x){x[2,]}
m1<-function(x){x[1,]}
for(i in 1:length(cut)){
  roc_rwr[i,]<-cutoff_rwr(p,ob,cut[i]) %>% New_ROC_PR(ll=l,top=nrow(p)) %>% mm }

Roc_rwr<-roc_rwr
for(i in 1:length(cut)){
  Roc_rwr[i,]<-cutoff_rwr(p,ob,cut[i]) %>% New_ROC_PR(ll=l,top=nrow(p)) %>% m1 }
######plot the graph ROC
plot(Roc_rwr,type="l",col=4,xlab="Recall",ylab = "Precision",main=" Precision VS Recall of RWR with differernt p-value")
plot(roc_rwr,type="l",col=2,xlab="False positive rate",ylab="True positive rate",main="RWR with different p-values")  
points(c(min(roc_rwr[,1]),max(roc_rwr[,1])),c(min(roc_rwr[,2]),max(roc_rwr[,2])),type="l",lty=3)









pa<-function(l,top){ nr<-length(unique(unlist(l)))
re<-sapply(l,function(x){length(unlist(x))})
re[which(re>10)]<-re[which(re>10)]*top/nr
return(re)}


