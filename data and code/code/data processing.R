####Data processing
##Get gene-disease from OUGene databases of publicly available data
# To load the gene-disease association
UG<-read.csv("D:/RWR/UG.csv") #  UG is the reactionship between disease and gene underexperssion
head(UG)
#how many genes are involved in 
length(unique(UG$Ensemble.ID)) #[1] 577

# To simplify the "UG" data.frame
UGs<-UG[,c(2,5)]
ma<-UGs[UGs$Ensemble.ID %in% rna_gene$ENSEMBLPROT,] 
nrow(ma)# only 34 matched result
#[1] 34 
OG<-read.csv("D:/RWR/OG.csv") #OG is the reactionship between disease and gene overexperssion

# To simplify the "OG" data.frame
OGs<-OG[,c(2,5)]
mb<-OGs[OGs$Ensemble.ID %in% rna_gene$ENSEMBLPROT,] 
#at this step, remarkably x<-OGs[ OGs$Ensemble.ID %in% unique(rna_gene$ENSEMBLPROT),] is the same as mb
# since therer are many repeat gene ID in "rna_gene$ENSEMBLPROT", but OGs[ unique(rna_gene$ENSEMBLPROT) %in% OGs$Ensemble.ID,] 
# is total different, the later one is wronng and it include worng OGs element
nrow(mb)
#[1] 4038
#######Due to the less match data of underexpression data, I decide to combind the underexpression and
##########overexpression gene-disease association
gene_disease<-unique(rbind(mb,ma)) #####have to unique() to remobe the repeat elements
nrow(gene_disease)
#[1] 1698
#To calculate how many gene involved in the gene-disease association data
length(unique(gene_disease$Ensemble.ID))
#[1] 262
length(unique(gene_disease$Disease))
#[1] 450
# this is to achieve the text-transform
gene_disease$Disease<-tolower(gene_disease$Disease)
head(gene_disease)
write_rds(gene_disease,"C:/Users/Administrator/Desktop/extend data/262gene_450disease_frame.rds")

#######################################################################################################

##### To process the RNA m6A-sites, which is obtained form the Xingyu
site_condition<-read_rds("C:/Users/Administrator/Desktop/extend data/sdm_filtered_testm.rds")
dim(site_condition)  # there are toal 28278 m6A-sites under 38 condictions
#[1]    38 28278

####To filter the genes which do not involved at both gene-disease anf gene-m6A site network
gene_disease<-read_rds("C:/Users/Administrator/Desktop/extend data/262gene_450disease_frame.rds")
rna_gene<-read_rds("C:/Users/Administrator/Desktop/extend data/1403gene_6576rna.rds")
rg<-rna_gene[rna_gene$ENSEMBLPROT %in% gene_disease$Ensemble.ID,]
nrow(na.omit(rg))
#[1] 1109
length(unique(rg$ENSEMBLPROT)) # it consists with the number of gene in "gene_disease" data.frame
#[1] 262
length(unique(rg$modName))
#[1] 1109
rna_condition<-site_condition[,rg$modName]
dim(rna_condition)
#[1]   38 1109
write_rds(rna_condition,"D:/RWR/extend data/1109m6A-sites_38conditions.rds")
##########################
# TO obtain RNA m6A-sites ce-expression network and these are server codes below
library(readr)
rd<-read_rds("D:/RWR/extend data/1109m6A-sites_38conditions.rds") #
dim(rd) 
#40 140574
r<-abs(cor(rd,method = "spearman")) #abs!!!
dim(r) [1] 1109 1109
###p-value
r2<-WGCNA::corPvalueFisher(r,nSamples = 38,twoSided = TRUE)# the 38 is the 40 coniditons, which means
# the original 38 groups data
r3<-multtest::mt.rawp2adjp(r2,proc="Bonferroni")#adjust p-value
adj_r<-r3$adjp[,2]# Round numbers
adjp_r<-matrix(adj_r[order(r3$index)],ncol=ncol(r2))
colnames(adjp_r)<-colnames(r2)
rownames(adjp_r)<-rownames(r2)
# To generate the 0/1 matrix
f_r<-r2<0.01 & adjp_r>=1
f_r1<-matrix(as.numeric(f_r),nrow = nrow(r2),ncol = ncol(r2))
colnames(f_r1)<-colnames(r2)
rownames(f_r1)<-rownames(r2)

#To add weight to each node
exchange<-function(x,y01){x[which(y01 == 0)] <- y01[which(y01== 0)]
return(x)}
rr<-exchange(r,f_r1) #[1] 1109 1109
write_rds(rr,"D:/RWR/extend data/1109rr_0weight_network.rds")
#######################################################################################################
######To build the gene-gene co-expression network
gd<-read_rds("D:/RWR/extend data/test_exp_f.rds")
dim(gd)
#[1]    38 18635

######to convert the gene ID
i<-AnnotationDbi::select(org.Hs.eg.db, keys = colnames(gd), columns = "ENSEMBLPROT", keytype = "SYMBOL")
nrow(i) #[1] 30831

ri<-na.omit(i)
nrow(ri)#it means that there are 13250 missing

ii<-ri[ri$SYMBOL %in% colnames(gd),]#which is the same as ncol(ri)
dim(na.omit(ii)) #[1] 17581    2

# To convert the gene ID
filter_gd<-gd[,ii$SYMBOL]
dim(filter_gd) #[1]    38 17581

colnames(filter_gd)<-ii$ENSEMBLPROT
write_rds(filter_gd,"D:/RWR/extend_data/17581gene_38condition.rds")
#In this netwotk, thare are 18653, thus I filter genes from 5580 genes are not enough, I have to
#filter again
# this is to extract the different gene IDs
#########################################################################################################
########################################################################################################
inf<- select(org.Hs.eg.db, 
             keys = colnames(gd),
             columns = c("ENSEMBL","ENSEMBLPROT"),
             keytype = "SYMBOL")
inf1<-inf[inf$ENSEMBLPROT %in% ii$ENSEMBLPROT,] #[1] 51051     3
rna_dna<-read_rds("site2gene.rds")
######## this is only 5749 site can match
length(na.omit(match(rna_dna$ENSEMBL,inf1$ENSEMBL))) #[1] 5749

#######
f_gr<-rna_dna[rna_dna$ENSEMBL %in% inf1$ENSEMBL,]
dim(f_gr) #[1] 5749    2

inf2<-merge(inf1,f_gr,by.x="ENSEMBL", by.y="ENSEMBL",all=T)
rinf2<-na.omit(inf2)[,c(3,4)]
length(unique(rinf2$modName)) #there are 5749 m6A-sites

length(unique(rinf2$ENSEMBLPROT)) #there are 3875 genes


# To calculate the number of gene that is related with m6A-sites
length(unique(rna_gene$ENSEMBLPROT)) #[1] 1403

length(unique(rna_gene$modName)) #[1] 6576

write_rds(rna_gene[,c(1,3)],"1403gene_6576rna.rds")
#########################################################################################################
########################################################################################################
# To load the gene-disease association
library(readr)
g17d<-read_rds("D:/RWR/extend_data/17581gene_38condition.rds")
UG<-read.csv("D:/RWR/UG.csv") #  UG is the reactionship between disease and gene underexperssion

#how many genes are involved in 
length(unique(UG$Ensemble.ID)) #[1] 577

# To simplify the "UG" data.frame
UGs<-UG[,c(2,5)]
ma<-UGs[UGs$Ensemble.ID %in% colnames(g17d),] 
nrow(ma)# only 34 matched result
#[1] 171 
# to verify the number of match 
length(na.omit(match(UGs$Ensemble.ID,colnames(g17d)))) #[1] 171

OG<-read.csv("D:/RWR/OG.csv") #OG is the reactionship between disease and gene overexperssion

# To simplify the "OG" data.frame
OGs<-OG[,c(2,5)]
mb<-OGs[OGs$Ensemble.ID %in% colnames(g17d),] 
nrow(mb) #[1] 14257

# to verify the number of match 
length(na.omit(match(OGs$Ensemble.ID,colnames(g17d)))) #[1] 14257

##########Due to the less match data of underexpression data, I decide to combind the underexpression and
##########overexpression gene-disease association
nrow(rbind(mb,ma)) #[1] 14428

gene_disease<-unique(rbind(mb,ma)) #####have to unique() to remobe the repeat elements
nrow(gene_disease) #[1] 7916

#To calculate how many gene involved in the gene-disease association data
length(unique(gene_disease$Ensemble.ID)) #[1] 1316

length(unique(gene_disease$Disease)) #[1] 836

# this is to achieve the text-transform
gene_disease$Disease<-tolower(gene_disease$Disease)

write_rds(gene_disease,"D:/RWR/extend_data/1316gene_836disease_frame.rds")

#######################################################################################################
##############due to the "desc2017.xml" is too big, the code below will run on the serve
#As we already get the related disease, next we can directly disease-disease association network
library(methods)
library(readr)
library(dplyr)
library(tidyr)
library(dendextend)
gd<-read_rds("1316gene_836disease_frame.rds")

#To load the disease ID database,remarkablely,there are four types mapping,
#mapping from cross reference (MFR); 2) mapping from synonyms (MFS); 3) mapping from inferring (MFI),Thus, in
#this sessction, we only use the mapping from synonyms (MFS) .
################################ example, why we dicede to use MFR & MFS
filter(mesh,mesh_id=="D011125")

filter(mesh,mesh_id=="D005736")

################################
mesh<-read.csv("mesh2do.csv")
mesh_RS<-filter(mesh,mapping !="MFI")
length(unique(mesh_RS$mesh_id)) #[1] 2936 there are total 2936 disease
mesh_f<-mesh_RS[,c(2,3)] #To samplify the data.frame
d_id<-merge(gd,mesh_f,by.x="Disease",by.y="do_term",all=T) %>% na.omit()

length(unique(d_id$mesh_id)) #after convering disease-id, there are 706 disease
#[1] 706
length(unique(d_id$Ensemble.ID))#[1] 1203

write_rds(d_id,"1203gene_706disease_aftermarchdiseaseID_dataframe.rds")
#######################################################################################################
#After discussing, our lab group have improve the way that construct the co-expression network,with that 
#in mind, we will use the way to build the gene co-expression network
#This is the fucntion that build the co-expression matrix
adjmake <- function(x, y, quant, pcut){
  cor <- cor(x, y, method = 'pearson')
  rawp <-WGCNA::corPvalueFisher(cor, 38) 
  #mt <- multtest::mt.rawp2adjp(rawp, proc = 'BH')
  #adj <- mt$adjp[,2]
  #adjp <- matrix(adj[order(mt$index)], ncol = ncol(cor))
  cor <- abs(cor)
  scc_cut <- IRanges::quantile(cor, quant)
  #cor[adjp > pcut | cor < scc_cut] <- 0
  cor[rawp > pcut | cor < scc_cut] <- 0
  cor[cor > 0] <- 1
  diag(cor) <- 0
  cor
}
adj_meth <- adjmake(x = new_testm, y = new_testm, quant = 0.80, pcut = 0.05)
#To build the 1316 genes's co-expression network 
library(readr)
library(dendextend)
gd<-read_rds("17581gene_38condition.rds") #there are 17581 gene expression level under 38 conditions
gene<-read_rds("1203gene_706disease_aftermarchdiseaseID_dataframe.rds")
un<-unique(gene$Ensemble.ID) %>% as.vector()
gdf<-gd[,un]
dim(gdf) #[1]   38 1203

#for self-network, we do not need 0/1, we need weigth, so
adjmake_w <- function(x, y, quant, pcut){
  cor <- cor(x, y, method = 'pearson')
  rawp <-WGCNA::corPvalueFisher(cor, 38) 
  mt <- multtest::mt.rawp2adjp(rawp, proc = 'BH') 
  p.adjust(rawp, method = "fdr", n = length(rawp)) 
  adj <- mt$adjp[,2]
  adjp <- matrix(adj[order(mt$index)], ncol = ncol(cor))
  cor <- abs(cor)
  scc_cut <- IRanges::quantile(cor, quant)
  cor[adjp > pcut | cor < scc_cut] <- 0
  diag(cor) <- 0
  cor
}
adj_meth <- adjmake_w(x = gdf, y = gdf, quant = 0.80, pcut = 0.05)

which(apply(adj_meth,2,sd)==0) #there are no 0 cloumn

write_rds(adj_meth,"1203gene_coexression_0weigth_network.rds")
#####################################################################################################################################################
#this is runniing on serve
xmldataframe <- xmlToDataFrame("desc2017.xml")
d_id<-read_rds("1203gene_706disease_aftermarchdiseaseID_dataframe.rds")
x<-select(xmldataframe,DescriptorUI,TreeNumberList,DescriptorName)
xm<-x[x$DescriptorUI %in% unique(d_id$mesh_id),]
nrow(xm)
#[1] 705
length(unique(xm$DescriptorUI))
write_rds(xm,"d_id_xm.rds")
xm<-read_rds("d_id_xm.rds")# it means that there are 705 diseaase

tree<-select(xm,TreeNumberList)
tr<-tree$TreeNumberList
n<-nrow(xm)
t<-vector("numeric",n)
for (i in 1:n ) {
  t[i]<- strsplit(tr[i],".",fixed = TRUE) 
  T<-as.list(t)
}
names(T)<-xm$DescriptorUI

# the fuvntion that can calculate the association between differernt diseases
dis<-function (x) { 
  N<-vector("numeric",n)
  U<-vector("numeric",n)
  for (i in 1:n) {
    N[i]<-length(na.omit(match(unlist(x),unlist(T[i]))))
    U[i]<-sum(length(unlist(x)),length(unlist(T[i])))
    D<-N/U
  } 
  return(D)
}
n<-length(T)
t<-matrix(0,n,n)
for (i in 1:n){
  t[,i]<-dis(T[i])
  colnames(t)<-xm$DescriptorUI
  rownames(t)<-xm$DescriptorUI
}

which(apply(t,2,sd)==0)
#named integer(0), there is no 0 cloumn
dim(t) #[1] 705 705
write_rds(t,"705disease-0weight-network.rds")
##########################################################################################################################
# to build the gene and all m6A-sites association network
library(readr)
if (!require("multtest")) {
  source("https://bioconductor.org/biocLite.R")
  install.packages("multtest",repos = "http://cran.us.r-project.org")
  
}
library(dendextend)
rna<-read_rds("sdm_filtered_testm.rds")
dim(rna)
gd<-read_rds("17581gene_38condition.rds")
gene<-read_rds("1203gene_706disease_aftermarchdiseaseID_dataframe.rds")
un<-unique(gene$Ensemble.ID) %>% as.vector()
gdf<-gd[,un]
adjmake <- function(x, y, quant, pcut){
  cor <- cor(x, y, method = 'pearson')
  rawp <-WGCNA::corPvalueFisher(cor, 38) 
  mt <- multtest::mt.rawp2adjp(rawp, proc = 'BH') 
  p.adjust(rawp, method = "fdr", n = length(rawp)) 
  adj <- mt$adjp[,2]
  adjp <- matrix(adj[order(mt$index)], ncol = ncol(cor))
  cor <- abs(cor)
  scc_cut <- IRanges::quantile(cor, quant)
  cor[adjp > pcut | cor < scc_cut] <- 0
  cor[cor > 0] <- 1
  diag(cor) <- 0
  cor
}
#there are problem, if we still keep the same parameter as the gene co-expression network, we will filter
# amounts of data. Hence, in order to preserve the data, we use the parameter below, after testing many
#different parameter
adj_meth <- adjmake(x = rna, y = gdf, quant = 0.80, pcut = 0.3)
write_rds(adj_meth,"grweigth.rds")
gr<-read_rds("0.8q_0.3p_gr.rds")
dim(gr) #[1] 28278  1203

#To check the 0 column
rzero<-which(apply(gr,1,sd)==0)
length(rzero) #[1] 25451

czero<-which(apply(gr,2,sd)==0)
length(czero) #[1]  123

#To remove the 0 column
gr_f<-gr[-rzero,-czero]
dim(gr_f)
#[1] 2827 1080
write_rds(gr_f,"2827site_1080gene_01matrix.rds")
############################################################################################################################################
#to build RNA m6A-sites co-expression
library(readr)
RNA<-read_rds("sdm_filtered_testm.rds")
dim(RNA)
#[1]    38 28278
gr<-read_rds("2827site_1080gene_01matrix.rds")
dim(gr) #[1] 2827 1080
site<-RNA[,rownames(gr)]
dim(site) #1]   38 2827
#for self-network, we do not need 0/1, we need weigth, sunique(gene$Ensemble.ID)o
adjmake_w <- function(x, y, quant, pcut){
  cor <- cor(x, y, method = 'pearson')
  rawp <-WGCNA::corPvalueFisher(cor, 38) 
  mt <- multtest::mt.rawp2adjp(rawp, proc = 'BH') 
  p.adjust(rawp, method = "fdr", n = length(rawp)) 
  adj <- mt$adjp[,2]
  adjp <- matrix(adj[order(mt$index)], ncol = ncol(cor))
  cor <- abs(cor)
  scc_cut <- IRanges::quantile(cor, quant)
  cor[adjp > pcut | cor < scc_cut] <- 0
  diag(cor) <- 0
  cor
}
adj_meth <- adjmake_w(x = site, y = site, quant = 0.80, pcut = 0.05)
dim(adj_meth)
#[1] 2827 2827
zero<-which(apply(adj_meth,1,sd)==0)
length(zero)
#[1]  62
rna<-adj_meth[-zero,-zero]
dim(rna) #[1]  2765 2765
write_rds(rna,"2765rna_coexpression_0weight.rds")
###############################################################################################################
#futher modifiy the data
library(readr)
library(dplyr)
rna<-read_rds("2765rna_coexpression_0weight.rds")
o_gr<-read_rds("2827site_1080gene_01matrix.rds")
rg<-o_gr[colnames(rna),] #[1] 2765 1080
write_rds(rg,"filter_2765site_1080ene_01weight.rds")
o_gd<-read_rds("1203gene_706disease_aftermarchdiseaseID_dataframe.rds")
head(o_gd)
rownames(o_gd)<-NULL #reorider
##############################################################
dim(unique(o_gd)) #[1] 6790    3
dim(o_gd) #[1] 6790    3
#some differernt disease names are combine into same mesh id, so we have to unique(tog)
tog<-select(o_gd,Ensemble.ID,mesh_id)
dim(tog) #[1] 6790    2
dim(unique(tog)) #[1] 6776    2
t<-table(unique(tog))
t1<-as.matrix.data.frame(t)
colnames(t1)<-colnames(t)
rownames(t1)<-rownames(t)
dim(t1)
#[1] 7238 4284
e<-unique(o_gd$Ensemble.ID)#[1] 1203
m<-unique(o_gd$mesh_id) #[1] 706
gd<-t1[e,m]
dim(gd) # [1] 1203  706
write_rds(gd,"706disease_1203gene_01mariex.rds")
gd<-read_rds("706disease_1203gene_01mariex.rds")
gr<-read_rds("2827site_1080gene_01matrix.rds") #[1] 2827 1080
dd<-read_rds("705disease-0weight-network.rds")
gdf<-gd[colnames(gr),colnames(dd)] 
dim(gdf)#[1] 1080  705
write_rds(gdf,"1080gne_705disease_01matrix.rds")
o_gg<-read_rds("1203gene_coexression_0weigth_network.rds")
gg<-o_gg[colnames(rg),colnames(rg)]
dim(gg)
write_rds(gg,"1080gene-coexpression_0weight.rds")
rr<-read_rds("2765rna_coexpression_0weight.rds")
########################################################
#Attention !!!! if remove the rna site, there are some gene is 0 clumn
f_gr<-gr[colnames(rr),]
dim(f_gr) #[1] 2765 1080
write_rds(f_gr,"filter2765site_1080gene_01matrix.rds")
