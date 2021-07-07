rm(list = ls())
library(tidyverse)
library(caper)
#library(dplyr)
library(phytools)

#testing binary traits (Sex, MSM) for association with phylogenetic tree
#read tree in newick format
gcTree <- read.newick("input/gubbins_customSNP_iqtree_contree.nwk")
gcTree <- multi2di(gcTree)
gcTree$edge.length[gcTree$edge.length==0]<-max(nodeHeights(gcTree))*1e-6
samples <- gcTree$tip.label

#remove internal node labels (bootstrap values) to avoid error for comparative.data function
#http://blog.phytools.org/2015/08/removing-node-labels-from-newick-string.html
text<-write.tree(gcTree)
text
strip.nodelabels<-function(text){
  obj<-strsplit(text,"")[[1]]
  cp<-grep(")",obj)
  csc<-c(grep(":",obj),length(obj))
  exc<-cbind(cp,sapply(cp,function(x,y) y[which(y>x)[1]],y=csc))
  exc<-exc[(exc[,2]-exc[,1])>1,]
  inc<-rep(TRUE,length(obj))
  if(nrow(exc)>0) for(i in 1:nrow(exc)) 
    inc[(exc[i,1]+1):(exc[i,2]-1)]<-FALSE
  paste(obj[inc],collapse="")
}
text2 <- strip.nodelabels(text)
gcTree2 <- read.tree(text=text2)

#metadata - traits
metadata <- read.table("input/metadata_TableS5.csv", header = T, sep = ",", comment.char = "")
metadata.sex <- metadata %>% 
  dplyr::select(Sample.ID, Sex) %>% 
  rename(Isolate = Sample.ID) %>% 
  mutate(Sex=ifelse(Sex== "U", NA, Sex))
metadata.MSM <- metadata %>% 
  dplyr::select(Sample.ID, MSM_final) %>% 
  rename(Isolate = Sample.ID) %>% 
  mutate(MSM = ifelse(MSM_final == "1", "yes", "no"))
metadata.F <- metadata %>% 
  dplyr::select(Sample.ID, MSM_final) %>% 
  rename(Isolate = Sample.ID) %>% 
  mutate(Fem = ifelse(MSM_final == "3", "yes", "no"))

gcTree.sex <- comparative.data(phy = gcTree2, data = metadata.sex, names.col = Isolate, na.omit = TRUE)
gcTree.MSM <- comparative.data(phy = gcTree2, data = metadata.MSM, names.col = Isolate, na.omit = TRUE)
gcTree.F <- comparative.data(phy = gcTree2, data = metadata.F, names.col = Isolate, na.omit = TRUE)

sex.d <- phylo.d(data = gcTree.sex, binvar = Sex, permut = 1000)
print(sex.d)
#p<0.05 (p=0.042)
#significant
#Estimated D :  0.8908957
##Probability of E(D) resulting from no (random) phylogenetic structure :  0.042
#robability of E(D) resulting from Brownian phylogenetic structure    :  0

MSM.d <- phylo.d(data = gcTree.MSM, binvar = MSM, permut = 1000)
print(MSM.d)
#Estimated D :  0.3523878
#Probability of E(D) resulting from no (random) phylogenetic structure :  0
#Probability of E(D) resulting from Brownian phylogenetic structure    :  0.001

F.d <- phylo.d(data = gcTree.F, binvar = Fem, permut = 1000)
print(MSM.d)
#Estimated D :  0.3523878
#Probability of E(D) resulting from no (random) phylogenetic structure :  0
#Probability of E(D) resulting from Brownian phylogenetic structure    :  0.001