rm(list = ls())
library(gplots)
library(RColorBrewer)
library(ape)
library(stats)
library(tidyverse)
library(reshape)
library(ggplot2)
library(ggpubr)
#####################################################
######## Pairwise SNP clustering - Figure 3 #########
#####################################################

### load DNA alignment
#Gubbins output 13,003 nucleotides (contains Ns) - tidied up with snp-sites to contain only ACGT sites (5,059 ACGT sites, 3,871 sites informative)
ali <- read.dna('input/gubbins_clean_variant_refremoved.fasta', format='fasta')
### Calculate pairwise distance btw SNPS
pw_dist <- dist.gene(ali)
#convert to dataframe
m <- as.matrix(pw_dist)
m2 <- melt(m)[melt(upper.tri(m))$value,]
names(m2) <- c("c1", "c2", "distance")
m2

#cluster analysis
cluster <- hclust(pw_dist) #perform hierarchical clustering
clus10 <- as.data.frame(cutree(cluster, h=10))  %>% tibble::rownames_to_column()
clus_results <- list(clus10) %>% reduce(left_join, by = "rowname")
colnames(clus_results)<- c("Isolate","SNP10")

#####################################################
######## Merge with metadata and analyse ############
#####################################################
metadata <- read.csv('input/metadata_TableS5.csv')
md <-  metadata %>% 
  select(Isolate, Sex, Age, DHB...Region, BodySite, MSM_final)
md %>% count(MSM_final)
#replace all MSM = 4 with NA, because 4 = unknown or F who do not have sex with males, but there are no F isolates in this category anyway
md2 <- md %>% 
  mutate(MSM_final=replace(MSM_final, MSM_final==4, NA)) %>% 
  as.data.frame()
df2 <- merge(clus_results, md2, by="Isolate")

#Summary counts for all clusters
SNP10_counts <- df2 %>% count(SNP10, sort=TRUE) 

#Look at percentage of MSM in individual clusters -
MSMperc <- df2 %>%
  select(MSM_final, SNP10 ) %>% 
  group_by(SNP10, MSM_final) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

#############################################################
######## Plotting results with MSM = NA included ############
############################################################
#### Horizontal bar plot cluster size
#filter all clusters with >= 9 isolates in cluster
SNP10cts <- SNP10_counts %>% 
  arrange(desc(n)) %>% 
  filter(n >= 9)
SNP10cts[, 'SNP10'] <- as.factor(SNP10cts[, 'SNP10'])

size <- ggplot(SNP10cts, aes(x=as.character(SNP10), y=n))+
  geom_bar(stat="identity") +
  xlab("Clusters")+
  ylab("Number of isolates")+
  theme_bw()+
  coord_flip()+
  theme(text = element_text(size=20),
        axis.text.x = element_text(hjust=1)) 

#### plot MSM percentage ####
#choose all clusters with 9 or more isolates for figure: cluster 6,4,43,1,30, 23,13,52
clstr8<- MSMperc %>% filter(SNP10 %in% c("6","4","43","1","30","23","13","52"))

#plotting order should be 6, 4, 43, 1, 30, 23, 13 52
MSM <- ggplot(clstr8, aes(x = as.character(SNP10), y = freq, fill = factor(MSM_final))) +
  geom_bar(stat = 'identity') +
  theme_bw()+
  #scale_y_continuous(labels = percent) +
  labs(x = '', 
       y = 'Relative percentage',
       fill = 'MSM_final')+
  coord_flip()+
  theme(legend.position = "none")+
  theme(text = element_text(size=20),
        axis.text.x = element_text(hjust=1))
  
#### plot age range in cluster ####
age<- df2 %>% filter(SNP10 %in% c("6","4","43","1","30","23","13","52"))
age_p <- ggplot(age, aes(x=as.character(SNP10), y=Age))+
  geom_jitter(stat="identity",width = 0.15, aes(color=factor(MSM_final))) +
  #labs(title="TITLE",x="Clusters", y="Age range", colour="Sexual orientation")+
  xlab("")+
  ylab("Age range")+
  theme_bw()+
  coord_flip()+
  scale_colour_discrete(name="Sexual orientation",
                      breaks=c("1", "2", "3", "NA"),
                      labels=c("MSM", "Heterosexual", "Female", "NA"))+
  theme(text = element_text(size=20),
        axis.text.x = element_text(hjust=1))

#Final graph - Figure 3
ggarrange(size, MSM, age_p,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)

#look at age range
ageMSM<- df2 %>% filter(SNP10 %in% c("1"))
median(ageMSM$Age, na.rm=TRUE)
#median 26, range 18 - 48

agehet<- df2 %>% filter(SNP10 %in% c("13","23","52"))
median(agehet$Age, na.rm=TRUE)
#median 26, range 17 - 38

ageMix <- df2 %>% filter(SNP10 %in% c("6","43","4","30"))
median(ageMix$Age, na.rm=TRUE)
#median 29, range 0 - 63

