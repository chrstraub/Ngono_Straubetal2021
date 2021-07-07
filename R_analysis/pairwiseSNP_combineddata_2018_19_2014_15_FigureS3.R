rm(list = ls())
library(gplots)
library(RColorBrewer)
library(ape)
library(stats)
library(tidyverse)

###############################################################################################
######## Pairwise SNP clustering for 2014/15 and 2018/19 data combined             ############
########                     counts + Figure S3                                    ############
###############################################################################################

### load DNA alignment
#Gubbins output 13,003 nucleotides (contains Ns) - tidied up with snp-sites to contain only ACGT sites (5,059 ACGT sites, 3,871 sites informative)
ali <- read.dna('input/gubbins_1815_40it_clean_variant_refremoved.fasta', format='fasta')

### Calculate pairwise distance btw SNPS
pw_dist <- dist.gene(ali)
head(pw_dist)

#cluster analysis
cluster <- hclust(pw_dist) #perform hierarchical clustering
clus10 <- as.data.frame(cutree(cluster, h=10))  %>% tibble::rownames_to_column()
clus_results <- list(clus10) %>% reduce(left_join, by = "rowname")
colnames(clus_results)<- c("Isolate", "SNP10")

metadata <- read.csv('input/metadata_1415_1819.csv')
df2 <- merge(clus_results, metadata, by="Isolate")

#Summary counts for all clusters
SNP10_counts <- clus_results%>% count(SNP10, sort=TRUE) 
#190 clusters - 94 singletons
#7 clusters that share isolates from 2014/15 and 2018/19
#2, 7, 15, 24, 92, 94, 107

#export cluster ids + associated isolates for highlighting in phylogenetic tree
isolates <- df2 %>% select(Isolate, SNP10, Year) %>%
filter(SNP10 %in% c("2","7","15","24","92","94","107"))

#relative percentage of year
year <- df2 %>%
  select(SNP10, Year ) %>% 
  group_by(SNP10, Year) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n)) %>% 
  filter(freq <1)

##### FIGURE S3 #####
ggplot(year, aes(x = as.character(SNP10), y = freq, fill = factor(Year))) +
  geom_bar(stat = 'identity') +
  theme_bw()+
  #scale_y_continuous(labels = percent) +
  labs(x = 'Cluster', 
       y = 'Relative percentage',
       fill = 'Year')+
  coord_flip()+
  theme(text = element_text(size=20),
        axis.text.x = element_text(hjust=1))
