rm(list = ls())
library(gplots)
library(RColorBrewer)
library(ape)
library(stats)
library(tidyverse)
library(ggpubr)

####################################################################################################
######### Pairwise SNP differences btw isolates per NG-MAST, NG-STAR and MLST ST           ######### 
####################################################################################################
## load typing data
typing_data <- read.csv('input/typing_data.csv')

### load DNA alignment
ali <- read.dna('input/gubbins_clean_variant_refremoved.fasta', format='fasta')
### Calculate pairwise distance btw SNPS
pw_dist <- dist.gene(ali)

#function to convert matrix to list
#from https://github.com/vmikk/metagMisc/tree/master/R
dist2list <- function (dist, tri=TRUE) {
  if (!class(dist) == "dist") { stop("Error: The input data must be a dist object.\n") }
  
  dat <- as.data.frame(as.matrix(dist))
  if (is.null(names(dat))) {
    rownames(dat) <- paste(1:nrow(dat))
  }
  value <- stack(dat)$values
  rnames <- rownames(dat)
  namecol <- expand.grid(rnames, rnames)
  colnames(namecol) <- c("col", "row")
  res <- data.frame(namecol, value)
  
  if(tri == TRUE){    # return only lower triangular part of dist
    res <- res[-which(upper.tri(as.matrix(dist), diag = T)), ]
  }
  
  return(res)
}

#create list
dist_matrix_longformat <- dist2list(pw_dist, tri=TRUE)

####### NGMAST ############3
NGMAST <-typing_data %>% select(Isolate, NGMAST) %>% group_by(NGMAST)

#rename columns for easier merging
colnames(dist_matrix_longformat)[1] <- "Isolate" 
colnames(dist_matrix_longformat)[2] <- "Isolate_row" 

#merge based on Isolate
ndf <- merge(dist_matrix_longformat, NGMAST, by="Isolate")
colnames(ndf)[1] <- "Isolate_col"
colnames(ndf)[2] <- "Isolate"
colnames(ndf)[3] <- "pairwiseSNP"
colnames(ndf)[4] <- "Isolate_col_NGMAST"
ndf2 <- merge(ndf, NGMAST, by="Isolate")
colnames(ndf2)[1] <- "Isolate_row"
colnames(ndf2)[5] <- "Isolate_row_NGMAST"

#if Isolate_row_NGMAST and isoalte_col_NGMAST match, then group and plot based on this
matches <- ndf2 %>% filter(ndf2$Isolate_col_NGMAST==ndf2$Isolate_row_NGMAST)
NGMASTp <-ggplot(matches, aes(x=as.character(Isolate_col_NGMAST), y=pairwiseSNP))+
  geom_jitter(width = 0.3, height = 0.1,alpha = .5,size=2)+
  theme_bw()+
  labs(x = 'NGMAST type', 
       y = 'pairwise SNPs')+
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle = 90,hjust=1))

####### NGSTAR ############3
NGSTAR <-typing_data %>% select(Isolate, ngSTar) %>% group_by(ngSTar)
#rename columns for easier merging
colnames(dist_matrix_longformat)[1] <- "Isolate" 
colnames(dist_matrix_longformat)[2] <- "Isolate_row" 

#merge based on Isolate
ndfSTAR <- merge(dist_matrix_longformat, NGSTAR, by="Isolate")
colnames(ndfSTAR)[1] <- "Isolate_col"
colnames(ndfSTAR)[2] <- "Isolate"
colnames(ndfSTAR)[3] <- "pairwiseSNP"
colnames(ndfSTAR)[4] <- "Isolate_col_NGSTAR"
ndf2STAR <- merge(ndfSTAR, NGSTAR, by="Isolate")
colnames(ndf2STAR)[1] <- "Isolate_row"
colnames(ndf2STAR)[5] <- "Isolate_row_NGSTAR"

#if Isolate_row_NGSTAR and isoalte_col_NGSTAR match, then group and plot based on this
matchesSTAR <- ndf2STAR %>% filter(ndf2STAR$Isolate_col_NGSTAR==ndf2STAR$Isolate_row_NGSTAR)
NGSTARp <- ggplot(matchesSTAR, aes(x=as.character(Isolate_col_NGSTAR), y=pairwiseSNP))+
  geom_jitter(width = 0.3, height = 0.1,alpha = .5,size=2,color="red")+
  theme_bw()+
  labs(x = 'NGSTAR type', 
       y = 'pairwise SNPs')+
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle = 90,hjust=1))

####### MLST ############3
MLST <-typing_data %>% select(Isolate, MLST) %>% group_by(MLST)
#rename columns for easier merging
colnames(dist_matrix_longformat)[1] <- "Isolate" 
colnames(dist_matrix_longformat)[2] <- "Isolate_row" 

#merge based on Isolate
ndfMLST <- merge(dist_matrix_longformat, MLST, by="Isolate")
colnames(ndfMLST)[1] <- "Isolate_col"
colnames(ndfMLST)[2] <- "Isolate"
colnames(ndfMLST)[3] <- "pairwiseSNP"
colnames(ndfMLST)[4] <- "Isolate_col_MLST"
ndfMLST2 <- merge(ndfMLST, MLST, by="Isolate")
colnames(ndfMLST2)[1] <- "Isolate_row"
colnames(ndfMLST2)[5] <- "Isolate_row_MLST"

#if Isolate_row_MLST and isoalte_col_MLST match, then group and plot based on this
matchesMLST <- ndfMLST2 %>% filter(ndfMLST2$Isolate_col_MLST==ndfMLST2$Isolate_row_MLST)

MLSTp <- ggplot(matchesMLST, aes(x=as.character(Isolate_col_MLST), y=pairwiseSNP))+
  geom_jitter(width = 0.3, height = 0.1,alpha = .5,size=2,color='blue')+
  theme_bw()+
  xlab("MLST type")+
  ylab("pairwise SNPs")  +
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle = 90, hjust=1))
  
############## Figure 4 ##############
ggarrange(MLSTp, NGSTARp, NGMASTp,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)
