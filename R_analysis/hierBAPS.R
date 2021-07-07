#Hierarchical clustering analysis
library(rhierbaps)
set.seed(1234)
#load data
snp.matrix <- load_fasta("input/gubbins_clean_variant.fasta")
snp.matrix2 <- load_fasta("input/gubbins_1815_40it_clean_variant.fasta")

hb.results <- hierBAPS(snp.matrix, max.depth = 2, n.pops = 65, n.cores=24, n.extra.rounds = Inf)
hb.results2 <- hierBAPS(snp.matrix2, max.depth = 2, n.pops = 100, n.cores=48, n.extra.rounds = Inf)
