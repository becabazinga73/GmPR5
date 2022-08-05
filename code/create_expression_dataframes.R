# Create a data frame of GmPR-5 gene expression on atlas and recently downloaded stress-related samples
# Jan 2021

library(tidyverse)

#----Create character vector of GmPR-5 genes----
pr5 <- read.csv(file = '/home/rebeca/Desktop/Projeto_PR-5/products/tables/table1.txt')[,1]

#----Create character vector of abiotic and biotic stress-related samples----
load(file = '/home/rebeca/Desktop/Projeto_PR-5/data/stress_samples_for_DGE_analysis.RData')
abiotic_samples <- Reduce(rbind, abiotic_samplelist)$BioSample
biotic_samples <- Reduce(rbind, biotic_samplelist)$BioSample
stress_samples <- c(abiotic_samples, biotic_samples)

#----Create data frame of raw counts for atlas samples----
load(file = '/home/rebeca/Desktop/Projeto_PR-5/data/atlas.rda')
atlas <- atlas[pr5, ]

 (genes Glyma.14G077267, Glyma.14G077333 não possuem amostras)

#----Confirm samples
print <- print(select(stressexp, everything(), matches("stress_samples")))

#----Load count matrix with new samples (on server)----
counts <- load(here("data", "new_samples.rda"))
counts <- as.data.frame(new_samples)
counts <- counts[counts$gene_id %in% pr5, ]
counts <- counts[, colnames(counts) %in% c("gene_id", stress_samples)]
counts <- counts[order(counts$gene_id), ]

(não tem o gene Glyma.12G031000, Glyma.14G077267 e Glyma.14G077333)

#----Create data frame of global expression: atlas + new samples----
globalexp <- merge(atlas, counts, by.x="row.names", by.y="gene_id", all=TRUE)
globalexp <- globalexp[-c(54,55), ]
rownames(globalexp) <- globalexp$Row.names
globalexp$Row.names <- NULL

#----Create data frame of stress-related expression----
stressexp <- globalexp[, colnames(globalexp) %in% stress_samples]

#----Save objects----
save(globalexp, stressexp, abiotic_samples, biotic_samples, stress_samples,
     file = here("data", "expression.RData"),
     compress = TRUE)