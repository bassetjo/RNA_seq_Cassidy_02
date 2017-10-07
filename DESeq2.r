source("http://bioconductor.org/biocLite.R")
biocLite("GenomeInfoDb")
library("DESeq2")

setwd("X:/SOM/Dermatology/Research/Research Division/Leachman_Cassidy_lab/John/RNASeq/RNA_seq_02")
load("fcts.RData")

coldata <- data.frame(Samples = fc1$targets, Condition =  c(2,1,2,1,2,1,2,1,2,1,1,2,1,2,2,1,2,1,2,1,2,1))

