#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library("DESeq2")
library("biomaRt")
#library("AnnotationDbi")
#library("org.Hs.eg.db")

source("Directories.txt")
##Directories.txt specifies only one directory used in this file, DEMasterDir
##DeMasterDir contains:
  ## DESeq2.r -- This file
  ## fcts.RData -- A file containing the .RData output of featureCounts
  
setwd(DEMasterDir)
load("fcts.RData")

###########
##Functions
###########
##The mapIds()function from AnnotationDbi using org.Hs.eg.db contained many missing symbols. getGeneNames() is used instead
##getGeneNames() takes a character vector of ensembl Gene IDs and attaches the corresponding gene name
##Any IDs which are not found in the bioMart database are attributed with "Not In Database" under gene name
##depends on library("bioMart)
getGeneNames <- function(geneIds){
  message("Running getGeneNames()")
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  message("Aquiring Gene Names from ensembl bioMart database")
  geneNames <- getBM(attributes = c("ensembl_gene_id","external_gene_name"), 
                     mart = ensembl, 
                     filters = "ensembl_gene_id", 
                     values = geneIds,
                     uniqueRows = FALSE,
                     verbose = FALSE)
  
  geneNamesNotFound <- subset(geneIds, !(geneIds %in% geneNames[,1]))
  
    if(length(geneNamesNotFound) > 0){
    message(paste("Warning:", length(geneNamesNotFound),"gene IDs not found in database"))
    geneNamesNotFound = data.frame(ensembl_gene_id = GeneNamesNotFound, 
                                   external_gene_name = rep("Not In Database", times = length(GeneNamesNotFound)) )
    geneNames <- rbind(geneNames, geneNamesNotFound)
  } 
  #check that geneNames matches geneIds
  inGeneIds = all(geneNames[,1] %in% geneIds)
  geneNames <- geneNames[match(geneIds, geneNames[,1]),]
  valuesMatched = all(geneIds == geneNames[,1], na.rm = TRUE)

    if(valuesMatched & inGeneIds){
    message("Order of input vector matches output \n done.")
    return(geneNames)
    } else{
    message(paste("All Values in geneNames are in GeneIds:", inGeneIds))
    message(paste("GeneIds is equal to gene names:", valuesMatched))
    return(NULL)
  }
}

############################################################
##Generate a condition vector from alphabetical designations
############################################################
##Treatment group is recorded as the last character in each sample name
des <- substring(fc1$targets, nchar(fc1$targets))
des <- replace(des, which(des =='A' | des == 'C'), "Control")
des <- replace(des, which(des =='B' | des == 'D'), "Treatment")

#####################################
##Create coldata perameter for DESeq2
#####################################

coldata <- data.frame(row.names= fc1$targets, Condition = des )
##check if caldata matches columns in counts table
all(rownames(coldata) %in% colnames(fc1$counts))
all(rownames(coldata) == colnames(fc1$counts))


#######################
##Generate DESeqDataSet
#######################
dds <- DESeqDataSetFromMatrix(countData = fc1$counts,
                              colData = coldata,
                              design = ~ Condition )

################################################
##Prefilter for low counts and set factor levels
################################################
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds$Condition <- factor(dds$Condition, levels=c("Control","Treatment"))


#########################
##Differential Expression
#########################
dds <- DESeq(dds)
res <- results(dds, contrast = c("Condition", "Treatment", "Control"))

####################################
##Filtering Results by BH Correction
####################################
resOrdered <- res[order(res$padj),]
resSig <- subset(resOrdered, padj < 0.1)

####################
##Annotating Results
####################

resSymbols<-getGeneNames(rownames(resSig))
resSig <- cbind(resSymbols, resSig)
all(rownames(resSig) == resSig$ensembl_gene_id)
resSig$ensembl_gene_id <- NULL


write.csv(resSig, file = 'DifferentialExpressionAnalysis.csv')
write.table(counts(dds)[rownames(resSig),], file = "deCounts.txt", sep = '\t', quote = FALSE, row.names = TRUE)

plotMA(res, ylim = c(-5,5))