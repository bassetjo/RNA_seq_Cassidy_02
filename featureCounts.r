#installation of Rsubread from biocLite
source("http://bioconductor.org/biocLite.R")
biocLite("Rsubread")
library(Rsubread)

setwd("/home/user1/Documents/R_projects/RNA_seq_01")

#####################
##Functions
#####################

##getfileNames returns a single comma separated string of paths to bam files in the given directory
getfileNames <- function(DIR){
	bamfiles <- system(paste('ls ', DIR, '/*bam', sep = ""), intern = TRUE)
	return(bamfiles)
	}

strFileNames <- getfileNames("/home/user1/Documents/R_projects/RNA_seq_01/Alignments")

######################
##RUN FEATURE COUNTS
######################
fc1 <- featureCounts(files = c(strFileNames),
	annot.ext = "/home/user1/Documents/R_projects/RNA_seq_01/Annotation_files/Homo_sapiens.GRCh38.89.chr.gtf",
	strandSpecific = 2, largestOverlap =TRUE, isGTFAnnotationFile = TRUE)

##remove all the irrelevent directory information from the collumn names in fc1 tables
fc1$targets <- regmatches(fc1$targets, regexpr("13708.*[^.bam]", fc1$targets))
colnames(fc1$counts) = fc1$targets
colnames(fc1$stat) = c("Status", fc1$targets)

##save the featurecounts output as an r object
save(fc1, file = "fcts.RData")

####################
##Supplemental Code
####################

##convert fc1 to object class data frame
##fcdf <- data.frame("Name" = rownames(fc1$counts), fc1$counts, check.rows = TRUE, check.names = TRUE)

##write individual text files for all the elements of the output of feature counts
##sapply(names(fc1), function (x) write.table(fc1[[x]], file=paste(x, "txt", sep =".")))

##write a text file specifically formatted for use in GSEA
##write.table(fcdf, file = "GSEA_COUNTS.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE, fileEncoding = "UTF-8")