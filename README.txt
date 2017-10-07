RNA_seq_02
Created 10/6/2017
The main difference between RNA_seq_01 and RNA_seq_02 is to change from using edgeR to using DESeq2. This choice was made largely as a way to learn about batch effects and controlling for multiple comparisons as per advice from Dr. Abhinav Nellore.

files:
DESeq2.r
	R script for performing analysis using DESeq2
featureCounts.r
	R script for generating read counts from .bam files using featureCounts
fcts.RData
	R data object containing the count data generated from featureCounts.r
