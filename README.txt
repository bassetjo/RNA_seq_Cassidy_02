RNA_seq_02
Created 10/6/2017

RNA_seq_02 is a collaborative RNA-seq analysis project performed by John Bassett, Reid Thompson, Abhinav Nellore, and Pamela Cassidy at Oregon Health & Science University.
The experiment and sequencing was performed by the Grossman laboratory at the University of Utah in collaboration with Pamela Cassidy.
The main difference between RNA_seq_01 and RNA_seq_02 is to change from using edgeR to using DESeq2. This choice was made largely as a way to learn about batch effects and controlling for multiple comparisons as per advice from Dr. Abhinav Nellore.

Methods:
	featureCounts
	DESeq2
	topGO

Files:

featureCounts.r
	*issue* Annotation file (AnnotationDir) is GRCH38.89 Compare with UTAH annotation file GRCh38.85. IS this acceptable?

	R script for generating read counts from .bam files using featureCounts

	Reference:
			Liao Y, Smyth GK and Shi W (2014). featureCounts: an efficient general purposeprogram for assigning sequence reads to genomic features. Bioinformatics,30(7):923-30.
			http://www.ncbi.nlm.nih.gov/pubmed/24227677


	Arguments used

	Input_files (files) : (.bam files) a list of bam files

	-a (annot.ext) : annotation file <string> directory of the annotation file. File which contains gene annotation information about the raw sequence.

	-s (isStrandSpecific) : Indicates if strand specific read counting should be performed. 0 (unstranded), 1 (stranded), 2 (reversely stranded)

	-T (nthreads) : <int> number of threads to using in computer processing. Between 1 and 32. 1 by default (only one thread used in JB analysis)

	--largestOverlap (largestOverlap) : if specified reads will be assigned to the target that has the largest number of overlapping bases.

	-F (isGTFAnnotationFile) : specify the format of the annotation file. by default C version is SAF; R version is GTF.

	Analysis performed by UTAH

		featureCounts -a /tomato/dev/data/Human/GRCh38.p7/Homo_sapiens.GRCh38.85.chr.gtf \
			-o 13708R_unique_counts.txt \
			-s 2 \
			-T 20 \
			--largestOverlap *bam
		*note* that in this analysis 20 threads are used for parallel processing this computation

	Analysis performed by John Bassett

		c1 <- featureCounts(files = c(strFileNames),
		annot.ext = AnnotationDir,
		strandSpecific = 2, 
		largestOverlap = TRUE, 
		isGTFAnnotationFile = TRUE)

	
DESeq2.r
	R script for performing differential expression analysis using DESeq2
	
	Reference:
			Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)

	running DESeq2.r will write 2 files to the working directory
		1. 'DifferentialExpressionAnalysis.csv'
			--a file containing the results of DESeq(dds) and annotated with HUGO gene symbols
		2. 'deCounts.txt'
			--a file containing the output of counts(dds) reduced to only differentially expressed genes with FDR < 0.1. Intended for GSEA. 
	

	Analysis Summary

		dds <- DESeqDataSetFromMatrix(countData = fc1$counts,
                              colData = coldata,
                              design = ~ Condition )
		dds <- DESeq(dds)
		res <- results(dds, contrast = c("Condition", "Treatment", "Control"))

			
topGO.r
	R script for performing enrichment analysis using topGO
	
	Reference:
			Adrian Alexa and Jorg Rahnenfuhrer (2016). topGO: Enrichment Analysis for Gene Ontology. R package version 2.30.0.
	
	This file is currently incomplete. Requires fixes and improvements that are documented within the file. Analysis and Statistical tests to be redone. 


fcts.RData
	R data object containing the output generated from featureCounts.r (not in repo)


DEAnalysis.RData
	R data object containing the differential expression results from DESeq2 (not in repo)
