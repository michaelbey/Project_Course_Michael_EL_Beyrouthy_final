# Project_Course_Michael_EL_Beyrouthy_final
This Github page contains the R code to run differential expression analysis of two publicly available datasets obtained from RNA-seq of different T cells in patients and healthy subjects

The data consists of the counts obtained after Illumina sequencing of different T cells in healthy subjects and patients

The pre-processing of the data was performed in Galaxy using thr following tools and options:

- FastQC Trimming (using Fastp) <br/>
- Alignment (using HISAT2) <br/>
- MultiQC Quantification (using FeatureCounts) <br/>
- Differential expression analysis (edgeR) <br/>

The code in R contains the following:
- Read in and format the data for downstream analysis <br/>
- Principle component analysis
- T test and adjustment using Benjamin-Hochberg method
- singnificant genes were extracted
- Heatmap was perfomed
- Volcano plot
- scatter plot of the most up and down regulated genes according to their log fold change (logFC)
- a heat map of the 10 most up and down regulated genes
- Box plot of the raw counts and significat genes obtained from the t test
- Density plot of raw counts and significant genes obtain fromm thee t test
