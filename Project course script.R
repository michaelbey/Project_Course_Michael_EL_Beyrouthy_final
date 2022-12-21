install.packages("ggplot2")
library(ggplot2)
install.packages("ggfortify")
library(ggfortify)
install.packages("dplyr")
library(dplyr)
install.packages("ggrepel")
library(ggrepel)
install.packages("pheatmap")
library(pheatmap)
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("edgeR")
library(edgeR)
BiocManager::install("affy")
library(affy)

table1 <- read.table("edgeR_normcounts.tabular", header = TRUE) #making a table out of the normalized data that we got out of Galaxy
table2 <- read.table("edgeR_Verboom-Roels.tabular", header = TRUE) 
annotate1 <- read.table("annotateMyIDs_edgeR_normcounts.tabular", header = TRUE) #Opening the annotation file as a data frame -> the file contain the names of the genes according to their access code -> making it easier to read
annotate2 <- read.table("annotateMyIDs_edgeR_Verboom-Roels.tabular", header = TRUE)

table1 <- na.omit(table1) #removing the non existant values to be able to apply math on the data frame
annotate1 <- na.omit(annotate1)
table2 <-  na.omit(table2)
annotate2 <- na.omit(annotate2)

merged <- merge(table1, annotate1, by.x = "GeneID", by.y = "ENTREZID") #merging the normalized data file and the annotation file so that each gene will have it's name in a comlumn on the side
merged2 <- merge(table2, annotate2, by.x = "GeneID", by.y = "ENTREZID")

#name the rows by GeneID on the merged file
rownames(merged) <- merged$SYMBOL #changing the row names into the gene name
rownames(merged2) <- merged2$SYMBOL

#remove the first column as the rows are named after GeneID -> removing the Gene ID number and the gene name as the rows are names after the gene name
merged <- merged[, -1] 
merged <- merged[, -9]
merged2 <- merged2[, -1]
merged2 <- merged2[, -6]


#--------------------------------------------------------------------------------------
#PCA

#before a PCA metadata file would be useful -  allows us to change the colors according to the groups
#metadata file  protocol
md <- data.frame(colnames(merged)) #making it so that the row name in the metadata is the same as the columns in the normalized and annotated data file -> so that we can classify the data ccording to weither its a control or a patient sample
colnames(md)[1] <- "SampleID" #naming the first column sample ID so that we can add wether its from a cancer patient or from control
md$group <- rep(c("control", "cancer"), times=c(4,4)) #added a column that was named group, allowing to classify the data as either cancer or contorl
md$Study <- rep(c("Roels", "Verboom"), times=c(4,4)) #we could also add a column and insert from which study the data was obtained
md <- md[, -3]


#PCA
PCA <- prcomp(t(merged)) #PCA file -> The t is so that we inverse the data instead of getting foor every gene we will get it for every sample accordng to the principle components of the gene
autoplot(PCA) #PCA plot

autoplot(PCA, data = md, col = "group", size = 4) #The autoplot function will make a polt based on the principal components

autoplot(PCA, data=md, colour="group", size=4, frame=TRUE, frame.type = 'norm')+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) #more complex PCA plot with more annotations


#PCA David's method
#Principal Component Analysis (PCA)
#------------------------------------------------------------------
# Apply PCA on the transposed data
output_pca_exercise <- prcomp(t(merged))

#-----------------------------------------------------------
#Scree plot
#-----------------------------------------------------------
# Calculate the variance of each PC
var_explained = output_pca_exercise$sdev^2 / sum(output_pca_exercise$sdev^2)
var_explained <- round(var_explained, digits = 2)

var_explained_df <- data.frame(PC = paste0("PC", 1:6),
                               var_explained = var_explained[1:6])

# Create the screen plot
ggplot(var_explained_df, aes(x = PC, y = var_explained)) +
  geom_col() + geom_bar(stat = "identity", fill = "royalblue4") +
  labs(title = "Scree Plot of Principal Components") + 
  ylab("Explained variance") + ylim(0, 1)

#------------------------------------------------------------------
#PCA plot
#----------------------------------------------------------------

# Extract PCs 1 and 2
output_pca_exercise <- as.data.frame(output_pca_exercise$x[, 1:2])

# Add metadata variables to the PCA output
output_pca_exercise$Group <- md$group
output_pca_exercise$SampleID <- md$SampleID

# Plot samples by their PC3 and PC4
ggplot(output_pca_exercise, aes_string(x = colnames(output_pca_exercise)[1], 
                                       y = colnames(output_pca_exercise)[2], 
                                       color = "SampleID", 
                                       shape = "Group")) +
  ggtitle("PCA Plot") +
  geom_point(size = 3, alpha = 1) +
  theme(text = element_text(size = 3)) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab(paste(colnames(output_pca_exercise)[1], var_explained[1], sep = ": ")) +
  ylab(paste(colnames(output_pca_exercise)[2], var_explained[2], sep = ": "))














#-----------------------------------------------------------
#T test

#-----------------------------------------------------------
#T test

#pvalue - empty file in order to add the p values to it
pvalue <-data.frame(p_value=rep(0,nrow(merged)))

#rowname in pvalue file the same as the in merged file
rownames(pvalue) <- rownames(merged)

#pvalue algorithm via t test
for (i in 1:nrow(merged)){pvalue$p_value[i] <-t.test(merged [i,1:4], 
                                                     merged [i,5:8], 
                                                     alternative = "two.sided")$p.value
}



#Benjamin test for 
pvalue$padjust <- p.adjust(as.numeric(unlist(pvalue)), method = "BH")


#LogFC
logFC <- (rowMeans(merged[1:13323,1:4]))-rowMeans((merged[1:13323,5:8]))
pvalue$LogFC <- logFC


#significant p value <0.05 according to the t test
sig_pval_t_test <- pvalue[pvalue$padjust<0.05,,drop=FALSE]

sig_pval_t_test$p_value = NULL
sig_pval_t_test$LogFC = NULL

pvalue$padjust<0.05 #just to view the p valaue in the console, quite useless

#steps to remove non significant genes from the merged file
# step 1 make a vector
names_sig_genes_t_test <-rownames(sig_pval_t_test)

#step 2 make a file containing only the significant genes
sig_genes_t_test <-subset(merged,rownames(merged) %in% names_sig_genes_t_test)

#before making a heat map you should make a data matrix
sig_genes_t_test <- data.matrix(sig_genes_t_test)

#then you make a heatmap
heatmap(sig_genes_t_test)



#--------------------------------------------------------------
#Volcano plot

#We need a data frame with all adjusted pvalues (padjusted) in a file called FDR_pvalue
FDR_pvalue <- data.frame(pvalue$padjust)

#You also need the Log2FC ------------ When applying this formula, it would be of note that you are claculating the rowmeans of the first 4 column - the rowmeans of the 4 second columns -> Not the rownmeans of the substraction
logFC <- (rowMeans(merged[1:13323,1:4]))-rowMeans((merged[1:13323,5:8]))
#now add the log2FC to the FDR_pvalue file
FDR_pvalue$logFC <- 0
FDR_pvalue$logFC <- logFC

pvalue$LogFC <- logFC

#now name the rows according tp the genes
row.names(FDR_pvalue) <- row.names(pvalue)

FDR_pvalue <- na.omit(FDR_pvalue)
is.na(FDR_pvalue)


# Standard volcano plot
p <- FDR_pvalue %>%
  ggplot(mapping = aes(x=logFC, y=-log10(pvalue.padjust)))+
  geom_point()+
  theme_minimal()

p

# Add a line to highlight the padjusted value of 0.05 and a threshold for the logFC. 

p2 <- p +
  geom_hline(yintercept = -log10(0.05), col="red") +
  geom_vline(xintercept= c(-2, 2), col="red")

p2


# Color the dots based on their differential expression
FDR_pvalue$diffexpressed <- "NO" # create a new column in the data frame and fill it with NO
FDR_pvalue$diffexpressed[FDR_pvalue$logFC > 2 & FDR_pvalue$pvalue.padjust < 0.05] <- "UP"
FDR_pvalue$diffexpressed[FDR_pvalue$logFC < -2 & FDR_pvalue$pvalue.padjust < 0.05] <- "DOWN"



#We then redo the volcano plot
p <- FDR_pvalue %>%
  ggplot(mapping=aes(x=logFC, y=-log10(pvalue.padjust), col=diffexpressed))+
  geom_point()+
  theme_minimal()

#and redo the second volcano plot based with the lines added
p2 <- p +
  geom_hline(yintercept = -log10(0.05), col="red")+
  geom_vline(xintercept= c(-2,2), col="red")


#here we created a file called mycolors where the gene are according to the up or downregulation and color
mycolors <- c("blue", "red", "black") # specifiy the colors you want to use

names(mycolors) <- c("DOWN", "UP", "NO") # add column names


#did a third volcano plot with the colors
p3 <- p2+
  scale_color_manual(values=mycolors)

p3

#----------------------------------------------------------
#Scatter plot of the most up and down regulated genes according to their logFC

sig_pvalue_FDR <- FDR_pvalue[FDR_pvalue$pvalue.padjust<0.05,,drop=FALSE] # keep only FDR p value<0.05

sig_pvalue_FDR <- sig_pvalue_FDR[order(sig_pvalue_FDR$logFC),,drop=FALSE] # order the genes based on logFC

top_genes <- sig_pvalue_FDR[c(1:10, 4229:4238),] # extract top 10 upregulated and top 10 downregulated

top_genes <- rownames(top_genes) # extract the names of the genes that change the most

top_genes_expression <- sig_genes_t_test[(rownames(sig_genes_t_test) %in% top_genes),] # subset the gene expression matrix

top_genes_expression <- data.frame(top_genes_expression) #an error was showing up until turning into a dataframe

top_genes_expression$ControlMean<-rowMeans(top_genes_expression[,1:4]) # add control mean 
top_genes_expression$PatientMean<-rowMeans(top_genes_expression[,5:8]) # add patient mean

head(top_genes_expression)

top_genes_expression$genename<-rownames(top_genes_expression) # add a column with genename to use as labels for the plot


ggplot(top_genes_expression, aes(x=ControlMean, y=PatientMean, label=genename))+
  geom_point()+
  geom_abline(intercept = 0)+
  theme_bw()+
  ggtitle(label="Top 10 up and downregulated genes")+
  aes(color=PatientMean)+
  geom_text_repel()+
  theme(plot.title = element_text(hjust = 0.5, size=20))

pheatmap(top_genes_expression[1:8])


#---------------------------------------------------------------------------------------------
# 2 graphs from david

#Box Plot
sig_genes_t_test_dataframe <- data.frame(sig_genes_t_test) #turn the result of the t test into a dataframe
sig_genes_t_test_dataframe[sig_genes_t_test_dataframe < 0] <- NA #remove negative counts as they are not allowed in dglists
sig_genes_t_test_dataframe2 <- na.omit(sig_genes_t_test_dataframe)


#make a dglist in order to do a boxplot
dglist1 <- DGEList(counts = sig_genes_t_test_dataframe2, genes = rownames(sig_genes_t_test_dataframe2))

boxplot(cpm(dglist1, log = TRUE), col = "indianred1", cex.axis = 0.7, 
        las = 1, ylab = "cpm of counts", main = "Counts")

#make a new file containing the raw counts
merged2 <- merged
merged2[merged2 < 0] <- NA #remove negative counts
merged2 <- na.omit(merged2)

#make a dglist
dglist2 <- DGEList(counts = merged2, genes = rownames(merged2))

boxplot(cpm(merged2, log = TRUE), col = "tan1", cex.axis = 0.7, 
        las = 1, ylab = "cpm of counts", main = "Counts")





#Density plot
#density plot of raw data dglist
plotDensity(cpm(dglist2, log = TRUE), col = "indianred1",
            xlab = "CPM of counts", ylab = "Density",
            main = "Raw counts")

#density plot of counts after t test
plotDensity(cpm(dglist1, log = TRUE), col = "tan1",
            xlab = "CPM of counts", ylab = "Density",
            main = "T test counts")


