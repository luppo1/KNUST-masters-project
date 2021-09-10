## This is HISAT2 aligner

# Load all required packages
library("DESeq2")
library("edgeR")
library("limma")
library("sva")
library("data.table")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("affy")
library("ggfortify")
library("Hmisc")
library("corrplot")
library("reshape2")
library("geneplotter")
library("VennDiagram")



# Set working directory
setwd("C:/Users/Albert Doughan/Desktop/My Project/R Datasets/Hisat/")
dir = "C:/Users/Albert Doughan/Desktop/My Project/R Datasets/Hisat/"

# Import metadata
metadatah = read.csv(file= "metadata.csv", header=T, sep = ",")
head(metadatah)

# Reading counts data from featureCounts
counts = read.csv(file = "hisat2_count.csv", header = T, sep = ",")
head(counts)

# Remove the Gene ID column
countdata <- counts[, -c(1)]

# Making "Geneid" column the rownames
rownames(countdata) <- counts[,1]
head(countdata)

# Check if the metadata and samples have the same names
table(colnames(countdata)==metadatah$SampleID)



##########################################
# Running differential expression analysis with DESeq2
##########################################

# Create the DESeqDataSet object from Matrix of counts and metadata
dds <- DESeqDataSetFromMatrix(countData = round(countdata), 
                              colData = metadatah, 
                              design = ~Condition)

nrow(dds) 


# Hidden batch effect detection and removal 
dds <- DESeq(dds)
dat1 <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat1) > 1
dat1 <- dat1[idx,]
mod <- model.matrix(~ as.factor(Condition), colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))

# Using the SVA package
svseq <- svaseq(dat1, mod, mod0, n.sv = 2)
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2  + Condition

# Visualize the batches
par(mfrow=c(2,1),oma = c(0, 0, 2, 0))
stripchart(svseq$sv[,1] ~ dds$Condition,vertical=TRUE,main="SV1")
abline(h=0)
stripchart(svseq$sv[,2] ~ dds$Condition,vertical=TRUE,main="SV2")
abline(h=0)

# Remove Genes with low counts
# Here we perform a minimal pre-filtering to keep only rows that have at least 10 reads total.
dds1 <- ddssva[rowSums(counts(ddssva)) > 100,]
nrow(dds1)

# Run DESeq function on the data to perform differential gene expression analysis
dds1 <- DESeq(dds1)
head(assay(dds1))

# Building out results table
res_table <- results(dds1)
summary(res_table)

# Working with alpha 0.05
res2 <- results(dds1, alpha=0.05)
summary(res2)

# How many adjusted p-values were less than 0.05?
sum(res2$padj < 0.05, na.rm=TRUE)

# We order our results table by the smallest p value:
res_small_p <- res2[order(res2$pvalue),]

# Select genes with p less than 0.05
res_sig <- subset(res_small_p, padj < 0.05)
dim(res_sig)

# Write final list to file
write.csv(as.data.frame(res_sig), "hisat_deseq_project.csv")


#####################
## DATA EXPLORATION
#####################
# Principal component analysis
PCAdata <- prcomp(t(assay(dds1)))
autoplot(PCAdata, data = metadatah, colour = "Condition",label = FALSE, size = 5)+
  theme_bw() +
  labs(colour="Condition")+
  theme(legend.title = element_text(size = 21),
        legend.text = element_text(size = 20),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        axis.text=element_text(size=21))


# Hierarchical clustering
clusters2 <- hclust(dist(t(assay(dds1))), method ="ward.D")
plot(clusters2, labels = FALSE)

# Density plot
plotDensity(assay(dds1), col=1:24,lwd=2,lty=1,xlab("Density"),ylab("Counts"))

# Data normalization
vst = vst(dds1, blind=FALSE)

# Visualize transformed data
par(mfrow=c(1,2))
plot(assay(dds1))
plot(assay(vst))

# Principal component analysis after normalization
PCAdata1 <- prcomp(t(assay(vst)))
autoplot(PCAdata1, data = metadatah, colour = "Condition",label = FALSE, size = 5)+
  theme_bw() +
  labs(colour="Condition")+
  theme(legend.title = element_text(size = 21),
        legend.text = element_text(size = 20),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        axis.text=element_text(size=21))

# Hierarchical clustering after normalization
clusters1 <- hclust(dist( t( assay(vst) ) ),method ="ward.D")
plot(clusters1, label = FALSE)

# Density plot after normalization
plotDensity(assay(vst), lwd=2,lty=1,xlab("Density"),ylab("Counts"), main = "Density plot")

# Heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(vst)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(vst$SampleID, vst$Condition, sep="-" )
colnames(sampleDistMatrix) <- paste(dds1$SampleID, dds1$Condition, sep="-")
colors <- colorRampPalette(rev(brewer.pal(9, "Reds")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, main = "Sample Distance Matrix ")


# We can use plotCounts fxn to compare the normalized counts
#between treated and control groups for our top 6 genes
par(mfrow=c(2,3))
plotCounts(dds1, gene="ENSG00000198648", intgroup="Condition")
plotCounts(dds1, gene="ENSG00000152642", intgroup="Condition")
plotCounts(dds1, gene="ENSG00000154165", intgroup="Condition")
plotCounts(dds1, gene="ENSG00000196549", intgroup="Condition")
plotCounts(dds1, gene="ENSG00000113916", intgroup="Condition")
plotCounts(dds1, gene="ENSG00000163453", intgroup="Condition")

## Spearman's correlation plots
## Compute correlation matrix
resu <- cor(assay(vst), method = "spearman", use = "complete.obs")
round(resu, 2)

## compute the significance levels
resa = rcorr(assay(vst), type = "spearman")

# Heatmap of sample-to-sample distances  
col<- colorRampPalette(c("blue", "white", "red"))(100)
heatmap(x = resu, col = col, symm = TRUE)

## MA-plot
plotMA(res_table, colSig = "red3", colNonSig = "gray32")






##########################################
# Running differential expression analysis with edgeR
##########################################

# Creating the edgeR gene list
y <- DGEList(counts=countdata,group=factor(metadatah$Condition))
dim(y)

# Remove low count genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y)

# After filtering, it is a good idea to reset the library sizes:
y$samples$lib.size <- colSums(y$counts)
y$samples

# Normalizing the data
y <- calcNormFactors(y)
y$samples

# Plot MDS
col=c(rep("black",50), rep("red",50))
plotMDS(y, labels = NULL, pch = 16, cex = 1, col = col)
legend("top", c("Normal","Disease"), pch = 16, col = c("black","red"))

# The Design matrix
group = y$samples$group
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

# Estimating the Dispersion
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion

# Plot the dispersion
plotBCV(y)

# Model the quasi-likelihood dispersions of variation within model
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)

# Testing for differential expression
# Set the contrast to be tested
DvsN <- makeContrasts(Disease-Normal, levels=design)
res <- glmQLFTest(fit, contrast=DvsN)
topTags(res)

# The total number of differentially expressed genes at FDR< 0:05 is:
is.de <- decideTestsDGE(res)
summary(is.de)

## Plot log-fold change against log-counts per million, with DE genes highlighted
plotMD(res, status = is.de, cex = .4)
abline(h=c(-1, 1), col="black")

# Order by FDR-corrected p-values
deg <- as.data.frame(topTags(res, n=Inf))
deg = as.data.frame(deg)
order_res <- deg[order(deg$FDR),]
dim(order_res)

# Select only genes with FDR-corrected p-values less that 0.05
sig_deg <- subset(order_res, FDR < 0.05)
dim(sig_deg)

# Write final list to file
write.csv(as.data.frame(sig_deg), "hisat_edger_project.csv")




##########################################
# Running differential expression analysis with limma+voom
##########################################

#Creating the gene list through edgeR
dge <- DGEList(countdata)
dim(dge)

# Transformations from the raw-scale
cpm <- cpm(dge)
lcpm <- cpm(dge, log=TRUE)
L <- mean(dge$samples$lib.size) * 1e-6
M <- median(dge$samples$lib.size) * 1e-6
c(L, M)

# Number of genes with 0 count in all samples 
table(rowSums(dge$counts==0)==100)

# Removing genes that are lowly expressed
keep.exprs <- filterByExpr(dge, group=metadatah$Condition)
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge)


# Plot density before and after filtering
cpm1 <- cpm(dge)
lcpm1 <- cpm(dge, log=TRUE)
lcpm.cutoff <- log2(10/M + 2/L)
plotDensity(lcpm, main = "A. Raw data", xlab = "Log CPM")
abline(v=lcpm.cutoff, lty=3)
plotDensity(lcpm1, main = "B. Filtered data", xlab = "Log CPM")
abline(v=lcpm.cutoff, lty=3)

# Data normalization
dge1 <- calcNormFactors(dge, method = "TMM")
dge1$samples$norm.factors

# Box plots before and after normalization
nsamples <- ncol(dge1)
col <- brewer.pal(12, "Paired")
un <- cpm(dge, log=TRUE)
n <- cpm(dge1, log=TRUE)
boxplot(un, main = "A: Unnormalized data", col=col)
boxplot(n, main = "B: Normalized data", col=col)

# Unsupervised clustering
col=c(rep("blue",50), rep("red",50))
plotMDS(dge1, labels = NULL, pch = 16, cex = 1, col = col)
legend("top", c("Normal","Disease"), pch = 16, col = c("blue","red"))

# Create a design matrix
design <- cbind("1"=1,"1vs2"=rep(c(1,2), each = nrow(metadatah)/2))

# Running the limma voom function
v <- voom(dge1, design, plot=TRUE, normalize="quantile")

# After this, the usual limma pipelines for differential expression is be applied.
fit <- lmFit(v, design)
fit <- eBayes(fit)
res <- topTable(fit, coef=ncol(design),number=Inf)
summary(decideTests(fit))
res_pvalue <- as.data.frame(subset(res, adj.P.Val < 0.05))
dim(res_pvalue)

# We order our results table by the smallest p value:
order_res <- res_pvalue[order(res_pvalue$adj.P.Val),]
dim(order_res)

# Display the top 6 most significant genes
topTreat(fit, coef=1, n=6)

# Write final list to file
write.csv(as.data.frame(order_res), "hisat_limma_project.csv")




## All HISAT common
des1 = read.csv("hisat_deseq_project.csv", header = T)
head(des1)  
edg1 = read.csv("hisat_edger_project.csv", header = T)
head(edg1)
lim1 = read.csv("hisat_limma_project.csv", header = T)
head(lim1)

inter1 = intersect(intersect(des1$X, edg1$X), lim1$X)
length(inter1)

write.csv(inter1, "hisat_common_project.csv")


## All Hisat common
x = list(des1$X, edg1$X, lim1$X)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = x,
  category.names = c("DESeq2" , "edgeR " , "limma"),
  filename = 'hisat_common.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#999999", "#E69F00", "#56B4E9"),
  
  # Numbers
  cex = .6,
  #fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)



