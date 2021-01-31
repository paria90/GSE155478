setwd("c:/users/volkan/Desktop/R related files/GSE155478")

library(Biobase)
library(limma)
library(reshape2)
library(ggplot2)
library(gplots)
library(GEOquery)
library(plyr)
library(pheatmap)
#1:
series <- "GSE155478"
gset <- getGEO(series,GSEMatrix = TRUE,AnnotGPL = TRUE,destdir = "data/")
length(gset)
class(gset)
names(gset)
platform <- "GSE155478"
if(length(gset)>1) {idx <- grep(platform,attr(gset,"names"))} else {idx <- 1}
gset <- gset[[idx]]
gr <- c(rep("ADRres",3),rep("ADRsen",3))
length(gr)

#2:
ex <- exprs(gset)
dim(ex)
max(ex)
min(ex)

#3_1:
pdf("results/boxplot.pdf",width = 6)
boxplot(ex)
dev.off()

#3_2:
pdf("results/CoreHeatmap.pdf")
pheatmap(cor(ex),labels_row = gr,labels_col = gr)
dev.off

#3_3:PCA for genes
ex.scale <- t(scale(t(ex),scale = FALSE))
pc <- prcomp(ex.scale)
pdf("results/pc_scale.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

#3_4:PCA for samples
pcr <- data.frame(pc$rotation[,1:3],Group=gr)
head(pcr)
pdf("results/pca_samples.pdf")
ggplot(pcr,aes(PC1,PC2,color=Group)) + geom_point(size=3)+theme_bw()
dev.off()

#4:
gr <- as.factor(gr)
gset$description <- gr
design <- model.matrix(~ description +0,gset)
colnames(design) <- levels(gr)
head(design)
fit <- lmFit(gset,design)
cont.matrix <- makeContrasts(ADRres-ADRsen,levels = design)
fit2 <- contrasts.fit(fit,cont.matrix) 
fit2 <- eBayes(fit2,0.01)
tT <- topTable(fit2,adjust.method = "fdr",sort.by = "B",number = Inf)
head(tT)
colnames(tT)
head(tT$GeneID)
tT <- subset(tT,select=c("GeneID","GeneSymbol","logFC","adj.P.Val"))
head(tT)

#5:
ADR.res <- subset(tT,logFC > 1 & adj.P.Val < 0.05)
ADR.res.genes <- unique(ADR.res$GeneSymbol)
write.table(ADR.res.genes,file = "results/ADRres-ADRsen-high.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)

ADR.sen <- subset(tT,logFC < -1 & adj.P.Val < 0.05)
ADR.sen.genes <- unique(ADR.sen$GeneSymbol)
write.table(ADR.sen.genes,file = "results/ADRres-ADRsen-low.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
