#!/usr/bin/env Rscript

## ----style, echo=FALSE, message=FALSE, warning=FALSE, results="asis"-------
library("BiocStyle")
library("knitr")
library("rmarkdown")
## options(width = 70) # Andrzej said this is not needed
opts_chunk$set(message = FALSE, error = FALSE, warning = FALSE,
               cache = FALSE, fig.width = 5, fig.height = 5)

## ----loadairway------------------------------------------------------------
library("airway")

## ----dir-------------------------------------------------------------------
indir <- system.file("extdata", package="airway", mustWork=TRUE)

## ----list.files------------------------------------------------------------
list.files(indir)

## ----sampleTable-----------------------------------------------------------

mydir=('.')
#csvfile <- file.path(indir, "sample_table.csv")
csvfile <- file.path(mydir, "sample_table.csv")
sampleTable <- read.csv(csvfile, row.names = 1)
sampleTable

## ----filenames-------------------------------------------------------------
#filenames <- file.path(indir, paste0(sampleTable$Run, "_subset.bam"))
filenames <- file.path(mydir, paste0(sampleTable$Run, ".s.22.bam"))
file.exists(filenames)

## ----Rsamtools-------------------------------------------------------------
library("Rsamtools")
bamfiles <- BamFileList(filenames, yieldSize=2000000)

## ----seqinfo---------------------------------------------------------------
seqinfo(bamfiles[1])

## ----genfeat---------------------------------------------------------------
library("GenomicFeatures")

## ----txdb------------------------------------------------------------------
#gtffile <- file.path(indir,"Homo_sapiens.GRCh37.75_subset.gtf")
#gtffile <- file.path(mydir,"Homo_sapiens.GRCh37.75.gtf")
gtffile <- file.path(mydir,"Homo_sapiens.GRCh37.75.22.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
txdb

## --------------------------------------------------------------------------
ebg <- exonsBy(txdb, by="gene")
ebg

## --------------------------------------------------------------------------
library("GenomicAlignments")
library("BiocParallel")

## --------------------------------------------------------------------------
register(SerialParam())

## --------------------------------------------------------------------------
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )

## ----sumexp, echo=FALSE----------------------------------------------------
par(mar=c(0,0,0,0))
plot(1,1,xlim=c(0,100),ylim=c(0,100),bty="n",
     type="n",xlab="",ylab="",xaxt="n",yaxt="n")
polygon(c(45,90,90,45),c(5,5,70,70),col="pink",border=NA)
polygon(c(45,90,90,45),c(68,68,70,70),col="pink3",border=NA)
text(67.5,40,"assay")
text(67.5,35,'e.g. "counts"')
polygon(c(10,40,40,10),c(5,5,70,70),col="skyblue",border=NA)
polygon(c(10,40,40,10),c(68,68,70,70),col="skyblue3",border=NA)
text(25,40,"rowRanges")
polygon(c(45,90,90,45),c(75,75,95,95),col="palegreen",border=NA)
polygon(c(45,47,47,45),c(75,75,95,95),col="palegreen3",border=NA)
text(67.5,85,"colData")

## --------------------------------------------------------------------------
se
dim(se)
assayNames(se)
head(assay(se), 3)
colSums(assay(se))

## --------------------------------------------------------------------------
rowRanges(se)

## --------------------------------------------------------------------------
str(metadata(rowRanges(se)))

## --------------------------------------------------------------------------
colData(se)

## --------------------------------------------------------------------------
colData(se) <- DataFrame(sampleTable)
colData(se)

## ----secell----------------------------------------------------------------
se$cell
se$dex

## ----sedex-----------------------------------------------------------------
library("magrittr")
se$dex %<>% relevel("untrt")
se$dex

## ----explaincmpass, eval = FALSE-------------------------------------------
#  se$dex <- relevel(se$dex, "untrt")

## --------------------------------------------------------------------------
#data("airway")
#se <- airway

## --------------------------------------------------------------------------
se$dex %<>% relevel("untrt")
se$dex

## --------------------------------------------------------------------------
round( colSums(assay(se)) / 1e6, 1 )

## --------------------------------------------------------------------------
colData(se)

## --------------------------------------------------------------------------
library("DESeq2")

## --------------------------------------------------------------------------
dds <- DESeqDataSet(se, design = ~ cell + dex)

## --------------------------------------------------------------------------
countdata <- assay(se)
head(countdata, 3)

## --------------------------------------------------------------------------
coldata <- colData(se)

## --------------------------------------------------------------------------
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~ cell + dex)

## --------------------------------------------------------------------------
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

## ----meanSdCts-------------------------------------------------------------
lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library("vsn")
meanSdPlot(cts, ranks = FALSE)

## ----meanSdLogCts----------------------------------------------------------
log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)

## ----rlog------------------------------------------------------------------
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

## ----vst-------------------------------------------------------------------
#vsd <- vst(dds, blind = FALSE)
#Error in vst(dds, blind = FALSE): less than 'nsub' rows,
#  it is recommended to use varianceStabilizingTransformation directly
#Traceback:
#
#1. vst(dds, blind = FALSE)
#2. stop("less than 'nsub' rows,\n  it is recommended to use varianceStabilizingTransformation directly")

vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

head(assay(vsd), 3)

## ----rldplot, fig.width = 6, fig.height = 2.5------------------------------
library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
  
colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

## --------------------------------------------------------------------------
sampleDists <- dist(t(assay(rld)))
sampleDists

## --------------------------------------------------------------------------
library("pheatmap")
library("RColorBrewer")

## ----distheatmap, fig.width = 6.1, fig.height = 4.5------------------------
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

## --------------------------------------------------------------------------
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))

## ----poisdistheatmap, fig.width = 6.1, fig.height = 4.5--------------------
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$dex, rld$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

## ----plotpca, fig.width=6, fig.height=4.5----------------------------------
plotPCA(rld, intgroup = c("dex", "cell"))

## --------------------------------------------------------------------------
pcaData <- plotPCA(rld, intgroup = c( "dex", "cell"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))

## ----ggplotpca, fig.width=6, fig.height=4.5--------------------------------
ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()

## ----mdsrlog, fig.width=6, fig.height=4.5----------------------------------
mds <- as.data.frame(colData(rld))  %>%
         cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed()

## ----mdspois, fig.width=6, fig.height=4.5----------------------------------
mdsPois <- as.data.frame(colData(dds)) %>%
   cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed()

## ----airwayDE--------------------------------------------------------------
dds <- DESeq(dds)

## --------------------------------------------------------------------------
res <- results(dds)
res

## --------------------------------------------------------------------------
res <- results(dds, contrast=c("dex","trt","untrt"))

## --------------------------------------------------------------------------
mcols(res, use.names = TRUE)

## --------------------------------------------------------------------------
summary(res)

## --------------------------------------------------------------------------
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)

## --------------------------------------------------------------------------
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

## --------------------------------------------------------------------------
results(dds, contrast = c("cell", "N061011", "N61311"))

## ----sumres----------------------------------------------------------------
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))

## --------------------------------------------------------------------------
sum(res$padj < 0.1, na.rm=TRUE)

## --------------------------------------------------------------------------
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])

## --------------------------------------------------------------------------
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

## ----plotcounts------------------------------------------------------------
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("dex"))

## ----ggplotcountsjitter, fig.width = 4, fig.height = 3---------------------
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex","cell"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = cell)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

## ----ggplotcountsgroup, fig.width = 4, fig.height = 3----------------------
ggplot(geneCounts, aes(x = dex, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()

## ----plotma----------------------------------------------------------------
res <- lfcShrink(dds, contrast=c("dex","trt","untrt"), res=res)
plotMA(res, ylim = c(-5, 5))

## ----plotmaNoShr-----------------------------------------------------------
res.noshr <- results(dds)
plotMA(res.noshr, ylim = c(-5, 5))

## ----plotmalabel-----------------------------------------------------------
plotMA(res, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

## ----histpvalue2-----------------------------------------------------------
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

## --------------------------------------------------------------------------
library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)

## ----genescluster----------------------------------------------------------
mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, c("cell","dex")])
pheatmap(mat, annotation_col = anno)

## ----sensitivityovermean, fig.width=6--------------------------------------
qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(resLFC1$pvalue, bins, function(p)
                          mean(p < .05, na.rm = TRUE))
barplot(fractionSig, xlab = "mean normalized count",
                     ylab = "fraction of small p values")

## --------------------------------------------------------------------------
library("AnnotationDbi")
library("org.Hs.eg.db")

## --------------------------------------------------------------------------
columns(org.Hs.eg.db)

## --------------------------------------------------------------------------
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

## --------------------------------------------------------------------------
resOrdered <- res[order(res$pvalue),]
head(resOrdered)

## ----eval=FALSE------------------------------------------------------------
#  resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
#  write.csv(resOrderedDF, file = "results.csv")

## ----eval=FALSE------------------------------------------------------------
#  library("ReportingTools")
#  htmlRep <- HTMLReport(shortName="report", title="My report",
#                        reportDirectory="./report")
#  publish(resOrderedDF, htmlRep)
#  url <- finish(htmlRep)
#  browseURL(url)

## --------------------------------------------------------------------------
resGR <- results(dds, lfcThreshold = 1, format = "GRanges")
resGR

## --------------------------------------------------------------------------
resGR$symbol <- mapIds(org.Hs.eg.db, names(resGR), "SYMBOL", "ENSEMBL")

## --------------------------------------------------------------------------
library("Gviz")

## --------------------------------------------------------------------------
window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)

## --------------------------------------------------------------------------
status <- factor(ifelse(resGRsub$padj < 0.1 & !is.na(resGRsub$padj),
                     "sig", "notsig"))

## ----gvizplot--------------------------------------------------------------
options(ucscChromosomeNames = FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name = "gene ranges", feature = status)
d <- DataTrack(resGRsub, data = "log2FoldChange", baseline = 0,
               type = "h", name = "log2 fold change", strand = "+")
plotTracks(list(g, d, a), groupAnnotation = "group",
           notsig = "grey", sig = "hotpink")

## --------------------------------------------------------------------------
library("sva")

## --------------------------------------------------------------------------
dat  <- counts(dds, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ dex, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)
svseq$sv

## ----svaplot---------------------------------------------------------------
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds$cell, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
 }

## --------------------------------------------------------------------------
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + dex

## ----svaDE-----------------------------------------------------------------
ddssva %<>% DESeq

## --------------------------------------------------------------------------
library("fission")
data("fission")
ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)

## ----fissionDE-------------------------------------------------------------
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ strain + minute)
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),], 4)

## ----fissioncounts, fig.width=6, fig.height=4.5----------------------------
fiss <- plotCounts(ddsTC, which.min(resTC$padj), 
                   intgroup = c("minute","strain"), returnData = TRUE)
ggplot(fiss,
  aes(x = as.numeric(minute), y = count, color = strain, group = strain)) + 
  geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10()

## --------------------------------------------------------------------------
resultsNames(ddsTC)
res30 <- results(ddsTC, name="strainmut.minute30", test="Wald")
res30[which.min(resTC$padj),]

## --------------------------------------------------------------------------
betas <- coef(ddsTC)
colnames(betas)

## ----fissionheatmap--------------------------------------------------------
topGenes <- head(order(resTC$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)

## --------------------------------------------------------------------------
devtools::session_info()

