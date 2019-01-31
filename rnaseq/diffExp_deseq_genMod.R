#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)
library( "biomaRt" )
library("pheatmap")
library("RColorBrewer")
library( "genefilter")


#######LOAD SAMPLES#####
directory <- "~/uproov_wd/rnaseq/readcounts/gen_modified/"


###ALL

sampleFiles <- grep("htseq",list.files(directory),value=TRUE)
sampleCondition <- sapply(sampleFiles, function(x){
    sample <- unlist(strsplit(x, '_'))[1]
    if (grepl("-P53", sample)){
        cond <- "P53"
    } else if (grepl("-R", sample)){
        cond <- "P53_RB1"
    } else {
        cond <- "WT"
    } 
    cond
})

sampleNames <- sub("(.*)(_htseq_counts.txt)", "\\1", sampleFiles)

sampleTable <- data.frame(sampleName  = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleCondition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)

#########FILTERING#######

#Remove genes with zero counts
dds <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 1, ]


ensembl = useMart( "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="grch37.ensembl.org")
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = rownames(dds),
                  mart = ensembl )
genemap <- genemap[match(rownames(dds),genemap$ensembl_gene_id),]

rownames(dds) <-genemap$hgnc_symbol
dds <- dds[!is.na(rownames(dds)),]
dds <- dds[rownames(dds) != "", ]
dds <- dds[!duplicated(rownames(dds)),]




#######ANALYSIS####
dds$condition <- factor(dds$condition, levels = c("WT","P53", "P53_RB1"))
dds <- DESeq(dds)

resultsNames(dds)


###Diff Exp

wt_p53 <- results(dds, contrast=c("condition", "WT", "P53"))
wt_RB1 <- results(dds, contrast=c("condition", "WT", "P53_RB1"))
p53_RB1 <- results(dds, contrast=c("condition", "P53", "P53_RB1"))


##VOLCANO 
 pretty_volcano <- function(resObject, title) {
     alpha <- 0.01 # Threshold on the p-value
     
     resObject$sig <- -log10(resObject$padj)
     resObject[is.infinite(resObject$sig),"sig"] <- 350
     # View(resultDESeq2[is.na(resultDESeq2$pvalue),])
     
     # Select genes with a defined p-value (DESeq2 assigns NA to some genes)
     genes.to.plot <- !is.na(resObject$padj)
     # sum(genes.to.plot)
     range(resObject[genes.to.plot, "log2FoldChange"])
     
     ## Volcano plot of adjusted p-values
     cols <- densCols(resObject$log2FoldChange, resObject$sig)
     cols[resObject$padj ==0] <- "purple"
     resObject$pch <- 19
     resObject$pch[resObject$padj ==0] <- 6
     plot(resObject$log2FoldChange, 
          resObject$sig, 
          col=cols, panel.first=grid(),
          main=title, 
          xlab="Effect size: log2(fold-change)",
          ylab="-log10(adjusted p-value)",
          pch=resObject$pch, cex=0.4)
     abline(v=0)
     abline(v=c(-1,1), col="brown")
     abline(h=-log10(alpha), col="brown")
     
     ## Plot the names of a reasonable number of genes, by selecting those begin not only significant but also having a strong effect size
     gn.selected <- abs(resObject$log2FoldChange) > 2 & resObject$padj < alpha 
     text(resObject$log2FoldChange[gn.selected],
          -log10(resObject$padj)[gn.selected],
          lab=rownames(resObject)[gn.selected ], cex=0.6)
     pval= 0.01
     mtext(paste("pval adj =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
     mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
 }

simple_volcano <- function(resObject, title){
    par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
    
    topT <- as.data.frame(resObject)
    
    #Adjusted P values (FDR Q values)
    with(topT, plot(log2FoldChange, -log10(padj), pch=20, main=title, cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
    
    with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

    #Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
    abline(v=0, col="black", lty=3, lwd=1.0)
    abline(v=-2, col="black", lty=4, lwd=2.0)
    abline(v=2, col="black", lty=4, lwd=2.0)
    abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
    mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1) 
    mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc), cex = 0.8, line = 0.5)
}
pdf("~/uproov_wd/rnaseq/genMod/diffexp_genmod_volcanos.pdf")
pretty_volcano(wt_p53, 'WT vs P53mut')
simple_volcano(wt_p53, 'WT vs P53mut')
pretty_volcano(wt_RB1, 'WT vs P53+RB1mut')
simple_volcano(wt_RB1, 'WT vs P53+RB1mut')
pretty_volcano(p53_RB1, 'P53mut vs P53+RB1mut')
simple_volcano(p53_RB1, 'WTmut vs P53+RB1mut')
dev.off()


##WRITE TO CSV
wt_p53_ord <- wt_p53[order(wt_p53$padj),]
wt_RB1_ord <- wt_RB1[order(wt_RB1$padj),]
p53_RB1_ord <- p53_RB1[order(p53_RB1$padj),]

write.csv(as.data.frame(wt_p53_ord), file="~/uproov_wd/rnaseq/genMod/wt_p53_diffExp.csv")
write.csv(as.data.frame(wt_RB1_ord), file="~/uproov_wd/rnaseq/genMod/wt_RB1_diffExp.csv")
write.csv(as.data.frame(p53_RB1_ord), file="~/uproov_wd/rnaseq/genMod/p53_RB1_diffExp.csv")

# ######DISTANCE PLOT######
# 
rlds <- rlog(dds, blind=FALSE)
rld <- assay(rlds)

colors <- colorRampPalette(brewer.pal(8, "RdBu"))(255)
revcolors <- rev(colors)


###VARIABLE GENES
rld_zero = rld[apply(rld, 1, function(row) all(row !=0 )),]
rld_sd = rld_zero[order(apply(rld_zero, 1, sd),decreasing = T),]


mil <- head(rld_sd, n = 1000)
milspecor <- cor(mil, use="complete.obs", method = "spearman")
milspedist<- dist(x=milspecor)
milclust=hclust(milspedist, method = "complete")
pdf("~/uproov_wd/rnaseq/genMod/spearman_distance_genMod_1000.pdf", width=14, height = 14)
pheatmap(as.matrix(milspecor),
         cluster_rows=milclust,
         cluster_cols=milclust,
         col=revcolors)
dev.off()

