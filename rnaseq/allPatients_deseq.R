#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)
library( "biomaRt" )
library("pheatmap")
library("RColorBrewer")
library( "genefilter")


#######LOAD SAMPLES#####
#DIRECTORY WITH ALL INDIVIDUAL HTSEQ READ COUNT FILES
directory <- "."

sampleFiles <- grep("htseq",list.files(directory),value=TRUE)

#I don't have power for real DE analysis, so no conditions

#COMMON
sampleNames <- sub("(.*)(_htseq_counts.txt)", "\\1", sampleFiles)

sampleTable <- data.frame(sampleName  = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ 1)
#63677 genes

#########FILTERING#######

#Remove genes with zero counts
dds <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 1, ]
#41966 genes


ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
#str(ensembl)
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = rownames(dds),
                  mart = ensembl )
genemap <- genemap[match(rownames(dds),genemap$ensembl_gene_id),]

rownames(dds) <-genemap$hgnc_symbol
dds <- dds[!is.na(rownames(dds)),]
#39597 genes
dds <- dds[rownames(dds) != "", ]
#27884 genes
dds <- dds[!duplicated(rownames(dds)),]
#27877 genes



#######ANALYSIS####

dds <- DESeq(dds)


res <- results(dds)
res
resultsNames(dds)
resOrdered <- res[order(res$padj),]

counts <- counts(dds, normalized=TRUE)
genemap2 <- genemap[match(rownames(counts),genemap$ensembl_gene_id),]
write.csv(counts, file="~/uproov_wd/rnaseq/normalizedCounts.csv")


# ######DISTANCE PLOT######
# 
rlds <- rlog(dds, blind=FALSE)
rld <- assay(rlds)
specor <- cor(rld, use="complete.obs", method = "spearman")
revcolors <- rev(colors)
spedist <- dist(x = t(specor), method = "e")
clust=hclust(spedist, method = "complete")
#ALL GENES
pdf("~/uproov_wd/rnaseq/spearman_distance.pdf", width=14, height = 14)
pheatmap(as.matrix(specor),
         cluster_rows=clust,
         cluster_cols=clust,
         col=colors)
dev.off()


### TOP 5000 VARIABLE GENES
rld_zero = rld[apply(rld, 1, function(row) all(row !=0 )),]
rld_sd = rld_zero[order(apply(rld_zero, 1, sd),decreasing = T),]

cincomil <- head(rld_sd, n = 5000)
cincomilspecor <- cor(cincomil, use="complete.obs", method = "spearman")
cincomilspedist<- dist(x=cincomilspecor)
cincomilclust=hclust(cincomilspedist, method = "complete")
pdf("~/uproov_wd/rnaseq/spearman_distance_5000.pdf", width=14, height = 14)
pheatmap(as.matrix(cincomilspecor),
         cluster_rows=cincomilclust,
         cluster_cols=cincomilclust,
         col=revcolors)
dev.off()



#EXAMPLE PLOT COUNTS FOR ONE GENE

# #TP53
# d <- plotCounts(dds, gene="ENSG00000141510", returnData = TRUE)
# library("ggplot2")
# ggplot(d, aes(x=condition, y=count)) + 
#     geom_point(position=position_jitter(w=0.1,h=0)) + 
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
# #scale_y_log10(breaks=c(25,100,400))

# ASCL2
# library("ggplot2")
# 
# df <- plotCounts(dds, gene="ASCL2", returnData = TRUE)
# 
# df$sample <- rownames(df)
# 
# df$gene <- "ASCL2"
# 
# #df <- rbind(d1, d2, d3)
# 
# #folr_subtype <- 
# 
# ggplot(df, aes(x=condition, y=count)) + geom_boxplot() + geom_jitter(width = 0.2) +  theme_bw() + 
#     theme(axis.text.x = element_text(angle=90)) + ggtitle("ASCL2")
# 
# # ggplot(df, aes(x=sample, y=count)) +
# #     geom_point(position=position_jitter(w=0.1,h=0)) +
# #     theme_bw() + facet_grid(. ~ condition)  +
# # theme(axis.text.x = element_text(angle = 90, hjust = 1),
# #       panel.grid.minor = element_blank())+
# #     ylab("Normalized Read Count") + xlab("Sample") 
# # 
# 
# 
# cnv <- read.table("~/uproov_wd/folr/ascl2/ascl2_cnv.txt", header = TRUE)
# library(reshape2)    
# cnv.melt <- melt(cnv, id.vars = c("SAMPLE"), variable.name = "gene", value.name = "cnv" )  
# colnames(cnv.melt)[1] <- "sample"
# df.org <- df
# 
# df.cnv <- merge(df.org, cnv.melt, by = c("sample", "gene"))
# df.cnv$cnv <- factor(df.cnv$cnv)
# 
# ggplot(df.cnv, aes(x=cnv, y=count)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(col= condition), width = 0.2, size = 3) +  
#     facet_grid(. ~ gene) + theme_classic() + 
#     theme(axis.text.x = element_text(angle=0))
# 
