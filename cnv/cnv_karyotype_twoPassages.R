
#####USE AS: 
#Rscript cnv_karyotype_twoPassages.R $BLOOD.bam_ratio.txt $ORG.bam_ratio.txt $ORG2.bam_ratio.txt "TITLE" $OUTFILE




#' read.data
#'
#' @description Reads the different files
#'
#' @param tumor file with copy number ratios of tumor
#' @param ref file with copy number ratios of reference
#' @param org file with copy number ratios of organoid
#' @keywords internal
#'
#'
read.data <- function(ref, tumor){

    r <- data.frame(read.table(ref, sep='\t', header=TRUE))
    r$type <-  "Reference"
    t <- data.frame(read.table(tumor, sep='\t', header=TRUE))
    t$type <- "Tumor"
    # o <- data.frame(read.table(org, sep='\t', header=TRUE))
    # o$type <- "Organoid"
    merged.data <- rbind(r, t)
return (merged.data)
}

#' get.medians
#' @keywords internal
get.medians <- function(df, type, nrChromsToplot, binsize, ploidy, chrom, i){
    grouping <- cut(df$Start, c(seq(0,max(df$Start),binsize), max(df$Start)), labels=seq(binsize/2, (max(df$Start)+(binsize/2)), binsize))

    medians <- data.frame(tapply(df$Ratio, grouping, median))

    colnames(medians) <- "Ratio"
    # Alternate chromosome banding
    if ((i%%2) == 0) {
        medians$Cols <- 'neutral1'
    } else {
        medians$Cols <- 'neutral2'
    }
    medians <- medians[!is.na(medians$Ratio),]


    medians$Chromosome <- chrom
    medians$GenomicPosition <- as.numeric(as.character(rownames(medians)))
    medians$CopyNumber <- medians$Ratio*ploidy

    medians$Cols[medians$CopyNumber>=ploidy+0.8] <- "gain"
    medians$Cols[medians$CopyNumber<=ploidy-0.8] <- "loss"
    medians$type <- type
    return(medians)
}

#' transform.data
#' @keywords internal
transform.data <- function(merged.data, nrChromsToplot, binsize, ploidy, maxLevelToPlot){
    merged.data <- subset(merged.data, Ratio != -1)
    merged.data <- merged.data[!is.na(merged.data$Ratio),]
    chromosomes <- mixedsort(as.character(unique((merged.data$Chromosome))))
    plotting.df <- data.frame()

    for (i in c(1:nrChromsToplot)) {
        chrom <- chromosomes[i]

        #Ref
        ref.medians <- get.medians(df=subset(merged.data, Chromosome==chrom & type=="Reference"), nrChromsToplot = 23, type = "Reference",
                                   binsize = binsize, ploidy = ploidy, chrom = chrom, i = i)
        #Tumor
        t.medians <- get.medians(df=subset(merged.data, Chromosome==chrom & type=="Tumor"), nrChromsToplot = 23, type = "Tumor",
                                 binsize = binsize, ploidy = ploidy, chrom = chrom, i = i)
        #Organoid
        # o.medians <- get.medians(df=subset(merged.data, Chromosome==chrom & type=="Organoid"), nrChromsToplot = 23, type = "Organoid",
        #                          binsize = binsize, ploidy = ploidy, chrom = chrom, i = i)
        # Merge data
        #print(c(chrom, nrow(ref.medians), nrow(t.medians), nrow(o.medians), nrow(plotting.df)))
        plotting.df <- rbind(plotting.df, ref.medians, t.medians)

    }
    plotting.df$CopyNumber[plotting.df$CopyNumber > maxLevelToPlot] <- maxLevelToPlot
    plotting.df$Chromosome <- factor(plotting.df$Chromosome, levels=chromosomes)
    return(plotting.df)
}


#' plot.data
#' @keywords internal
plot.data <- function(plotting.df, maxLevelToPlot, title){
    p <- ggplot(plotting.df,
                aes(GenomicPosition, CopyNumber, group=Chromosome)) +
        geom_point(aes(colour=type), size = 0.8)  +
        ylim(0,maxLevelToPlot) +
        facet_wrap(~Chromosome, nrow=1, scales="free_x") + scale_colour_manual(limits = c("Reference", "Tumor"), values=c("orchid2", "#418042"))
    # + scale_shape_manual(values=c(3, 17, 19))
    # + scale_colour_manual(name="event", values=myCols)
    p <- p + theme(
        axis.line.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),

        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="grey80", linetype=2, size=0.5),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        plot.background=element_blank()) +
        ggtitle(title)
    return(p)
}


#' run.karyotipe
#' Plots a karyotype from the ratio files of the reference, the tumor, and the organoid
#' @export

run.karyotype <- function(ref, tumor, nrChromsToplot = 24, binsize = 1000000, ploidy = 2, maxLevelToPlot=4, title = "Title"){

    merged.data <- read.data(ref = ref, tumor = tumor)
    plotting.df <- transform.data(merged.data = merged.data, nrChromsToplot = nrChromsToplot, binsize = binsize,
                                  ploidy = ploidy, maxLevelToPlot =  maxLevelToPlot)
    plot.object <- plot.data(plotting.df, maxLevelToPlot, title)
    plot.data(plotting.df, maxLevelToPlot, title)
    return(plot.object)
}

library(gtools)
library(ggplot2)
library(RColorBrewer)

args <- commandArgs()
referenceFile <- args[6]
tumorFile <- args[7]
title <- args[8]
outFile <- args[9]
pl <- run.karyotype(ref = referenceFile, tumor = tumorFile, nrChromsToplot = 23, title = title)

pdf(file=outFile, width=14, height=4, pointsize=6, useDingbats=FALSE)
print(pl)
dev.off()
