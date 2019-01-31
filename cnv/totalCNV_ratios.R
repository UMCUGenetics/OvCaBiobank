library(reshape2)
library(ggplot2)
library(gtools)
library(ggrepel)
library(grid)
library(gridExtra)


#Get proper sample names from the files
nameFromPath <- function(path)
{
    name.list <- unlist(strsplit(path, '/'))
    name.tmp <- name.list[length(name.list)]
    name.tmp <- unlist(strsplit(name.tmp, '\\.'))[1]
    name <- unlist(strsplit(name.tmp, '_dedup'))[1]
    splitname <- paste(unlist(strsplit(name, 'T|O'))[2], unlist(strsplit(name, '[0-9]'))[1])
    splitname
}

#DIRECTORY WITH RATIO FILES OUT OF FREEC
CNV_DIR="/home/cog/jvalleinclan/uproov_wd/cnv/revision/ratios_final/"
setwd(CNV_DIR)
#Files
cnv_files <- list.files(path=CNV_DIR, pattern = "^(T|O).*_ratio\\.txt$", full.names = T) 
#Get sample names
cnv_names <- lapply(cnv_files, nameFromPath)



binsize=1000000
nrChromsToplot = 23
maxLevelToPlot = 4

get.medians <- function(df, nrChromsToplot, binsize, ploidy, chrom, i){
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
    
    return(medians)
}

read_ratio <- function(filename){
    table <- read.table(filename, header=TRUE, sep='\t', stringsAsFactors = FALSE)
    table <- subset(table, Ratio != -1)
    table <- table[!is.na(table$Ratio),]
    chromosomes <- mixedsort(as.character(unique((table$Chromosome))))
    sample <- tail(strsplit(head(strsplit(filename, ".bam")[[1]], n=1), "/")[[1]], n = 1)
    # print(sample)
    
    if (sample == "O13" || sample == "O13-S") {
        ploidy <- 3
    } else {
        ploidy <- 2
    }
    
    plotting.df <- data.frame()
    for (i in c(1:nrChromsToplot)) {
        chrom <- chromosomes[i]
        medians <- get.medians(df=subset(table, Chromosome==chrom), nrChromsToplot = 23,
                               binsize = binsize, ploidy = ploidy, chrom = chrom, i = i)
        plotting.df <- rbind(plotting.df, medians)

    }
    plotting.df$CopyNumber[plotting.df$CopyNumber > maxLevelToPlot] <- maxLevelToPlot
    plotting.df$Chromosome <- factor(plotting.df$Chromosome, levels=chromosomes)
    return(plotting.df)

}

#Read all files at once
cnv.list <- lapply(cnv_files, read_ratio)

names(cnv.list) <- cnv_names
#We only want one data frame with all the data


###################BASIC DF#######################
cnv.df <- melt(cnv.list, id.vars=c('Ratio','Cols', 'Chromosome', 'GenomicPosition', 'CopyNumber'))
cnv.df[is.na(cnv.df)] <- 0
colnames(cnv.df) <- c('Ratio','Cols', 'chr', 'GenomicPosition', 'CopyNumber', 'rawsample')
cnv.df$start <- cnv.df$GenomicPosition - binsize/2
cnv.df$end <- cnv.df$GenomicPosition + binsize/2

####################FILTERS#########################

#REGIONS FROM CIRCOS KARYOTYPE. THIS ARE CENTROMERES AND p-arm of acrocentric chromosomes
#Centromeres
exclude <- read.table("~/uproov_wd/cnv/regions_to_exclude.Rbased.cent.txt", header = T, sep='\t', stringsAsFactors = F)
exclude[is.na(exclude)] <- 0
row.names(exclude) <- exclude$chr
excludefunc <- function(x){
    c <- x["chr"]
    r <- (x["start"] <= exclude[c, "end"] & x["end"] >= exclude[c, "start"])
    r
}

exclude.vector <- apply(cnv.df, 1, excludefunc)

cnv.filtered <- cnv.df
cnv.filtered$CopyNumber[exclude.vector] <- 2
#cnv.filtered$type[exclude.vector] <- "normal"
cnv.filtered <- subset(cnv.filtered, cnv.filtered$chr!="Y")

#Stalk 
exclude.stalk <- read.table("~/uproov_wd/cnv/regions_to_exclude.Rbased_stalk.txt", header = T, sep='\t', stringsAsFactors = F)
exclude.stalk[is.na(exclude.stalk)] <- 0
row.names(exclude.stalk) <- exclude.stalk$chr
excludeStalkfunc <- function(x){
    c <- x["chr"]
    r <- (x["start"] <= exclude.stalk[c, "end"] & x["end"] >= exclude.stalk[c, "start"])
    r
}

exclude.stalk.vector <- apply(cnv.df, 1, excludeStalkfunc)

cnv.filtered$CopyNumber[exclude.talk.vector] <- 2

###############################PLOT##########################

cnv.plot <- cnv.filtered

sampledict <- c("MC-1.1", "MC-1.2", "LGS-1.1", "LGS-1.2", "LGS-1.3", "LGS-1.4", "HGS-3.1", "HGS-3.1.S", "HGS-3.2", "HGS-3.2.S", "MBT-1", 
                "LGS-2.2", "HGS-6", "HGS-6.S", "MBT-2.1", "HGS-2", "HGS-2.S", "HGS-4.1", "HGS-4.3", "MBT-2.2", "HGS-4.5", 
                "HGS-1", "HGS-1.S", "HGS-1-R1", "HGS-1-R2", "HGS-1-R2.S", "HGS-1-R3",
                "HGS-19", "HGS-5.1", "MC-2.1", "MC-2.2", "SBT-4.1", "SBT-4.2",
                "SBT-3.1", "SBT-3.2", "LGS-3.1", "LGS-3.2", "HGS-13.3", "HGS-13.4", "HGS-13.4", "LGS-5.2", "LGS-5.4", "OSE-H", 
                "L-OC135","CCC-1", "END-1", "END-1", "END-1")
names(sampledict) <- c("1", "2", "4", "5", "7", "8", "9", "9-S", "10", "10-S","11", 
                       "12", "13", "13-S", "14", "15", "15-S", "16", "17", "18", "21", 
                       "25", "25-S", "37", "61", "61-S", "74", 
                       "27", "29", "40", "41", "52", "53", 
                       "55", "56", "58", "59", "65", "66", "66-B","68", "70", "-OSE10", 
                       "-L-OC135", "75", "76", "76-A", "76-B")


cnv.plot$label <- apply(cnv.plot, 1, function(x){
    if (grepl("-S", x["rawsample"])){
        sample_raw <- unlist(strsplit(x["rawsample"], ' '))[1]
        sample <- unlist(strsplit(x["rawsample"], '-'))[1]
        type <- paste(unlist(strsplit(x["rawsample"], ' '))[2],"-S", sep="")
        
        
    } else if(grepl("-A", x["rawsample"])) {
        sample_raw <- unlist(strsplit(x["rawsample"], ' '))[1]
        sample <- unlist(strsplit(x["rawsample"], '-'))[1]
        type <- paste(unlist(strsplit(x["rawsample"], ' '))[2],"-A", sep="")
        } else if(grepl("-B", x["rawsample"])) {
            sample_raw <- unlist(strsplit(x["rawsample"], ' '))[1]
            sample <- unlist(strsplit(x["rawsample"], '-'))[1]
            type <- paste(unlist(strsplit(x["rawsample"], ' '))[2],"-B", sep="")
        } else {
        sample <- unlist(strsplit(x["rawsample"], ' '))[1]
        type <- unlist(strsplit(x["rawsample"], ' '))[2]
    } 
    
    patient <- sampledict[sample]
    
    label <- paste(patient,type, sep='.')
    label
})

cnv.plot$patient <- apply(cnv.plot, 1, function(x){
    pat <- unlist(strsplit(unlist(x["label"])[1], '(\\.)'))[1]
    pat
})

cnv.plot <- cnv.plot[cnv.plot$patient != 'NA',]
cnv.plot$CopyNumber[cnv.plot$CopyNumber > 4] <- 4

chromosomes <- mixedsort(as.character(unique((cnv.plot$chr))))
cnv.plot$chr <- factor(cnv.plot$chr, levels=chromosomes)
cnv.plot$patient <- factor(cnv.plot$patient, 
                           levels=c("HGS-1", "HGS-1-R1", "HGS-1-R2", "HGS-1-R3",
                                    "HGS-2", "HGS-3", "HGS-4", "HGS-5", "HGS-6",
                                    "HGS-13","HGS-19",
                                    "LGS-1", "LGS-2","LGS-3","LGS-5",
                                    "CCC-1", "END-1","MC-1", "MC-2","MBT-1", "MBT-2",
                                    "SBT-3","SBT-4"))
cnv.plot[cnv.plot$label=="HGS-13.4.T-B", "label"] <- "HGS-13.4.T"


#Extract Legend 
g_legend<-function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    return(legend)} 





p <- ggplot(cnv.plot)

cnv_type <- p + theme_bw() + 
    geom_segment(aes(x=start, xend=end, y=label, yend=label, color=CopyNumber), size = 3) +
    #scale_color_distiller(palette=rev("RdBu")) +
    scale_color_gradient2(low = "midnightblue",
                          midpoint = 2,
                          mid = "white",
                          high = "orangered3",
                          space="Lab",
                          breaks=c(0,1,2,3,4),
                          labels=c(0,1,2,3,"4 or more")) +
    facet_grid(patient ~ chr, scales = "free", space = "free_y" )+#, space = "free") + 
    theme(panel.spacing=unit(.01, "lines"),
          panel.border = element_rect(color = "black", fill = NA, size = .2), 
          axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y =element_text(size=8, hjust=0,vjust=0.5, face = "bold"),
          strip.text = element_text(size = 8, angle= 0, hjust = 0.5, vjust=0.5, face="bold"),
          strip.text.y = element_text(size = 8, angle= 0, hjust = 0.5, vjust=0.5, face="bold"),
          strip.placement = "outside",
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_line(size=1, linetype="dashed"),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          legend.position="top",
          legend.title = element_text(face="bold"), 
          panel.background=element_blank(),
          plot.background=element_blank(),
          legend.key = element_rect(fill = "white"),
          strip.background = element_rect(color = "black", fill = NA, size = .5)
    ) +
    theme(
        legend.key = element_rect(fill = "white"),
        legend.position = "top",
        legend.box = "vertical",
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        legend.direction='vertical',
        legend.box.just = "left",
        legend.spacing.y = unit(-0.55, "cm"),
        legend.background = element_rect(fill="transparent")
    ) +
    labs(title="",x = "", y="")

cnv_type
cnv_type_legend <- g_legend(cnv_type) 
cnv_type

pdf("~/uproov_wd/PLOTS/revision/cnv_ratios.pdf", width=16, height=10, useDingbats=FALSE, pointsize=6)
cnv_type
dev.off()
pdf("~/uproov_wd/PLOTS/revision/cnv_ratios_legend.pdf", pointsize=10)
cnv_type_legend
dev.off()
