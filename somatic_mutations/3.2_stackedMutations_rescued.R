
library(gtools)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

directory <- "."


filenames <- list.files(directory, recursive = FALSE, pattern = '*.shares.af.0.05.txt',full.names = TRUE)
df.list <- lapply(filenames, read.table, sep = '\t', header = TRUE, row.names = 1, stringsAsFactors=FALSE)



#Rescue we convert the 2s into 1s
df.rescue.list <- lapply(df.list, function(x) {
    tmp <- x
    
    #!!!!THIS ASSUMES THAT ORGANOID COMES FIRST AND TUMOR SECOND CAREFUL
    #Try to change with regexp or something 
    names(tmp) <- c("org", "tumor")
    tmp$type <- apply(tmp, 1, function(y){
        if(y["org"] == 1 && y["tumor"] == 1){
            type <- "shared"
        } else if (y["org"] == 1 && y["tumor"] == 0) {
            type <- "org"
        } else if (y["org"] == 0 && y["tumor"] == 1) {
            type <- "tumor"
        } else if (y["org"] == 2 && y["tumor"] == 1) {
            type <- "orescued"
        } else if (y["org"] == 1 && y["tumor"] == 2) {
            type <- "trescued"
        } else {
            ###YOU CANNOT HAVE A 0 in one and a 2 in the other, can you?
            print ("ERROR")
            print (c(y["org"],y["tumor"]))
        }
    })
    names(tmp) <- c(names(x), "type")
    tmp
})

df.rescue.percents <- lapply(df.rescue.list, function(x) {
    tmp <- table(x["type"])
    tmp$sample <- unlist(strsplit(colnames(x)[1], '(T|O)'))[2]
    #print(names(tmp))
    if (!("orescued"  %in% names(tmp))){
        tmp$orescued <- 0
    }
    if (!("trescued"  %in% names(tmp))){
        tmp$trescued <- 0
    }
    if (!("tumor"  %in% names(tmp))){
        tmp$tumor <- 0
    }
    if (!("org"  %in% names(tmp))){
        tmp$orescued <- 0
    }
    if (!("shared"  %in% names(tmp))){
        tmp$orescued <- 0
    }
    #Order elements to avoid problems when binding
    tmp <- tmp[c("orescued", "trescued", "tumor", "org", "shared", "sample")]
    tmp
})

#Merge the different matrices into one
df.rescue <- (as.data.frame(t(do.call("cbind", df.rescue.percents ))))
row.names(df.rescue) <- df.rescue$sample



sampledict <- c("MC-1.1", "MC-1.2", "LGS-1.1", "LGS-1.2", "LGS-1.3", "LGS-1.4", "HGS-3.1", "HGS-3.1.S", "HGS-3.2", "HGS-3.2.S", "MBT-1", 
                "LGS-2.2", "HGS-6", "HGS-6.S", "MBT-2.1", "HGS-2", "HGS-2.S", "HGS-4.1", "HGS-4.3", "MBT-2.2", "HGS-4.5", 
                "HGS-1", "HGS-1.S", "HGS-1-R1", "HGS-1-R2", "HGS-1-R2.S", "HGS-1-R3",
                "HGS-19", "HGS-5.1", "MC-2.1", "MC-2.2", "SBT-4.1", "SBT-4.2","HGS-15", 
                "SBT-3.1", "SBT-3.2", "LGS-3.1", "LGS-3.2", "HGS-13.3", "HGS-13.4", "LGS-5.2", "LGS-5.4", "OSE-H", "L-OC135", "HGS-10",
                "CCC-1", "END-1", "END-1.A", "END-1.B")
names(sampledict) <- c("") #Tumor & organoid numbers removed for anonimity



df.rescue$label <- apply(df.rescue, 1, function(x){
   
    #FOR SOME REASON SAMPLE IS A LIST
    sample.raw <- unlist(x["sample"])[1]
    sample <- sampledict[sample.raw]
    sample
    
})

df.rescue$patient <- apply(df.rescue, 1, function(x){
    pat <- unlist(strsplit(unlist(x["label"])[1], '(\\.)'))[1]
    pat
})


#After 5476362 tries to melt and modify the plotting df, gonna do it manually
df.rescue.plot.raw <- data.frame(sample=rep(df.rescue$label, each=2), 
                             patient=rep(df.rescue$patient, each=2),
                             type=rep(c("tumor", "organoid"), length(df.rescue$label)),
                             shared=rep(as.numeric(as.character(unlist(df.rescue$shared))), each=2),
                             unique=c(rbind(as.numeric(as.character(unlist(df.rescue$tumor))), unlist(df.rescue$org))),
                             #Trescued should be in the organoid and of course viceversa
                             #Therefore order inverted
                             rescued=c(rbind(as.numeric(as.character(unlist(df.rescue$orescued))),
                                              unlist(df.rescue$trescued)))
                             )
df.rescue.plot.raw$total <- apply(df.rescue.plot.raw, 1, function(x){
    as.numeric(x["shared"]) + as.numeric(x["unique"]) + as.numeric(x["rescued"])
})

df.rescue.plot.raw$plotSample <- apply(df.rescue.plot.raw, 1, function(x){
    if (x["type"] == "tumor"){
        tmp = "T"
    } else if (x["type"] == "organoid") {
        tmp = "O"
    }
    paste(x["sample"], tmp, sep = ".")
})

df.rescue.plot.raw <- df.rescue.plot.raw[!grepl("S\\.T", df.rescue.plot.raw$plotSample),]

df.shared.raw <- df.rescue.plot.raw
df.shared.raw$shared <- apply(df.rescue.plot.raw, 1, function(x){
    as.numeric(x["shared"]) + as.numeric(x["rescued"])
})

#df.rescue.plot.raw$shared <- NULL
df.shared.raw$rescued <- NULL


###SHARED
df.shared.plot <- melt(df.shared.raw, id.vars = c("sample", "patient", "total", "type", "plotSample"), 
                       variable.name = "typeMut")


df.shared.plot$typeMut <- factor(df.shared.plot$typeMut, levels=rev(c("shared", "unique")))

df.shared.plot$patient <- factor(df.shared.plot$patient, 
                                 levels=c("HGS-1", "HGS-1-R1", "HGS-1-R2", "HGS-1-R3",
                                          "HGS-2", "HGS-3", "HGS-4", "HGS-5", "HGS-6",
                                          "HGS-13","HGS-19",
                                          "LGS-1", "LGS-2","LGS-3","LGS-5",
                                          "CCC-1", "END-1","MC-1", "MC-2","MBT-1", "MBT-2",
                                          "SBT-3","SBT-4"))
df.shared.plot <- df.shared.plot[!is.na(df.shared.plot$patient),]

shared <- ggplot(df.shared.plot, aes(x = plotSample, y = value, fill=typeMut)) +
    geom_bar(position="stack", stat="identity", width=.9) + theme_bw()

shared.plot <- shared +
    facet_grid(patient ~ ., scales = "free", space = "free") +
    coord_flip() +
    theme(panel.grid = element_blank(),
          axis.line.x = element_line(colour = "black"),
          panel.spacing.y=unit(.3, "lines"),
          panel.border = element_rect(color = "black", fill = NA, size = .5),
          strip.text.y = element_text(size = 8, angle= 0, hjust = 0, vjust=0.5, face="bold"),
          axis.text.y =element_text(size=8, hjust=0,vjust=0.5, face = "bold"),
          axis.text.x =element_text(size=8, vjust=0.2, face="bold"),
          strip.background = element_blank()) +
    scale_fill_manual(name="Type", values = c("#e79f00", "#0072B2"), limits=c("shared", "unique"), labels = c("shared", "unique")) +
    ylab("# SNVs") + xlab("Sample") +
    scale_y_continuous(breaks = seq(0,24000,2000), labels = seq(0,24000, 2000), expand = c(0.01, 0))

shared.plot

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}

slegend <- g_legend(shared.plot)

shared.plot <- shared.plot + theme(legend.position = 'none')


pdf(file =  "~/uproov_wd/PLOTS/revision/sharedSNV_noBlood_new.pdf", width=10, height = 16)
shared.plot
dev.off()
ggsave(filename = "~/uproov_wd/PLOTS/revision/sharedSNV_legend.pdf", plot = slegend, device = pdf,
       width=8, height = 5)


