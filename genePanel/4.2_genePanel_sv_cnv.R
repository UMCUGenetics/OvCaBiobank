library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
filename = args[1]
outfile = args[2]

df <- read.table(filename, header = FALSE, sep = '\t', stringsAsFactors = FALSE)

colnames(df) <- c("patient", "type", "gene", "chr", "snv1", "snv1Gen", "snv2", "snv2Gen", "sv", "svGen", "cn", "cnType")


df.second <- df #legacy issue
df.second$cn[df.second$cn > 4] <- 4

df.second$realPatient <- apply(df.second, 1, function(x){
    label = x["patient"]
    p = unlist(strsplit(label, "\\."))[1]
    p
})

df.second$realLabel <- apply(df.second, 1, function(x){
    label = x["patient"]
    type = x["type"]
    realLabel = paste(label, type, sep = "-")
    realLabel
})

df.second$stroke <- (df.second$cnType != 'normal' | df.second$snv1 != '.')

df.second$sv <- factor(df.second$sv)

df.second$cn <- factor(df.second$cn)

df.second$stroke <- factor(df.second$stroke)

dfSize <- c(6,9,12)
names(dfSize) <- c('loss', 'normal', 'gain')
df.second$size <- dfSize[df.second$cnType]
df.second$size <- factor(df.second$size)



####WITH DIVIDE SQUARES

snv2plot <- c("chromosome_number_variation", 'exon_loss_variant', 'frameshift_variant', 'stop_gained', 'stop_lost', 'start_lost', 
             'splice_acceptor_variant', 'splice_donor_variant', 'rare_amino_acid_variant', 'missense_variant', 'inframe_insertion',
             'disruptive_inframe_insertion', 'inframe_deletion', 'disruptive_inframe_deletion', '5_prime_UTR_truncation+exon_loss_variant',
             '3_prime_UTR_truncation+exon_loss', 'splice_branch_variant', 'splice_region_variant', 'splice_branch_variant', 
             'initiator_codon_variant', 'initiator_codon_variant+non_canonical_start_codon', '.')
######CAREFUL IN CASE THERE ARE NOT ANY OTHER, I KNOW IN MY DATA THERE IS NOT ANYTHING ELSE I WANT TO PLOT###
snv2plot <- c('.', "missense_variant","frameshift_variant","splice_region_variant" ,
              "stop_gained", "inframe_insertion", "inframe_deletion")


mixedPalette <- c('white', '#e79f00', 'deeppink3',  'darkslategray', '#018571', '#0072B2',  "khaki3", "orangered3", "midnightblue")
names(mixedPalette) <- c('.', "missense_variant","frameshift_variant","splice_region_variant" ,
"stop_gained", "inframe_insertion","inframe_deletion", 'gain', 'loss')

df.second[df.second$cnType == "normal", "cnType"] <- '.'
df.second[!(df.second$snv1 %in% snv2plot), "snv1"] <- '.'


 df.second$realPatient <- factor(df.second$realPatient,
                             levels=c("HGS-1", "HGS-1-R1", "HGS-1-R2", "HGS-1-R3",
                                      "HGS-2", "HGS-3", "HGS-4", "HGS-5", "HGS-6",
                                      "HGS-13","HGS-19",
                                      "LGS-1", "LGS-2","LGS-3","LGS-5",
                                      "CCC-1", "END-1","MC-1", "MC-2","MBT-1", "MBT-2",
                                      "SBT-3","SBT-4"))
 df.filtered <- df.second[!is.na(df.second$realPatient),]
q <-  ggplot(df.filtered, aes(x=gene, y=realLabel))
p <- q +
    #SNV
    #geom_point(data = df.filtered[(df.filtered$snv1 %in% snv2plot),], aes(col=snv1), shape = "\u25E3", size = 7) +
    geom_point(data = df.filtered[df.filtered$cnType != '.' | df.filtered$snv1 != '.', ], shape = 0, size = 5, 
               col = "black",stroke = 1, show.legend = FALSE) + 
    geom_point(data = df.filtered, aes(col=snv1), shape = "\u25E3", size = 5,stroke = 1) +
    
    geom_point(data = df.filtered, aes(col=cnType), shape = "\u25E5", size = 5,stroke = 1)+
    
    scale_colour_manual(values = mixedPalette, limits = names(mixedPalette),  
                       labels=c("", "Missense", "Frameshift","Splice region","Stop gained",
                                "Inframe insertion", "Inframe Deletion", 'Gain', "Loss")) +
    
    #PRETTINESS
    facet_grid(realPatient ~ ., scales = "free", space = "free") +
    theme(panel.spacing.y=unit(.3, "lines"),
          panel.border = element_rect(color = "black", fill = NA, size = 1), 
          #strip.background = element_blank(),
          strip.background = element_rect(color = "black", size =.5, fill = NA),
          axis.line.x=element_blank(),
          axis.text.x=element_text(size = 8, angle = 90, vjust = .5, hjust = 1, face = "bold"),
          axis.text.y=element_text(size = 8, vjust = .5, hjust = 1, face = "bold"),
          axis.ticks.x=element_blank(),
          legend.key = element_rect(fill = "white"),
          legend.position = "none",
          legend.box = "vertical",
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.background=element_blank(),
          legend.title = element_text(face="bold"), 
          #panel.border=element_rect(size = .2, fill = "black"),
          plot.background=element_blank(),
          strip.text.y = element_text(size = 10, angle= 0, hjust = 0.5, vjust=0.5, face = "bold")
    ) 
p
outfile = "~/uproov_wd/genePanel_revision/noblood_filtered/genePanel_revision.svg"
ggsave(outfile, plot=p, device=svg, width=20, height=20)
