library(tidyverse)
library(RColorBrewer)
#plot stuff for qc
reads_num <- read.csv("data/track_reads.csv")

reads_num <- reads_num %>% mutate(Human=Host-Human,Host=Fastp-Host,Fastp=Raw-Fastp) 
#reads_num[,4:7] <- reads_num[,4:7] / rowSums(reads_num[,4:7])
reads_num <- reads_num %>% pivot_longer(cols=Raw:Human)
reads_num$sample <- factor(reads_num$sample)
reads_num$type <- factor(reads_num$type)
reads_num$species <- factor(reads_num$species)
colour <- brewer.pal(4,"Set1")[c(1,3,4,2)]
#Figure 2
svg("qc_plot.svg", width=14, height=8)
ggplot(reads_num, aes(fill=name,y=value,x=sample)) + facet_wrap(type~species, scales = "free_x", drop=T) +
    geom_bar(position="stack", stat="identity") + theme_light() + coord_cartesian(ylim = c(4e7, 6.7e7)) +
    theme(strip.background =element_rect(fill="grey80"),strip.text = element_text(colour = 'black')) + 
    scale_fill_manual(values=colour, name="Reads Processing", labels=c("Low Quality","Host contamination","Human contamination","Pass QC")) + 
    ylab("Number of reads") + xlab("Sample Number") + 
    theme(legend.text=element_text(size=12), legend.title=element_text(size=14), axis.text=element_text(size=13), axis.title=element_text(size=16,face="bold"))
dev.off()

#kraken
kraken_num <- read.csv("kraken2_initial.csv", stringsAsFactors = T)
kraken_num$sample <- factor(kraken_num$sample)
#svg("kraken_mapping.svg", width=14, height=7)
p1 <- ggplot(kraken_num, aes(y=percent_mapped,x=sample)) + facet_wrap(type~species, scales="free_x") +
    geom_bar(fill=as.vector(pals::cols25(1)), stat="identity") + theme_light() + coord_cartesian(ylim = c(0,100)) +
    theme(strip.background =element_rect(fill="grey80"),strip.text = element_text(colour = 'black'), panel.grid.major.x = element_blank()) + 
    #scale_fill_manual(values=as.vector(pals::cols25(1)), name="Reads", labels ="Mapped") + 
    ylab("Percentage Mapped") + xlab("Sample Number")

#dev.off()

#bracken_plotting
df <- read.csv("data/bracken_combined.tsv", sep="\t", row.names = "taxa")
rownames(df) <- gsub(" ", "_", rownames(df))
df <- t(df)
df <- df[,colSums(df==0) < 25]
tmp <- as.data.frame(df/rowSums(df))
means <- colMeans(tmp)
keep <- names(means[!means < 0.001])
tmp <- tmp[,which(colnames(tmp) %in% keep)]
tmp <- tmp[,order(apply(tmp, 2, median, na.rm = T), decreasing = T)]
tmp$Low_counts <- 1-rowSums(tmp)
tmp$type <- substr(rownames(tmp), 1,1)
tmp$sample <- factor(sub(".", "", rownames(tmp)), 
                     levels= c("1_SC","2_SL","3_KM","4_SL","5_SC","6_SC","7_SC",
                               "8_KM","9_KM","10_SL","11_KM","12_SL","13_KM","14_SL","15_SC"))

tmp <- reshape2::melt(tmp)
#svg("bracken_prop.svg", width=14, height=7)
p2 <- ggplot(tmp, aes(x=type,y=value,fill=variable)) + geom_bar(position="fill", stat="identity") + facet_grid(~sample) +
    scale_fill_manual(values=as.vector(pals::cols25(15))) + 
    labs(y= "Proportion", x = "Sample") +
    guides(fill=guide_legend(ncol=1, title = "Phylum")) + theme_minimal() + 
    theme(strip.background = element_rect(colour="grey80", fill="grey80"))
#dev.off()

#join figs
library(cowplot)
#p1 = kraken, p2 = bracken
#Figure 3
svg("kraken_combined.svg", width=12, height=15)
plot_grid(p1,p2,ncol=1, align = "v", axis = "lr", labels = "AUTO", label_size = 18, rel_heights = c(1,1.7))
dev.off()
