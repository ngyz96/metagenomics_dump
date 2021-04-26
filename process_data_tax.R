#read in data, rows as variables, columns as observations
df <- read.csv("TAX.tsv", sep="\t") #choose from EC, EGGNOG, GTDB, INTERPRO2GO, SEED or TAX
rownames(df) <- df$NUM
df <- df[,-c(1,2)]
df <- as.data.frame(t(df))

#reading in 4M reads
bulk_4m <- read.csv("bulk4m_TAX.tsv", sep="\t")
root_4m <- read.csv("root4m_TAX.tsv", sep="\t")
num_col <- dim(bulk_4m)[2]
#sum up the counts
for (i in seq(4, num_col,2)) {
  bulk_4m[,i] <- bulk_4m[,i-1] + bulk_4m[,i]
  root_4m[,i] <- root_4m[,i-1] + root_4m[,i]
}
bulk_4m <- bulk_4m[,c(1,2,seq(4, num_col,2))]
root_4m <- root_4m[,c(1,2,seq(4, num_col,2))]
rownames(bulk_4m) <- bulk_4m$NUM
bulk_4m <- bulk_4m[,-c(1,2)]
bulk_4m <- as.data.frame(t(bulk_4m))
rownames(root_4m) <- root_4m$NUM
root_4m <- root_4m[,-c(1,2)]
root_4m <- as.data.frame(t(root_4m))
#some magic to combine them (i hope) Spoilers: it DID!
df$site <- rownames(df)
bulk_4m$site <- rownames(bulk_4m)
root_4m$site <- rownames(root_4m)
combined <- rbind(reshape2::melt(df,id="site"),
                  reshape2::melt(bulk_4m,id="site"),reshape2::melt(root_4m,id="site"))
combined <- reshape2::dcast(data=combined,formula=site ~ variable,fun.aggregate=sum)
rownames(combined) <- combined$site
df <- combined[, !(colnames(combined) %in% "site")]

# Some code to plot volin plots of count distribution
# library(tidyr)
# library(magrittr)
# library(ggplot2)
# test <- as.data.frame(df) 
# test$type <- c(rep("B", 15), rep("R", 15))
# test <- test %>% tibble::rownames_to_column(var="ID") %>% pivot_longer(c(-ID,-type))
# test$value <- sapply(test$value, "+", 1e-9)
# ggplot(test, aes(x=ID, y=value, fill = type)) + geom_violin() + 
#   scale_y_continuous(trans='log10') #+ geom_jitter()

# Read in NCBI taxid to scientific name mapping
taxid_map <- read.table("taxid_map.tsv", sep="\t", quote="")
# generated from NCBI names.dmp file with the following cmd line
# grep "scientific name" names.dmp | cut -d "|" -f1,2 | sed 's/\t|//;s/[\t]*$//;s/ /_/g' > taxid_map.tsv
colnames(taxid_map) <- c("taxid", "sci_name")
length(unique(taxid_map$taxid)) == length(unique(taxid_map$sci_name)) # looks like there are repeated names...
taxid_map$sci_name <- make.names(taxid_map$sci_name, unique = T, allow_ = T) #add .# behind duplicated value, where #=count
length(unique(taxid_map$taxid)) == length(unique(taxid_map$sci_name)) #check that it worked!

#collapse to a specific taxonomic level (PLEASE DONT RUN)
library(taxonomizr)
taxids <- as.numeric(colnames(df))
tmp <- sapply(taxids, getTaxonomy,'accessionTaxa.sql')
tmp <- t(tmp) #make same shape as taxid_map
rownames(tmp) <- taxids
colnames(tmp) <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
tmp <- gsub(" ", "_", tmp) #make same convention as taxid_map
saveRDS(tmp, "collapse_taxid_map.rds")
collapse_taxid_map <- readRDS("collapse_taxid_map.rds")
all(rownames(collapse_taxid_map) == colnames(df))
sci_names <- character()
headers <- colnames(df)
for (i in seq(1:ncol(df))) {
  id <- headers[i]
  curr_vec <- collapse_taxid_map[id,]
  if (!is.na(curr_vec[5])) {      #change to desired levels, phylum == 2 here
    sci_names[i] <- curr_vec[5]
  } else {
    curr_vec<-curr_vec[!is.na(curr_vec)]
    sci_names[i] <- ifelse(length(curr_vec) < 0, "root", "missing")
                           #taxid_map$sci_name[which(taxid_map$taxid == as.numeric(id))])
  }
}
colnames(df) <- sci_names

#collapse columns of same names
uniq_names <- unique(colnames(df))
df_temp <- data.frame(matrix(nrow=nrow(df), ncol=length(uniq_names)))
rownames(df_temp) <- rownames(df)
for (i in seq(1:length(uniq_names))){
  name <- uniq_names[i]
  df_temp[,i] <- rowSums(df[,which(colnames(df)==name), drop = FALSE])
  colnames(df_temp)[i] <- name
}
df <- df_temp[,which(!(colnames(df_temp) == "missing"))] #to remove parent taxonomic assignment
#df <- df_temp # leave taxonomic assignment in

# #map colnames of data to scientific names (old)
# colnames(df) <- taxid_map$sci_name[match(colnames(df), taxid_map$taxid)]

#get some metadata from sample names
names <- rownames(df)
names <- strsplit(names, split="_")
species <- as.factor(sapply(names, "[[", 2))
type <- sapply(names,"[[",1)
type <- as.factor(substr(type, 1,1))
index <- sapply(names,"[[",1)
index <- as.numeric(gsub("^\\D", "", index))
batch <- as.factor(
    ifelse(index==1|index==2|index==3,
           "1",
           ifelse(index==4|index==5|index==6,
                  "2",
                  ifelse(index==8|index==9,
                         "3",
                         ifelse(index==10|index==11|index==12,
                                "4",
                                ifelse(index==7|index==13|index==14|index==15,
                                       "5", "NA")))))
)
family <- ifelse(substring(species,1,1) == "S", "Dipterocarpaceae", "Fabaceae")
metadata <- data.frame(species,family,type,batch,index)
rownames(metadata) <- rownames(df)
rownames(df) == paste0(metadata$type,metadata$index,"_",metadata$species) #check if data corresponds

#save files 
save(df,metadata, file="5mdata.RData")

#Start from here 
load(file="5mdata.RData")

#check percent of zeroes (sparse matrix)
per_zeroes <- sum(df==0) / (dim(df)[1]*dim(df)[2])
per_zeroes
hist(log10(1 + t(df)), 100)
barplot(rowSums(df))
library(RColorBrewer)
colours <- sample(colorRampPalette(brewer.pal(8, "Set3"))(dim(df)[2]))
barplot(as.matrix(t(df/rowSums(df))), col=colours)
#ggplot ver 
library(ggplot2)
tmp <- df/rowSums(df)
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
svg("diamond_prop.svg", width=16, height=8)
ggplot(tmp, aes(x=type,y=value,fill=variable)) + geom_bar(position="fill", stat="identity") + facet_grid(~sample) +
    scale_fill_manual(values=as.vector(pals::cols25(24))) + labs(y= "Proportion", x = "Sample") +
    guides(fill=guide_legend(ncol=1, title = "Phylum")) + theme_minimal() + 
    theme(strip.background = element_rect(colour="grey80", fill="grey80"))
dev.off()

#compositional normalisation
library(zCompositions)
library(ALDEx2)
library(CoDaSeq)
df <- df[,colSums(df==0) < 25]
#alde.clr <- aldex.clr(t(df), conds=metadata$type, mc.samples = 128, denom="all", verbose=FALSE, useMC=FALSE)
df.nozero <- cmultRepl(t(df), label=0, method="CZM", output="p-counts")
df.nozero.clr <- codaSeq.clr(df.nozero, samples.by.row=F)
#another clr
library(compositions)
df.clr <- acomp(df, warn.na = T)

pdf(file="test.pdf", width=50, height=100)
# plot(df.clr)
# barplot(df.clr)
# boxplot(df.clr)
# boxplot(rcomp(df.clr))
pie(df.clr)
barplot(mean(df.clr))
dev.off()
#EDA with PCA
library(PCAtools)
#res_pca <- prcomp(df, center = T, scale = T) #scale = F for sparse matrix (old ver)
res_pca <- pca(df.nozero.clr, center = T, scale = F, transposed = F, metadata = metadata)
screeplot(res_pca) #only N- PCs due to centering
biplot(res_pca, showLoadings = T, colby = "species", shape = "type", ntopLoadings =10)
pairsplot(res_pca, components = getComponents(res_pca, seq_len(5)), colby = "family", shape = "type", 
          pointSize = 2, plotaxes = F, lab=rownames(df),labSize = 1.4)
plotloadings(res_pca, labSize = 2, shapeSizeRange=c(5,5))
# basic_plot <- fviz_pca_ind(res_pca)
# ggplot(cbind(basic_plot$data,metadata),aes(x=x,y=y,col=species,shape=type, label=name)) + 
#     geom_point(size=3)+theme_light() + geom_label_repel()

#multivariate analysis
library(vegan)
library(ggord)
library(ggvegan)
#ord.dca <- decorana(t(df.nozero.clr))
#summary(ord.dca)
ord.dbrda <- dbrda(t(df.nozero.clr)~dist_mat, distance="euclidean")
ord.rda <- rda(t(df.nozero.clr) ~ metadata$species + metadata$type + Condition(metadata$batch))
ord.rda2 <- rda(t(df.nozero.clr) ~ metadata$family + metadata$type) 
ord.rda3 <- rda(t(df.nozero.clr) ~ metadata$family + metadata$species + metadata$type) 
#ggord(ord.rda)
allscores <- scores(ord.rda, display = "sites")
allscores <- cbind(allscores, metadata)
centroids <- as.data.frame(scores(ord.rda, display="cn"))
speciesscore <- as.data.frame(scores(ord.rda, display = "species"))
speciesscore$dist <- sqrt(speciesscore$RDA1^2 + speciesscore$RDA2^2)
ggplot() + 
    geom_point(aes(x=RDA1,y=RDA2,shape=type,col=species), data=allscores, size=3) + scale_shape_manual(values=c(16,17))+geom_label_repel(aes(x=RDA1,y=RDA2,label=rownames(allscores)),data = allscores)+
    geom_point(aes(RDA1,RDA2), data=centroids, size=3, shape=23, fill="grey25", col="grey39", alpha = 0.7) + 
    geom_text_repel(aes(RDA1,RDA2, label=row.names(centroids)), data=centroids) +
    #geom_point(aes(RDA1,RDA2), data=speciesscore,size=2, shape=3, alpha=0.7,col="red") +
    #geom_text(aes(RDA1,RDA2, label=row.names(speciesscore)), data=speciesscore, size = 2.3, alpha = 0.6) +
      theme_light()
#perm test
anova(ord.rda, by = "margin", permutations = how(nperm=100))
anova(ord.rda2, ord.rda3, permutations = how(nperm=10000))
anova(ord.rda3, by="term", permutations = how(nperm=100000))


###INTERPRO2GO SECTION
#interpro2go mapping
interpro2go_mapping <- read.table("INTERPRO2GO_MAP_CLEANED.tsv", 
                                  col.names = c("IPR", "IPRTERM", "GO", "GOTERM"),
                                  header = F, sep = "\t",  quote="")
#read in data, rows as variables, columns as observations
df <- read.csv("INTERPRO2GO.tsv", sep="\t") #choose from EC, EGGNOG, GTDB, INTERPRO2GO, SEED or TAX
rownames(df) <- df$NUM
df <- df[,-c(1,2)]
df <- as.data.frame(t(df))

#reading in 4M reads
bulk_4m <- read.csv("bulk4m_INTERPRO2GO.tsv", sep="\t")
root_4m <- read.csv("root4m_INTERPRO2GO.tsv", sep="\t")
num_col <- dim(bulk_4m)[2]
#sum up the counts
for (i in seq(4, num_col,2)) {
    bulk_4m[,i] <- bulk_4m[,i-1] + bulk_4m[,i]
    root_4m[,i] <- root_4m[,i-1] + root_4m[,i]
}
bulk_4m <- bulk_4m[,c(1,2,seq(4, num_col,2))]
root_4m <- root_4m[,c(1,2,seq(4, num_col,2))]
rownames(bulk_4m) <- bulk_4m$NUM
bulk_4m <- bulk_4m[,-c(1,2)]
bulk_4m <- as.data.frame(t(bulk_4m))
rownames(root_4m) <- root_4m$NUM
root_4m <- root_4m[,-c(1,2)]
root_4m <- as.data.frame(t(root_4m))
#some magic to combine them (i hope) Spoilers: it DID!
df$site <- rownames(df)
bulk_4m$site <- rownames(bulk_4m)
root_4m$site <- rownames(root_4m)
combined <- rbind(reshape2::melt(df,id="site"),
                  reshape2::melt(bulk_4m,id="site"),reshape2::melt(root_4m,id="site"))
combined <- reshape2::dcast(data=combined,formula=site ~ variable,fun.aggregate=sum)
rownames(combined) <- combined$site
df <- combined[, !(colnames(combined) %in% "site")]

#save interpro
save(df, metadata, interpro2go_mapping, file="5M_interpro.RData")

#start here
load("5M_interpro.RData")

dir.create("goatools_data")
df <- t(df)
### for loop since idk how else to convert (NOT DIS)
for (i in 1:ncol(df)){
    tmp <- df[,i]
    tmp <- tmp[complete.cases(ifelse(tmp > 0, 1, NA))]
    write.table(as.numeric(names(tmp)),file=paste0("goatools_data/",colnames(df)[i],".txt"), 
                row.names=F, col.names=F, sep = "")
}

#get background
tmp <- rowSums(df)
tmp <- tmp[complete.cases(ifelse(tmp > 0, tmp, NA))]
write.table(as.numeric(names(tmp)),file="goatools_data/background.txt", 
            row.names=F, col.names=F, sep = "")

# goatools for top differentiating IPS 
sd_rda1 <- sd(speciesscore$RDA1)
mean_rda1 <- mean(speciesscore$RDA1)
ips_rda1_pos <- rownames(speciesscore[which(speciesscore$RDA1 > (mean_rda1+1*sd_rda1)),])
ips_rda1_neg <- rownames(speciesscore[which(speciesscore$RDA1 < -(mean_rda1+1*sd_rda1)),])
sd_rda2 <- sd(speciesscore$RDA2)
mean_rda2 <- mean(speciesscore$RDA2)
ips_rda2_pos <- rownames(speciesscore[which(speciesscore$RDA2 > (mean_rda2+1*sd_rda2)),])
ips_rda2_neg <- rownames(speciesscore[which(speciesscore$RDA2 < -(mean_rda2+1*sd_rda2)),])
write.table(as.numeric(ips_rda1_pos), file="goatools_data/ips_rda1_pos.txt", row.names=F, col.names=F, sep="")
write.table(as.numeric(ips_rda1_neg), file="goatools_data/ips_rda1_neg.txt", row.names=F, col.names=F, sep="")
write.table(as.numeric(ips_rda2_pos), file="goatools_data/ips_rda2_pos.txt", row.names=F, col.names=F, sep="")
write.table(as.numeric(ips_rda2_neg), file="goatools_data/ips_rda2_neg.txt", row.names=F, col.names=F, sep="")

#get euclidean distance
load("bukit_timah_big_trees.stem1.rdata")
metadata <- read.csv("metadata.tsv", sep="\t")
#somehow i lost dis will get back at some later time
bukit_timah_big_trees.stem1$tag <- as.numeric(bukit_timah_big_trees.stem1$tag)
trees <- bukit_timah_big_trees.stem1[bukit_timah_big_trees.stem1$tag %in% metadata$TAG_ID,]
trees <- trees[,c("tag", "gx", "gy")]
samples <- metadata$INDEX[match(trees$tag, metadata$TAG_ID)]
trees <- rbind(trees,trees)
rownames(trees) <- c(paste0("B", samples),paste0("R", samples))
trees <- trees[order(row.names(trees)), ]
dist_mat <- as.matrix(dist(trees[,c(2,3)]))
heatmap(dist_mat)


#ANCOM
library(ANCOMBC)
library(phyloseq)
# #create phyloseq object 
# tax_mat <- otu_table(df, taxa_are_rows = F) #pre collapse df
# tax_tab <- tax_table(collapse_taxid_map)
# samp_data <- sample_data(metadata)
# ps <- phyloseq(tax_mat, tax_tab, samp_data)
# saveRDS(ps, file="5mdata_phyloseq.rds")
# phylum phyloseq
tax_mat <- otu_table(df, taxa_are_rows = F)
samp_data <- sample_data(metadata)
ps <- phyloseq(tax_mat, samp_data)
saveRDS(ps, file="5mdata_phylum_phyloseq_ips.rds")
# load the phyloseq object
ps <- readRDS("5mdata_phylum_phyloseq.rds")
ancom_out <- ancombc(ps, formula="type*species", conserve=T, group="species", global = T)
library(DESeq2)
ds2 <- phyloseq_to_deseq2(ps, design = ~ batch + species + type)
res <- DESeq(ds2)
EnhancedVolcano(deseq2_res, lab = rownames(deseq2_res), x = 'log2FoldChange',y = 'pvalue')
ancom_out$res$diff_abn
source("ANCOM-2.1/scripts/ancom_v2.1.R")
res <- ANCOM(df, metadata, struc_zero = NULL, main_var="type", p_adj_method = "BH", 
alpha = 0.05)

