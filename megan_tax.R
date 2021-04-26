library(ggplot2)
library(cowplot)
library(zCompositions)
library(compositions)
library(PCAtools)
library(vegan)
library(ALDEx2)


#read in data
df <- read.csv("data/rma_TAX.tsv", sep="\t")
rownames(df) <- df$NUM
df <- df[,-c(1,2)]
df <- as.data.frame(t(df))
metadata <- read.csv("data/metadata.tsv", sep="\t")
metadata$BATCH <- as.factor(metadata$BATCH)
rownames(metadata) <- paste0(metadata$SAMPLE_ID, "_",metadata$SPECIES)

# get total counts for plotting
dia_num <- as.data.frame(rowSums(df)) / 100000
colnames(dia_num) <- "percent_mapped"
dia_num$species <- sapply(strsplit(rownames(dia_num), split="_"), "[[", 2)
dia_num$type <- ifelse(substr(sapply(strsplit(rownames(dia_num), split="_"), "[[", 1),1,1) == "B",
                       "Bulk", "Root")
dia_num$sample <- gsub("^.", "", sapply(strsplit(rownames(dia_num), split="_"), "[[", 1))
p1 <-ggplot(dia_num, aes(y=percent_mapped,x=sample)) + facet_wrap(type~species, scales="free_x") +
    geom_bar(fill=as.vector(pals::cols25(1)), stat="identity") + theme_light() + coord_cartesian(ylim = c(0,100)) +
    theme(strip.background =element_rect(fill="grey80"),strip.text = element_text(colour = 'black'), panel.grid.major.x = element_blank()) + 
    #scale_fill_manual(values=as.vector(pals::cols25(1)), name="Reads", labels ="Mapped") + 
    ylab("Percentage Mapped") + xlab("Sample Number")

#collapse to a specific taxonomic level (PLEASE DONT RUN UNLESS NEEDED)
library(taxonomizr)
taxids <- as.numeric(colnames(df))
tmp <- sapply(taxids, getTaxonomy,'accessionTaxa.sql')
tmp <- t(tmp) #make same shape as taxid_map
rownames(tmp) <- taxids
colnames(tmp) <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
tmp <- gsub(" ", "_", tmp) #make same convention as taxid_map
saveRDS(tmp, "read/collapse_taxid_map.rds")
collapse_taxid_map <- readRDS("data/collapse_taxid_map.rds")
all(rownames(collapse_taxid_map) == colnames(df))
sci_names <- character()
headers <- colnames(df)
for (i in seq(1:ncol(df))) {
    id <- headers[i]
    curr_vec <- collapse_taxid_map[id,]
    if (!is.na(curr_vec[2])) {      #change to desired levels, phylum = 2 here
        sci_names[i] <- curr_vec[2]
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
df <- df_temp[,which(!(colnames(df_temp) == "missing"))] # to remove parent taxonomic assignment
#df <- df_temp # leave taxonomic assignment in

#save files 
save(df,metadata,p1, file="data/5m_rma_data.RData")

#Start from here 
load(file="data/5m_rma_data.RData")

#check percent of zeroes (sparse matrix)
per_zeroes <- sum(df==0) / (dim(df)[1]*dim(df)[2])
per_zeroes
hist(log10(1 + t(df)), 100)
barplot(rowSums(df))

#diamond proportion plot
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
#svg("diamond_prop.svg", width=16, height=8)
p2 <- ggplot(tmp, aes(x=type,y=value,fill=variable)) + geom_bar(position="fill", stat="identity") + facet_grid(~sample) +
    scale_fill_manual(values=as.vector(pals::cols25(24))) + labs(y= "Proportion", x = "Sample") +
    guides(fill=guide_legend(ncol=1, title = "Phylum")) + theme_minimal() + 
    theme(strip.background = element_rect(colour="grey80", fill="grey80"))
#dev.off()
#figure 4
svg("diamond_combined.svg", width=12, height=15)
plot_grid(p1,p2,ncol=1, align = "v", axis = "lr", labels = "AUTO", label_size = 18, rel_heights = c(1,1.7))
dev.off()

#calculate alpha diversity first
Shannon <- diversity(x=df, index = "shannon")
Simpson <- diversity(x=df, index = "simpson")
richness <- specnumber(df, MARGIN = 1)
Pielou <- Shannon/log(richness)
indices <- as.data.frame(cbind(Shannon, Simpson, Pielou))
indices$type <- substr(rownames(indices), 1,1)
indices$species <- sapply(strsplit(rownames(indices), "_"), "[[", 2)
indices_plot <- tidyr::pivot_longer(indices, cols =c("Shannon", "Simpson", "Pielou"))
indices_plot$name <- factor(indices_plot$name, levels = c("Shannon", "Simpson", "Pielou"))

svg("a_diversity.svg", width=12, height=5)
alpha <- ggplot(indices_plot) + facet_wrap(~name, scale="free") + 
    geom_boxplot(aes(x=type, y=value, fill = species)) + 
    #scale_fill_manual(values=as.vector(pals::cols25(3))) + 
    theme_light() + theme(strip.background = element_rect( fill="grey80"),strip.text = element_text(colour = 'black'))
dev.off()

summary(aov(Shannon ~ type, indices))
summary(aov(Simpson ~ type, indices))
summary(aov(Pielou ~ type, indices))


#clr preprocessing
df <- df[,colSums(df==0) < 25]
df.nozero <- cmultRepl(t(df), label=0, method="CZM", output="p-counts")
df.clr <- clr(df.nozero)
df.clr <- t(as.data.frame(df.clr))

#pca
res_pca <- PCAtools::pca(df.diff, center = T, scale = T, transposed = T, metadata = metadata.root)
screeplot(res_pca) #only N- PCs due to centering
PCAtools::biplot(res_pca, showLoadings = T, colby = "BATCH", shape = "TYPE", ntopLoadings=3,xlim=c(-6,6),ylim=c(-6,6), ellipse = F) + 
    theme_light()
pairsplot(res_pca, components = getComponents(res_pca, seq_len(5)), colby = "SPECIES", shape = "TYPE", 
          pointSize = 2, plotaxes = F, lab=rownames(df),labSize = 1.4) 
plotloadings(res_pca, labSize = 2, shapeSizeRange=c(5,5))

#custom pca
xidx <- order(abs(res_pca$loadings[, "PC1"]), decreasing = TRUE)
yidx <- order(abs(res_pca$loadings[, "PC2"]), decreasing = TRUE)
vars <- unique(c(rownames(res_pca$loadings)[xidx][seq_len(4)], 
                 rownames(res_pca$loadings)[yidx][seq_len(4)]))
r <- min((max(res_pca$rotated[, "PC1"]) - min(res_pca$rotated[, "PC1"])/(max(res_pca$loadings[, 
            "PC1"]) - min(res_pca$loadings[,"PC1"]))), (max(res_pca$rotated[, "PC2"]) - min(res_pca$rotated[, 
            "PC2"])/(max(res_pca$loadings[, "PC2"]) - min(res_pca$loadings[, "PC2"]))))

pca.plot <- ggplot() + 
    geom_point(aes(x=res_pca$rotated[, "PC1"],y=res_pca$rotated[, "PC2"],shape=res_pca$metadata$TYPE,col=res_pca$metadata$SPECIES), size=3) + 
    geom_text_repel(aes(x=res_pca$rotated[, "PC1"],y=res_pca$rotated[, "PC2"],label=res_pca$yvars))+ 
    geom_segment(data = res_pca$loadings[vars,], aes(x = 0, y = 0, xend = res_pca$loadings[vars,"PC1"] * r * 1, 
                yend = res_pca$loadings[vars,"PC2"] * r * 1), arrow = arrow(length = unit(0.3, "cm"))) +
    geom_label_repel(aes(res_pca$loadings[vars,"PC1"] * r * 1,res_pca$loadings[vars,"PC2"] * r * 1,label=vars, alpha = 0.6)) + 
    xlab(paste0("PC1, ", round(res_pca$variance["PC1"], digits = 2), "% variation")) +
    ylab(paste0("PC2, ", round(res_pca$variance["PC2"], digits = 2), "% variation")) +
    theme_light() + theme(legend.position = "none")

#split to bulk and root
df.clr.bulk <- df.clr[which(substr(rownames(df.clr),1,1) == "B"),]
metadata.bulk <- metadata[which(substr(rownames(metadata),1,1) == "B"),]
df.clr.root <- df.clr[which(substr(rownames(df.clr),1,1) == "R"),]
metadata.root <- metadata[which(substr(rownames(metadata),1,1) == "R"),]

#get euclidean distance
load("data/bukit_timah_big_trees.stem1.rdata")
bukit_timah_big_trees.stem1$tag <- as.numeric(bukit_timah_big_trees.stem1$tag)
trees <- bukit_timah_big_trees.stem1[bukit_timah_big_trees.stem1$tag %in% metadata$TAG_ID,]
trees <- trees[,c("tag", "gx", "gy")]
trees[,c("gx","gy")] <- trees[,c("gx","gy")] - rep(apply(trees[,c("gx","gy")], 2, min),each=15)
trees[,c("gx","gy")] <- trees[,c("gx","gy")] / rep(apply(trees[,c("gx","gy")], 2, max),each=15)
samples <- metadata.root$INDEX[match(trees$tag, metadata$TAG_ID)]
rownames(trees) <- c(paste0("R", samples))#,paste0("R", samples))
trees <- trees[order(row.names(trees)), ]
metadata.root[,c("gx", "gy")] <- trees[,c("gx","gy")]
metadata[,c("gx", "gy")] <- rbind(trees[,c("gx","gy")],trees[,c("gx","gy")])

#normal rda
ord.rda <- rda(df.clr ~ SPECIES+TYPE + Condition(gx+gy), metadata)
ord.rda2 <- rda(df.clr ~ SPECIES+TYPE, metadata)
anova(ord.rda, permutations = how(nperm=100000)) # check model sig
anova(ord.rda2, ord.rda, permutations = how(nperm=100000)) # check condition sig
anova(ord.rda, by="term", permutations = how(nperm=100000)) # get term sig
anova(ord.rda, by="axis", permutations = how(nperm=100000)) # get axis sig
#anova(ord.rda, by="margin", permutations = how(nperm=100000))

#root rda
# Root ~ species
ord.rda.root.sp <- rda(df.clr.root ~ SPECIES + Condition(gx+gy), metadata.root) 
ord.rda2.root.sp <- rda(df.clr.root ~ SPECIES, metadata.root)
anova(ord.rda.root.sp, permutations = how(nperm=100000)) # check model sig
anova(ord.rda2.root.sp, ord.rda.root.sp, permutations = how(nperm=100000)) # check condition sig
anova(ord.rda.root.sp, by="term", permutations = how(nperm=100000)) # get term sig
anova(ord.rda.root.sp, by="axis", permutations = how(nperm=100000)) # get axis sig
# Root ~ family
ord.rda.root.f <- rda(df.clr.root ~ FAMILY + Condition(gx+gy), metadata.root)
ord.rda2.root.f <- rda(df.clr.root ~ FAMILY, metadata.root)
anova(ord.rda.root.f, permutations = how(nperm=100000)) # check model sig
anova(ord.rda2.root.f, ord.rda.root.f, permutations = how(nperm=100000)) # check condition sig
anova(ord.rda.root.f, by="term", permutations = how(nperm=100000)) # get term sig
anova(ord.rda.root.f, by="axis", permutations = how(nperm=100000)) # get axis sig
# Root ~ Habitat
ord.rda.root.h <- rda(df.clr.root ~ HABITAT + Condition(gx+gy), metadata.root) 
ord.rda2.root.h <- rda(df.clr.root ~ HABITAT, metadata.root)
anova(ord.rda.root.h, permutations = how(nperm=100000)) # check model sig
anova(ord.rda2.root.h, ord.rda.root.h, permutations = how(nperm=100000)) # check condition sig
anova(ord.rda.root.h, by="term", permutations = how(nperm=100000)) # get term sig
anova(ord.rda.root.h, by="axis", permutations = how(nperm=100000)) # get axis sig

#rda plotting for bulk/root + species
allscores <- scores(ord.rda, display = "sites")
allscores <- cbind(allscores, metadata.root)
centroids <- as.data.frame(scores(ord.rda, display="cn"))
speciesscore <- as.data.frame(scores(ord.rda, display = "species"))
speciesscore$dist <- sqrt(speciesscore$RDA1^2 + speciesscore$RDA2^2)
speciesscore_top <- speciesscore[order(speciesscore$dist, decreasing=T),][1:8,]
rda.summary <- summary(ord.rda)
rda.plot <- ggplot() + 
    geom_point(aes(x=RDA1,y=RDA2,shape=TYPE,col=SPECIES), data=allscores, size=3) + 
    labs(shape = "Type", col = "Species") + 
    geom_text_repel(aes(x=RDA1,y=RDA2,label=rownames(allscores)),data = allscores) + 
    geom_point(aes(RDA1,RDA2), data=centroids, shape=23, fill="grey25", col="grey39", alpha = 0.7, size = 3) + 
    geom_text_repel(aes(RDA1,RDA2, label=c("KM", "SC", "SL", "Bulk", "Root")), data=centroids) +
    geom_segment(aes(x=0,y=0,xend=RDA1,yend=RDA2), data=speciesscore_top,arrow = arrow(length = unit(0.3, "cm"))) +
    geom_label_repel(aes(RDA1,RDA2, label=row.names(speciesscore_top)), data=speciesscore_top, alpha = 0.6) +
    xlab(paste0("RDA1, ", round(rda.summary$cont$importance[2,1]*100, digits = 2), "% variation")) +
    ylab(paste0("RDA2, ", round(rda.summary$cont$importance[2,2]*100, digits = 2), "% variation")) +
    theme_light()
rda.plot2 <- rda.plot+ theme(legend.position = "none")
legend <- get_legend(rda.plot)
plots <- plot_grid(pca.plot, rda.plot2, align = "v", axis = "lr", ncol=1, labels="AUTO",label_size = 16)
#Figure 5
svg("ordi_combined.svg", width=12, height=15)
plot_grid(plots,legend,rel_widths = c(1,0.1))
dev.off()

#rda plotting for Root ~ Species
allscores <- scores(ord.rda.root.sp, display = "sites")
allscores <- cbind(allscores, metadata.root)
centroids <- as.data.frame(scores(ord.rda.root.sp, display="cn"))
speciesscore <- as.data.frame(scores(ord.rda.root.sp, display = "species"))
speciesscore$dist <- sqrt(speciesscore$RDA1^2 + speciesscore$RDA2^2)
speciesscore_top <- speciesscore[order(speciesscore$dist, decreasing=T),][1:6,]
speciesscore_top$id <- rownames(speciesscore_top)
rda.summary <- summary(ord.rda.root.sp)
rda.root.sp.plot <- ggplot() + 
    geom_point(aes(x=RDA1,y=RDA2,col=SPECIES), data=allscores, size=3) + 
    labs(col = "Species") + 
    geom_text_repel(aes(x=RDA1,y=RDA2,label=rownames(allscores)),data = allscores) + 
    geom_point(aes(RDA1,RDA2), data=centroids, shape=23, fill="grey25", col="grey39", alpha = 0.7, size = 3) + 
    geom_text_repel(aes(RDA1,RDA2, label=c("KM", "SC", "SL")), data=centroids) +
    geom_segment(aes(x=0,y=0,xend=RDA1,yend=RDA2), data=speciesscore_top,arrow = arrow(length = unit(0.3, "cm"))) +
    geom_label_repel(aes(RDA1,RDA2, label=id), data=speciesscore_top, alpha = 0.6) +
    xlab(paste0("RDA1, ", round(rda.summary$cont$importance[2,1]*100, digits = 2), "% variation")) +
    ylab(paste0("RDA2, ", round(rda.summary$cont$importance[2,2]*100, digits = 2), "% variation")) +
    theme_light()

#rda plotting for Root ~ Family
allscores <- scores(ord.rda.root.f, display = "sites")
allscores <- cbind(allscores, metadata.root)
centroids <- as.data.frame(scores(ord.rda.root.f, display="cn"))
speciesscore <- as.data.frame(scores(ord.rda.root.f, display = "species"))
speciesscore$dist <- sqrt(speciesscore$RDA1^2 + speciesscore$PC1^2)
speciesscore_top2 <- speciesscore[order(speciesscore$dist, decreasing=T),][1:6,]
speciesscore_top2$id <- rownames(speciesscore_top2)
rda.summary <- summary(ord.rda.root.f)
rda.root.f.plot <- ggplot() + 
    geom_point(aes(x=RDA1,y=PC1,col=FAMILY), data=allscores, size=3) + 
    labs(col = "Family") + 
    geom_text_repel(aes(x=RDA1,y=PC1,label=rownames(allscores)),data = allscores) + 
    geom_point(aes(RDA1,PC1), data=centroids, shape=23, fill="grey25", col="grey39", alpha = 0.7, size = 3) + 
    geom_text_repel(aes(RDA1,PC1, label=c("Dipterocarpaceae", "Fabaceae")), data=centroids) +
    geom_segment(aes(x=0,y=0,xend=RDA1,yend=PC1), data=speciesscore_top2,arrow = arrow(length = unit(0.3, "cm"))) +
    geom_label_repel(aes(RDA1,PC1, label=id), data=speciesscore_top2, alpha = 0.6) +
    xlab(paste0("RDA1, ", round(rda.summary$cont$importance[2,1]*100, digits = 2), "% variation")) +
    ylab(paste0("PC1, ", round(rda.summary$cont$importance[2,2]*100, digits = 2), "% variation")) +
    theme_light()

#rda plotting for Root ~ Habitat
allscores <- scores(ord.rda.root.h, display = "sites")
allscores <- cbind(allscores, metadata.root)
centroids <- as.data.frame(scores(ord.rda.root.h, display="cn"))
speciesscore <- as.data.frame(scores(ord.rda.root.h, display = "species"))
speciesscore$dist <- sqrt(speciesscore$RDA1^2 + speciesscore$PC1^2)
speciesscore_top3 <- speciesscore[order(speciesscore$dist, decreasing=T),][1:6,]
speciesscore_top3$id <- rownames(speciesscore_top3)
rda.summary <- summary(ord.rda.root.h)
rda.root.h.plot <- ggplot() + 
    geom_point(aes(x=RDA1,y=PC1,col=HABITAT), data=allscores, size=3) + 
    labs(col = "Habitat") + 
    geom_text_repel(aes(x=RDA1,y=PC1,label=rownames(allscores)),data = allscores) + 
    geom_point(aes(RDA1,PC1), data=centroids, shape=23, fill="grey25", col="grey39", alpha = 0.7, size = 3) + 
    geom_text_repel(aes(RDA1,PC1, label=c("Highland", "Lowland")), data=centroids) +
    geom_segment(aes(x=0,y=0,xend=RDA1,yend=PC1), data=speciesscore_top3,arrow = arrow(length = unit(0.3, "cm"))) +
    geom_label_repel(aes(RDA1,PC1, label=id), data=speciesscore_top3, alpha = 0.6) +
    xlab(paste0("RDA1, ", round(rda.summary$cont$importance[2,1]*100, digits = 2), "% variation")) +
    ylab(paste0("PC1, ", round(rda.summary$cont$importance[2,2]*100, digits = 2), "% variation")) +
    theme_light()

#Figure 6
svg("rda_root.svg", width=12, height=21)
plot_grid(rda.root.sp.plot,rda.root.f.plot,rda.root.h.plot,
          ncol = 1, align = "v", axis = "lr", labels="AUTO", label_size = 16)
dev.off()

aldex_res <- aldex(t(df), conditions=metadata$TYPE, mc.samples = 128, test="t", effect=T,denom = "all")
aldex.plot(aldex_res, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
top_ip2g <- aldex_res[aldex_res$wi.eBH < 0.05,]

#ANCOMBC 
library(ANCOMBC)
library(phyloseq)
tax_mat <- otu_table(df, taxa_are_rows = F)
samp_data <- sample_data(metadata)
ps <- phyloseq(tax_mat, samp_data)
ancom_out <- ancombc(ps, formula="TYPE*SPECIES", conserve=T, group="SPECIES", global=F)
ancom_out$res$diff_abn


#propr
library(ALDEx2)
library(propr)
alde.clr <- aldex.clr(t(df), conds=metadata$TYPE, mc.samples = 128, denom="all", verbose=FALSE, useMC=FALSE)
df.rho <- aldex2propr(alde.clr)
diag(df.rho@matrix) <- 0
rownames(df.rho@matrix) <- colnames(df.rho@counts)
colnames(df.rho@matrix) <- colnames(df.rho@counts)
d.rho <- df.rho

d.hi.rho <- d.rho[">", 0.85]
# get the positions of the OTU pairs in the matrix
pairs <- arrayInd(d.hi.rho@pairs,dim(d.hi.rho@matrix))

# replace indices with OTU ids
pairs.OTU <- apply(pairs, 2, function(x){ colnames(d.rho@counts)[x] })

coeff <- matrix(data=NA, nrow=nrow(pairs), ncol=3)
for(i in 1:nrow(pairs)){
    coeff[i,c(1,2)] <- round(as.numeric(lm(d.rho@logratio[,pairs[i,1]]~d.rho@logratio[,pairs[i,2]])$coefficients), 2)
    coeff[i,3] <- round(as.numeric(cor(d.rho@logratio[,pairs[i,1]] , d.rho@logratio[,pairs[i,2]])),2)
}
pairs.info <- cbind(pairs.OTU[,1], df[pairs.OTU[,1],"genus"], pairs.OTU[,2], df[pairs.OTU[,2],"genus"], round(d.rho@matrix[d.hi.rho@pairs],2),coeff)

colnames(pairs.info) <- c("OTU1","OTU2", "E(rho)", "intercept", "slope", "corr")

# igraph: convert the connections into a graphical object
g <- igraph::graph.data.frame(pairs.OTU, directed=FALSE)

# igraph: find the clusters
g.clust <- igraph::clusters(g)

# make a table to examine the cluster membership by hand
g.df <- data.frame(Systematic.name=V(g)$name, cluster=g.clust$membership,
                   cluster.size=g.clust$csize[g.clust$membership])

# generate a set of clusters larger than some size # minimum is 2 (obviously)
big <- g.df[which(g.df$cluster.size >= 3),]
colnames(big) <- colnames(g.df)
big$genus <- t.df[rownames(big),"genus"]

#ANCOM2
source("ancom_v2.1.R")
metadata$sample <- rownames(metadata)
prepro <- feature_table_pre_process(t(df), meta_data=metadata, sample_var="sample", group_var="TYPE", lib_cut=1, neg_lb=F)
feature_table <- prepro$feature_table # Preprocessed feature table
meta_data <- prepro$meta_data # Preprocessed metadata
struc_zero <- prepro$structure_zeros # Structural zero info
res = ANCOM(feature_table, meta_data, struc_zero, main_var = "TYPE")

n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig = res$fig +  
    geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
    geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
              size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig  
