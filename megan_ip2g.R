library(ggplot2)
library(cowplot)
library(zCompositions)
library(compositions)
library(PCAtools)
library(vegan)

#read in data
df <- read.csv("data/rma_INTERPRO2GO.tsv", sep="\t")
rownames(df) <- df$NUM
df <- df[,-c(1,2)]
df <- as.data.frame(t(df))
metadata <- read.csv("data/metadata.tsv", sep="\t")
rownames(metadata) <- paste0(metadata$SAMPLE_ID, "_",metadata$SPECIES)
metadata$BATCH <- as.factor(metadata$BATCH)

#save files 
save(df,metadata, file="data/5m_rma_ip2g.RData")

#Start from here 
load(file="data/5m_rma_ip2g.RData")

#clr preprocessing
df <- df[,colSums(df==0) < 25]
df.nozero <- cmultRepl(t(df), label=0, method="CZM", output="p-counts")
df.clr <- clr(df.nozero)
df.clr <- t(as.data.frame(df.clr))

#pca
res_pca <- pca(df.clr, center = F, scale = T, transposed = T, metadata = metadata)
# screeplot(res_pca) #only N- PCs due to centering
# biplot(res_pca, showLoadings = T, colby = "BATCH", shape = "TYPE", ntopLoadings=3,xlim=c(-6,6),ylim=c(-6,6), ellipse = F) + 
#     theme_light()
# pairsplot(res_pca, components = getComponents(res_pca, seq_len(5)), colby = "SPECIES", shape = "TYPE", 
#           pointSize = 2, plotaxes = F, lab=rownames(df),labSize = 1.4) 
# plotloadings(res_pca, labSize = 2, shapeSizeRange=c(5,5))

#custom pca
xidx <- order(abs(res_pca$loadings[, "PC1"]), decreasing = TRUE)
yidx <- order(abs(res_pca$loadings[, "PC2"]), decreasing = TRUE)
vars <- unique(c(rownames(res_pca$loadings)[xidx][seq_len(2)], 
                 rownames(res_pca$loadings)[yidx][seq_len(2)]))
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
ord.rda <- rda(df.clr ~ SPECIES + TYPE + Condition(gx+gy), metadata)
ord.rda2 <- rda(df.clr ~ TYPE, metadata)
anova(ord.rda, permutations = how(nperm=5000)) # check model sig
anova(ord.rda2, ord.rda, permutations = how(nperm=5000)) # check condition sig
anova(ord.rda, by="term", permutations = how(nperm=5000)) # get term sig
anova(ord.rda, by="axis", permutations = how(nperm=5000)) # get axis sig
#anova(ord.rda, by="margin", permutations = how(nperm=50000))

#root rda
# Root ~ species
ord.rda.root.sp <- rda(df.clr.root ~ SPECIES + Condition(gx+gy), metadata.root) 
ord.rda2.root.sp <- rda(df.clr.root ~ SPECIES, metadata.root)
anova(ord.rda.root.sp, permutations = how(nperm=5000)) # check model sig
anova(ord.rda2.root.sp, ord.rda.root.sp, permutations = how(nperm=5000)) # check condition sig
anova(ord.rda.root.sp, by="term", permutations = how(nperm=5000)) # get term sig
anova(ord.rda.root.sp, by="axis", permutations = how(nperm=5000)) # get axis sig
# Root ~ family
ord.rda.root.f <- rda(df.clr.root ~ FAMILY + Condition(gx+gy), metadata.root)
ord.rda2.root.f <- rda(df.clr.root ~ FAMILY, metadata.root)
anova(ord.rda.root.f, permutations = how(nperm=5000)) # check model sig
anova(ord.rda2.root.f, ord.rda.root.f, permutations = how(nperm=5000)) # check condition sig
anova(ord.rda.root.f, by="term", permutations = how(nperm=5000)) # get term sig
anova(ord.rda.root.f, by="axis", permutations = how(nperm=5000)) # get axis sig
# Root ~ Habitat
ord.rda.root.h <- rda(df.clr.root ~ HABITAT + Condition(gx+gy), metadata.root) 
ord.rda2.root.h <- rda(df.clr.root ~ HABITAT, metadata.root)
anova(ord.rda.root.h, permutations = how(nperm=5000)) # check model sig
anova(ord.rda2.root.h, ord.rda.root.h, permutations = how(nperm=5000)) # check condition sig
anova(ord.rda.root.h, by="term", permutations = how(nperm=5000)) # get term sig
anova(ord.rda.root.h, by="axis", permutations = how(nperm=5000)) # get axis sig


allscores <- scores(ord.rda, display = "sites")
allscores <- cbind(allscores, metadata)
centroids <- as.data.frame(scores(ord.rda, display="cn"))
speciesscore <- as.data.frame(scores(ord.rda, display = "species"))
speciesscore$dist <- sqrt(speciesscore$RDA1^2 + speciesscore$RDA2^2)
speciesscore_top <- speciesscore[order(speciesscore$dist, decreasing=T),][1:6,] # > 0.6 is shown
rda.summary <- summary(ord.rda)
rda.plot <- ggplot() + 
    geom_point(aes(x=RDA1,y=RDA2,shape=TYPE,col=SPECIES), data=allscores, size=3) + 
    labs(shape = "Type", col = "Species") + 
    geom_text_repel(aes(x=RDA1,y=RDA2,label=rownames(allscores)),data = allscores)+ 
    geom_point(aes(RDA1,RDA2), data=centroids, shape=23, fill="grey25", col="grey39", alpha = 0.7, size = 3) + 
    geom_text_repel(aes(RDA1,RDA2, label=c("KM", "SC", "SL", "Bulk", "Root")), data=centroids) +
    geom_segment(aes(x=0,y=0,xend=RDA1,yend=RDA2), data=speciesscore_top,arrow = arrow(length = unit(0.3, "cm"))) +
    geom_label_repel(aes(RDA1,RDA2, label=row.names(speciesscore_top)), data=speciesscore_top, alpha = 0.6) +
    xlab(paste0("RDA1, ", round(rda.summary$cont$importance[2,1]*100, digits = 2), "% variation")) +
    ylab(paste0("RDA2, ", round(rda.summary$cont$importance[2,2]*100, digits = 2), "% variation")) +
    theme_light()
#cowplot to combine
rda.plot2 <- rda.plot+ theme(legend.position = "none")
legend <- get_legend(rda.plot)
plots <- plot_grid(pca.plot, rda.plot2, align = "v", axis = "lr", ncol=1, labels="AUTO",label_size = 16)
svg("ordi_ip2g_combined.svg", width=12, height=15)
plot_grid(plots,legend,rel_widths = c(1,0.1))
dev.off()

#rda plotting for Root ~ Species
allscores <- scores(ord.rda.root.sp, display = "sites")
allscores <- cbind(allscores, metadata.root)
centroids <- as.data.frame(scores(ord.rda.root.sp, display="cn"))
speciesscore <- as.data.frame(scores(ord.rda.root.sp, display = "species"))
speciesscore$dist <- sqrt(speciesscore$RDA1^2 + speciesscore$RDA2^2)
speciesscore_top <- speciesscore[order(speciesscore$dist, decreasing=T),][1:8,]
rda.summary <- summary(ord.rda.root.sp)
rda.root.sp.plot <- ggplot() + 
    geom_point(aes(x=RDA1,y=RDA2,col=SPECIES), data=allscores, size=3) + 
    labs(col = "Species") + 
    geom_text_repel(aes(x=RDA1,y=RDA2,label=rownames(allscores)),data = allscores) + 
    geom_point(aes(RDA1,RDA2), data=centroids, shape=23, fill="grey25", col="grey39", alpha = 0.7, size = 3) + 
    geom_text_repel(aes(RDA1,RDA2, label=c("KM", "SC", "SL")), data=centroids) +
    geom_segment(aes(x=0,y=0,xend=RDA1,yend=RDA2), data=speciesscore_top,arrow = arrow(length = unit(0.3, "cm"))) +
    geom_label_repel(aes(RDA1,RDA2, label=row.names(speciesscore_top)), data=speciesscore_top, alpha = 0.6) +
    xlab(paste0("RDA1, ", round(rda.summary$cont$importance[2,1]*100, digits = 2), "% variation")) +
    ylab(paste0("RDA2, ", round(rda.summary$cont$importance[2,2]*100, digits = 2), "% variation")) +
    theme_light()

#rda plotting for Root ~ Family
allscores <- scores(ord.rda.root.f, display = "sites")
allscores <- cbind(allscores, metadata.root)
centroids <- as.data.frame(scores(ord.rda.root.f, display="cn"))
speciesscore <- as.data.frame(scores(ord.rda.root.f, display = "species"))
speciesscore$dist <- sqrt(speciesscore$RDA1^2 + speciesscore$PC1^2)
speciesscore_top <- speciesscore[order(speciesscore$dist, decreasing=T),][1:8,]
rda.summary <- summary(ord.rda.root.f)
rda.root.f.plot <- ggplot() + 
    geom_point(aes(x=RDA1,y=PC1,col=FAMILY), data=allscores, size=3) + 
    labs(col = "Family") + 
    geom_text_repel(aes(x=RDA1,y=PC1,label=rownames(allscores)),data = allscores) + 
    geom_point(aes(RDA1,PC1), data=centroids, shape=23, fill="grey25", col="grey39", alpha = 0.7, size = 3) + 
    geom_text_repel(aes(RDA1,PC1, label=c("Dipterocarpaceae", "Fabaceae")), data=centroids) +
    geom_segment(aes(x=0,y=0,xend=RDA1,yend=PC1), data=speciesscore_top,arrow = arrow(length = unit(0.3, "cm"))) +
    geom_label_repel(aes(RDA1,PC1, label=row.names(speciesscore_top)), data=speciesscore_top, alpha = 0.6) +
    xlab(paste0("RDA1, ", round(rda.summary$cont$importance[2,1]*100, digits = 2), "% variation")) +
    ylab(paste0("PC1, ", round(rda.summary$cont$importance[2,2]*100, digits = 2), "% variation")) +
    theme_light()

#rda plotting for Root ~ Habitat
allscores <- scores(ord.rda.root.h, display = "sites")
allscores <- cbind(allscores, metadata.root)
centroids <- as.data.frame(scores(ord.rda.root.h, display="cn"))
speciesscore <- as.data.frame(scores(ord.rda.root.h, display = "species"))
speciesscore$dist <- sqrt(speciesscore$RDA1^2 + speciesscore$PC1^2)
speciesscore_top <- speciesscore[order(speciesscore$dist, decreasing=T),][1:8,]
rda.summary <- summary(ord.rda.root.h)
rda.root.h.plot <- ggplot() + 
    geom_point(aes(x=RDA1,y=PC1,col=HABITAT), data=allscores, size=3) + 
    labs(col = "Habitat") + 
    geom_text_repel(aes(x=RDA1,y=PC1,label=rownames(allscores)),data = allscores) + 
    geom_point(aes(RDA1,PC1), data=centroids, shape=23, fill="grey25", col="grey39", alpha = 0.7, size = 3) + 
    geom_text_repel(aes(RDA1,PC1, label=c("Highland", "Lowland")), data=centroids) +
    geom_segment(aes(x=0,y=0,xend=RDA1,yend=PC1), data=speciesscore_top,arrow = arrow(length = unit(0.3, "cm"))) +
    geom_label_repel(aes(RDA1,PC1, label=row.names(speciesscore_top)), data=speciesscore_top, alpha = 0.6) +
    xlab(paste0("RDA1, ", round(rda.summary$cont$importance[2,1]*100, digits = 2), "% variation")) +
    ylab(paste0("PC1, ", round(rda.summary$cont$importance[2,2]*100, digits = 2), "% variation")) +
    theme_light()

#Figure 6
svg("rda_root_ip2g.svg", width=12, height=21)
plot_grid(rda.root.sp.plot,rda.root.f.plot,rda.root.h.plot,
          ncol = 1, align = "v", axis = "lr", labels="AUTO", label_size = 16)
dev.off()

#ANCOMBC 
library(ANCOMBC)
library(phyloseq)
tax_mat <- otu_table(df, taxa_are_rows = F)
samp_data <- sample_data(metadata)
ps <- phyloseq(tax_mat, samp_data)
ancom_out <- ancombc(ps, formula="TYPE", conserve=T, group="SPECIES", global=F)
ancom_out$res$diff_abn

#aldex
library(ALDEx2)
#library(propr)
aldex_res <- aldex(t(df), conditions=metadata$TYPE, mc.samples = 128, test="t", effect=T,denom = "all", )
df.root <- df[which(substr(rownames(df),1,1) == "R"),]
aldex_res.root <- aldex(t(df.root), conditions=metadata.root$SPECIES, mc.samples = 128, test="kw",denom = "all")
aldex.plot(aldex_res, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
top_ip2g_down <- rownames(aldex_res[aldex_res$wi.eBH < 0.05 & aldex_res$effect < -1,])
top_ip2g_top <- rownames(aldex_res[aldex_res$wi.eBH < 0.05 & aldex_res$effect > 1,])
top_ip2g.root <- rownames(aldex_res[aldex_res.root$kw.eBH < 0.05,])

write.table(stringr::str_pad(top_ip2g_down, 6, pad="0"), file="ip2g_down.txt", quote=F, row.names = F, col.names = F)
write.table(stringr::str_pad(top_ip2g_top, 6, pad="0"), file="ip2g_up.txt", quote=F, row.names = F, col.names = F)
write.table(stringr::str_pad(colnames(df), 6, pad="0"), file="background.txt", quote=F, row.names = F, col.names = F)

#ANCOM2
source("ancom_v2.1.R")
metadata$sample <- rownames(metadata)
metadata.root <- metadata[which(substr(rownames(metadata),1,1) == "R"),]
df.root <- df[which(substr(rownames(metadata),1,1) == "R"),]
prepro <- feature_table_pre_process(t(df.root), meta_data=metadata.root, sample_var= "sample", group_var="SPECIES", 
                                   lib_cut=1, neg_lb=F)
feature_table <- prepro$feature_table # Preprocessed feature table
meta_data <- prepro$meta_data # Preprocessed metadata
struc_zero <- prepro$structure_zeros # Structural zero info
res <- ANCOM(feature_table, meta_data, struc_zero, main_var="SPECIES")
