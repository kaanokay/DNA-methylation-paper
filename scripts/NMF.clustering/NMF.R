### Nonnegative Matrix Factorization (NMF) classification of EM genes relying on DNAm of mouse protein coding genes


library(bsseq)
library(NMF)
library(factoextra)
library(dendextend)
library(tidyverse)
library(pheatmap)


### Upload BSseq object

BSseq.obj <- readRDS("BSSeq.object.smoothed.no.any.filtering.rds")

### Remove outliers

filtered_data <- subset(pData(BSseq.obj), !(metadata == "Wild_type.Batch.13.07.2023" |
                                                metadata == "Ctr1.Batch.09.07.2023"))


BSseq.obj <- BSseq.obj[, rownames(filtered_data)]



dim(pData(BSseq.obj))


### Remove CpGs with low coverage


cov <- getCoverage(BSseq.obj)
keepLoci.ex <- which(rowSums(cov[, BSseq.obj$metadata.2 == "EM_gene"] >= 3) >= 36 &
                       rowSums(cov[, BSseq.obj$metadata.2 == "Control"] >= 3) >= 10)



length(keepLoci.ex)



BSseq.obj <- BSseq.obj[keepLoci.ex, ]


### Load promoter regions of all mouse protein coding genes


load("promoters.of.all.mouse.protein.coding.genes.rda")


### Get smoothed methylation of promoters of all mouse-protein coding genes (i.e., 20482 genes in total)


methylation.values.for.all <- getMeth(BSseq.obj, type = "smooth", what = "perRegion", regions = promoters.of.all.mouse.protein.coding.genes)


### Assign NCBI gene ids to rows


rownames(methylation.values.for.all) <- promoters.of.all.mouse.protein.coding.genes$gene_id


### Summarize the number of NA values in each row


na_count_per_row <- rowSums(is.na(methylation.values.for.all))


### Identify rows that have one or more NA values


rows_with_na <- which(na_count_per_row > 0)
print(rows_with_na)


### Remove rows with NA


methylation.values.for.all.genes.filtered <- na.omit(methylation.values.for.all)


### Convert object into data frame


methylation.values.for.all.genes.filtered <- as.data.frame(methylation.values.for.all.genes.filtered)


### Merge samples by averaging DNAm values


Dnmt1.samples <- methylation.values.for.all.genes.filtered[, grepl("Dnmt1", names(methylation.values.for.all.genes.filtered))]
Tet1.samples <- methylation.values.for.all.genes.filtered[, grepl("Tet1", names(methylation.values.for.all.genes.filtered))]
Chd1.samples <- methylation.values.for.all.genes.filtered[, grepl("Chd1", names(methylation.values.for.all.genes.filtered))]
Chd4.samples <- methylation.values.for.all.genes.filtered[, grepl("Chd4", names(methylation.values.for.all.genes.filtered))]
Smarca1.samples <- methylation.values.for.all.genes.filtered[, grepl("Smarca1", names(methylation.values.for.all.genes.filtered))]
Kmt2a.samples <- methylation.values.for.all.genes.filtered[, grepl("Kmt2a", names(methylation.values.for.all.genes.filtered))]
Rere.samples <- methylation.values.for.all.genes.filtered[, grepl("Rere", names(methylation.values.for.all.genes.filtered))]
Srcap.samples <- methylation.values.for.all.genes.filtered[, grepl("Srcap", names(methylation.values.for.all.genes.filtered))]
Kdm3b.samples <- methylation.values.for.all.genes.filtered[, grepl("Kdm3b", names(methylation.values.for.all.genes.filtered))]
Ctrl.samples <- methylation.values.for.all.genes.filtered[, grepl("Ctr", names(methylation.values.for.all.genes.filtered))]


### Check whether all samples subset


dim(Dnmt1.samples)
dim(Tet1.samples)
dim(Chd1.samples)
dim(Chd4.samples)
dim(Smarca1.samples)
dim(Kmt2a.samples)
dim(Rere.samples)
dim(Srcap.samples)
dim(Kdm3b.samples)
dim(Ctrl.samples)


### Average DNAm for each promoter region across all samples


Dnmt1.samples.merged <- as.data.frame(rowSums(Dnmt1.samples)/ncol(Dnmt1.samples))
Tet1.samples.merged <- as.data.frame(rowSums(Tet1.samples)/ncol(Tet1.samples))
Chd1.samples.merged <- as.data.frame(rowSums(Chd1.samples)/ncol(Chd1.samples))
Chd4.samples.merged <- as.data.frame(rowSums(Chd4.samples)/ncol(Chd4.samples))
Smarca1.samples.merged <- as.data.frame(rowSums(Smarca1.samples)/ncol(Smarca1.samples))
Kmt2a.samples.merged <- as.data.frame(rowSums(Kmt2a.samples)/ncol(Kmt2a.samples))
Rere.samples.merged <- as.data.frame(rowSums(Rere.samples)/ncol(Rere.samples))
Srcap.samples.merged <- as.data.frame(rowSums(Srcap.samples)/ncol(Srcap.samples))
Kdm3b.samples.merged <- as.data.frame(rowSums(Kdm3b.samples)/ncol(Kdm3b.samples))
Ctrl.samples.merged <- as.data.frame(rowSums(Ctrl.samples)/ncol(Ctrl.samples))


### Check whether samples merged into one


dim(Dnmt1.samples.merged)
dim(Tet1.samples.merged)
dim(Chd1.samples.merged)
dim(Chd4.samples.merged)
dim(Smarca1.samples.merged)
dim(Kmt2a.samples.merged)
dim(Rere.samples.merged)
dim(Srcap.samples.merged)
dim(Kdm3b.samples.merged)
dim(Ctrl.samples.merged)


### Check whether DNAm values are ranging from 0 to 1


range(Dnmt1.samples.merged)
range(Tet1.samples.merged)
range(Chd1.samples.merged)
range(Chd4.samples.merged)
range(Smarca1.samples.merged)
range(Kmt2a.samples.merged)
range(Rere.samples.merged)
range(Srcap.samples.merged)
range(Kdm3b.samples.merged)
range(Ctrl.samples.merged)


### Remove individual sample that have replicates (i.e., Dnmt1, Tet1 etc.)


methylation.values.for.all.genes.filtered <- methylation.values.for.all.genes.filtered[, !grepl("Dnmt1", names(methylation.values.for.all.genes.filtered))]
methylation.values.for.all.genes.filtered <- methylation.values.for.all.genes.filtered[, !grepl("Tet1", names(methylation.values.for.all.genes.filtered))]
methylation.values.for.all.genes.filtered <- methylation.values.for.all.genes.filtered[, !grepl("Chd1", names(methylation.values.for.all.genes.filtered))]
methylation.values.for.all.genes.filtered <- methylation.values.for.all.genes.filtered[, !grepl("Chd4", names(methylation.values.for.all.genes.filtered))]
methylation.values.for.all.genes.filtered <- methylation.values.for.all.genes.filtered[, !grepl("Smarca1", names(methylation.values.for.all.genes.filtered))]
methylation.values.for.all.genes.filtered <- methylation.values.for.all.genes.filtered[, !grepl("Kmt2a", names(methylation.values.for.all.genes.filtered))]
methylation.values.for.all.genes.filtered <- methylation.values.for.all.genes.filtered[, !grepl("Rere", names(methylation.values.for.all.genes.filtered))]
methylation.values.for.all.genes.filtered <- methylation.values.for.all.genes.filtered[, !grepl("Srcap", names(methylation.values.for.all.genes.filtered))]
methylation.values.for.all.genes.filtered <- methylation.values.for.all.genes.filtered[, !grepl("Kdm3b", names(methylation.values.for.all.genes.filtered))]
methylation.values.for.all.genes.filtered <- methylation.values.for.all.genes.filtered[, !grepl("Ctr", names(methylation.values.for.all.genes.filtered))]


### Add samples that were avaraged on their DNAm


methylation.values.for.all.genes.filtered <- cbind(methylation.values.for.all.genes.filtered, Dnmt1.samples.merged)
methylation.values.for.all.genes.filtered <- cbind(methylation.values.for.all.genes.filtered, Tet1.samples.merged)
methylation.values.for.all.genes.filtered <- cbind(methylation.values.for.all.genes.filtered, Chd1.samples.merged)
methylation.values.for.all.genes.filtered <- cbind(methylation.values.for.all.genes.filtered, Chd4.samples.merged)
methylation.values.for.all.genes.filtered <- cbind(methylation.values.for.all.genes.filtered, Smarca1.samples.merged)
methylation.values.for.all.genes.filtered <- cbind(methylation.values.for.all.genes.filtered, Kmt2a.samples.merged)
methylation.values.for.all.genes.filtered <- cbind(methylation.values.for.all.genes.filtered, Rere.samples.merged)
methylation.values.for.all.genes.filtered <- cbind(methylation.values.for.all.genes.filtered, Srcap.samples.merged)
methylation.values.for.all.genes.filtered <- cbind(methylation.values.for.all.genes.filtered, Kdm3b.samples.merged)
methylation.values.for.all.genes.filtered <- cbind(methylation.values.for.all.genes.filtered, Ctrl.samples.merged)


### Check how many samples left in data, should be 38 (37 EM KO plus controls)


ncol(methylation.values.for.all.genes.filtered)


### Rename columns


colnames(methylation.values.for.all.genes.filtered) <- c("Kdm5c", "Kdm5b", "Kdm2b", "Kdm1a", "Hdac8", "Kdm6a", "Hdac6", "Ezh2",
 "Crebbp", "Phf8", "Kmt2d", "Kat6b", "Asxl3", "Prdm8", "Kmt2c", "Kmt2b", "Kat6a", "Ehmt1", "Dnmt3b", "Chd8", "Chd5", "Atrx",
  "Ash1l", "Taf1", "Setd5", "Rai1", "Mecp2", "Bptf","Dnmt1", "Tet1", "Chd1","Chd4","Smarca1","Kmt2a", "Rere", "Srcap", "Kdm3b" ,"Ctrl")


### Use random initialization (ensure reproducibility by setting a seed)

set.seed(123)

### Run NMF; specify how many clusters in data


rank <- 3 # generate three cluster for promoters and samples in total


### Run NMF by seed argument for reproducibility


res.nmf <- nmf(methylation.values.for.all.genes.filtered, rank, method = "brunet", seed = "random")


### The nmfSeed() function in the R package for Non-negative Matrix Factorization (NMF) is used to set a starting seed for the factorization process. 
### This function provides a controlled initialization of the 𝑊 and 𝐻 matrices, which can be essential for achieving reproducible results when running NMF analyses.

### Extract the W and H matrices


W <- basis(res.nmf) # cluster of promoters
H <- coef(res.nmf) # cluster of samples, to access the coefficient matrix from the NMF result



colnames(W) <- c("cluster1", "cluster2", "cluster3") # change colnames as cluster names
rownames(H) <- c("cluster1", "cluster2", "cluster3") # change rownames as cluster names



write.csv(W, "W.matrix.csv", quote = F, row.names = F)
write.csv(H, "H.matrix.csv", quote = F, row.names = F)


### basis() function which respectively extract and set the matrix of basis components of an NMF model (i.e., the first matrix factor, W = basis matrix).

### coef() function respectively extract and set the coefficient matrix of an NMF model (i.e. the second matrix factor, H = coefficient matrix).

### V matrix is data matrix whil is being decomposed to W and H matrices.

### Gene promoters can be clustered based on W matrix, samples can be clustered based on H matrix.

### After running NMF, you can extract the membership of each sample to the clusters

### The membership of each sample to a particular cluster can be determined by the highest coefficient in the corresponding column of the coefficient matrix from the NMF result.

### Determine the cluster membership for each sample by finding the index of the max value in each column


cluster_of_samples <- apply(H, 2, which.max)


### Counts of samples in each cluster


table(cluster_of_samples)


### See which cluster each sample belongs to


which(cluster_of_samples == "1")
which(cluster_of_samples == "2")
which(cluster_of_samples == "3")


### visualize samples as dendogram plot


hc <- hclust(dist(t(H)), method = "ward.D2")


### Create dendogram plot


pdf("dendogram.plot.pdf")

fviz_dend(hc, k = 4, # Cut in four groups
          cex = 1, # label size
          k_colors = c("#1b9e77", "#d95f02", "#7570b3", "#a6761d"),
          color_labels_by_k = TRUE, # color labels by groups
          ggtheme = theme_classic() # Change theme
          )
dev.off()


### Draw a coefficient heatmap of NMF (H matrix) to see membership of each sample to each cluster


pdf("heatmap.of.H.matrix.pdf")

pheatmap(H, show_rownames=T, color=colorRampPalette(c("navy", "gray90", "red"))(50), cellheight= 15, cellwidth = 10)

dev.off()


### Find cluster membership of promoters

clusters.of.promoters <- apply(W, 1, which.max)


### See distribution of promoters in clusters

table(clusters.of.promoters)

which(clusters.of.promoters == "1")
which(clusters.of.promoters == "2")
which(clusters.of.promoters == "3")


### visualize groups as heatmap plot

pdf("heatmap.of.W.matrix.pdf")
pheatmap(W, show_rownames=F, color=colorRampPalette(c("navy", "gray90", "red"))(50))
dev.off()


### Subset clusters for promoters

first.cluster <- clusters.of.promoters[clusters.of.promoters == "1"]
second.cluster <- clusters.of.promoters[clusters.of.promoters == "2"]
third.cluster <- clusters.of.promoters[clusters.of.promoters == "3"]

### Get gene symbols of genes from NCBI genes

NCBI.genes <- fread("gene_info.gz")

first.cluster.gene.symbols <- subset(NCBI.genes, GeneID %in% names(first.cluster))
second.cluster.gene.symbols <- subset(NCBI.genes, GeneID %in% names(second.cluster))
third.cluster.gene.symbols <- subset(NCBI.genes, GeneID %in% names(third.cluster))