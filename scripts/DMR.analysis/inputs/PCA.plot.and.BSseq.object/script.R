# Loading libraries

library(bsseq)
library(ggforce)
library(BiocParallel)
library(parallel)
library(data.table)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)

# Read bedMethly files from modkit

EM_KO <- fread("EM_KO.bed", nThread = 32)

colnames(EM_KO) <- c("chrom", "start", "end", "modified_base_code", "score", "strand", "tstart", "tend", "color", "Nvalid_cov", "Nmod_/_Nvalid_cov", "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall")

EM_KO <- subset(EM_KO, EM_KO$modified_base_code == "m")

# Generate BSSeq object

EM_KO.BSSeq.object <- BSseq(chr = as.character(EM_KO$chrom), pos = EM_KO$start +1, M = matrix(EM_KO$Nmod), Cov = matrix(EM_KO$Nmod + EM_KO$Ncanonical), sampleNames = c("EM_KO"))

sampleNames(EM_KO.BSSeq.object) <- "EM_KO"

#

control <- fread("control.bed", nThread = 32)

colnames(control) <- c("chrom", "start", "end", "modified_base_code", "score", "strand", "tstart", "tend", "color", "Nvalid_cov", "Nmod_/_Nvalid_cov", "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall")

control <- subset(control, control$modified_base_code == "m")

# Generate BSSeq object

control.BSSeq.object <- BSseq(chr = as.character(control$chrom), pos = control$start +1, M = matrix(control$Nmod), Cov = matrix(control$Nmod + control$Ncanonical), sampleNames = c("control"))

sampleNames(control.BSSeq.object) <- "control"

# Merge BSseq objects

merged_BSseq <- combine(
EM_KO.BSSeq.object,
control.BSSeq.object)

merged_BSseq <- collapseBSseq(merged_BSseq, group = c("EM_KO", "Control"))

nrow(merged_BSseq)

Samples <- as.character(c("EM_KO", "Control"))

# sample description with metadata

pData(merged_BSseq)$Samples <- Samples
pData(merged_BSseq)

# Group description

Groups <- as.character(c("EM_KO", "Control"))

pData(merged_BSseq)$Groups <- Groups
pData(merged_BSseq)

# Save merged BSseq object

saveRDS(merged_BSseq, "BSseq.object.unsmoothed.and.unphased.rds")

# --- Filtering number of CpGs with 0 coverage in all samples start ---

merged_BSseq.cov <- getCoverage(merged_BSseq)
keep <- which(rowSums(merged_BSseq.cov) !=0)
merged_BSseq_filtered <- merged_BSseq[keep,]
nrow(merged_BSseq_filtered)

# Write a custom function for PCA plot

PCA <- function(
  BSseq.obj,
  genome,
  tilewidth,
  cut.last.tile.in.chrom,
  CpG,
  nudge_x,
  nudge_y,
  size
)

{
  
  # Tile the UCSC mm10 mouse genome
  
  mm10.mouse.genome.1kb.tiles <- GenomicRanges::tileGenome(seqinfo(genome), tilewidth = tilewidth, cut.last.tile.in.chrom = cut.last.tile.in.chrom)
  
  # Filtering number of CpGs with 0 coverage across all samples
  
  coverage <- getCoverage(BSseq.obj)
  keep <- which(rowSums(coverage) !=0)
  BSseq.obj.filtered <- BSseq.obj[keep,]

  # Find which 1kb tile has more than a number of CpGs (default is 3)
  
  one.kb.tiled.genome.more.number.of.CpGs <- mm10.mouse.genome.1kb.tiles[which(GenomicRanges::countOverlaps(mm10.mouse.genome.1kb.tiles, BSseq.obj.filtered) >= CpG), ]

  # Find which CpG in BSseq.object overlap with 1kb tiled coordinates with more than number of CpGs and then keep those CpGs in the BSseq.obj
  
  BSseq.obj.CpGs.within.1kb.tilled.genome <- BSseq.obj.filtered[as.data.frame(GenomicRanges::findOverlaps(BSseq.obj.filtered, one.kb.tiled.genome.more.number.of.CpGs))$queryHits, ]
  
  # Get raw methylation values
  
  BSseq.obj.methylation <- getMeth(BSseq.obj.CpGs.within.1kb.tilled.genome, type = "raw")
  
  # Filter tiles with NA methylation values
  
  BSseq.obj.methylation.filtered <- BSseq.obj.methylation[rowSums(is.na(BSseq.obj.methylation)) == 0, ]
  
  # Perform PCA
  pca_data <- prcomp(t(BSseq.obj.methylation.filtered))
  pca_data_perc <- round(100 * pca_data$sdev^2 / sum(pca_data$sdev^2), 1)
  df_pca_data <- data.frame(
    PC1 = pca_data$x[, 1],
    PC2 = pca_data$x[, 2])
    
    # Write down PCA coordinates

  write.table(df_pca_data, "PCA_coordinates.txt")
  
  pdf("PCA.plot.pdf")  
  
  # Plot PCA
  p <- ggplot(df_pca_data, aes(PC1, PC2, label = row.names(df_pca_data))) +
    geom_point() +
    labs(
      x = paste0("PC1 (", pca_data_perc[1], ")"),
      y = paste0("PC2 (", pca_data_perc[2], ")")
    ) +
    geom_text(nudge_x = nudge_x, nudge_y = nudge_y, size = size) +
    theme_classic()
  
  print(p)
  
  dev.off()
}

# Draw PCA plot

PCA(BSseq.obj = merged_BSseq_filtered, genome = Mmusculus, tilewidth=1000, cut.last.tile.in.chrom= TRUE, CpG = 3, nudge_x = 2.5, nudge_y = 2.5, size = 3)