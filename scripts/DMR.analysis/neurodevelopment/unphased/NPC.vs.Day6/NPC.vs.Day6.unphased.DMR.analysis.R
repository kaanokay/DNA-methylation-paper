# NPC vs Day6 unphased through dmrseq

# Upload BSseq obj

library(locfit)
library(dmrseq)
library(BiocParallel)
library(rtracklayer)

neuro.BSseq <- readRDS("merged_BSseq_unsmoothed.rds")

# Keep NPCs and Day 6 neurons

neuro.BSseq <- neuro.BSseq[, !grepl("D2", pData(neuro.BSseq)$Samples)]

# Rename samples

neuro.BSseq$Samples <- c("2.Day6", "1.NPC", "2.Day6", "1.NPC", "2.Day6", "1.NPC")

# Remove chrX and chrY

neuro.BSseq <- subset(neuro.BSseq,
                          !seqnames %in% c("chrX", "chrY"))

# Coverage filtering

cov <- getCoverage(neuro.BSseq)
keepLoci.ex <- which(rowSums(cov[, neuro.BSseq$Samples == "1.NPC"] > 0) >= 3 &
                       rowSums(cov[, neuro.BSseq$Samples == "2.Day6"] > 0) >= 3)

neuro.BSseq <- neuro.BSseq[keepLoci.ex, ]

length(keepLoci.ex)

# Assign colors to each group

neuro.BSseq$col <- c("#7570b3", "#d95f02", "#7570b3", "#d95f02", "#7570b3", "#d95f02")

# DMR analysis

register(BPPARAM = MulticoreParam(workers = 40))

regions <- dmrseq(bs = neuro.BSseq, testCovariate = "Samples", maxGap = 300, BPPARAM = MulticoreParam(workers = 40), minNumRegion = 3)

# get genome-wide significant and insignificant DMRs

genome.wide.significant.DMRs <- regions[which(regions$qval < 0.05)]

genome.wide.insignificant.DMRs <- regions[which(regions$qval > 0.05)]

# Plot DMRs

genome.wide.significant.DMRs.first.1000 <- genome.wide.significant.DMRs[1:1000]

genome.wide.insignificant.DMRs.first.1000 <- genome.wide.insignificant.DMRs[1:1000]

# To plot many regions in dmrseq pipeline:

pdf("genome.wide.significant.DMRs.first.1000.pdf")

for (i in seq_along(genome.wide.significant.DMRs.first.1000)) {
  plotDMRs(
    BSseq = neuro.BSseq, 
    regions = genome.wide.significant.DMRs.first.1000[i],  # Plot one region at a time
    testCovariate = "Samples",
    extend = 5000
  )              
}

dev.off()

pdf("genome.wide.insignificant.DMRs.first.1000.pdf")

for (i in seq_along(genome.wide.insignificant.DMRs.first.1000)) {
  plotDMRs(
    BSseq = neuro.BSseq,
    regions = genome.wide.insignificant.DMRs.first.1000[i],  # Plot one region at a time
    testCovariate = "Samples",
    extend = 5000
  )              
}

dev.off()