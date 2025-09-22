# Kmt2a unphased through dmrseq

# Upload BSseq obj

library(locfit)
library(dmrseq)
library(BiocParallel)
library(rtracklayer)

Kmt2a.BSseq <- readRDS("merged.BSseq.object.without.any.filtering.unsmoothed.rds")

# Remove chrX and chrY

Kmt2a.BSseq <- subset(Kmt2a.BSseq,
                          !seqnames %in% c("chrX", "chrY"))

# Coverage filtering

cov <- getCoverage(Kmt2a.BSseq)
keepLoci.ex <- which(rowSums(cov[, Kmt2a.BSseq$Groups == "Kmt2a"] > 0) >= 2 &
                       rowSums(cov[, Kmt2a.BSseq$Groups == "Control"] > 0) >= 4)

Kmt2a.BSseq <- Kmt2a.BSseq[keepLoci.ex, ]

length(keepLoci.ex)

# Assign colors to each group

Kmt2a.BSseq$col <- c("#d95f02", "#d95f02", "#d95f02", "#d95f02", "#7570b3", "#7570b3")

# DMR analysis

register(BPPARAM = MulticoreParam(workers = 40))

regions <- dmrseq(bs = Kmt2a.BSseq, testCovariate = "Groups", maxGap = 300, BPPARAM = MulticoreParam(workers = 40), minNumRegion = 3)

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
    BSseq = Kmt2a.BSseq, 
    regions = genome.wide.significant.DMRs.first.1000[i],  # Plot one region at a time
    testCovariate = "Groups",
    extend = 5000
  )              
}

dev.off()

pdf("genome.wide.insignificant.DMRs.first.1000.pdf")

for (i in seq_along(genome.wide.insignificant.DMRs.first.1000)) {
  plotDMRs(
    BSseq = Kmt2a.BSseq,
    regions = genome.wide.insignificant.DMRs.first.1000[i],  # Plot one region at a time
    testCovariate = "Groups",
    extend = 5000
  )              
}

dev.off()