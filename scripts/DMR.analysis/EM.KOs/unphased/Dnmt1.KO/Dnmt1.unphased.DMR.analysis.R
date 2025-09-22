# Dnmt1 unphased through dmrseq

# Upload BSseq obj

library(locfit)
library(dmrseq)
library(BiocParallel)
library(rtracklayer)

Dnmt1.BSseq <- readRDS("merged.BSseq.object.without.any.filtering.unsmoothed.rds")

# Remove chrX and chrY

Dnmt1.BSseq <- subset(Dnmt1.BSseq,
                          !seqnames %in% c("chrX", "chrY"))

# Coverage filtering

cov <- getCoverage(Dnmt1.BSseq)
keepLoci.ex <- which(rowSums(cov[, Dnmt1.BSseq$Groups == "Dnmt1"] > 0) >= 3 &
                       rowSums(cov[, Dnmt1.BSseq$Groups == "Control"] > 0) >= 3)

Dnmt1.BSseq <- Dnmt1.BSseq[keepLoci.ex, ]

length(keepLoci.ex)

# Assign colors to each group

Dnmt1.BSseq$col <- c("#d95f02", "#d95f02", "#d95f02", "#7570b3", "#7570b3", "#7570b3")

# DMR analysis

register(BPPARAM = MulticoreParam(workers = 40))

regions <- dmrseq(bs = Dnmt1.BSseq, testCovariate = "Groups", maxGap = 300, BPPARAM = MulticoreParam(workers = 40), minNumRegion = 3)

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
    BSseq = Dnmt1.BSseq, 
    regions = genome.wide.significant.DMRs.first.1000[i],  # Plot one region at a time
    testCovariate = "Groups",
    extend = 5000
  )              
}

dev.off()


pdf("genome.wide.insignificant.DMRs.first.1000.pdf")

for (i in seq_along(genome.wide.insignificant.DMRs.first.1000)) {
  plotDMRs(
    BSseq = Dnmt1.BSseq, 
    regions = genome.wide.insignificant.DMRs.first.1000[i],  # Plot one region at a time
    testCovariate = "Groups",
    extend = 5000
  )              
}

dev.off()