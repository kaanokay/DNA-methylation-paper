library(bsseq)
library(BiocParallel)

BSseq.obj <- readRDS("BSseq.object.unsmoothed.and.unphased.rds")

BSseq.obj.smoothed <- BSmooth(
    BSseq = BSseq.obj, 
    BPPARAM = MulticoreParam(workers = 40), 
    verbose = TRUE)

saveRDS(BSseq.obj.smoothed, "BSseq.object.smoothed.and.unphased.rds")