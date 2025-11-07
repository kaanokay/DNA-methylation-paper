library(biomaRt)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(bsseq)
library(data.table)

# Get canonical transcripts of protein coding mouse genes

mart <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")

BM.info <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id",
                              "external_gene_name", "transcript_is_canonical",
                              "gene_biotype"), mart = mart)

canonical.transcripts.of.protein.coding.genes <- subset(BM.info, BM.info$transcript_is_canonical == "1" &
                                                          BM.info$gene_biotype == "protein_coding")

# UCSC mm10

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Finding intergenic regions
# Source: https://stackoverflow.com/questions/29253412/finding-intergenic-regions

# Get transcripts from UCSC with ENSEMBL ids

transcipts.from.UCSC <- transcripts(txdb)

# Clean transcript ids

tx_names_clean <- lapply(transcipts.from.UCSC$tx_name, function(x) sub("\\..*", "", x))

transcipts.from.UCSC$tx_name_clean <- tx_names_clean

# subset intronic parts to canonical and protein-coding genes

is_canonical <- vapply(
  transcipts.from.UCSC$tx_name_clean,
  function(x) any(x %in% canonical.transcripts.of.protein.coding.genes$ensembl_transcript_id),
  logical(1)
)

# Now subset

canonical.and.protein.coding.transcripts <- transcipts.from.UCSC[is_canonical]

# Keep only conventional chromosomes:

canonical.and.protein.coding.transcripts <- subset(canonical.and.protein.coding.transcripts, seqnames(canonical.and.protein.coding.transcripts)
                                                   %in% c("chr1", "chr2", "chr3", "chr4",
                                                           "chr5", "chr6", "chr7", "chr8",
                                                           "chr9", "chr10", "chr11", "chr12",
                                                           "chr13","chr14", "chr15", "chr16",
                                                           "chr17","chr18", "chr19","chrX",
                                                           "chrY", "chrM"))


# Merge overlapping regions

canonical.and.protein.coding.transcripts.2 <- reduce(canonical.and.protein.coding.transcripts, ignore.strand = TRUE)
intergenic <- gaps(canonical.and.protein.coding.transcripts.2)
intergenic.2 <- intergenic[strand(intergenic) == "*"]

# Save object

saveRDS(intergenic.2, "intergenic.regions.in.Musmusculus.mm10.coordinates.rds")

# Bsseq obj (CpG coordinates in this BSseq object should be converted to B6 coordinates from intermediate coordinates)

BSseq.obj <- readRDS("BSseq.obj.rds")

print(BSseq.obj) # BSseq object
length(BSseq.obj) # Total number of CpGs

# CpGs were taken from BSseq.obj above to convert intermediate coordinates to B6 coordinates

full.set.of.CpGs <- fread("all.CpGs.from.haplotype.specific.analysis.in.B6.coordinates.bed")
full.set.of.CpGs <- GRanges(full.set.of.CpGs$V1, IRanges(start = full.set.of.CpGs$V2, end = full.set.of.CpGs$V3))
length(full.set.of.CpGs)

# Get genome-wide significant DMRs

DMRs <- fread("genome.wide.significant.DMRs.in.B6.coordinates.bed")
DMRs <- GRanges(DMRs$V1, IRanges(start = DMRs$V2, end = DMRs$V3))
length(DMRs) # Number of genome-wide significant DMRs

# Find overlapping CpGs in intergenic and DMRs

overlapping.CpGs.in.intergenic <- findOverlaps(full.set.of.CpGs, intergenic.2)
overlapping.CpGs.in.sig.DMRs <- findOverlaps(full.set.of.CpGs, DMRs)
BSseq.overlapping.CpGs.in.intergenic <- full.set.of.CpGs[queryHits(overlapping.CpGs.in.intergenic), ]
BSseq.overlapping.CpGs.in.DMRs <- full.set.of.CpGs[queryHits(overlapping.CpGs.in.sig.DMRs), ]

# Calculate number of CpGs for Fisher's exact test

number.of.CpGs.in.both <- length(findOverlaps(BSseq.overlapping.CpGs.in.DMRs, BSseq.overlapping.CpGs.in.intergenic)) # number of CpGs in both significant DMRs and intergenic
number_of_CpGs.in.DMRs.not.in.intergenic <- length(setdiff(BSseq.overlapping.CpGs.in.DMRs, BSseq.overlapping.CpGs.in.intergenic)) # number of CpGs in significant DMRs but not in intergenic
number_of_CpGs_in.intergenic.not_in.DMRs <- length(setdiff(BSseq.overlapping.CpGs.in.intergenic, BSseq.overlapping.CpGs.in.DMRs)) # number of CpGs in intergenic but not in significant DMRs
number_of_CpGs_in_neither_intergenic_nor_DMRs <- length(full.set.of.CpGs) - (number.of.CpGs.in.both + number_of_CpGs.in.DMRs.not.in.intergenic + number_of_CpGs_in.intergenic.not_in.DMRs) # Number of CpGs that is in neither DMRs nor intergenic:

# Fisher exact test to see whether CpGs in DMrs are enriched in gene intergenic

contingency_table <- data.frame(
  "CpGs.in.intergenic" = c(number.of.CpGs.in.both, number_of_CpGs_in.intergenic.not_in.DMRs),
  "CpGs.in.outside.of.intergenic" = c(number_of_CpGs.in.DMRs.not.in.intergenic, number_of_CpGs_in_neither_intergenic_nor_DMRs),
  row.names = c("CpGs.in.DMRs", "CpGs.in.outside.of.DMRs"),
  stringsAsFactors = FALSE
)

print(contingency_table)

contingency_table.2 <- fisher.test(contingency_table)

print(contingency_table.2)

result_df <- data.frame(
  p_value = contingency_table.2$p.value,
  odds_ratio = contingency_table.2$estimate,
  conf_low = contingency_table.2$conf.int[1],
  conf_high = contingency_table.2$conf.int[2]
)

write.table(result_df, "fisher.exact.test.results.txt", quote = F, row.names = F)