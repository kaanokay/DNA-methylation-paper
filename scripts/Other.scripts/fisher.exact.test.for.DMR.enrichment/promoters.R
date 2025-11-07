# Are CpGs in DMRs enriched in promoters of canonical transcripts?

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

# Get promoter regions of all transcripts of mouse genes

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Get promoter regions
# upstream = 2000 bp, downstream = 500 bp (adjust as needed)
promoters_mm10 <- promoters(txdb, upstream = 2000, downstream = 500)

# remove version extensions from transcript IDs

names(promoters_mm10) <- sub("\\..*", "", names(promoters_mm10))

# promoters of canonical transcripts of protein-coding genes

canonical.promoters <- subset(promoters_mm10, names(promoters_mm10) %in% 
                                canonical.transcripts.of.protein.coding.genes$ensembl_transcript_id)


# Keep canonical chromosomes                                

canonical.promoters <- subset(canonical.promoters, seqnames(canonical.promoters)
                                                   %in% c("chr1", "chr2", "chr3", "chr4",
                                                           "chr5", "chr6", "chr7", "chr8",
                                                           "chr9", "chr10", "chr11", "chr12",
                                                           "chr13","chr14", "chr15", "chr16",
                                                           "chr17","chr18", "chr19","chrX",
                                                           "chrY", "chrM"))

# Create lookup vector from transcript ID to gene name
tx2gene <- setNames(canonical.transcripts.of.protein.coding.genes$external_gene_name,
                    canonical.transcripts.of.protein.coding.genes$ensembl_transcript_id)

# Add gene name to mcols
mcols(canonical.promoters)$gene_name <- tx2gene[names(canonical.promoters)]

# Save object

saveRDS(canonical.promoters, "promoters.of.canonical.transcripts.in.Musmusculus.mm10.coordinates.rds")

# Merge overlapping promoters

canonical.promoters.merged <- reduce(canonical.promoters)

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

# Find overlapping CpGs in promoters and DMRs

overlapping.CpGs.in.promoters <- findOverlaps(full.set.of.CpGs, canonical.promoters.merged)
overlapping.CpGs.in.sig.DMRs <- findOverlaps(full.set.of.CpGs, DMRs)
BSseq.overlapping.CpGs.in.promoters <- full.set.of.CpGs[queryHits(overlapping.CpGs.in.promoters), ]
BSseq.overlapping.CpGs.in.DMRs <- full.set.of.CpGs[queryHits(overlapping.CpGs.in.sig.DMRs), ]

# Calculate number of CpGs for Fisher's exact test

number.of.CpGs.in.both <- length(findOverlaps(BSseq.overlapping.CpGs.in.DMRs, BSseq.overlapping.CpGs.in.promoters)) # number of CpGs in both significant DMRs and promoters
number_of_CpGs.in.DMRs.not.in.promoters <- length(setdiff(BSseq.overlapping.CpGs.in.DMRs, BSseq.overlapping.CpGs.in.promoters)) # number of CpGs in significant DMRs but not in promoters
number_of_CpGs_in.promoters.not_in.DMRs <- length(setdiff(BSseq.overlapping.CpGs.in.promoters, BSseq.overlapping.CpGs.in.DMRs)) # number of CpGs in promoters but not in significant DMRs
number_of_CpGs_in_neither_promoters_nor_DMRs <- length(full.set.of.CpGs) - (number.of.CpGs.in.both + number_of_CpGs.in.DMRs.not.in.promoters + number_of_CpGs_in.promoters.not_in.DMRs) # Number of CpGs that is in neither DMRs nor promoters:

# Fisher exact test to see whether CpGs in DMrs are enriched in gene promoters

contingency_table <- data.frame(
  "CpGs.in.promoters" = c(number.of.CpGs.in.both, number_of_CpGs_in.promoters.not_in.DMRs),
  "CpGs.in.outside.of.promoters" = c(number_of_CpGs.in.DMRs.not.in.promoters, number_of_CpGs_in_neither_promoters_nor_DMRs),
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