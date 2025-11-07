# Are CpGs in DMRs enriched in exons of canonical transcripts?

library(biomaRt)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(bsseq)
library(org.Mm.eg.db)
library(data.table)

# Get latest GRCm38 data from ENSEMBL

mart <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", version = 102)

# 

exon_data <- getBM(
  attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name",
                 "chromosome_name", "exon_chrom_start", "exon_chrom_end",
                 "strand", "rank"),
  mart = mart
)

# Keep only conventional chromosomes:

exon_data_2 <- subset(exon_data, exon_data$chromosome_name %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9",
                                                                  "10", "11", "12", "13", "14", "15", "16", "17",
                                                                  "18", "19","X", "Y", "MT"))
# Update chromosome names:

exon_data_2$chromosome_name <- ifelse(
  exon_data_2$chromosome_name == "MT",
  "chrM",
  paste0("chr", exon_data_2$chromosome_name)
)

# To get canonical transcripts from GRCm39, remove version = 102

mart.2 <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")

BM.info <- getBM(attributes=c("ensembl_gene_id", "ensembl_transcript_id",
                              "external_gene_name", "transcript_is_canonical",
                              "gene_biotype"), mart = mart.2)

canonical.transcripts.of.protein.coding.genes <- subset(BM.info, BM.info$transcript_is_canonical == "1" &
                                                          BM.info$gene_biotype == "protein_coding")

# promoters of canonical transcripts of protein-coding genes

canonical.exons <- subset(exon_data_2, exon_data_2$ensembl_transcript_id %in% 
                                canonical.transcripts.of.protein.coding.genes$ensembl_transcript_id)

# Recreate exon object as GRanges

canonical.exons.2 <- GRanges(canonical.exons$chromosome_name,
                             IRanges(start = canonical.exons$exon_chrom_start,
                             end = canonical.exons$exon_chrom_end),
                             ensembl_transcript_id = canonical.exons$ensembl_transcript_id,
                             ensembl_gene_id = canonical.exons$ensembl_gene_id,
                             external_gene_name = canonical.exons$external_gene_name,
                             strand = canonical.exons$strand,
                             nth.exon = canonical.exons$rank)

# Save obj

saveRDS(canonical.exons.2, "GRCm38.ENSEMBL.mouse.exon.coordinates.of.canonical.transcripts.rds")

# Merge overlapping exons

canonical.exons.merged <- reduce(canonical.exons.2)

# Bsseq obj (CpG coordinates in this BSseq object should be converted to B6 coordinates from intermediate coordinates)

BSseq.obj <- readRDS("BSseq.obj.rds")

print(BSseq.obj) # BSseq object
length(BSseq.obj) # Total number of CpGs

# CpGs were taken from BSseq.obj above to convert intermediate coordinates to B6 coordinates

full.set.of.CpGs <- fread("all.CpGs.in.B6.coordinates.bed")
full.set.of.CpGs <- GRanges(full.set.of.CpGs$V1, IRanges(start = full.set.of.CpGs$V2, end = full.set.of.CpGs$V3))
length(full.set.of.CpGs)

# Get genome-wide significant DMRs

DMRs <- fread("genome.wide.significant.DMRs.in.B6.coordinates.bed")
DMRs <- GRanges(DMRs$V1, IRanges(start = DMRs$V2, end = DMRs$V3))
length(DMRs) # Number of genome-wide significant DMRs

# Find overlapping CpGs in exons and DMRs

overlapping.CpGs.in.exons <- findOverlaps(full.set.of.CpGs, canonical.exons.merged)
overlapping.CpGs.in.sig.DMRs <- findOverlaps(full.set.of.CpGs, DMRs)
BSseq.overlapping.CpGs.in.exons <- full.set.of.CpGs[queryHits(overlapping.CpGs.in.exons), ]
BSseq.overlapping.CpGs.in.DMRs <- full.set.of.CpGs[queryHits(overlapping.CpGs.in.sig.DMRs), ]

# Calculate number of CpGs for Fisher's exact test

number.of.CpGs.in.both <- length(findOverlaps(BSseq.overlapping.CpGs.in.DMRs, BSseq.overlapping.CpGs.in.exons)) # number of CpGs in both significant DMRs and exons
number_of_CpGs.in.DMRs.not.in.exons <- length(setdiff(BSseq.overlapping.CpGs.in.DMRs, BSseq.overlapping.CpGs.in.exons)) # number of CpGs in significant DMRs but not in exons
number_of_CpGs_in.exons.not_in.DMRs <- length(setdiff(BSseq.overlapping.CpGs.in.exons, BSseq.overlapping.CpGs.in.DMRs)) # number of CpGs in exons but not in significant DMRs
number_of_CpGs_in_neither_exons_nor_DMRs <- length(full.set.of.CpGs) - (number.of.CpGs.in.both + number_of_CpGs.in.DMRs.not.in.exons + number_of_CpGs_in.exons.not_in.DMRs) # Number of CpGs that is in neither DMRs nor exons:

# Fisher exact test to see whether CpGs in DMrs are enriched in gene exons

contingency_table <- data.frame(
  "CpGs.in.exons" = c(number.of.CpGs.in.both, number_of_CpGs_in.exons.not_in.DMRs),
  "CpGs.in.outside.of.exons" = c(number_of_CpGs.in.DMRs.not.in.exons, number_of_CpGs_in_neither_exons_nor_DMRs),
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