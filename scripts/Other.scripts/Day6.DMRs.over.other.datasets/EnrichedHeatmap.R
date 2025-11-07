# Visualization unphased NPC vs Day6 DMRs over unphased Dnmt1 samples
# source: https://github.com/hansenlab/BrainEpigenomeNN/blob/7bd8f1ecadb2a7bddf4882c4d80bc125482e41df/nonCG/analysis/EDA_heatmaps.R#L156

library(viridis)
library(rtracklayer)
library(bsseq)
library(EnrichedHeatmap)
library(data.table)
library(Cairo)
library(matrixStats)
library(circlize)
library(scales)
library(BiocParallel)
library(parallel)
library(magick)

# Set a working directory

# Color for heatmap

# col <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

col <- viridis(100)

# Get genome-wide significant Day 6 DMRs

NPC.vs.Day6.unphased.DMRs <- readRDS("all.DMRs.rds")
NPC.vs.Day6.unphased.DMRs <- as.data.frame(NPC.vs.Day6.unphased.DMRs)
NPC.vs.Day6.unphased.DMRs.significant <- NPC.vs.Day6.unphased.DMRs[NPC.vs.Day6.unphased.DMRs$qval < 0.05 & (NPC.vs.Day6.unphased.DMRs$width/NPC.vs.Day6.unphased.DMRs$L) < 120 & NPC.vs.Day6.unphased.DMRs$L > 5, ]
NPC.vs.Day6.unphased.DMRs.significant <- GRanges(NPC.vs.Day6.unphased.DMRs.significant$seqnames, IRanges(start =NPC.vs.Day6.unphased.DMRs.significant$start, end =NPC.vs.Day6.unphased.DMRs.significant$end))
length(NPC.vs.Day6.unphased.DMRs.significant) # Number of genome-wide significant DMRs

NPC.vs.Day6.DMRs <- sort(NPC.vs.Day6.unphased.DMRs.significant)

# Dnmt1 data

Dnmt1.BSseq <- readRDS("BSseq.object.rds")

# neurodevelopment data

neurodevelopment <- readRDS("BSseq.obj.rds")

# Kmt2a data

Kmt2a <- readRDS("BSseq.object.rds")

# Combine BSseq objects

combined.BSseq <- combine(Dnmt1.BSseq, neurodevelopment, Kmt2a)

# Keep CpGs with > 0 coverage in all samples to keep same set of CpGs for each sample

cov_matrix <- getCoverage(combined.BSseq, type = "Cov")
keep_CpGs <- rowSums(cov_matrix > 0) == ncol(cov_matrix)
length(keep_CpGs)
filtered_BSseq <- combined.BSseq[keep_CpGs, ]

filtered_BSseq <- BSmooth(
    BSseq = filtered_BSseq,
    BPPARAM = MulticoreParam(workers = 10),
    verbose = TRUE)

# Separate samples

Dnmt1.controls <- filtered_BSseq[, c(1:3)]
Dnmt1.KOs <- filtered_BSseq[, c(4:6)]
mNPCs <- filtered_BSseq[, c(9, 12, 15)]
Day2.neurons <- filtered_BSseq[, c(7,10,13)]
Day6.neurons <- filtered_BSseq[, c(8,11,14)]
Kmt2a.controls <- filtered_BSseq[, c(16:19)]
Kmt2a.KOs <- filtered_BSseq[, c(20:21)]

# Calculate average methylation

avg_methylation <- rowMeans(getMeth(Dnmt1.controls, type="smooth", what="perBase"), na.rm = TRUE)
avg_methylation.2 <- rowMeans(getMeth(Dnmt1.KOs, type="smooth", what="perBase"), na.rm = TRUE)
avg_methylation.3 <- rowMeans(getMeth(mNPCs, type="smooth", what="perBase"), na.rm = TRUE)
avg_methylation.4 <- rowMeans(getMeth(Day2.neurons, type="smooth", what="perBase"), na.rm = TRUE)
avg_methylation.5 <- rowMeans(getMeth(Day6.neurons, type="smooth", what="perBase"), na.rm = TRUE)
avg_methylation.6 <- rowMeans(getMeth(Kmt2a.controls, type="smooth", what="perBase"), na.rm = TRUE)
avg_methylation.7 <- rowMeans(getMeth(Kmt2a.KOs, type="smooth", what="perBase"), na.rm = TRUE)

# Create GRanges object

Dnmt1_control_GRanges <- granges(Dnmt1.controls)
Dnmt1_GRanges <- granges(Dnmt1.KOs)
mNPCs_GRanges <- granges(mNPCs)
Day2.neurons_GRanges <- granges(Day2.neurons)
Day6.neurons_GRanges <- granges(Day6.neurons)
Kmt2a.controls_GRanges <- granges(Kmt2a.controls)
Kmt2a.KOs_GRanges <- granges(Kmt2a.KOs)

# Add average methylation to each CpG

Dnmt1_control_GRanges$avg_methylation <- avg_methylation
Dnmt1_GRanges$avg_methylation.2 <- avg_methylation.2
mNPCs_GRanges$avg_methylation.3 <- avg_methylation.3
Day2.neurons_GRanges$avg_methylation.4 <- avg_methylation.4
Day6.neurons_GRanges$avg_methylation.5 <- avg_methylation.5
Kmt2a.controls_GRanges$avg_methylation.6 <- avg_methylation.6
Kmt2a.KOs_GRanges$avg_methylation.7 <- avg_methylation.7

# Create methylation matrix

Dnmt1_controls_avg_methylation_matrix <- normalizeToMatrix(
  Dnmt1_control_GRanges, NPC.vs.Day6.DMRs,
  value_column = "avg_methylation", 
  mean_mode = "absolute", 
  w = 100, 
  background = NA,
  smooth = TRUE,
  extend = 5000,
  target_ratio = 0.35
)

Dnmt1_avg_methylation_matrix <- normalizeToMatrix(
  Dnmt1_GRanges, NPC.vs.Day6.DMRs,
  value_column = "avg_methylation.2", 
  mean_mode = "absolute", 
  w = 100, 
  background = NA,
  smooth = TRUE,
  extend = 5000,
  target_ratio = 0.35
)

mNPCs_avg_methylation_matrix <- normalizeToMatrix(
  mNPCs_GRanges, NPC.vs.Day6.DMRs,
  value_column = "avg_methylation.3", 
  mean_mode = "absolute", 
  w = 100, 
  background = NA,
  smooth = TRUE,
  extend = 5000,
  target_ratio = 0.35
)

Day2.neurons_avg_methylation_matrix <- normalizeToMatrix(
  Day2.neurons_GRanges, NPC.vs.Day6.DMRs,
  value_column = "avg_methylation.4",
  mean_mode = "absolute",
  w = 100, 
  background = NA,
  smooth = TRUE,
  extend = 5000,
  target_ratio = 0.35
)

Day6.neurons_avg_methylation_matrix <- normalizeToMatrix(
  Day6.neurons_GRanges, NPC.vs.Day6.DMRs,
  value_column = "avg_methylation.5",
  mean_mode = "absolute",
  w = 100, 
  background = NA,
  smooth = TRUE,
  extend = 5000,
  target_ratio = 0.35
)

Kmt2a_controls_avg_methylation_matrix <- normalizeToMatrix(
  Kmt2a.controls_GRanges, NPC.vs.Day6.DMRs,
  value_column = "avg_methylation.6",
  mean_mode = "absolute",
  w = 100, 
  background = NA,
  smooth = TRUE,
  extend = 5000,
  target_ratio = 0.35
)

Kmt2a_KOs_avg_methylation_matrix <- normalizeToMatrix(
  Kmt2a.KOs_GRanges, NPC.vs.Day6.DMRs,
  value_column = "avg_methylation.7",
  mean_mode = "absolute",
  w = 100, 
  background = NA,
  smooth = TRUE,
  extend = 5000,
  target_ratio = 0.35
)

# Generate enrich heatmap

Dnmt1.controls.EH <- EnrichedHeatmap(
                mat = Dnmt1_controls_avg_methylation_matrix,
                col = col,
                pos_line = FALSE,
                column_title_gp = gpar(fontsize = 15),
                axis_name_gp = gpar(fontsize = 6),
                column_title = "NPC vs Day 6 DMRs over Dnmt1 controls",
                width = 1,
                axis_name = c("-5000", "start", "end", "5000"),
                use_raster = TRUE,
                raster_device = "CairoPNG", top_annotation = HeatmapAnnotation(
                  enriched = anno_enriched(
                    ylim = c(0, 1),
                    axis_param = list(
                      at = c(0, 0.5, 1),
                      labels = c("0", "0.5", "1"),
                      side = "right",
                      facing = "outside"
                    ))))


Dnmt1.KOs.EH <- EnrichedHeatmap(
                mat = Dnmt1_avg_methylation_matrix,
                col = col,
                pos_line = FALSE,
                column_title_gp = gpar(fontsize = 15),
                axis_name_gp = gpar(fontsize = 6),
                column_title = "NPC vs Day 6 DMRs over Dnmt1 KOs",
                width = 1,
                axis_name = c("-5000", "start", "end", "5000"),
                use_raster = TRUE,
                raster_device = "CairoPNG", top_annotation = HeatmapAnnotation(
                  enriched = anno_enriched(
                    ylim = c(0, 1),
                    axis_param = list(
                      at = c(0, 0.5, 1),
                      labels = c("0", "0.5", "1"),
                      side = "right",
                      facing = "outside"
                    ))))

mNPCs.EH <- EnrichedHeatmap(
                mat = mNPCs_avg_methylation_matrix,
                col = col,
                pos_line = FALSE,
                column_title_gp = gpar(fontsize = 15),
                axis_name_gp = gpar(fontsize = 6),
                column_title = "NPC vs Day 6 DMRs over mNPCs",
                width = 1,
                axis_name = c("-5000", "start", "end", "5000"),
                use_raster = TRUE,
                raster_device = "CairoPNG", top_annotation = HeatmapAnnotation(
                  enriched = anno_enriched(
                    ylim = c(0, 1),
                    axis_param = list(
                      at = c(0, 0.5, 1),
                      labels = c("0", "0.5", "1"),
                      side = "right",
                      facing = "outside"
                    ))))

Day2.neurons.EH <- EnrichedHeatmap(
                mat = Day2.neurons_avg_methylation_matrix,
                col = col,
                pos_line = FALSE,
                column_title_gp = gpar(fontsize = 15),
                axis_name_gp = gpar(fontsize = 6),
                column_title = "NPC vs Day 6 DMRs over Day 2 neurons",
                width = 1,
                axis_name = c("-5000", "start", "end", "5000"),
                use_raster = TRUE,
                raster_device = "CairoPNG", top_annotation = HeatmapAnnotation(
                  enriched = anno_enriched(
                    ylim = c(0, 1),
                    axis_param = list(
                      at = c(0, 0.5, 1),
                      labels = c("0", "0.5", "1"),
                      side = "right",
                      facing = "outside"
                    ))))                                        


Day6.neurons.EH <- EnrichedHeatmap(
                mat = Day6.neurons_avg_methylation_matrix,
                col = col,
                pos_line = FALSE,
                column_title_gp = gpar(fontsize = 15),
                axis_name_gp = gpar(fontsize = 6),
                column_title = "NPC vs Day 6 DMRs over Day 6 neurons",
                width = 1,
                axis_name = c("-5000", "start", "end", "5000"),
                use_raster = TRUE,
                raster_device = "CairoPNG", top_annotation = HeatmapAnnotation(
                  enriched = anno_enriched(
                    ylim = c(0, 1),
                    axis_param = list(
                      at = c(0, 0.5, 1),
                      labels = c("0", "0.5", "1"),
                      side = "right",
                      facing = "outside"
                    ))))

Kmt2a.controls.EH <- EnrichedHeatmap(
                mat = Kmt2a_controls_avg_methylation_matrix,
                col = col,
                pos_line = FALSE,
                column_title_gp = gpar(fontsize = 15),
                axis_name_gp = gpar(fontsize = 6),
                column_title = "NPC vs Day 6 DMRs over Kmt2a controls",
                width = 1,
                axis_name = c("-5000", "start", "end", "5000"),
                use_raster = TRUE,
                raster_device = "CairoPNG", top_annotation = HeatmapAnnotation(
                  enriched = anno_enriched(
                    ylim = c(0, 1),
                    axis_param = list(
                      at = c(0, 0.5, 1),
                      labels = c("0", "0.5", "1"),
                      side = "right",
                      facing = "outside"
                    ))))                    

Kmt2a.KOs.EH <- EnrichedHeatmap(
                mat = Kmt2a_KOs_avg_methylation_matrix,
                col = col,
                pos_line = FALSE,
                column_title_gp = gpar(fontsize = 15),
                axis_name_gp = gpar(fontsize = 6),
                column_title = "NPC vs Day 6 DMRs over Kmt2a KOs",
                width = 1,
                axis_name = c("-5000", "start", "end", "5000"),
                use_raster = TRUE,
                raster_device = "CairoPNG", top_annotation = HeatmapAnnotation(
                  enriched = anno_enriched(
                    ylim = c(0, 1),
                    axis_param = list(
                      at = c(0, 0.5, 1),
                      labels = c("0", "0.5", "1"),
                      side = "right",
                      facing = "outside"
                    ))))

# Merge all data into one pdf

pdf("NPC.vs.Day6.unphased.significant.DMRs.on.other.samples.same.set.of.CpGs.with.viridis.colored.pdf", width = 20)

ht_list <- mNPCs.EH + Day2.neurons.EH + Day6.neurons.EH + Kmt2a.controls.EH + Kmt2a.KOs.EH + Dnmt1.controls.EH + Dnmt1.KOs.EH

draw(ht_list, ht_gap = unit(c(8, 8, 8, 8, 8, 8), "mm")) # draw() function put space between heatmaps

# Close the PDF device
dev.off()