Scripts for figures
================
KatrínM
2025-11-11

# Install packages

``` r
# List all packages you need
packages <- c(
  "tidyverse", "ggplot2", "dplyr", "readr", "tibble", "grid",
  "ggalt", "stringr", "ggrepel", "universalmotif", "VennDiagram", "patchwork",
  "biomaRt", "GenomicFeatures", "AnnotationDbi","TxDb.Mmusculus.UCSC.mm10.knownGene", 
  "pheatmap", "RColorBrewer", "data.table", "colorspace", "viridis"
)


# Install missing ones (CRAN first)
to_install <- packages[!packages %in% installed.packages()[,"Package"]]
if (length(to_install) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(to_install, ask = FALSE, update = FALSE)
}

# Load them all
invisible(lapply(packages, library, character.only = TRUE))
```

## Figure 2D

``` r
#load the data
top30KDnmt1_DMRs_vs_matchedbg <- read.delim("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Juan_RNAseq/ame_top30K_Dnmt1_DMRs_vs_newbackground.tsv", sep="\t",header=TRUE)

top30KDnmt1_DMRs_vs_matchedbg <- top30KDnmt1_DMRs_vs_matchedbg %>%
  mutate(
    cleaned_motif = str_to_title(str_extract(motif_ID, "^[^_]+"))
  )

#Create the enrichment ratio by dividing the %positive / %negative
top30KDnmt1_DMRs_vs_matchedbg <- top30KDnmt1_DMRs_vs_matchedbg %>%
  mutate(enrichment_ratio = X.TP / X.FP)

#Calculate the odds ratios 
top30KDnmt1_DMRs_vs_matchedbg <- top30KDnmt1_DMRs_vs_matchedbg %>%
  mutate(odds_ratio = (X.TP/100) / (1 - X.TP/100) / ((X.FP/100) / (1 - X.FP/100)))

#Identify the top10 motifs to label
DMRtop10_motifs <- top30KDnmt1_DMRs_vs_matchedbg %>%
  arrange(adj_p.value) %>%   # Ascending order = smallest E-values
  dplyr::slice(1:10)

# Add Epas1 from cleaned_motif (filtering by name)
Epas1_row <- top30KDnmt1_DMRs_vs_matchedbg %>%
  filter(cleaned_motif == "Epas1")

# Combine
DMRtop10_motifs <- bind_rows(DMRtop10_motifs, Epas1_row)

#Now plot Figure 2D
ggplot(top30KDnmt1_DMRs_vs_matchedbg, aes(x = log2(odds_ratio), y = -log10(adj_p.value), color = TP)) +
  geom_point(size=4, alpha = 0.8) +  # Scatter points
  scale_color_gradient(low = "gray80", high = "purple3") +
  labs(
    title = "Motif enrichment of top 30K Dnmt1-DMRs",
    color = "#Positive",
    #size = "Overlap (%)",
    x = "Log2(OR)",
    y = "-log10(adj.pval)"
  ) +
  coord_cartesian(xlim = c(0, 3)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),           # Remove all grid lines
    axis.line = element_line(color = "black", linewidth = 1), # Black axis lines (both x and y)
    axis.ticks = element_line(color = "black", linewidth = 1) # Black ticks
  ) +
  # Add labels to the top 10 points
  geom_text_repel(data = DMRtop10_motifs, aes(label = cleaned_motif), 
                  size = 5, color = "black", 
                  max.overlaps = Inf, box.padding = 0.5)
```

![](Figure_codes_files/figure-gfm/Figure_2D-1.png)<!-- -->

``` r
#ggsave("top30K_Dnmt1_DMRmotifs.pdf", width = 8, height = 5, dpi = 300, device = "pdf")
```

## Supplementary Figure 2A

``` r
#1. Find the Dnmt1-interactors from BioGrid
f <- "C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Juan_RNAseq/Dnmt1_BioGrid_Interactors/BIOGRID-GENE-108123-4.4.249.tab3.txt"

df <- read_tsv(f, show_col_types = FALSE)

target_symbols <- c("DNMT1","Dnmt1")
target_entrez  <- c(1786) 

hitA <- (df$`Official Symbol Interactor A` %in% target_symbols) |
        (if ("Entrez Gene Interactor A" %in% names(df)) df$`Entrez Gene Interactor A` %in% target_entrez else FALSE)

hitB <- (df$`Official Symbol Interactor B` %in% target_symbols) |
        (if ("Entrez Gene Interactor B" %in% names(df)) df$`Entrez Gene Interactor B` %in% target_entrez else FALSE)

Dnmt1_interactors <- ifelse(hitA, df$`Official Symbol Interactor B`,
                 ifelse(hitB, df$`Official Symbol Interactor A`, NA))

Dnmt1_interactors <- Dnmt1_interactors[!is.na(Dnmt1_interactors)] |> unique() |> sort()
#writeLines(Dnmt1_interactors, "DNMT1_Dnmt1_interactors.txt")
#cat(length(Dnmt1_interactors), "Dnmt1_interactors written to DNMT1_Dnmt1_interactors.txt\n")

Dnmt1_interactors <- tolower(Dnmt1_interactors)
Dnmt1_interactors <- paste0(toupper(substr(Dnmt1_interactors, 1, 1)),
                      substr(Dnmt1_interactors, 2, nchar(Dnmt1_interactors)))

#2. #Download the Hocomoco v11 core motifs from mouse, import and use as universe
motifs <- read_meme("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Juan_RNAseq/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme")

# ensure we always have a list of universalmotif objects
if (inherits(motifs, "universalmotif")) motifs <- list(motifs)

# extract motif names
raw_names <- vapply(motifs, function(m) methods::slot(m, "name"), character(1))

# keep everything before first "_"
base_names <- sub("_.*", "", raw_names)

# capitalize first letter only
Hocomoco_motifs <- paste0(toupper(substr(base_names, 1, 1)),
                       tolower(substr(base_names, 2, nchar(base_names))))


#3. Plot the BioGrid interactor overlap with the enriched Hocomoco motifs
motif_overlap <- intersect(Hocomoco_motifs, Dnmt1_interactors)

top30KDnmt1_DMRs_vs_matchedbg <- top30KDnmt1_DMRs_vs_matchedbg %>%
  mutate(is_overlap = cleaned_motif %in% motif_overlap)

ggplot(top30KDnmt1_DMRs_vs_matchedbg, aes(x = log2(odds_ratio), y = -log10(adj_p.value), color = TP)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_gradient(low = "gray80", high = "purple3") +
  labs(
    title = "Motif enrichment of top 30K Dnmt1-DMRs",
    color = "#Positive",
    x = "log(OR)",
    y = "-log10(adj.pval)"
  ) +
  coord_cartesian(xlim = c(0, 3)) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1)) +
  geom_text_repel(
    data = subset(top30KDnmt1_DMRs_vs_matchedbg, is_overlap),
    aes(label = cleaned_motif),
    size = 5, color = "black",
    max.overlaps = Inf, box.padding = 0.5, force = 5
  )
```

![](Figure_codes_files/figure-gfm/Supplementary_Figure_2A-1.png)<!-- -->

``` r
#ggsave("Dnmt1_top30Kmotifs_interactors.pdf", device = "pdf", dpi = 300, width = 8, height = 5 )
```

## Load in RNAseq data for rest of plotting

``` r
load("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Juan_RNAseq/Kmt2a_filtered_CPM_level_count_matrix.1.rda")
Kmt2a_CPM <- filtered_expression_data

load("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Juan_RNAseq/Dnmt1_filtered_CPM_level_count_matrix.rda")
Dnmt1_CPM <- filtered_expression_data

Kmt2a_DESeq2 <- read.table("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Juan_RNAseq/Differential_expression_analysis_results_of_Kmt2a_knockouts_without_padj_filtering.csv", sep=";",row.names=1, header = TRUE)

Dnmt1_DESeq2 <- read.table("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Juan_RNAseq/Differential_expression_analysis_results_of_Dnmt1_knockouts_without_padj_filtering.csv", sep=",",row.names=1,header = TRUE)

#Make into one file
Kmt2a_results <- na.omit(merge(as.data.frame(Kmt2a_CPM), Kmt2a_DESeq2, by = 'row.names', all =TRUE))
Kmt2a_results <- dplyr::rename(Kmt2a_results, GeneID = Row.names)
Kmt2a_results <- Kmt2a_results[,c(1,10:15,2:9)]

#Make into one file
Dnmt1_results <- na.omit(merge(as.data.frame(Dnmt1_CPM), Dnmt1_DESeq2, by = 'row.names', all =TRUE))
Dnmt1_results <- dplyr::rename(Dnmt1_results, GeneID = Row.names)
Dnmt1_results <- Dnmt1_results[,c(1,10:15,2:9)]
```

## Figure 2E

``` r
#Load in HOMER annotated Dnmt1-DMR list
Dnmt1_HOMER <- read.delim("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Juan_RNAseq/DEG_DMR_analysis/unphased.Dnmt1.significant.DMRs.HOMER.annotation.txt", sep = "\t",header = TRUE)

## Unique lists
genes_with_HOMER_dmr <- unique(Dnmt1_HOMER$Entrez.ID)
detected_genes <- unique(Dnmt1_results$GeneID)

# Restrict DMR list to detected genes only
dmr_HOMER_detected <- intersect(genes_with_HOMER_dmr, detected_genes)

# split into DEGs vs non-DEGs:
Dnmt1_DEGs <- subset(Dnmt1_results, padj < 0.05)$GeneID

#Now find the overlap
deg_dmr_list <- intersect(dmr_HOMER_detected, Dnmt1_DEGs)

# Add plotting info to Dnmt1_results
Dnmt1_results$group <- ifelse(
  Dnmt1_results$padj > 0.05, "nonsign",
  ifelse(Dnmt1_results$GeneID %in% deg_dmr_list, "dmr", "sign")
)
Dnmt1_results$group <- factor(Dnmt1_results$group, levels = c("nonsign", "sign", "dmr"))

# Define custom group colors
my_colors <- c(
  nonsign = "black",
  sign = "#204664",
  dmr = "purple3"
)

# Subset data for layered plotting
nonsign_df <- subset(Dnmt1_results, group == "nonsign")
sign_df    <- subset(Dnmt1_results, group == "sign")
dmr_df    <- subset(Dnmt1_results, group == "dmr")

# Plot in order: nonsign ➝ sign ➝ gsx1
ggplot() +
  geom_point(data = nonsign_df,
             aes(x = log2FoldChange, y = -log10(padj), color = group),
             size = 2, alpha = 0.5) +
  
  geom_point(data = sign_df,
             aes(x = log2FoldChange, y = -log10(padj), color = group),
             size = 2, alpha = 0.8) +
  
  geom_point(data = dmr_df,
             aes(x = log2FoldChange, y = -log10(padj), color = group),
             size = 2, alpha = 0.5) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +  # padj=0.05

  scale_color_manual(values = my_colors) +
  
  labs(
    title = "sgDnmt1 vs sgCtrl",
    subtitle = "DEGs with hypomethylated regions (padj < 0.05)",
    x = "log2(Fold Change)",
    y = expression(-log[10](padj)),
    color = "Group"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 1),
    axis.ticks = element_line(color = "black", size = 1)
  )
```

![](Figure_codes_files/figure-gfm/Figure_2E-1.png)<!-- -->

``` r
#ggsave("Dnmt1_DMRs_volcanoplot.pdf", width = 7, height = 5,dpi = 300, device = "pdf")
```

## Find overlaps in datasets

``` r
# 1. Merge both datasets by GeneID
sharedD_Kgenes <- merge(
  Dnmt1_results, Kmt2a_results, 
  by = "GeneID", 
  suffixes = c(".Dnmt1", ".Kmt2a")
)

# 2. Genes significant in either dataset
bothsignDEGs <- sharedD_Kgenes[sharedD_Kgenes$padj.Dnmt1 < 0.05 | sharedD_Kgenes$padj.Kmt2a < 0.05, ]

# 3. Genes significant in both datasets
sharedsignDEGs <- sharedD_Kgenes[sharedD_Kgenes$padj.Dnmt1 < 0.05 & sharedD_Kgenes$padj.Kmt2a < 0.05, ]
```

## Plot Figure 3A

``` r
venn_diagram_shared_DEGs <- venn.diagram(
  x = list(
    "Dnmt1 DEGs" = subset(Dnmt1_results, padj < 0.05)$GeneID,
    "Kmt2a DEGs" = subset(Kmt2a_results, padj < 0.05)$GeneID
  ),
  filename = NULL,          # return a grob
  category.names = c("Dnmt1-KO DEGs", "Kmt2a-KO DEGs"),
  print.mode = "raw",
  main = "Overlap",
  main.cex = 1.5,
  disable.logging = TRUE,
  euler.d = TRUE,
  scaled = TRUE,
  hyper.test = TRUE,
  fill = c("orange1", "darkorchid3"),
  alpha = 0.5,
  cat.pos = c(180, 170),
  lwd = 2,
  cex = 1.4,
  cat.cex = 0.8,
  cat.dist = 0.07
)
```

    ## INFO [2025-11-17 14:34:43] $x
    ## INFO [2025-11-17 14:34:43] list(`Dnmt1 DEGs` = subset(Dnmt1_results, padj < 0.05)$GeneID, 
    ## INFO [2025-11-17 14:34:43]     `Kmt2a DEGs` = subset(Kmt2a_results, padj < 0.05)$GeneID)
    ## INFO [2025-11-17 14:34:43] 
    ## INFO [2025-11-17 14:34:43] $filename
    ## INFO [2025-11-17 14:34:43] NULL
    ## INFO [2025-11-17 14:34:43] 
    ## INFO [2025-11-17 14:34:43] $category.names
    ## INFO [2025-11-17 14:34:43] c("Dnmt1-KO DEGs", "Kmt2a-KO DEGs")
    ## INFO [2025-11-17 14:34:43] 
    ## INFO [2025-11-17 14:34:43] $print.mode
    ## INFO [2025-11-17 14:34:43] [1] "raw"
    ## INFO [2025-11-17 14:34:43] 
    ## INFO [2025-11-17 14:34:43] $main
    ## INFO [2025-11-17 14:34:43] [1] "Overlap"
    ## INFO [2025-11-17 14:34:43] 
    ## INFO [2025-11-17 14:34:43] $main.cex
    ## INFO [2025-11-17 14:34:43] [1] 1.5
    ## INFO [2025-11-17 14:34:43] 
    ## INFO [2025-11-17 14:34:43] $disable.logging
    ## INFO [2025-11-17 14:34:43] [1] TRUE
    ## INFO [2025-11-17 14:34:43] 
    ## INFO [2025-11-17 14:34:43] $euler.d
    ## INFO [2025-11-17 14:34:43] [1] TRUE
    ## INFO [2025-11-17 14:34:43] 
    ## INFO [2025-11-17 14:34:43] $scaled
    ## INFO [2025-11-17 14:34:43] [1] TRUE
    ## INFO [2025-11-17 14:34:43] 
    ## INFO [2025-11-17 14:34:43] $hyper.test
    ## INFO [2025-11-17 14:34:43] [1] TRUE
    ## INFO [2025-11-17 14:34:43] 
    ## INFO [2025-11-17 14:34:43] $fill
    ## INFO [2025-11-17 14:34:43] c("orange1", "darkorchid3")
    ## INFO [2025-11-17 14:34:43] 
    ## INFO [2025-11-17 14:34:43] $alpha
    ## INFO [2025-11-17 14:34:43] [1] 0.5
    ## INFO [2025-11-17 14:34:43] 
    ## INFO [2025-11-17 14:34:43] $cat.pos
    ## INFO [2025-11-17 14:34:43] c(180, 170)
    ## INFO [2025-11-17 14:34:43] 
    ## INFO [2025-11-17 14:34:43] $lwd
    ## INFO [2025-11-17 14:34:43] [1] 2
    ## INFO [2025-11-17 14:34:43] 
    ## INFO [2025-11-17 14:34:43] $cex
    ## INFO [2025-11-17 14:34:43] [1] 1.4
    ## INFO [2025-11-17 14:34:43] 
    ## INFO [2025-11-17 14:34:43] $cat.cex
    ## INFO [2025-11-17 14:34:43] [1] 0.8
    ## INFO [2025-11-17 14:34:43] 
    ## INFO [2025-11-17 14:34:43] $cat.dist
    ## INFO [2025-11-17 14:34:43] [1] 0.07
    ## INFO [2025-11-17 14:34:43]

``` r
grid::grid.newpage()
grid::grid.draw(venn_diagram_shared_DEGs)
```

![](Figure_codes_files/figure-gfm/Figure%203_A-1.png)<!-- -->

## Plot Figure 3B

``` r
ggplot(bothsignDEGs, aes(x = log2FoldChange.Kmt2a, y = log2FoldChange.Dnmt1)) +
  geom_point(size = 2.5, alpha = 0.8, color = "gray60") +  # background points
  # Highlight shared significant DEGs
  geom_point(data = sharedsignDEGs,
             aes(x = log2FoldChange.Kmt2a, y = log2FoldChange.Dnmt1),
             shape = 21,
             fill = "gold2",           # Inner color
             color = "gold3",          # Border color
             stroke = 0.7, 
             size = 3, 
             alpha = 0.9) +
   # Linear regression line for shared DEGs
  geom_smooth(data = sharedsignDEGs,
              aes(x = log2FoldChange.Kmt2a, y = log2FoldChange.Dnmt1),
              method = "lm", se = FALSE, color = "gold4", linetype = "dashed", size=0.5) +
  xlim(-2, 2) +
  ylim(-2, 3) +
  labs(
    title = "Comparison of significant genes from Kmt2a and Dnmt1 knockouts",
    x = "log2 Fold Change (Kmt2a)",
    y = "log2 Fold Change (Dnmt1)"
  ) +
  theme_minimal(base_size = 14) 
```

![](Figure_codes_files/figure-gfm/Figure%203B-1.png)<!-- -->

``` r
#ggsave("sharedSign_bothSign.svg", width = 6, height = 5, dpi = 300, device = "svg")
```

## Webgestalt analysis

Now export the shared significant DEGs and load them into
webgestalt.org.  
Choose the “Mus Musculus” genome.  
Method of interest: “Over-representation”.  
“Gene Ontology”  
Functional database: “Biological Process”  
Reference set: “genome protein-coding”  
In “Advanced Parameters” only change “Significance level” to FDR 0.05  
keep the rest on default

## Figure 3C

``` r
#load the data
GO_Bio_UP <- read.delim("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Juan_RNAseq/Webgestalt_GO_Biologicalprocess_SharedKmt2a-Dnmt1/enrichment_results_wg_result1750331313_ShareUP.txt", sep="\t",header=TRUE)

#generate an extra column, that is the % overlap of the genes in each GO term
GO_Bio_UP <- GO_Bio_UP %>%
  mutate(`%overlap` = overlap / size)
GO_Bio_UP$`%overlap` <- 100 * GO_Bio_UP$overlap / GO_Bio_UP$size

#then filter the GO terms with 10 or more genes:
GO_Bio_UP_filtered <- GO_Bio_UP %>%
  filter(size >= 10)

#Filter the top13 GO terms based on FDR values
top13_UP <- GO_Bio_UP_filtered[order(GO_Bio_UP_filtered$FDR), ][1:13, ]

#Change "description" to factor and order based on enrichment ratio
top13_UP$description <- factor(top13_UP$description, levels = top13_UP$description[order(top13_UP$enrichmentRatio)])

#Plot
ggplot(top13_UP, aes(x = enrichmentRatio, y = description, color = -log10(FDR))) +
  geom_point(size = 6) + 
  #scale_color_viridis_c(option = "D") +
  scale_color_gradient(low = "gray80", high = "purple3") +
  #scale_size_continuous(range = c(2, 10)) +
  labs(
    title = "GO-term enrichment analysis",
    x = "Enrichment ratio",
    y = "GO Term",
    color = "-Log10 FDR",
    size = "Gene overlap"
  ) +
  theme_minimal(base_size = 14) +
 theme(
    #panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16),
    
    # Adjust legend to take up less space
    legend.key.size = unit(0.5, "cm"),  # Smaller legend keys
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.text = element_text(size = 10),  # Adjust legend text size
    legend.position = "right",  # Keep legend on right, but compact
    
    # Adjust plot margins to avoid excess empty space
    plot.margin = margin(10, 10, 10, 10) # Adjust the aspect ratio to make the plot wider
)
```

![](Figure_codes_files/figure-gfm/Figure_3C-1.png)<!-- -->

``` r
#Save the plot
#ggsave("GO_Bio_top13UP.svg", width = 8, height = 5, dpi = 300, device = "svg")
```

## Figure 3D

``` r
GO_Bio_Down <- read.delim("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Juan_RNAseq/Webgestalt_GO_Biologicalprocess_SharedKmt2a-Dnmt1/enrichment_results_wg_result1750331313_ShareDown.txt", sep="\t",header=TRUE)

#generate an extra column, that is the % overlap of the genes in each GO term
GO_Bio_Down <- GO_Bio_Down %>%
  mutate(`%overlap` = overlap / size)
GO_Bio_Down$`%overlap` <- 100 * GO_Bio_Down$overlap / GO_Bio_Down$size

#then filter the GO terms with 10 or more genes:
GO_Bio_Down_filtered <- GO_Bio_Down %>%
  filter(size >= 10)

#Find the top enriched terms
top13_Down <- GO_Bio_Down_filtered[order(GO_Bio_Down_filtered$FDR), ][1:13, ]
top13_Down$description <- factor(top13_Down$description, levels = top13_Down$description[order(top13_Down$enrichmentRatio)])

#Plot the results
ggplot(top13_Down, aes(x = enrichmentRatio, y = description, color = -log10(FDR))) +
  geom_point(size = 6) + 
  #scale_color_viridis_c(option = "D") +
  scale_color_gradient(low = "gray80", high = "purple3") +
  #scale_size_continuous(range = c(2, 10)) +
  labs(
    title = "GO-term enrichment analysis",
    x = "Enrichment ratio",
    y = "GO Term",
    color = "-Log10 FDR",
    size = "Gene overlap"
  ) +
  theme_minimal(base_size = 14) +
 theme(
    #panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16),
    
    # Adjust legend to take up less space
    legend.key.size = unit(0.5, "cm"),  # Smaller legend keys
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.text = element_text(size = 10),  # Adjust legend text size
    legend.position = "right"  # Keep legend on right, but compact
)
```

![](Figure_codes_files/figure-gfm/Figure_3D-1.png)<!-- -->

``` r
#ggsave("GO_Bio_top13Down.svg", width = 6, height = 5, dpi = 300, device = "svg")
```

## Supplementary Figure 3G

``` r
#list the EM genes from epigeneticmachinery.org
EM_genes <- c("Aire","Akap1","Alg13","Ash1l","Asxl1","Asxl2","Asxl3","Atad2","Atad2b","Atrx","Bahcc1","Bahd1","Baz1a","Baz1b","Baz2a","Baz2b","Bptf","Brd1","Brd2","Brd3","Brd4","Brd7","Brd8","Brd9","Brdt","Brpf1","Brpf3","Brwd1","Brwd3","C14orf169","Cbx1","Cbx2","Cbx3","Cbx4","Cbx5","Cbx6","Cbx7","Cbx8","Cdy1","Cdy2a","Cdyl","Cdyl2","Cecr2","Chd1","Chd2","Chd3","Chd4","Chd5","Chd6","Chd7","Chd8","Chd9","Crebbp","Cxxc1","Cxxc4","Cxxc5","Dido1","Dnmt1","Dnmt3a","Dnmt3b","Dnmt3l","Dpf1","Dpf2","Dpf3","Eed","Ehmt1","Ehmt2","Ep300","Ep400","Ezh1","Ezh2","Fbxl19","G2e3","Glyr1","Hdac1","Hdac10","Hdac11","Hdac2","Hdac3","Hdac4","Hdac5","Hdac6","Hdac7","Hdac8","Hdac9","Hdgf","Hdgfl1","Hdgfrp2","Hdgfrp3","Hif1an","Hr","Hspbap1","Ing1","Ing2","Ing3","Ing4","Ing5","Ino80","Ints12","Jade1","Jade2","Jade3","Jarid2","Jmjd1c","Jmjd4","Jmjd6","Jmjd7","Jmjd8","Kat2a","Kat2b","Kat5","Kat6a","Kat6b","Kat7","Kat8","Kdm1a","Kdm1b","Kdm2a","Kdm2b","Kdm3a","Kdm3b","Kdm4a","Kdm4b","Kdm4c","Kdm4d","Kdm4e","Kdm5a","Kdm5b","Kdm5c","Kdm5d","Kdm6a","Kdm6b","Kdm7a","Kdm8","Neurodev","Kmt2b","Kmt2c","Kmt2d","Kmt2e","Setd8","Kmt5b","Kmt5c","L3mbtl1","L3mbtl2","L3mbtl3","L3mbtl4","Lbr","Mbd1","Mbd2","Mbd3","Mbd4","Mbd5","Mbd6","Mbtd1","Mecp2","Mina","Mllt10","Mllt6","Morc1","Morc2","Morc3","Morc4","Mphosph8","Msh6","Msl3","Mta1","Mta2","Mta3","Mtf2","Mum1","Mum1l1","Nsd1","Orc1","Pbrm1","Phf1","Phf10","Phf11","Phf12","Phf13","Phf14","Phf19","Phf2","Phf20","Phf20l1","Phf21a","Phf21b","Phf23","Phf3","Phf6","Phf7","Phf8","Phip","Phrf1","Prdm1","Prdm10","Prdm11","Prdm12","Prdm13","Prdm14","Prdm15","Prdm16","Prdm2","Prdm4","Prdm5","Prdm6","Prdm7","Prdm8","Prdm9","Psip1","Pwwp2a","Pwwp2b","Pygo1","Pygo2","Rag2","Rai1","Rere","Rnf17","Rsf1","Scmh1","Scml2","Setd1a","Setd1b","Setd2","Setd3","Setd4","Setd5","Setd6","Setd7","Setd9","Setdb1","Setdb2","Setmar","Sfmbt1","Sfmbt2","Shprh","Sirt1","Sirt2","Sirt3","Sirt4","Sirt5","Sirt6","Sirt7","Smarca1","Smarca2","Smarca4","Smarca5","Smn1","Smndc1","Smyd1","Smyd2","Smyd3","Smyd4","Smyd5","Snd1","Sp110","Sp140","Sp140l","Srcap","Stk31","Suv39h1","Suv39h2","Taf1","Taf1l","Taf3","Tcf19","Tcf20","Tdrd1","Tdrd10","Tdrd12","Tdrd15","Tdrd3","Tdrd5","Tdrd6","Tdrd7","Tdrd9","Tdrkh","Tet1","Tet2","Tet3","Tnrc18","Trim24","Trim28","Trim33","Trim66","Tyw5","Ubr7","Uhrf1","Uhrf2","Uty","Whsc1","Whsc1l1","Zcwpw1","Zcwpw2","Zmynd11","Zmynd8")

bothsign_EMgenes <- bothsignDEGs %>% filter(GeneID %in% EM_genes)
sharedsign_EMgenes <- sharedsignDEGs %>% filter(GeneID %in% EM_genes)

Interesting_EMs <- c("Ezh2","Dnmt1","Kmt2a","Uhrf1","Trim28","Kat2b","Setdb1")
interesting_EM <- bothsign_EMgenes %>% filter(GeneID %in% Interesting_EMs)

ggplot(bothsign_EMgenes, aes(x = log2FoldChange.Kmt2a, y = log2FoldChange.Dnmt1)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  geom_point(aes(color = "Either signifcant EM genes", fill = "Either signifcant EM genes", shape = "Either signifcant EM genes"), size = 2.5, alpha = 0.5) +
  geom_point(data = sharedsign_EMgenes,
             aes(color = "Shared significant EMs", fill = "Shared significant EMs", shape = "Shared significant EMs"),
             size = 3, alpha = 0.5, stroke = 0.7) +
  geom_text_repel(
  data = interesting_EM,
  aes(label = GeneID),
  color = "black", size = 4, max.overlaps = 20,
  show.legend = FALSE ,          # keeps labels out of the legend
  box.padding = 0.5) +
  scale_shape_manual(name = "Group",
    values = c("Either signifcant EM genes" = 16, "Shared significant EMs" = 21)) +
  scale_color_manual(name = "Group",
    values = c("Either signifcant EM genes" = "gray60", "Shared significant EMs" = "gold3")) +
  scale_fill_manual(name = "Group",
    values = c("Either signifcant EM genes" = NA, "Shared significant EMs" = "gold2")) +
  xlim(-0.5, 0.5) + ylim(-1.5, 1) +
  labs(title = "Comparison of EM genes from Kmt2a and Dnmt1 knockouts",
       x = "log2 Fold Change (Kmt2a)", y = "log2 Fold Change (Dnmt1)") +
  theme_minimal(base_size = 14)
```

![](Figure_codes_files/figure-gfm/Supplementary_Figure_3G-1.png)<!-- -->

``` r
#ggsave("shared_EMgenes.pdf", width = 7, height = 4, dpi = 300, device = "pdf")
```

## Motif enrichment analysis of shared DEG promoters vs all protein-coding promoters

To perform motif enrichment on shared DEG promoters, we first define
promoter regions.  
We select canonical protein-coding transcripts from Ensembl, then use a
TxDb transcript database for mm10 to extract promoter coordinates (2 kb
upstream, 500 bp downstream of the TSS).  
We then split promoters into those belonging to shared DEGs and the
remaining protein-coding promoters, which serve as the background.

``` r
ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl") #Ensembl release 102 uses GRCm38

# Query all canonical + protein-coding transcripts
canonical_protein_coding_tx <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "ensembl_transcript_id",
    "transcript_is_canonical",
    "transcript_biotype"
  ),
  mart = ensembl
)

# Get all transcripts from TxDb
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
txs <- transcripts(txdb)

# Strip version numbers from tx_name
txs$tx_id_clean <- sub("\\.\\d+$", "", txs$tx_name)

# Filter for canonical protein-coding
canonical_protein_coding_tx <- canonical_protein_coding_tx %>%
  filter(transcript_is_canonical == 1, transcript_biotype == "protein_coding")

# Subset to your canonical Ensembl transcript IDs
txs_canonical_protein_coding <- txs[txs$tx_id_clean %in% canonical_protein_coding_tx$ensembl_transcript_id]

# Define promoter regions (default: 2000 bp upstream and 500 bp downstream of TSS)
promoters_canonical_proteincoding <- promoters(txs_canonical_protein_coding, upstream = 2000, downstream = 500)

promoters_canonical_proteincoding <- trim(promoters_canonical_proteincoding)

# Now find the shared significant promoters
sharedsig_genes <- sharedsignDEGs$GeneID

# Subset promoters based on the shared significant genes
canonical_tx_sigshared <- canonical_protein_coding_tx %>%
  filter(external_gene_name %in% sharedsig_genes)

# Subset to your canonical Ensembl transcript IDs
txs_canonical_sigshared <- txs[txs$tx_id_clean %in% canonical_tx_sigshared$ensembl_transcript_id]

# Extract the promoter regions:
# Define promoter regions (default: 2000 bp upstream and 500 bp downstream of TSS)
promoters_canonical_sigshared <- promoters(txs_canonical_sigshared, upstream = 2000, downstream = 500)

promoters_canonical_sigshared <- trim(promoters_canonical_sigshared)

#export(promoters_canonical_sigshared, con = "shared_significant_canonical_promoters_mm10.bed", format = "BED")

#and then the appropriate background:
promoters_sharedbackground <- setdiff(promoters_canonical_proteincoding, promoters_canonical_sigshared)

#export(promoters_sharedbackground, con = "Proteincoding_canonical_promoters_wosharedDEGs_mm10.bed", format = "BED")
```

## AME of promoters

Now open <https://meme-suite.org/meme/tools/ame>  
Load in the “shared_significant_canonical_promoters_mm10.bed” as the
input BED file  
Load in the “Proteincoding_canonical_promoters_wosharedDEGs_mm10.bed” as
the control BED files  
be sure to choose the mm10 genome for both

In the motif selection choose:  
MOUSE (mus musculus) DNA  
HOCOMOCO Mouse (v11 CORE)  
Sequence scoring method: Average Odds score  
Motif enrichment test: Fisher’s exact test  
Export the AME file

## Figure 3I

``` r
sharedMotifs_vs_promoters <- read.delim("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Juan_RNAseq/ame_sharedDEGs_vsallpromoters-woDEGs.tsv", sep="\t",header=TRUE)

sharedMotifs_vs_promoters <- sharedMotifs_vs_promoters %>%
  mutate(
    cleaned_motif = str_to_title(str_extract(motif_ID, "^[^_]+"))
  )

sharedMotifs_vs_promoters <- sharedMotifs_vs_promoters %>%
  mutate(odds_ratio = (X.TP/100) / (1 - X.TP/100) / ((X.FP/100) / (1 - X.FP/100)))

sharedtop10_motifs <- sharedMotifs_vs_promoters %>%
  arrange(adj_p.value) %>%   # Ascending order = smallest E-values
  dplyr::slice(1:10)

# Add Egr1 and binding partners from cleaned_motif (filtering by name)

Extra_rows <- c("Myc", "Mycn", "Nfib")
Extra <- sharedMotifs_vs_promoters %>% filter(cleaned_motif %in% Extra_rows)

# Combine
sharedtop10_motifs <- bind_rows(sharedtop10_motifs, Extra)


ggplot(sharedMotifs_vs_promoters, aes(x = log2(odds_ratio), y = -log10(adj_p.value), color = TP)) +
  geom_point(size=4, alpha = 0.8) +  # Scatter points
  scale_color_gradient(low = "gray80", high = "purple3") +
  labs(
    title = "Motif enrichment of top 30K Dnmt1-DMRs",
    color = "#Positive",
    #size = "Overlap (%)",
    x = "Log2(OR)",
    y = "-log10(adj.pval)"
  ) +
  coord_cartesian(xlim = c(0, 3)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),           # Remove all grid lines
    axis.line = element_line(color = "black", size = 1), # Black axis lines (both x and y)
    axis.ticks = element_line(color = "black", size = 1), # Black ticks
    #axis.text = element_text(color = "black"),  # Black axis text
    #axis.title = element_text(color = "black")#, # Black axis titles
    #plot.title = element_text(hjust = 0.5, size = 16) # Center the plot title
  ) +
  # Add labels to the top 10 points
  geom_text_repel(data = sharedtop10_motifs, aes(label = cleaned_motif), 
                  size = 5, color = "black", 
                  max.overlaps = Inf, box.padding = 0.5)
```

![](Figure_codes_files/figure-gfm/Figure_3I-1.png)<!-- -->

``` r
#ggsave("shared_DEGmotifs.pdf", width = 6, height = 5, dpi = 300, device = "pdf")
```

## GSEA analysis

This was done exactly as previously described in Rafnsdottir et
al. 2024, Cell Reports, DOI: 10.1016/j.celrep.2024.114554  
GSEA was used to identify significantly enriched gene sets within our
RNAseq data. For this, we generated ranked gene lists comprising all
detected genes assigned to a rank value. The rank value was calculated
as log10(q-value) × (sign of fold change). The genes were then ranked
from the highest rank value to the lowest. These lists were generated in
R:

``` r
# First rank the Dnmt1 dataset
Dnmt1.preranked.score <- Dnmt1_DESeq2
Dnmt1.preranked.score.up <- filter(Dnmt1.preranked.score, log2FoldChange >0)
Dnmt1.preranked.score.up$score <- log10(1/Dnmt1.preranked.score.up$padj)

Dnmt1.preranked.score.down <- filter(Dnmt1.preranked.score, log2FoldChange < 0)
Dnmt1.preranked.score.down$score <- -log10(1/Dnmt1.preranked.score.down$padj)

Dnmt1.preranked.score.final <- rbind(Dnmt1.preranked.score.up, Dnmt1.preranked.score.down)
Dnmt1.preranked.score.final <- data.frame(
  gene  = rownames(Dnmt1.preranked.score.final),
  score = Dnmt1.preranked.score.final[, 7],
  row.names = NULL
)
Dnmt1.preranked.score.final <- Dnmt1.preranked.score.final[order(-Dnmt1.preranked.score.final$score),]
#write.csv(Dnmt1.preranked.score.final, "Dnmt1_preranked_score.csv")

#Now do the same for Kmt2a dataset:
Kmt2a.preranked.score <- Kmt2a_DESeq2
Kmt2a.preranked.score.up <- filter(Kmt2a.preranked.score, log2FoldChange >0)
Kmt2a.preranked.score.up$score <- log10(1/Kmt2a.preranked.score.up$padj)

Kmt2a.preranked.score.down <- filter(Kmt2a.preranked.score, log2FoldChange < 0)
Kmt2a.preranked.score.down$score <- -log10(1/Kmt2a.preranked.score.down$padj)

Kmt2a.preranked.score.final <- rbind(Kmt2a.preranked.score.up, Kmt2a.preranked.score.down)
Kmt2a.preranked.score.final <- data.frame(
  gene  = rownames(Kmt2a.preranked.score.final),
  score = Kmt2a.preranked.score.final[, 7],
  row.names = NULL
)
Kmt2a.preranked.score.final <- Kmt2a.preranked.score.final[order(-Kmt2a.preranked.score.final$score),]
#write.csv(Kmt2a.preranked.score.final, "Kmt2a_preranked_score.csv")
```

This file was then opened in excel, headers of all columns were deleted
as well as all columns except for gene names and score. The file was
saved as a tab-delimited .txt rile with the .rnk extension. This was
then loaded into GSEA preranked program (v. 4.3.2). We used the
following settings:

Gene set:
ftp.broadinstitute.org://pub/gsea/msigdb/mouse/gene_sets/mh.all.v2024.1.Mm.symbols.gmt  
Collapse/remap to gene symbols: No_Collapse  
Chip platform:
ftp.broadinstitute.org://pub/gsea/msigdb/mouse/annotations/Mouse_Ensembl_Gene_ID_MSigDB.v2024.1.Mm.chip  
Number of Permutations: 1000  
Enrichment statistic: classic  
Max size: 5000  
Min size: 15  
Normalisation mode: meandiv

## Supplemental Figure 3A

``` r
Dnmt1_GSEA <- read.delim("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Juan_RNAseq/GSEA/Dnmt1_Deseq2_2.GseaPreranked.1759761114674/gsea_report_for_na_neg_1759761114674.tsv", sep="\t",header = TRUE)

Dnmt1_sign_GSEA <- Dnmt1_GSEA[Dnmt1_GSEA$FDR.q.val<0.05,] %>% arrange(desc(NES)) %>%
  mutate(NAME = factor(NAME, levels = NAME))

#FDRs are so small that they are represented as 0. 
min_pos <- min(Dnmt1_sign_GSEA$FDR.q.val[Dnmt1_sign_GSEA$FDR.q.val > 0], na.rm = TRUE)
eps <- min_pos / 10
Dnmt1_sign_GSEA <- Dnmt1_sign_GSEA %>% mutate(logFDR = -log10(pmax(FDR.q.val, eps)))



ggplot(Dnmt1_sign_GSEA, aes(x = (-1)*NES, y = NAME, color = -log10(FDR.q.val))) +
 geom_point(aes(size = SIZE, color = logFDR))  + 
  #scale_color_viridis_c(option = "D") +
  scale_color_gradient(low = "gray80", high = "purple3") +
  #scale_size_continuous(range = c(2, 10)) +
  labs(
    title = "Dnmt1 GSEA",
    x = "Normalized Enrichment Score (NES)",
    y = "HALLMARK term",
    color = "-Log10 FDR"
    ) +
  theme_minimal(base_size = 14) +
 theme(
    #panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16),
    
    # Adjust legend to take up less space
    legend.key.size = unit(0.5, "cm"),  # Smaller legend keys
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.text = element_text(size = 10),  # Adjust legend text size
    legend.position = "right"  # Keep legend on right, but compact
    
    # Adjust plot margins to avoid excess empty space
    #plot.margin = margin(10, 10, 10, 10) # Adjust the aspect ratio to make the plot wider
)
```

![](Figure_codes_files/figure-gfm/Supplemental_Figure_3A-1.png)<!-- -->

``` r
#ggsave("Dnmt1_GSEA.svg", width = 7, height = 4, dpi = 300, device = "svg")
```

## Supplemental Figure 3B

``` r
Kmt2a_GSEA <- read.delim("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Juan_RNAseq/GSEA/Kmt2a_Deseq2_2.GseaPreranked.1759761134160/gsea_report_for_na_neg_1759761134160.tsv", sep="\t",header = TRUE)

Kmt2a_sign_GSEA <- Kmt2a_GSEA[Kmt2a_GSEA$FDR.q.val<0.05,] %>% arrange(desc(NES)) %>%
  mutate(NAME = factor(NAME, levels = NAME))

#FDRs are so small that they are represented as 0. 
min_pos <- min(Kmt2a_sign_GSEA$FDR.q.val[Kmt2a_sign_GSEA$FDR.q.val > 0], na.rm = TRUE)
eps <- min_pos / 10
Kmt2a_sign_GSEA <- Kmt2a_sign_GSEA %>% mutate(logFDR = -log10(pmax(FDR.q.val, eps)))



ggplot(Kmt2a_sign_GSEA, aes(x = (-1)*NES, y = NAME, color = -log10(FDR.q.val))) +
 geom_point(aes(size = SIZE, color = logFDR))  + 
  #scale_color_viridis_c(option = "D") +
  scale_color_gradient(low = "gray80", high = "purple3") +
  #scale_size_continuous(range = c(2, 10)) +
  labs(
    title = "Kmt2a GSEA",
    x = "Normalized Enrichment Score (NES)",
    y = "HALLMARK term",
    color = "-Log10 FDR"
    ) +
  theme_minimal(base_size = 14) +
 theme(
    #panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16),
    
    # Adjust legend to take up less space
    legend.key.size = unit(0.5, "cm"),  # Smaller legend keys
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.text = element_text(size = 10),  # Adjust legend text size
    legend.position = "right"  # Keep legend on right, but compact
    
    # Adjust plot margins to avoid excess empty space
    #plot.margin = margin(10, 10, 10, 10) # Adjust the aspect ratio to make the plot wider
)
```

![](Figure_codes_files/figure-gfm/Supplemental_Figure_3B-1.png)<!-- -->

``` r
#ggsave("Kmt2a_GSEA.svg", width = 7, height = 4, dpi = 300, device = "svg")
```

## Find overlaps

Now find Kmt2a BioGrid Interactors and what is shared between Dnmt1 and
Kmt2a

``` r
f <- "C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Juan_RNAseq/Kmt2a_BioGrid_Interactors/BIOGRID-GENE-110443-4.4.249.tab3.txt"

df <- read_tsv(f, show_col_types = FALSE)

target_symbols <- c("KMT2A","Kmt2a")
target_entrez  <- c(1786, 13433) # human, mouse (use what you need)

hitA <- (df$`Official Symbol Interactor A` %in% target_symbols) |
        (if ("Entrez Gene Interactor A" %in% names(df)) df$`Entrez Gene Interactor A` %in% target_entrez else FALSE)

hitB <- (df$`Official Symbol Interactor B` %in% target_symbols) |
        (if ("Entrez Gene Interactor B" %in% names(df)) df$`Entrez Gene Interactor B` %in% target_entrez else FALSE)

Kmt2a_interactors <- ifelse(hitA, df$`Official Symbol Interactor B`,
                 ifelse(hitB, df$`Official Symbol Interactor A`, NA))

Kmt2a_interactors <- Kmt2a_interactors[!is.na(Kmt2a_interactors)] |> unique() |> sort()
#writeLines(Kmt2a_interactors, "Kmt2a_interactors.txt")
#cat(length(Kmt2a_interactors), "interactors written to Kmt2a_interactors.txt\n")

Kmt2a_interactors <- tolower(Kmt2a_interactors)
Kmt2a_interactors <- paste0(toupper(substr(Kmt2a_interactors, 1, 1)),
                      substr(Kmt2a_interactors, 2, nchar(Kmt2a_interactors)))

#Find the overlaps:
ol_interactors <- intersect(Dnmt1_interactors, Kmt2a_interactors)

#Find overlaps between interactors and motifs to plot
olinteractors_vs_motifs <- intersect(sharedMotifs_vs_promoters$cleaned_motif, ol_interactors)
```

## Fisher’s exact test

To perform Fisher’s test in the overlap in interactors, we decided to
use all interactions reported for humans, that are expressed in the
RNAseq, as Universe

``` r
f <- "C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Juan_RNAseq/All_BioGrid_Interactions/BIOGRID-ORGANISM-Homo_sapiens-4.4.249.tab3.txt"

df <- read_tsv(f, show_col_types = FALSE)

universe <- df |>
  dplyr::select(`Official Symbol Interactor A`, `Official Symbol Interactor B`) |>
  pivot_longer(everything(), names_to = "role", values_to = "symbol") |>
  transmute(symbol = trimws(symbol)) |>
  filter(!is.na(symbol), symbol != "") |>
  distinct(symbol) |>
  arrange(symbol)

# Capitalize first letter only
universe$symbol <- paste0(
  toupper(substr(universe$symbol, 1, 1)),
  tolower(substr(universe$symbol, 2, nchar(universe$symbol)))
)

universe <- unique(universe["symbol"])

# universe_symbols: all BioGRID-observed symbols (>=1 interaction) for your organism
# expressed_symbols: genes expressed/detectable in your system
U <- intersect(universe$symbol, sharedD_Kgenes$GeneID)

A_set <- Kmt2a_interactors       # as in earlier code
B_set <- Dnmt1_interactors

# restrict sets to the expression-aware universe
A_set <- intersect(A_set, U)
B_set <- intersect(B_set, U)

a <- length(intersect(A_set, B_set))
b <- length(setdiff(A_set, B_set))
c <- length(setdiff(B_set, A_set))
d <- length(U) - (a + b + c)

fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  matrix(c(a, b, c, d), nrow = 2)
    ## p-value < 2.2e-16
    ## alternative hypothesis: true odds ratio is greater than 1
    ## 95 percent confidence interval:
    ##  13.63268      Inf
    ## sample estimates:
    ## odds ratio 
    ##   18.50512

## Supplementary Figure 3H

``` r
#Now plot only the shared interactors that are expressed
venn_diagram_Dnmt1_imprinting <- venn.diagram(
    x = list("Dnmt1 interactors" = B_set,
             "Kmt2a interactors" = A_set), 
    filename = NULL,  # Set to NULL to plot to the current device
    category.names = c("Dnmt1 interactors" , 
                       "Kmt2a interactors"),
    print.mode = c("raw"),
    main = "Overlap",
    main.cex = 1.5,
    disable.logging = TRUE,
    euler.d = TRUE,
    scaled = TRUE,
    hyper.test = TRUE,
    fill = c("orange1", "darkorchid3"),
    alpha = 0.50,
    cat.pos = c(180,170),
    lwd = 2,
    cex = 1.4,
    cat.cex = 0.8,
    cat.dist = 0.07
)
```

    ## INFO [2025-11-17 14:36:08] $x
    ## INFO [2025-11-17 14:36:08] list(`Dnmt1 interactors` = B_set, `Kmt2a interactors` = A_set)
    ## INFO [2025-11-17 14:36:08] 
    ## INFO [2025-11-17 14:36:08] $filename
    ## INFO [2025-11-17 14:36:08] NULL
    ## INFO [2025-11-17 14:36:08] 
    ## INFO [2025-11-17 14:36:08] $category.names
    ## INFO [2025-11-17 14:36:08] c("Dnmt1 interactors", "Kmt2a interactors")
    ## INFO [2025-11-17 14:36:08] 
    ## INFO [2025-11-17 14:36:08] $print.mode
    ## INFO [2025-11-17 14:36:08] c("raw")
    ## INFO [2025-11-17 14:36:08] 
    ## INFO [2025-11-17 14:36:08] $main
    ## INFO [2025-11-17 14:36:08] [1] "Overlap"
    ## INFO [2025-11-17 14:36:08] 
    ## INFO [2025-11-17 14:36:08] $main.cex
    ## INFO [2025-11-17 14:36:08] [1] 1.5
    ## INFO [2025-11-17 14:36:08] 
    ## INFO [2025-11-17 14:36:08] $disable.logging
    ## INFO [2025-11-17 14:36:08] [1] TRUE
    ## INFO [2025-11-17 14:36:08] 
    ## INFO [2025-11-17 14:36:08] $euler.d
    ## INFO [2025-11-17 14:36:08] [1] TRUE
    ## INFO [2025-11-17 14:36:08] 
    ## INFO [2025-11-17 14:36:08] $scaled
    ## INFO [2025-11-17 14:36:08] [1] TRUE
    ## INFO [2025-11-17 14:36:08] 
    ## INFO [2025-11-17 14:36:08] $hyper.test
    ## INFO [2025-11-17 14:36:08] [1] TRUE
    ## INFO [2025-11-17 14:36:08] 
    ## INFO [2025-11-17 14:36:08] $fill
    ## INFO [2025-11-17 14:36:08] c("orange1", "darkorchid3")
    ## INFO [2025-11-17 14:36:08] 
    ## INFO [2025-11-17 14:36:08] $alpha
    ## INFO [2025-11-17 14:36:08] [1] 0.5
    ## INFO [2025-11-17 14:36:08] 
    ## INFO [2025-11-17 14:36:08] $cat.pos
    ## INFO [2025-11-17 14:36:08] c(180, 170)
    ## INFO [2025-11-17 14:36:08] 
    ## INFO [2025-11-17 14:36:08] $lwd
    ## INFO [2025-11-17 14:36:08] [1] 2
    ## INFO [2025-11-17 14:36:08] 
    ## INFO [2025-11-17 14:36:08] $cex
    ## INFO [2025-11-17 14:36:08] [1] 1.4
    ## INFO [2025-11-17 14:36:08] 
    ## INFO [2025-11-17 14:36:08] $cat.cex
    ## INFO [2025-11-17 14:36:08] [1] 0.8
    ## INFO [2025-11-17 14:36:08] 
    ## INFO [2025-11-17 14:36:08] $cat.dist
    ## INFO [2025-11-17 14:36:08] [1] 0.07
    ## INFO [2025-11-17 14:36:08]

``` r
grid::grid.newpage()
grid::grid.draw(venn_diagram_Dnmt1_imprinting)
```

![](Figure_codes_files/figure-gfm/Supplementary_Figure_3H-1.png)<!-- -->

# Neurodevelopmental RNAseq datasets

## load in the data

``` r
Day2_results <- read.table("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Neurodev_analysis/NPCvsDay2.results/NPCvsDay2.comparison.all.genes.csv", sep = ",", header=TRUE)
colnames(Day2_results)[1]<- "GeneID"

Day6_results <- read.table("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Neurodev_analysis/NPCvsDay6.results/NPCvsDay6.comparison.all.genes.csv", sep = ",", header=TRUE)
colnames(Day6_results)[1]<- "GeneID"
```

## Webgestalt analysis

Now export the significant DEGs for either dataset and load them into
webgestalt.org  
Choose the Mus Musculus genome  
Method of interest: Over-representation  
Gene Ontology  
Functional database: Biological Process  
Reference set: genome protein-coding  
In “Advanced Parameters” only change “Significance level” to FDR 0.05  
keep the rest on default

## Supplementary Figure 5C

``` r
GO_D2_Bio_UP <- read.delim("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Neurodev_analysis/Webgestalt_GO_Biological_NeurodevDay2/enrichment_results_wg_result1750595464_Neurod2up.txt", sep="\t",header=TRUE)

#generate an extra column, that is the % overlap of the genes in each GO term
GO_D2_Bio_UP <- GO_D2_Bio_UP %>%
  mutate(`%overlap` = overlap / size)
GO_D2_Bio_UP$`%overlap` <- 100 * GO_D2_Bio_UP$overlap / GO_D2_Bio_UP$size

#then filter the GO terms with 10 or more genes:
GO_D2_Bio_UP_filtered <- GO_D2_Bio_UP %>%
  filter(size >= 10)

#Find the top enriched terms
top13_D2_UP <- GO_D2_Bio_UP_filtered[order(GO_D2_Bio_UP_filtered$FDR), ][1:13, ]
top13_D2_UP$description <- factor(top13_D2_UP$description, levels = top13_D2_UP$description[order(top13_D2_UP$enrichmentRatio)])

#Then plot them
ggplot(top13_D2_UP, aes(x = enrichmentRatio, y = description, fill = -log10(FDR))) +
  geom_col() +
  scale_fill_gradient(low = "gray80", high = "#204664", name = "-log10(FDR)") +
  labs(x = "Enrichment Ratio", y = "Pathway", title = "Pathway Enrichment") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black", size = 1)
        )
```

![](Figure_codes_files/figure-gfm/Supplementary_Figure_5C-1.png)<!-- -->

``` r
#ggsave("GO_D2_Bio_top13UP_bars.svg", width = 8, height = 5, dpi = 300, device = "svg")
```

## Supplementary Figure 5D

``` r
GO_D6_Bio_UP <- read.delim("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/Neurodev_analysis/NPCvsDay6.results/Webgestalt_GO_BiolProc_NeurodevDay6/enrichment_results_wg_result1750598038_Neurod6up.txt", sep="\t",header=TRUE)

#generate an extra column, that is the % overlap of the genes in each GO term
GO_D6_Bio_UP <- GO_D6_Bio_UP %>%
  mutate(`%overlap` = overlap / size)
GO_D6_Bio_UP$`%overlap` <- 100 * GO_D6_Bio_UP$overlap / GO_D6_Bio_UP$size

#then filter the GO terms with 10 or more genes:
GO_D6_Bio_UP_filtered <- GO_D6_Bio_UP %>%
  filter(size >= 10)

#Find the top enriched terms
top13_D6_UP <- GO_D6_Bio_UP_filtered[order(GO_D6_Bio_UP_filtered$FDR), ][1:13, ]
top13_D6_UP$description <- factor(top13_D6_UP$description, levels = top13_D6_UP$description[order(top13_D6_UP$enrichmentRatio)])

#Plot these:
ggplot(top13_D6_UP, aes(x = enrichmentRatio, y = description, fill = -log10(FDR))) +
  geom_col() +
  scale_fill_gradient(low = "gray80", high = "#204664", name = "-log10(FDR)") +
  labs(x = "Enrichment Ratio", y = "Pathway", title = "Pathway Enrichment") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(color = "black", size = 1)
        )
```

![](Figure_codes_files/figure-gfm/Supplementary_Figure_5D-1.png)<!-- -->

``` r
#ggsave("GO_D6_Bio_top13UP_bars.svg", width = 8, height = 5, dpi = 300, device = "svg")
```

## Plot all the overlaps with Dnmt1 and Kmt2a KO datasets

### Plot Dnmt1_Day2 overlap

``` r
# 1. Merge both datasets by GeneID
sharedD_D2genes <- merge(
  Dnmt1_results, Day2_results, 
  by = "GeneID", 
  suffixes = c(".Dnmt1", ".Day2")
)

# 2. Genes significant in either dataset
eitherD_D2signDEGs <- sharedD_D2genes[sharedD_D2genes$padj.Dnmt1 < 0.05 | sharedD_D2genes$padj.Day2 < 0.05, ]

# 3. Genes significant in both datasets
sharedD_D2signDEGs <- sharedD_D2genes[sharedD_D2genes$padj.Dnmt1 < 0.05 & sharedD_D2genes$padj.Day2 < 0.05, ]

p1 <- ggplot(eitherD_D2signDEGs, aes(x = log2FoldChange.Day2, y = log2FoldChange.Dnmt1)) +
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  
  geom_point(size = 1, stroke = 0.3, alpha = 0.8, color = "gray60") +  # background points
 # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  
  # Highlight shared significant DEGs
  geom_point(data = sharedD_D2signDEGs,
             aes(x = log2FoldChange.Day2, y = log2FoldChange.Dnmt1),
             shape = 21,
             fill = "#FFC75F",           # Inner color
             color = "#CDAD00",           # Border color
             stroke = 0.3, 
             size = 1, 
             alpha = 0.9) +
   # Linear regression line for shared DEGs
  geom_smooth(data = sharedD_D2signDEGs,
              aes(x = log2FoldChange.Day2, y = log2FoldChange.Dnmt1),
              method = "lm", se = FALSE, color = "gold4", fullrange = TRUE, size=0.5) +

  xlim(-10, 10) +
  ylim(-2, 2) +
  labs(
    title = "Dnmt1-KO vs Day2",
    x = "log2FoldChange(Day2)",
    y = "log2FoldChange(Dnmt1)"
  ) +
  theme_minimal(base_size = 12)


#ggsave("sharedDnmt1-Day2_scatterplot.pdf", width = 6, height = 5, dpi = 300, device = "pdf")
```

### Plot Dnmt1_Day6 overlap

``` r
# 1. Merge both datasets by GeneID
sharedD_D6genes <- merge(
  Dnmt1_results, Day6_results, 
  by = "GeneID", 
  suffixes = c(".Dnmt1", ".Day6")
)

# 2. Genes significant in either dataset
eitherD_D6signDEGs <- sharedD_D6genes[sharedD_D6genes$padj.Dnmt1 < 0.05 | sharedD_D6genes$padj.Day6 < 0.05, ]

# 3. Genes significant in both datasets
sharedD_D6signDEGs <- sharedD_D6genes[sharedD_D6genes$padj.Dnmt1 < 0.05 & sharedD_D6genes$padj.Day6 < 0.05, ]

p2 <- ggplot(eitherD_D6signDEGs, aes(x = log2FoldChange.Day6, y = log2FoldChange.Dnmt1)) +

  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  
  geom_point(size = 1, stroke = 0.3, alpha = 0.8, color = "gray60") +  # background points
 # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  
  # Highlight shared significant DEGs
  geom_point(data = sharedD_D6signDEGs,
             aes(x = log2FoldChange.Day6, y = log2FoldChange.Dnmt1),
             shape = 21,
             fill = "#FFC75F",           # Inner color
             color = "#CDAD00",           # Border color
             stroke = 0.3, 
             size = 1, 
             alpha = 0.9) +
   # Linear regression line for shared DEGs
  geom_smooth(data = sharedD_D6signDEGs,
              aes(x = log2FoldChange.Day6, y = log2FoldChange.Dnmt1),
              method = "lm", se = FALSE, color = "gold4", fullrange = TRUE, size=0.5) +

  xlim(-10, 10) +
  ylim(-2, 2) +
  labs(
    title = "Dnmt1-KO vs Day6",
    x = "log2FoldChange(Day6)",
    y = "log2FoldChange(Dnmt1)"
  ) +
  theme_minimal(base_size = 12)

#ggsave("sharedDnmt1-Day6_scatterplot.pdf", width = 6, height = 5, dpi = 300, device = "pdf")
```

### Plot Kmt2a_Day2 overlap

``` r
# 1. Merge both datasets by GeneID
sharedK_D2genes <- merge(
  Kmt2a_results, Day2_results, 
  by = "GeneID", 
  suffixes = c(".Kmt2a", ".Day2")
)

# 2. Genes significant in either dataset
eitherK_D2signDEGs <- sharedK_D2genes[sharedK_D2genes$padj.Kmt2a < 0.05 | sharedK_D2genes$padj.Day2 < 0.05, ]

# 3. Genes significant in both datasets
sharedK_D2signDEGs <- sharedK_D2genes[sharedK_D2genes$padj.Kmt2a < 0.05 & sharedK_D2genes$padj.Day2 < 0.05, ]

p3 <- ggplot(eitherK_D2signDEGs, aes(x = log2FoldChange.Day2, y = log2FoldChange.Kmt2a)) +

    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  
    geom_point(size = 1, stroke = 0.3, alpha = 0.8, color = "gray60") +  # background points
 # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  
  # Highlight shared significant DEGs
  geom_point(data = sharedK_D2signDEGs,
             aes(x = log2FoldChange.Day2, y = log2FoldChange.Kmt2a),
             shape = 21,
             fill = "#FFC75F",           # Inner color
             color = "#CDAD00",           # Border color
             stroke = 0.3, 
             size = 1, 
             alpha = 0.9) +
   # Linear regression line for shared DEGs
  geom_smooth(data = sharedK_D2signDEGs,
              aes(x = log2FoldChange.Day2, y = log2FoldChange.Kmt2a),
              method = "lm", se = FALSE, color = "gold4", fullrange = TRUE, size=0.5) +

  xlim(-10, 10) +
  ylim(-2, 2) +
  labs(
    title = "Kmt2a-KO vs Day2",
    x = "log2FoldChange(Day2)",
    y = "log2FoldChange(Kmt2a)"
  ) +
  theme_minimal(base_size = 12)

#ggsave("sharedKmt2a-Day2_scatterplot.pdf", width = 6, height = 5, dpi = 300, device = "pdf")
```

### Plot Kmt2a_Day6 overlap

``` r
# 1. Merge both datasets by GeneID
sharedK_D6genes <- merge(
  Kmt2a_results, Day6_results, 
  by = "GeneID", 
  suffixes = c(".Kmt2a", ".Day6")
)

# 2. Genes significant in either dataset
eitherK_D6signDEGs <- sharedK_D6genes[sharedK_D6genes$padj.Kmt2a < 0.05 | sharedK_D6genes$padj.Day6 < 0.05, ]

# 3. Genes significant in both datasets
sharedK_D6signDEGs <- sharedK_D6genes[sharedK_D6genes$padj.Kmt2a < 0.05 & sharedK_D6genes$padj.Day6 < 0.05, ]

p4 <- ggplot(eitherK_D6signDEGs, aes(x = log2FoldChange.Day6, y = log2FoldChange.Kmt2a)) +
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  
  geom_point(size = 1, stroke = 0.3, alpha = 0.8, color = "gray60") +  # background points
 # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  
  # Highlight shared significant DEGs
  geom_point(data = sharedK_D6signDEGs,
             aes(x = log2FoldChange.Day6, y = log2FoldChange.Kmt2a),
             shape = 21,
             fill = "#FFC75F",           # Inner color
             color = "#CDAD00",          # Border color
             stroke = 0.3, 
             size = 1, 
             alpha = 0.9) +
   # Linear regression line for shared DEGs
  geom_smooth(data = sharedK_D6signDEGs,
              aes(x = log2FoldChange.Day6, y = log2FoldChange.Kmt2a),
              method = "lm", se = FALSE, color = "gold4", fullrange = TRUE, size=0.5) +

  xlim(-10, 10) +
  ylim(-2, 2) +
  labs(
    title = "Kmt2a-KO vs Day6",
    x = "log2FoldChange(Day6)",
    y = "log2FoldChange(Kmt2a)"
  ) +
  theme_minimal(base_size = 12)

#ggsave("sharedKmt2a-Day6_scatterplot.pdf", width = 6, height = 5, dpi = 300, device = "pdf")
```

## Figure 4C-F

Now all the plots are combined together

``` r
#Now try to plot them all together
(p1 | p2) / (p3 | p4)
```

![](Figure_codes_files/figure-gfm/all_overlaps_together-1.png)<!-- -->

``` r
combined_plot <- (p1 | p2) / (p3 | p4)

#ggsave("scatter_figure_overlaps.pdf", combined_plot,
#       width = 8.27, units = "in")
```

## Figure 5A

``` r
H.matrix <- fread("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/NMF analysis/H.matrix.csv")

colnames(H.matrix) <- c("Kdm5c", "Kdm5b", "Kdm2b", "Kdm1a", "Hdac8", "Kdm6a", "Hdac6", "Ezh2", "Crebbp", "Phf8",
                        "Kmt2d", "Kat6b", "Asxl3", "Prdm8", "Kmt2c", "Kmt2b", "Kat6a", "Ehmt1", "Dnmt3b", "Chd8",
                        "Chd5", "Atrx", "Ash1l", "Taf1", "Setd5", "Rai1", "Mecp2", "Bptf", "Chd3", "Setd1a",
                        "Setd1b", "Setd4", "Kmt2e", "Lbr", "Msl3", "Tcf20", "Tdrd3", "Dnmt1", "Tet1", "Chd1",
                        "Chd4", "Smarca1", "Kmt2a", "Rere", "Srcap", "Kdm3b", "Ctrl")

Dnmt1 <- openxlsx::read.xlsx("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/NMF analysis/EM.genes.and.DNAm.gene.interactions.xlsx", rowNames = TRUE, sheet = 1)
Dnmt3a <- openxlsx::read.xlsx("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/NMF analysis/EM.genes.and.DNAm.gene.interactions.xlsx", rowNames = TRUE, sheet = 2)
Dnmt3b <- openxlsx::read.xlsx("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/NMF analysis/EM.genes.and.DNAm.gene.interactions.xlsx", rowNames = TRUE, sheet = 3)
Tet1 <- openxlsx::read.xlsx("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/NMF analysis/EM.genes.and.DNAm.gene.interactions.xlsx", rowNames = TRUE, sheet = 4)
Tet2 <- openxlsx::read.xlsx("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/NMF analysis/EM.genes.and.DNAm.gene.interactions.xlsx", rowNames = TRUE, sheet = 5)
Tet3 <- openxlsx::read.xlsx("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/NMF analysis/EM.genes.and.DNAm.gene.interactions.xlsx", rowNames = TRUE, sheet = 6)
all.DNAm.machinery.genes <- openxlsx::read.xlsx("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/NMF analysis/EM.genes.and.DNAm.gene.interactions.xlsx", rowNames = TRUE, sheet = 7)


# Combine annotations into a single data frame
combined_annotation <- cbind(Dnmt1, Dnmt3a, Dnmt3b,
                             Tet1, Tet2, Tet3,
                             all.DNAm.machinery.genes)
colnames(combined_annotation) <- gsub("^interaction.with.", "", colnames(combined_annotation))

# Original colors
original_color <- "#bcbddc"

# Darken the color
darker_color <- darken(original_color, amount = 0.20)  # Darken by 20%

Dnmt1_colour = list(Dnmt1 = c(Yes = "#287D8E", No = "#bdbdbd"))
Dnmt3a_colour = list(Dnmt3a = c(Yes = "#287D8E", No = "#bdbdbd"))
Dnmt3b_colour = list(Dnmt3b = c(Yes = "#287D8E", No = "#bdbdbd"))
Tet1_colour = list(Tet1 = c(Yes = "#287D8E", No = "#bdbdbd"))
Tet2_colour = list(Tet2 = c(Yes = "#287D8E", No = "#bdbdbd"))
Tet3_colour = list(Tet3 = c(Yes = "#287D8E", No = "#bdbdbd"))
all_DNAm.machinery.colour = list(all.enzymatic.DNAm.machinery.genes = c(Yes = "#287D8E", No = "#bdbdbd"))


# Combine color schemes
annotation_colors <- c(Dnmt1_colour, Dnmt3a_colour, Dnmt3b_colour,
                       Tet1_colour, Tet2_colour, Tet3_colour,
                       all_DNAm.machinery.colour)


# Adjust scale for heatmap

breaks_seq <- seq(0, 0.4, length.out = 50)
#my_colors <- colorRampPalette(brewer.pal(9, "Purples"))(50)

my_colors <- viridis(50, option = "viridis")

#pdf("NMF_plot.pdf", width = 20, height = 12)

pheatmap(
  H.matrix,
  show_rownames = T,
  annotation_col = combined_annotation,
  annotation_colors = annotation_colors,
  color=my_colors,
  cellwidth = 10, cellheight = 15,
  breaks = breaks_seq,
  border_color = "NA"
)
```

![](Figure_codes_files/figure-gfm/Figure_5A-1.png)<!-- -->

``` r
#dev.off()
```

## Figure 5B

``` r
GO_results <- read.table(
  "C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/NMF analysis/GO_Biological_Process_2023_table.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  fill = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)
GO_results$Term <- gsub("\\(GO:\\d+\\)", "", GO_results$Term)
GO_results$Term <- trimws(GO_results$Term)  # optional: remove trailing spaces

# Create scatterplot
GO_results$log_pval <- -log10(GO_results$'Adjusted P-value')
GO_results$log_OR <- log2(GO_results$'Odds Ratio')
GO_results$Gene_Count <- as.numeric(sub("/.*", "", GO_results$Overlap))
GO_results$Total_Genes <- as.numeric(sub(".*/", "", GO_results$Overlap))
GO_results$Gene_ratio <- GO_results$Gene_Count / GO_results$Total_Genes
GO_filtered <- GO_results[GO_results$Total_Genes >= 10, ]

GO_sign <- GO_filtered[GO_filtered$'Adjusted P-value'<0.01,]

GO_sign_sorted <- GO_sign %>%  arrange(log_OR)

# Convert GO_term to a factor with levels ordered according to Combined.Score
GO_sign_sorted$Term <- factor(GO_sign_sorted$Term, levels = GO_sign_sorted$Term)

top13 <- GO_sign_sorted[order(GO_sign_sorted$log_pval,decreasing = TRUE), ][1:13, ]

ggplot(top13, aes(x = log_OR, y = Term, color = log_pval, size = Gene_ratio)) +
  geom_point() + 
  #scale_color_viridis_c(option = "D") +
  scale_color_gradient(low = "gray80", high = "purple3") +
  #scale_size_continuous(range = c(2, 10)) +
  labs(
    title = "GO-term enrichment analysis",
    x = "Log2 Odds Ratio (Log2OR)",
    y = "GO Term",
    color = "-Log10 p-value",
    size = "Gene ratio"
  ) +
  theme_minimal(base_size = 14) +
 theme(
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16),
    
    # Adjust legend to take up less space
    legend.key.size = unit(0.5, "cm"),  # Smaller legend keys
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.text = element_text(size = 10),  # Adjust legend text size
    legend.position = "right",  # Keep legend on right, but compact
    
    # Adjust plot margins to avoid excess empty space
    plot.margin = margin(10, 10, 10, 10) # Adjust the aspect ratio to make the plot wider
)
```

![](Figure_codes_files/figure-gfm/Figure_5B-1.png)<!-- -->

``` r
#ggsave("NMF_GOresults4.pdf", width = 8, height = 5, dpi = 300, device = "pdf")
```

## Supplementary Figure 6D

``` r
#Load in EnrichR results:

GO_processresults <- read.table("C:/Users/kmoller/OneDrive - Menntaský/Documents/Post-doc application/Kaan/NMF analysis/BioPlanet_2019_table.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  fill = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)
GO_processresults$Term <- gsub("\\(GO:\\d+\\)", "", GO_processresults$Term)
GO_processresults$Term <- trimws(GO_processresults$Term)  # optional: remove trailing spaces

# Create scatterplot
GO_processresults$log_pval <- -log10(GO_processresults$'Adjusted P-value')
GO_processresults$log_OR <- log2(GO_processresults$'Odds Ratio')
GO_processresults$Gene_Count <- as.numeric(sub("/.*", "", GO_processresults$Overlap))
GO_processresults$Total_Genes <- as.numeric(sub(".*/", "", GO_processresults$Overlap))
GO_processresults$Gene_ratio <- (GO_processresults$Gene_Count / GO_processresults$Total_Genes)
GO_process_filtered <- GO_processresults[GO_processresults$Total_Genes >= 10, ]

GO_process_sign <- GO_process_filtered[GO_process_filtered$'Adjusted P-value'<0.01,]

GO_sign_sorted2 <- GO_process_sign %>%  arrange(log_OR)

# Convert GO_term to a factor with levels ordered according to Combined.Score
GO_sign_sorted2$Term <- factor(GO_sign_sorted2$Term, levels = GO_sign_sorted2$Term)

top10 <- GO_sign_sorted2[order(GO_sign_sorted2$log_pval,decreasing = TRUE), ][1:10, ]

ggplot(top10, aes(x = log_OR, y = Term, color = log_pval, size = Gene_ratio)) +
  geom_point() + 
  #scale_color_viridis_c(option = "D") +
  scale_color_gradient(low = "gray80", high = "purple3") +
  #scale_size_continuous(range = c(2, 10)) +
  labs(
    title = "GO-term enrichment analysis",
    x = "Log2 Odds Ratio (Log2OR)",
    y = "GO Term",
    color = "-Log10 p-value",
    size = "Gene ratio"
  ) +
  theme_minimal(base_size = 14) +
 theme(
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16),
    
    # Adjust legend to take up less space
    legend.key.size = unit(0.5, "cm"),  # Smaller legend keys
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.text = element_text(size = 10),  # Adjust legend text size
    legend.position = "right",  # Keep legend on right, but compact
    
    # Adjust plot margins to avoid excess empty space
    plot.margin = margin(10, 10, 10, 10) # Adjust the aspect ratio to make the plot wider
)
```

![](Figure_codes_files/figure-gfm/Supplementary_Figure_6D-1.png)<!-- -->

``` r
#ggsave("NMF_GOprocess.pdf", width = 8, height = 5, dpi = 300, device = "pdf")
```

## Print session info

``` r
sessionInfo()
```

    ## R version 4.5.0 (2025-04-11 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 10 x64 (build 19044)
    ## 
    ## Matrix products: default
    ##   LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] LC_COLLATE=Icelandic_Iceland.utf8  LC_CTYPE=Icelandic_Iceland.utf8   
    ## [3] LC_MONETARY=Icelandic_Iceland.utf8 LC_NUMERIC=C                      
    ## [5] LC_TIME=Icelandic_Iceland.utf8    
    ## 
    ## time zone: Atlantic/Reykjavik
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats4    grid      stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] viridis_0.6.5                            
    ##  [2] viridisLite_0.4.2                        
    ##  [3] colorspace_2.1-2                         
    ##  [4] data.table_1.17.8                        
    ##  [5] RColorBrewer_1.1-3                       
    ##  [6] pheatmap_1.0.13                          
    ##  [7] TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0
    ##  [8] GenomicFeatures_1.60.0                   
    ##  [9] AnnotationDbi_1.70.0                     
    ## [10] Biobase_2.68.0                           
    ## [11] GenomicRanges_1.60.0                     
    ## [12] GenomeInfoDb_1.44.3                      
    ## [13] IRanges_2.42.0                           
    ## [14] S4Vectors_0.46.0                         
    ## [15] BiocGenerics_0.54.1                      
    ## [16] generics_0.1.4                           
    ## [17] biomaRt_2.64.0                           
    ## [18] patchwork_1.3.2                          
    ## [19] VennDiagram_1.7.3                        
    ## [20] futile.logger_1.4.3                      
    ## [21] universalmotif_1.26.2                    
    ## [22] ggrepel_0.9.6                            
    ## [23] ggalt_0.4.0                              
    ## [24] lubridate_1.9.4                          
    ## [25] forcats_1.0.1                            
    ## [26] stringr_1.6.0                            
    ## [27] dplyr_1.1.4                              
    ## [28] purrr_1.2.0                              
    ## [29] readr_2.1.5                              
    ## [30] tidyr_1.3.1                              
    ## [31] tibble_3.3.0                             
    ## [32] ggplot2_4.0.0                            
    ## [33] tidyverse_2.0.0                          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rstudioapi_0.17.1           jsonlite_2.0.0             
    ##  [3] magrittr_2.0.4              farver_2.1.2               
    ##  [5] rmarkdown_2.30              BiocIO_1.18.0              
    ##  [7] vctrs_0.6.5                 memoise_2.0.1              
    ##  [9] Rsamtools_2.24.0            RCurl_1.98-1.17            
    ## [11] htmltools_0.5.8.1           S4Arrays_1.8.1             
    ## [13] progress_1.2.3              lambda.r_1.2.4             
    ## [15] curl_6.3.0                  SparseArray_1.8.0          
    ## [17] proj4_1.0-15                KernSmooth_2.23-26         
    ## [19] httr2_1.2.1                 futile.options_1.0.1       
    ## [21] cachem_1.1.0                GenomicAlignments_1.44.0   
    ## [23] lifecycle_1.0.4             pkgconfig_2.0.3            
    ## [25] Matrix_1.7-3                R6_2.6.1                   
    ## [27] fastmap_1.2.0               GenomeInfoDbData_1.2.14    
    ## [29] MatrixGenerics_1.20.0       digest_0.6.37              
    ## [31] RSQLite_2.4.3               filelock_1.0.3             
    ## [33] labeling_0.4.3              timechange_0.3.0           
    ## [35] httr_1.4.7                  abind_1.4-8                
    ## [37] mgcv_1.9-4                  compiler_4.5.0             
    ## [39] bit64_4.6.0-1               withr_3.0.2                
    ## [41] S7_0.2.0                    BiocParallel_1.42.1        
    ## [43] DBI_1.2.3                   Rttf2pt1_1.3.14            
    ## [45] maps_3.4.3                  MASS_7.3-65                
    ## [47] rappdirs_0.3.3              DelayedArray_0.34.1        
    ## [49] rjson_0.2.23                tools_4.5.0                
    ## [51] zip_2.3.3                   extrafontdb_1.1            
    ## [53] glue_1.8.0                  restfulr_0.0.15            
    ## [55] nlme_3.1-168                gtable_0.3.6               
    ## [57] tzdb_0.5.0                  hms_1.1.4                  
    ## [59] xml2_1.4.1                  XVector_0.48.0             
    ## [61] pillar_1.11.1               vroom_1.6.6                
    ## [63] splines_4.5.0               BiocFileCache_2.16.2       
    ## [65] lattice_0.22-7              rtracklayer_1.68.0         
    ## [67] bit_4.6.0                   tidyselect_1.2.1           
    ## [69] Biostrings_2.76.0           knitr_1.50                 
    ## [71] gridExtra_2.3               SummarizedExperiment_1.38.1
    ## [73] xfun_0.52                   matrixStats_1.5.0          
    ## [75] stringi_1.8.7               UCSC.utils_1.4.0           
    ## [77] yaml_2.3.10                 evaluate_1.0.5             
    ## [79] codetools_0.2-20            extrafont_0.20             
    ## [81] cli_3.6.5                   ash_1.0-15                 
    ## [83] dichromat_2.0-0.1           Rcpp_1.1.0                 
    ## [85] dbplyr_2.5.1                png_0.1-8                  
    ## [87] XML_3.99-0.18               parallel_4.5.0             
    ## [89] blob_1.2.4                  prettyunits_1.2.0          
    ## [91] bitops_1.0-9                scales_1.4.0               
    ## [93] openxlsx_4.2.8.1            crayon_1.5.3               
    ## [95] rlang_1.1.6                 KEGGREST_1.48.1            
    ## [97] formatR_1.14
