# Dowloading expression and methylation data from GEO

# Analysis of Dnmt1 knockout expression data

Dnmt1 knockout short-read RNA sequencing data is deposit at GEO under accession id GSE308011.

For instance, raw gene-level counts of Dnmt1_1 sample is in GSM9235598_Dnmt1_1.raw.counts.txt.gz file. First five rows of that data look like this:

| Gene symbol | Dnmt1_1 |
|----------|----------|
| Zc3h14    | 2144   |
| Troap    | 320   |
| Clca2    | 0   |
| 1810013L24Rik    | 1803   |
| Srsf7    | 2790   |

### Once other samples are downloaded, a raw count matrix can be generated in R, where rownames are gene symbols but columns are samples.

```{R}
library("DESeq2")
library("edgeR")
```

```{R}
files <- list.files(path = "/path/to/raw.counts", pattern = "\\.txt$", full.names = TRUE)
```

```{R}
tables <- lapply(files, function(f) {
  read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
})
```

```{R}
merged.raw.counts <- Reduce(function(x, y) merge(x, y, by = "Gene.symbol"), tables)
```

```{R}
head(merged.raw.counts)
```

| Gene.symbol   | Dnmt1_1 | Dnmt1_2 | Dnmt1_3 | Dnmt1_4 | Ctr3_1 | Ctr3_2 | Ctr3_3 | Ctr3_4 |
|---------------|---------|---------|---------|---------|--------|--------|--------|--------|
| 0610005C13Rik |      14 |      17 |      16 |      11 |     15 |     19 |     20 |     14 |
| 0610009B22Rik |     291 |     353 |     280 |     291 |    382 |    397 |    392 |    328 |
| 0610009E02Rik |      17 |      24 |      24 |      16 |     21 |     22 |     21 |     11 |
| 0610009L18Rik |      91 |      95 |     105 |     107 |    103 |    121 |    100 |     77 |
| 0610010F05Rik |    1155 |    1205 |    1122 |    1125 |   1444 |   1426 |   1776 |   1148 |
| 0610010K14Rik |     817 |    1133 |     913 |    1045 |   1106 |   1619 |   1421 |   1056 |

### Calculate library size of each sample

```{R}
library_size <- colSums(merged.raw.counts)
```

### Create a meta data annotating samples in the expression matrix respect to column order of samples

```{R}
group <- factor(c("Dnmt1", "Dnmt1", "Dnmt1", "Dnmt1", "Control", "Control", "Control", "Control"))
group <- relevel(group, "Control")
```

### Do CPM normalization by taking account library size in the expression matrix through cpm() function of the R package edgeR

```{R}
CPM_normalized_expression_values <- cpm(merged.raw.counts, lib.size = library_size, group = group)
```

### Keep genes with > 1 CPM

```{R}
CPM_cutoff <- 1
filtered_expression_data <- CPM_normalized_expression_values[apply(CPM_normalized_expression_values,1,function(x){max(x)}) > CPM_cutoff,]
```

### Obtain raw expression values of genes with > 1 CPM for DESeq2 DE analysis

```{R}
filtered_expression_data_2 <- subset(merged.raw.counts, rownames(merged.raw.counts) %in% rownames(filtered_expression_data))
```

### Sample-trait annotation

```{R}
sampleInfo <- data.frame(group, row.names = colnames(filtered_expression_data_2))
```

### Differential expression analysis

```{R}
dds <- DESeqDataSetFromMatrix(countData = filtered_expression_data_2, colData = sampleInfo, design = ~ group)
dds$group = relevel(dds$group, "Control")
dds <- DESeq(dds)
DE_results <- results(dds)
```

### Get differentially expressed genes with < 0.05 padj

```{R}
DEGs <- subset(DE_results, padj < 0.05)
```
# Analysis of Dnmt1 knockout (unphased) methylation data

```{R}
# loading libraries
library(locfit)
library(dmrseq)
library(BiocParallel)
library(rtracklayer)
```

Dnmt1 knockout long-read WGS sequencing data is deposit at GEO under accession id GSE308669.

| chr   | start   | end     | modified_base_code | score | strand | tstart  | tend    | color   | Nvalid_cov  | Nmod_/_Nvalid_cov  | Nmod | Ncanonical | Nother_mod | Ndelete | Nfail | Ndiff | Nnocall |
|-------|---------|---------|------|-----|-----|---------|---------|---------|----|--------|----|----|----|----|----|----|----|
| chr10 | 3100021 | 3100022 | h    | 10  | .   | 3100021 | 3100022 | 255,0,0 | 10 | 0.00   | 0  | 3  | 7  | 0  | 0  | 0  | 0  |
| chr10 | 3100021 | 3100022 | m    | 10  | .   | 3100021 | 3100022 | 255,0,0 | 10 | 70.00  | 7  | 3  | 0  | 0  | 0  | 0  | 0  |
| chr10 | 3100160 | 3100161 | h    | 9   | .   | 3100160 | 3100161 | 255,0,0 | 9  | 11.11  | 1  | 2  | 6  | 0  | 1  | 0  | 0  |
| chr10 | 3100160 | 3100161 | m    | 9   | .   | 3100160 | 3100161 | 255,0,0 | 9  | 66.67  | 6  | 2  | 1  | 0  | 1  | 0  | 0  |
| chr10 | 3100241 | 3100242 | h    | 11  | .   | 3100241 | 3100242 | 255,0,0 | 11 | 0.00   | 0  | 1  | 10 | 0  | 0  | 0  | 0  |
| chr10 | 3100241 | 3100242 | m    | 11  | .   | 3100241 | 3100242 | 255,0,0 | 11 | 90.91  | 10 | 1  | 0  | 0  | 0  | 0  | 0  |
| chr10 | 3100272 | 3100273 | h    | 9   | .   | 3100272 | 3100273 | 255,0,0 | 9  | 0.00   | 0  | 0  | 9  | 0  | 1  | 0  | 1  |
| chr10 | 3100272 | 3100273 | m    | 9   | .   | 3100272 | 3100273 | 255,0,0 | 9  | 100.00 | 9  | 0  | 0  | 0  | 1  | 0  | 1  |
| chr10 | 3100440 | 3100441 | h    | 9   | .   | 3100440 | 3100441 | 255,0,0 | 9  | 11.11  | 1  | 4  | 4  | 0  | 2  | 0  | 0  |
| chr10 | 3100440 | 3100441 | m    | 9   | .   | 3100440 | 3100441 | 255,0,0 | 9  | 44.44  | 4  | 4  | 1  | 0  | 2  | 0  | 0  |

```{R}
# Generation of BSseq object for given sample
Dnmt1 <- fread("Dnmt1.bed", nThread = 32)
colnames(Dnmt1) <- c("chrom", "start", "end", "modified_base_code", "score", "strand", "tstart", "tend", "color", "Nvalid_cov", "Nmod_/_Nvalid_cov", "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall")
Dnmt1 <- subset(Dnmt1, Dnmt1$modified_base_code == "m")
Dnmt1.BSSeq.object <- BSseq(chr = as.character(Dnmt1$chrom), pos = Dnmt1$start +1, M = matrix(Dnmt1$Nmod), Cov = matrix(Dnmt1$Nmod + Dnmt1$Ncanonical), sampleNames = c("Dnmt1"))
sampleNames(Dnmt1.BSSeq.object) <- "Dnmt1"
```




