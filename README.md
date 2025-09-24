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
# Analysis of Dnmt1 knockout methylation data



