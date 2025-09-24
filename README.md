# Dowloading expression and methylation data from GEO

Command-line tool for estimation of knockout efficiency in long-read whole genome sequencing. 

This tool is able to process multiple BAM files for multiple sgRNAs to detect CRISPR events around cut site of Cas9 in target region.

### Dependencies
In order to run samCRISPR the following software components and packages are required:
- SAMtools >=1.21
- R computing environment >=4.4.2
- Linux environment

### Installation
The samCRISPR script can be directly downloaded from this repository:
```{bash}
wget https://github.com/kaanokay/samCRISPR/blob/main/script/samCRISPR.sh
```

### Change the execution permissions:
```{bash}
chmod u+x samCRISPR.sh
```

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

### Once other samples are downloaded, a raw count matrix can be generated in R:

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


# Analysis of Dnmt1 knockout methylation data



