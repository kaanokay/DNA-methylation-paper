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

# Analysis of Dnmt1 knockout methylation data
