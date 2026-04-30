# Keto Bioinformatics Project: GSE248297

This RStudio project analyzes liver RNA-seq raw counts from GEO accession **GSE248297** to study how long-term ketogenic diet changes hepatic gene expression and pathways in mice.

## Folder Layout

- `data/raw`: downloaded GEO files
- `data/processed`: cleaned metadata and normalized expression objects
- `scripts`: R analysis scripts
- `results/tables`: differential expression and pathway tables
- `results/figures`: PCA, volcano, heatmap, and pathway figures
- `docs`: topic selection and data-source notes

## Run

Open `keto_bioinformatics_gse248297.Rproj` in RStudio, then run:

```r
source("scripts/01_deseq2_pathway_analysis.R")
```

The script uses DESeq2 for differential expression and clusterProfiler/fgsea-style workflows for GO and Hallmark pathway analysis.
