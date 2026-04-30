# Data Sources

## Primary Dataset

- GEO accession: GSE248297
- Title: Liver RNA-seq from male C57Bl/6J on a ketogenic diet
- Organism: Mus musculus
- Platform: Illumina NovaSeq 6000
- Public release: 2025-08-18
- GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248297
- Raw count file: `GSE248297_20462R_ALF_raw_counts.txt.gz`
- Metadata file: `GSE248297_series_matrix.txt.gz`

## Design Notes From GEO

- Male C57BL/6J mice were placed on one of three diets for 9 months:
  - 10% low-fat control diet
  - 60% high-fat diet
  - 89.9% fat ketogenic diet
- Tissue: liver
- One control sample, `20462X1_Liver_Control_ALF_replicate1`, was marked as contaminated and excluded from analysis.

## Analysis Focus

Primary contrast:

- Ketogenic diet vs low-fat control diet

Contextual contrast:

- High-fat diet vs low-fat control diet

The main output is differential expression plus GO and Hallmark pathway enrichment/GSEA where possible.
