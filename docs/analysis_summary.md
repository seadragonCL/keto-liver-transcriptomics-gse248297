# Analysis Summary

## Topic

**Long-term ketogenic diet and liver metabolic pathway rewiring**

Dataset: GEO **GSE248297**, liver RNA-seq from male C57BL/6J mice after 9 months on low-fat control, high-fat, or ketogenic diet.

## Workflow

- Downloaded raw count matrix and GEO series metadata into `data/raw`.
- Excluded contaminated control sample `20462X1`, following the GEO note.
- Ran DESeq2 with design `~ diet`.
- Primary contrast: ketogenic diet vs low-fat control.
- Added mouse gene annotations with `org.Mm.eg.db`.
- Ran GO biological-process over-representation analysis for higher/lower keto genes.
- Ran ranked GO GSEA and MSigDB Hallmark GSEA.

## Key Results

- Samples downloaded: 15
- Samples analyzed after excluding contaminated sample: 14
- Genes in raw count matrix: 56,980
- Genes after low-count filter: 15,063
- Keto vs control FDR < 0.05 genes: 2,432
- Keto vs control FDR < 0.05 and abs(log2FC) >= 0.58 genes: 1,656
- Higher in keto by threshold: 998 mapped Entrez genes
- Lower in keto by threshold: 538 mapped Entrez genes

## Top Differential Genes

Higher in ketogenic diet:

- `Cidec`, `Vnn1`, `Gpat3`, `Osbpl3`, `Cyp4a14`, `Aldh3a2`, `Cd36`, `Pex11a`

Lower in ketogenic diet:

- `Gpi1`, `Pklr`, `Adgrf1`, `Gm2a`, `Enpep`, `Acly`, `Fabp5`

## Pathway Interpretation

The strongest signal is a hepatic shift toward lipid utilization:

- GO over-representation among genes higher in keto: fatty acid metabolic process, fatty acid beta-oxidation, fatty acid catabolic process, fatty acid oxidation, lipid oxidation.
- Hallmark GSEA higher in keto: fatty acid metabolism, adipogenesis, peroxisome.
- GO GSEA higher in keto: fatty acid beta-oxidation and fatty acid oxidation.

The lower-in-keto side is consistent with reduced carbohydrate/glycolytic and biosynthetic activity:

- Lower genes include `Gpi1` and `Pklr`, both tied to glucose/glycolytic metabolism.
- GO over-representation among genes lower in keto includes carbohydrate, monosaccharide, hexose, organic acid biosynthetic, and carboxylic acid biosynthetic processes.
- Hallmark GSEA lower in keto includes unfolded protein response, MTORC1 signaling, and MYC targets.

## Main Takeaway

This dataset supports a clear transcriptomic signature of chronic ketogenic feeding in mouse liver: increased fatty-acid oxidation/peroxisomal lipid metabolism and reduced carbohydrate/biosynthetic expression programs. This makes the project a good basis for a bioinformatics research question such as:

**"How does long-term ketogenic diet remodel hepatic lipid oxidation and nutrient-sensing pathways at the transcriptomic level?"**

## Main Outputs

- Differential expression: `results/tables/deseq2_keto_vs_control.csv`
- Hallmark GSEA: `results/tables/gsea_hallmark_keto_vs_control.csv`
- GO GSEA: `results/tables/gsea_go_bp_keto_vs_control.csv`
- GO ORA higher/lower in keto: `results/tables/go_bp_enrichment_higher_in_keto.csv`, `results/tables/go_bp_enrichment_lower_in_keto.csv`
- Figures: PCA, volcano, heatmap, GO dotplots, Hallmark barplot in `results/figures`
