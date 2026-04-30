# Research Article Draft for Genes & Nutrition

## Title Page

**Long-term ketogenic diet reprograms the mouse liver transcriptome toward fatty acid oxidation and peroxisomal metabolism: a public RNA-seq reanalysis**

**Short title:** Keto diet and liver transcriptomic reprogramming

**Authors:**  
Yixuan Liang1*; [Author 2]2; [Author 3]3

**Affiliations:**  
1 Department of Neuroscience, Johns Hopkins University, Baltimore, MD, USA  
2 [Department], [Institution], [City], [Country]  
3 [Department], [Institution], [City], [Country]

**Corresponding author:**  
Yixuan Liang, [email address]

## Abstract

### Background

Ketogenic diets are increasingly studied as nutritional interventions for neurological, psychiatric, and metabolic conditions, yet their long-term effects on tissue-specific gene regulation remain incompletely characterized. Because the liver coordinates ketogenesis, fatty acid oxidation, lipid handling, and nutrient sensing, hepatic transcriptomic data can clarify molecular adaptations to chronic ketogenic feeding. We reanalyzed a public mouse liver RNA-seq dataset to identify differential gene expression and pathway-level signatures associated with long-term ketogenic diet exposure.

### Results

Raw count data from GEO accession GSE248297 were analyzed for male C57BL/6J mice fed control, high-fat, or ketogenic diets for 9 months. One contaminated control sample flagged by the data submitters was excluded, leaving 14 samples for analysis. After low-count filtering, 15,063 genes were tested with DESeq2. In the ketogenic diet versus control contrast, 2,432 genes were differentially expressed at false discovery rate (FDR) < 0.05. Applying an additional absolute log2 fold-change threshold of 0.58 identified 1,656 genes, including 1,085 higher and 571 lower in ketogenic diet liver. Genes higher in ketogenic diet included Cidec, Vnn1, Gpat3, Cyp4a14, Cd36, and Pex11a, whereas genes lower in ketogenic diet included Gpi1, Pklr, Acly, and Fabp5. Gene Ontology enrichment among genes higher in ketogenic diet highlighted fatty acid metabolic process, fatty acid beta-oxidation, fatty acid catabolic process, fatty acid oxidation, and lipid oxidation. Hallmark gene set enrichment analysis showed positive enrichment of fatty acid metabolism, peroxisome, adipogenesis, bile acid metabolism, and xenobiotic metabolism, and negative enrichment of unfolded protein response, MTORC1 signaling, PI3K-AKT-MTOR signaling, and MYC target programs.

### Conclusions

Long-term ketogenic diet exposure is associated with a strong hepatic transcriptomic shift toward lipid catabolism, fatty acid oxidation, and peroxisomal metabolism, together with lower expression of carbohydrate and biosynthetic programs. These findings provide a reproducible pathway-level framework for understanding hepatic molecular adaptation to chronic ketogenic feeding and for generating mechanistic hypotheses relevant to nutritional ketosis.

**Keywords:** ketogenic diet; nutrigenomics; RNA-seq; liver; fatty acid oxidation; peroxisome; DESeq2; pathway analysis; mouse

## Background

Ketogenic diets are high-fat, very-low-carbohydrate dietary patterns that promote nutritional ketosis and alter systemic substrate availability. Historically used in refractory epilepsy, ketogenic diets are now being investigated in obesity, type 2 diabetes, neurological disorders, and psychiatric or behavioral conditions [1,2]. Recent clinical pilot work, including a dietitian-supported ketogenic intervention for binge-eating disorder, suggests feasibility and potential symptom improvement in selected populations, while also emphasizing the need for controlled studies and mechanistic understanding [3]. The biological effects of ketogenic diets are therefore relevant not only to weight regulation, but also to nutrient sensing, intermediary metabolism, inflammation, and tissue-specific gene regulation.

The liver is central to adaptation during carbohydrate restriction because it coordinates fatty acid uptake, beta-oxidation, ketone body production, lipoprotein metabolism, and endocrine responses to fasting-like nutrient states. Earlier mouse work showed that ketogenic feeding can induce a metabolic state distinct from both chow feeding and high-fat/high-sucrose feeding, with marked changes in hepatic gene expression and lipid metabolism [4]. However, the molecular effects of chronic ketogenic feeding may vary with diet composition, duration, age, sex, tissue, and analytic strategy. Publicly available RNA-seq datasets now allow independent reanalysis using current statistical and pathway-enrichment approaches.

The present study reanalyzed GEO accession GSE248297, a public liver RNA-seq dataset from male C57BL/6J mice maintained for 9 months on low-fat control, high-fat, or ketogenic diets [5]. We focused on the ketogenic diet versus control contrast and used differential expression, Gene Ontology enrichment, and Hallmark gene set enrichment analysis to define hepatic pathways associated with long-term ketogenic feeding.

## Methods

### Study design and data source

This study is a secondary bioinformatic reanalysis of publicly available RNA-seq data deposited in the Gene Expression Omnibus (GEO), an archive for functional genomics datasets [6]. The primary dataset was GEO accession GSE248297, titled "Liver RNA-seq from male C57Bl/6J on a ketogenic diet" [5]. According to the GEO record, male C57BL/6J mice were placed on one of three diets for 9 months: a 10% fat low-fat control diet, a 60% fat high-fat diet, or an 89.9% fat ketogenic diet. Liver tissue was collected at sacrifice, RNA was extracted, and sequencing was performed using the Illumina NovaSeq 6000 platform. The GEO record notes that sample 20462X1_Liver_Control_ALF_replicate1 was contaminated and excluded from analysis; we followed this exclusion.

### Data acquisition and preprocessing

The raw count matrix file `GSE248297_20462R_ALF_raw_counts.txt.gz` and the GEO series matrix metadata file were downloaded from GEO. The raw matrix contained 56,980 Ensembl gene identifiers across 15 samples. After excluding the contaminated control sample, 14 samples remained: 4 control, 5 high-fat diet, and 5 ketogenic diet samples. Genes were retained for differential expression testing if they had at least 10 counts in at least 3 samples, leaving 15,063 genes.

### Differential expression analysis

Differential expression was performed in R version 4.4.1 using DESeq2 [7]. A DESeqDataSet was constructed from the raw count matrix and sample metadata using the design formula `~ diet`. Size-factor estimation, dispersion estimation, model fitting, and Wald testing were performed with DESeq2 defaults. The primary contrast was ketogenic diet versus control diet. A contextual high-fat diet versus control diet contrast was also computed but was not the primary focus of this manuscript. Genes were considered differentially expressed at FDR < 0.05 after Benjamini-Hochberg adjustment [8]. For a more interpretable effect-size subset, genes with FDR < 0.05 and absolute log2 fold-change >= 0.58 were also reported, corresponding approximately to a 1.5-fold expression difference.

### Gene annotation

Mouse Ensembl gene identifiers were annotated to gene symbols, Entrez Gene identifiers, and gene names using `org.Mm.eg.db` through the Bioconductor `AnnotationDbi` interface. Where one Ensembl identifier mapped to multiple annotations, the first returned mapping was used for descriptive reporting. Genes without Entrez identifiers were retained in differential expression tables but omitted from Entrez-based pathway analyses.

### Pathway and gene set analysis

Gene Ontology (GO) biological-process over-representation analysis was performed with clusterProfiler [9] using Entrez Gene identifiers for genes higher or lower in ketogenic diet at FDR < 0.05 and absolute log2 fold-change >= 0.58. The background universe was defined as all tested genes with Entrez identifiers. GO terms were interpreted using the Gene Ontology knowledgebase [10].

Ranked gene set enrichment analysis (GSEA) was performed using the DESeq2 Wald statistic as the ranking metric. GO biological-process GSEA was run with clusterProfiler. Hallmark GSEA used mouse-native MSigDB Hallmark gene sets accessed with `msigdbr` and analyzed with clusterProfiler's GSEA implementation. GSEA follows the principle of testing whether predefined gene sets are non-randomly distributed toward the top or bottom of a ranked gene list [11]. Hallmark gene sets were used because they reduce redundancy and summarize coherent biological states and processes [12].

### Visualization and reproducibility

Variance-stabilized counts from DESeq2 were used for principal component analysis (PCA) and heatmap visualization. Volcano plots, PCA plots, GO dotplots, and Hallmark GSEA barplots were generated in R using `ggplot2`, `ggrepel`, and `pheatmap`. The analysis script, processed tables, and generated figures are organized in a reproducible RStudio project. Before submission, the code and derived results should be deposited in a public repository such as GitHub and archived with a DOI through Zenodo.

## Results

### Sample inclusion and global expression structure

The downloaded GEO count matrix included 15 liver RNA-seq samples. Following the GEO submitters' annotation, one contaminated control sample was excluded, leaving 14 samples for analysis. After filtering low-count genes, 15,063 genes were retained. PCA of variance-stabilized expression values showed that the first two principal components explained 37.8% and 18.9% of the expression variance, respectively, and provided an initial visualization of diet-associated transcriptomic structure (Fig. 1).

### Long-term ketogenic diet alters hepatic gene expression

In the ketogenic diet versus control diet contrast, 2,432 genes were differentially expressed at FDR < 0.05. Applying an additional absolute log2 fold-change threshold of 0.58 yielded 1,656 genes, including 1,085 genes higher and 571 genes lower in ketogenic diet liver. The strongest genes higher in ketogenic diet included Cidec, Vnn1, Gpat3, Osbpl3, Cyp4a14, Aldh3a2, Cd36, and Pex11a. Several of these genes are consistent with hepatic lipid handling, fatty acid metabolism, peroxisome biology, or oxidative metabolism. Genes lower in ketogenic diet included Gpi1, Pklr, Acly, and Fabp5, consistent with reduced expression of selected carbohydrate metabolism and lipogenic/biosynthetic components (Fig. 2; Table 1).

**Table 1. Top differential genes in ketogenic diet versus control liver**

| Direction | Gene symbol | log2 fold-change | FDR | Gene name |
|---|---:|---:|---:|---|
| Higher in keto | Cidec | 4.94 | 7.14e-45 | cell death-inducing DFFA-like effector c |
| Higher in keto | Vnn1 | 3.88 | 1.34e-41 | vanin 1 |
| Higher in keto | Gpat3 | 3.97 | 2.38e-40 | glycerol-3-phosphate acyltransferase 3 |
| Higher in keto | Osbpl3 | 3.37 | 5.23e-34 | oxysterol binding protein-like 3 |
| Higher in keto | Cyp4a14 | 8.01 | 5.74e-32 | cytochrome P450 family 4 subfamily a polypeptide 14 |
| Higher in keto | Cd36 | 2.96 | 1.96e-26 | CD36 molecule |
| Higher in keto | Pex11a | 1.99 | 3.16e-25 | peroxisomal biogenesis factor 11 alpha |
| Lower in keto | Gpi1 | -1.02 | 1.26e-35 | glucose-6-phosphate isomerase 1 |
| Lower in keto | Pklr | -2.22 | 1.56e-32 | pyruvate kinase liver and red blood cell |
| Lower in keto | Acly | -2.77 | 1.50e-16 | ATP citrate lyase |
| Lower in keto | Fabp5 | -2.90 | 7.84e-16 | fatty acid binding protein 5 |

### Genes higher in ketogenic diet are enriched for fatty acid oxidation and lipid catabolism

GO biological-process over-representation analysis of genes higher in ketogenic diet revealed strong enrichment for lipid and fatty acid metabolic processes. Top enriched terms included fatty acid metabolic process, fatty acid beta-oxidation, fatty acid catabolic process, monocarboxylic acid catabolic process, fatty acid oxidation, long-chain fatty acid metabolic process, lipid oxidation, and lipid catabolic process (Table 2; Fig. 3). These results indicate that chronic ketogenic feeding increases hepatic transcriptional programs associated with lipid utilization and oxidative substrate metabolism.

**Table 2. Selected GO biological-process terms enriched among genes higher in ketogenic diet**

| GO ID | Description | FDR | Count | Gene ratio |
|---|---|---:|---:|---|
| GO:0006631 | fatty acid metabolic process | 2.75e-14 | 82 | 82/957 |
| GO:0006635 | fatty acid beta-oxidation | 3.30e-10 | 28 | 28/957 |
| GO:0009062 | fatty acid catabolic process | 3.82e-10 | 33 | 33/957 |
| GO:0072329 | monocarboxylic acid catabolic process | 1.99e-09 | 36 | 36/957 |
| GO:0019395 | fatty acid oxidation | 2.00e-09 | 33 | 33/957 |
| GO:0001676 | long-chain fatty acid metabolic process | 2.00e-09 | 32 | 32/957 |
| GO:0034440 | lipid oxidation | 2.72e-09 | 33 | 33/957 |
| GO:0016042 | lipid catabolic process | 1.96e-08 | 54 | 54/957 |

### Genes lower in ketogenic diet are enriched for carbohydrate and biosynthetic processes

Genes lower in ketogenic diet were enriched for sulfur compound metabolic process, carboxylic acid biosynthetic process, organic acid biosynthetic process, monosaccharide metabolic process, carbohydrate metabolic process, hexose metabolic process, steroid metabolic process, and xenobiotic metabolic process (Fig. 3). The presence of Gpi1 and Pklr among the most significantly lower genes supports a reduction in selected glycolytic or carbohydrate-handling expression programs under chronic ketogenic feeding.

### Ranked gene set analysis confirms lipid oxidative remodeling

Ranked GO GSEA confirmed pathway-level enrichment of fatty acid beta-oxidation and fatty acid oxidation among genes higher in ketogenic diet, whereas translation- and ribosome-related GO terms showed negative enrichment. Hallmark GSEA further identified positive enrichment of fatty acid metabolism, epithelial-mesenchymal transition, peroxisome, adipogenesis, TNFA signaling via NFKB, bile acid metabolism, interferon responses, and xenobiotic metabolism. Hallmarks with negative enrichment included unfolded protein response, MTORC1 signaling, PI3K-AKT-MTOR signaling, MYC targets, and hypoxia (Fig. 4). The strongest Hallmark term was fatty acid metabolism (normalized enrichment score [NES] = 2.85, FDR = 2.5e-09), followed by peroxisome (NES = 2.19, FDR = 1.43e-06) and adipogenesis (NES = 1.98, FDR = 4.26e-07).

## Discussion

This reanalysis identifies a clear hepatic transcriptomic signature of long-term ketogenic diet exposure in male mice. The dominant pattern is increased expression of lipid catabolic and oxidative programs, including fatty acid beta-oxidation, fatty acid oxidation, peroxisomal biology, and broader fatty acid metabolism. In parallel, several carbohydrate and biosynthetic expression programs were lower in ketogenic diet liver. These findings align with the expected metabolic logic of carbohydrate restriction: reduced reliance on glucose-derived carbon and increased hepatic capacity for fatty acid utilization and lipid-derived energy production.

The gene-level results are biologically coherent. Cyp4a14 and Cd36 are linked to fatty acid metabolism and lipid uptake/handling; Pex11a is involved in peroxisome biology; Gpat3 and Cidec reflect lipid storage or remodeling programs; and Aldh3a2 participates in fatty aldehyde metabolism. Conversely, lower expression of Gpi1 and Pklr is consistent with reduced glycolytic flux or carbohydrate-processing transcriptional programs, while lower Acly may reflect altered cytosolic acetyl-CoA production from citrate and reduced lipogenic biosynthesis. Together, these single-gene findings and pathway-level results support the conclusion that long-term ketogenic feeding produces a coordinated hepatic metabolic switch.

The findings extend earlier animal work showing that ketogenic diets induce a metabolic state distinct from conventional high-fat feeding [4]. They also provide a mechanistic bridge to the growing clinical literature on ketogenic interventions. Clinical pilot studies, including recent work in binge-eating disorder [3], are primarily designed to test feasibility, acceptability, symptoms, safety, or metabolic outcomes; they do not directly reveal tissue-level transcriptional mechanisms. The present liver transcriptomic analysis therefore adds molecular context rather than clinical efficacy evidence. It suggests that long-term nutritional ketosis is accompanied by durable hepatic pathway remodeling that should be considered when interpreting systemic metabolic effects of ketogenic interventions.

The negative enrichment of MTORC1 signaling, PI3K-AKT-MTOR signaling, MYC targets, and unfolded protein response is also notable. Reduced mTOR-related and MYC-related signatures may reflect altered nutrient sensing and biosynthetic demand under carbohydrate restriction. The lower unfolded protein response signal contrasts with some settings in which high-fat or metabolic stress diets induce endoplasmic reticulum stress, suggesting that ketogenic diet effects may not be equivalent to high-fat diet effects despite high fat content. This distinction is important for nutrigenomics, where macronutrient composition, carbohydrate availability, fat type, protein content, feeding duration, and energy balance can produce divergent molecular phenotypes.

Several limitations should be emphasized. First, this was a secondary analysis of a public dataset and cannot establish causality beyond the original experimental design. Second, only male mice were included, limiting generalizability to females. Third, the analysis used bulk liver RNA-seq, so cell-type-specific effects cannot be resolved. Fourth, transcript abundance does not necessarily indicate protein abundance, enzyme activity, metabolite flux, or ketone production. Fifth, the dataset included a long-term 9-month exposure; the results may not generalize to short-term ketogenic interventions. Finally, pathway enrichment depends on gene annotation and gene-set definitions, and should be interpreted as hypothesis-generating.

Future work should integrate transcriptomics with hepatic lipidomics, metabolomics, ketone measurements, insulin/glucose phenotypes, and histology. Comparing ketogenic diet with multiple high-fat non-ketogenic diets would help distinguish effects of carbohydrate restriction from high fat exposure. Additional validation across independent liver datasets, time courses, sexes, and tissues would strengthen inference. In clinical translational contexts, paired blood transcriptomic, metabolomic, and endocrine data may help determine whether accessible biomarkers reflect hepatic nutrient-sensing and lipid oxidation pathways.

## Conclusions

Long-term ketogenic diet exposure in male mouse liver is associated with broad transcriptomic remodeling characterized by increased fatty acid oxidation, lipid catabolism, and peroxisomal metabolism, together with lower carbohydrate and biosynthetic programs. These results provide a reproducible nutrigenomic framework for studying hepatic adaptation to chronic ketogenic feeding and generate testable hypotheses for future multi-omics and translational studies.

## List of abbreviations

BED: binge-eating disorder  
FDR: false discovery rate  
GEO: Gene Expression Omnibus  
GO: Gene Ontology  
GSEA: gene set enrichment analysis  
KD: ketogenic diet  
NES: normalized enrichment score  
PCA: principal component analysis  
RNA-seq: RNA sequencing  
VST: variance-stabilizing transformation

## Declarations

### Ethics approval and consent to participate

Not applicable. This study reanalyzed publicly available, aggregate mouse RNA-seq data from GEO and did not involve new animal experiments, human participants, human tissue, or identifiable human data.

### Consent for publication

Not applicable.

### Availability of data and materials

The dataset supporting the conclusions of this article is available in the Gene Expression Omnibus repository under accession GSE248297: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248297. The raw count matrix is available from the GEO supplementary files. Analysis code and derived result tables should be deposited in a public repository before submission, for example GitHub with a Zenodo archived DOI: [repository URL to be added].

### Competing interests

The authors declare that they have no competing interests. [Revise as needed before submission.]

### Funding

This secondary analysis received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors. [Revise as needed before submission.]

### Authors' contributions

YL conceived the study, curated the dataset, performed the bioinformatic analysis, interpreted the results, and drafted the manuscript. [Author 2 initials] contributed to biological interpretation and manuscript revision. [Author 3 initials] contributed to analysis review and manuscript revision. All authors read and approved the final manuscript. [Revise as needed before submission.]

### Acknowledgements

The authors thank the investigators who generated and deposited the GSE248297 dataset in GEO, enabling reproducible secondary analysis.

### Authors' information

Not applicable.

### Use of generative artificial intelligence tools

[If applicable at submission: During manuscript preparation, the authors used OpenAI ChatGPT/Codex to assist with drafting, code organization, and language refinement. The authors reviewed and verified all scientific content, analyses, citations, and interpretations, and accept full responsibility for the final manuscript.]

## Figure Legends

**Fig. 1 PCA of variance-stabilized liver RNA-seq expression.** Principal component analysis of variance-stabilized counts after excluding one contaminated control sample. PC1 and PC2 explained 37.8% and 18.9% of total variance, respectively.

**Fig. 2 Differential expression in ketogenic diet versus control liver.** Volcano plot showing log2 fold-change and adjusted P value for the ketogenic diet versus control contrast. Dashed lines indicate absolute log2 fold-change of 0.58 and FDR of 0.05.

**Fig. 3 GO biological-process enrichment of genes higher and lower in ketogenic diet.** Dotplots show selected enriched Gene Ontology biological-process terms from over-representation analysis among genes higher or lower in ketogenic diet liver.

**Fig. 4 Hallmark gene set enrichment analysis.** Barplot of selected MSigDB Hallmark terms ranked by normalized enrichment score. Positive scores indicate enrichment among genes higher in ketogenic diet; negative scores indicate enrichment among genes lower in ketogenic diet. Fill color indicates FDR.

## Additional Files

**Additional file 1. Differential expression results.** CSV file containing DESeq2 results for ketogenic diet versus control diet, including Ensembl identifiers, gene symbols, log2 fold-change, P values, adjusted P values, Entrez identifiers, and gene names.

**Additional file 2. Pathway enrichment results.** CSV files containing GO over-representation analysis, ranked GO GSEA, and Hallmark GSEA results.

**Additional file 3. Reproducible R analysis script.** R script used to download/process metadata, run differential expression, generate plots, and perform pathway analyses.

## References

1. Paoli A, Rubini A, Volek JS, Grimaldi KA. Beyond weight loss: a review of the therapeutic uses of very-low-carbohydrate (ketogenic) diets. Eur J Clin Nutr. 2013;67:789-96. doi:10.1038/ejcn.2013.116.

2. Rho JM. Mechanisms of ketogenic diet action. Epilepsia. 2017;58 Suppl 1:43-58. doi:10.1111/epi.13999.

3. Liang Y, Taylor A, Haas V, et al. A ketogenic diet for the management of binge-eating disorder: a pilot study. Eat Weight Disord. 2026. doi:10.1007/s40519-026-01843-7.

4. Kennedy AR, Pissios P, Otu H, Roberson R, Xue B, Asakura K, et al. A high-fat, ketogenic diet induces a unique metabolic state in mice. Am J Physiol Endocrinol Metab. 2007;292:E1724-39. doi:10.1152/ajpendo.00717.2006.

5. Gallop M, Chaix A. Liver RNA-seq from male C57Bl/6J on a ketogenic diet. Gene Expression Omnibus. 2025. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248297. Accessed 30 Apr 2026.

6. Barrett T, Wilhite SE, Ledoux P, Evangelista C, Kim IF, Tomashevsky M, et al. NCBI GEO: archive for functional genomics data sets--update. Nucleic Acids Res. 2013;41:D991-5. doi:10.1093/nar/gks1193.

7. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15:550. doi:10.1186/s13059-014-0550-8.

8. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57:289-300.

9. Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012;16:284-7. doi:10.1089/omi.2011.0118.

10. The Gene Ontology Consortium, Aleksander SA, Balhoff J, Carbon S, Cherry JM, Drabkin HJ, et al. The Gene Ontology knowledgebase in 2023. Genetics. 2023;224:iyad031. doi:10.1093/genetics/iyad031.

11. Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A. 2005;102:15545-50. doi:10.1073/pnas.0506580102.

12. Liberzon A, Birger C, Thorvaldsdottir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database Hallmark gene set collection. Cell Syst. 2015;1:417-25. doi:10.1016/j.cels.2015.12.004.
