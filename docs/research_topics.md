# Bioinformatics + Ketogenic Diet: Candidate Topics

1. **Long-term ketogenic diet and liver metabolic pathway rewiring**
   - Molecular data: liver RNA-seq from mice on control, high-fat, or ketogenic diet.
   - Why it is strong: direct diet intervention, raw counts available, enough replicates for DESeq2, and liver is central to ketogenesis and lipid metabolism.
   - Main question: Which hepatic pathways are shifted by chronic ketogenic diet relative to low-fat control diet?

2. **Ketogenic diet, epilepsy, and hippocampal astrocyte/neuroinflammatory pathways**
   - Molecular data: mouse hippocampal RNA-seq in wild type and Kcna1 knockout epilepsy model with standard or ketogenic diet.
   - Why it is strong: disease-relevant and mechanistically interesting, especially astrocyte/neuroinflammation pathways.
   - Tradeoff: the genotype-by-diet design is a little more complex for a first-pass analysis.

3. **Ketogenic diet therapy and the gut-brain axis in pediatric epilepsy**
   - Molecular data: metagenomics/metabolomics plus mouse brain transcriptomics from fecal microbiota transfer experiments.
   - Why it is strong: translational and multi-omics.
   - Tradeoff: integrating microbiome, metabolomics, and transcriptomics needs a larger workflow.

4. **Modified ketogenic diet and epigenomic change in adult epilepsy**
   - Molecular data: DNA methylation array data before/after dietary treatment.
   - Why it is strong: human clinical samples and paired design.
   - Tradeoff: methylation pathway interpretation is less direct than RNA-seq differential expression.

5. **Ketogenic diet and exercise tolerance in skeletal muscle**
   - Molecular data: mouse gastrocnemius transcriptome.
   - Why it is strong: connects diet to performance physiology.
   - Tradeoff: less central to ketogenesis than liver and less clinically mature than epilepsy.

## Chosen Topic

I chose **long-term ketogenic diet and liver metabolic pathway rewiring** using GEO accession **GSE248297**. It is a recent public RNA-seq dataset with raw counts, clear sample metadata, and a direct control-vs-ketogenic diet contrast. It is also well suited for pathway analysis because liver gene expression should capture fatty-acid oxidation, lipid handling, mitochondrial metabolism, and stress-response programs.
