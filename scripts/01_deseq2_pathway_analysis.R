initial_project_dir <- normalizePath(getwd(), mustWork = TRUE)
if (basename(initial_project_dir) == "scripts") {
  initial_project_dir <- normalizePath(file.path(initial_project_dir, ".."), mustWork = TRUE)
}
dir.create(file.path(initial_project_dir, ".cache"), recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_USER_CACHE_DIR = file.path(initial_project_dir, ".cache"))

suppressPackageStartupMessages({
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(msigdbr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(readr)
  library(tibble)
  library(stringr)
})

project_dir <- initial_project_dir
if (!dir.exists(file.path(project_dir, "data"))) {
  stop("Run this script from the RStudio project root or from the scripts directory.")
}

raw_counts_file <- file.path(project_dir, "data/raw/GSE248297_20462R_ALF_raw_counts.txt.gz")
processed_dir <- file.path(project_dir, "data/processed")
tables_dir <- file.path(project_dir, "results/tables")
figures_dir <- file.path(project_dir, "results/figures")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading raw counts...")
counts_tbl <- read_tsv(raw_counts_file, show_col_types = FALSE)
counts_mat <- counts_tbl %>%
  column_to_rownames("geneid") %>%
  as.matrix()
storage.mode(counts_mat) <- "integer"

sample_info <- tibble(
  sample_id = colnames(counts_mat),
  geo_accession = c(
    "GSM7911127", "GSM7911128", "GSM7911129", "GSM7911130", "GSM7911131",
    "GSM7911132", "GSM7911133", "GSM7911134", "GSM7911135", "GSM7911136",
    "GSM7911137", "GSM7911138", "GSM7911139", "GSM7911140", "GSM7911141"
  ),
  title = c(
    "20462X1_Liver_Control_ALF_replicate1",
    "20462X2_Liver_Control_ALF_replicate2",
    "20462X3_Liver_HFD_ALF_replicate1",
    "20462X4_Liver_HFD_ALF_replicate2",
    "20462X8_Liver_Keto_ALF_replicate1",
    "20462X9_Liver_Keto_ALF_replicate2",
    "20462X13_Liver_HFD_ALF_replicate3",
    "20462X14_Liver_HFD_ALF_replicate4",
    "20462X15_Liver_HFD_ALF_replicate5",
    "20462X18_Liver_Control_ALF_replicate3",
    "20462X19_Liver_Control_ALF_replicate4",
    "20462X20_Liver_Control_ALF_replicate5",
    "20462X23_Liver_Keto_ALF_replicate3",
    "20462X24_Liver_Keto_ALF_replicate4",
    "20462X25_Liver_Keto_ALF_replicate5"
  )
) %>%
  mutate(
    diet = case_when(
      str_detect(title, "Control") ~ "Control",
      str_detect(title, "HFD") ~ "HFD",
      str_detect(title, "Keto") ~ "Keto",
      TRUE ~ NA_character_
    ),
    replicate = str_extract(title, "replicate\\d+"),
    contaminated = sample_id == "20462X1",
    diet = factor(diet, levels = c("Control", "HFD", "Keto"))
  )

write_csv(sample_info, file.path(processed_dir, "sample_metadata.csv"))

analysis_samples <- sample_info %>%
  filter(!contaminated)

counts_use <- counts_mat[, analysis_samples$sample_id, drop = FALSE]
keep_genes <- rowSums(counts_use >= 10) >= 3
counts_use <- counts_use[keep_genes, ]

message("Running DESeq2...")
dds <- DESeqDataSetFromMatrix(
  countData = counts_use,
  colData = as.data.frame(column_to_rownames(analysis_samples, "sample_id")),
  design = ~ diet
)
dds <- DESeq(dds)

vst_mat <- assay(vst(dds, blind = FALSE))
saveRDS(dds, file.path(processed_dir, "dds_gse248297.rds"))
saveRDS(vst_mat, file.path(processed_dir, "vst_matrix_gse248297.rds"))

annotate_results <- function(res) {
  res_tbl <- as.data.frame(res) %>%
    rownames_to_column("ensembl_gene_id") %>%
    mutate(
      symbol = mapIds(org.Mm.eg.db, keys = ensembl_gene_id, keytype = "ENSEMBL",
                      column = "SYMBOL", multiVals = "first"),
      entrez_id = mapIds(org.Mm.eg.db, keys = ensembl_gene_id, keytype = "ENSEMBL",
                         column = "ENTREZID", multiVals = "first"),
      gene_name = mapIds(org.Mm.eg.db, keys = ensembl_gene_id, keytype = "ENSEMBL",
                         column = "GENENAME", multiVals = "first")
    ) %>%
    arrange(padj, pvalue)
  res_tbl
}

keto_vs_control <- annotate_results(results(dds, contrast = c("diet", "Keto", "Control")))
hfd_vs_control <- annotate_results(results(dds, contrast = c("diet", "HFD", "Control")))

write_csv(keto_vs_control, file.path(tables_dir, "deseq2_keto_vs_control.csv"))
write_csv(hfd_vs_control, file.path(tables_dir, "deseq2_hfd_vs_control.csv"))

normalized_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_gene_id") %>%
  mutate(symbol = mapIds(org.Mm.eg.db, keys = ensembl_gene_id, keytype = "ENSEMBL",
                         column = "SYMBOL", multiVals = "first"))
write_csv(normalized_counts, file.path(processed_dir, "normalized_counts_deseq2.csv"))

message("Creating PCA plot...")
pca <- prcomp(t(vst_mat))
pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("sample_id") %>%
  left_join(analysis_samples, by = "sample_id")
percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = diet, label = sample_id)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) +
  labs(
    title = "GSE248297 liver RNA-seq PCA",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Diet"
  ) +
  theme_minimal(base_size = 12)
ggsave(file.path(figures_dir, "pca_vst_deseq2.png"), p_pca, width = 7, height = 5, dpi = 300, bg = "white")

message("Creating volcano plot...")
volcano_df <- keto_vs_control %>%
  mutate(
    neg_log10_padj = -log10(padj),
    label = if_else(row_number() <= 15 & !is.na(symbol), symbol, NA_character_),
    direction = case_when(
      padj < 0.05 & log2FoldChange >= 0.58 ~ "Higher in keto",
      padj < 0.05 & log2FoldChange <= -0.58 ~ "Lower in keto",
      TRUE ~ "Not significant"
    )
  )

p_volcano <- ggplot(volcano_df, aes(log2FoldChange, neg_log10_padj, color = direction)) +
  geom_point(alpha = 0.7, size = 1.2) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey55") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey55") +
  ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf, na.rm = TRUE) +
  scale_color_manual(values = c("Higher in keto" = "#B23A48", "Lower in keto" = "#287C8E", "Not significant" = "grey70")) +
  labs(
    title = "Ketogenic diet vs control liver differential expression",
    x = "log2 fold-change",
    y = "-log10 adjusted p-value",
    color = NULL
  ) +
  theme_minimal(base_size = 12)
ggsave(file.path(figures_dir, "volcano_keto_vs_control.png"), p_volcano, width = 7, height = 5, dpi = 300, bg = "white")

message("Creating top-gene heatmap...")
top_gene_ids <- keto_vs_control %>%
  filter(!is.na(padj)) %>%
  slice_min(padj, n = 40) %>%
  pull(ensembl_gene_id)
heat_mat <- vst_mat[top_gene_ids, , drop = FALSE]
heat_labels <- keto_vs_control %>%
  filter(ensembl_gene_id %in% top_gene_ids) %>%
  mutate(label = if_else(is.na(symbol), ensembl_gene_id, symbol)) %>%
  select(ensembl_gene_id, label)
rownames(heat_mat) <- heat_labels$label[match(rownames(heat_mat), heat_labels$ensembl_gene_id)]
annotation_col <- analysis_samples %>%
  select(sample_id, diet) %>%
  column_to_rownames("sample_id")
png(file.path(figures_dir, "heatmap_top40_keto_vs_control.png"), width = 1800, height = 1600, res = 220, bg = "white")
pheatmap(
  heat_mat,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = TRUE,
  fontsize_row = 6,
  main = "Top 40 DE genes: keto vs control"
)
dev.off()

message("Running GO over-representation analysis...")
sig_up_entrez <- keto_vs_control %>%
  filter(!is.na(padj), padj < 0.05, log2FoldChange >= 0.58, !is.na(entrez_id)) %>%
  pull(entrez_id) %>%
  unique()
sig_down_entrez <- keto_vs_control %>%
  filter(!is.na(padj), padj < 0.05, log2FoldChange <= -0.58, !is.na(entrez_id)) %>%
  pull(entrez_id) %>%
  unique()
background_entrez <- keto_vs_control %>%
  filter(!is.na(entrez_id)) %>%
  pull(entrez_id) %>%
  unique()

run_enrich_go <- function(genes, label) {
  if (length(genes) < 10) {
    return(NULL)
  }
  ego <- enrichGO(
    gene = genes,
    universe = background_entrez,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  ego_df <- as.data.frame(ego)
  write_csv(ego_df, file.path(tables_dir, paste0("go_bp_enrichment_", label, ".csv")))
  if (nrow(ego_df) > 0) {
    p <- dotplot(ego, showCategory = min(15, nrow(ego_df))) +
      ggtitle(paste("GO biological process:", label, "in keto vs control"))
    ggsave(file.path(figures_dir, paste0("go_bp_dotplot_", label, ".png")), p, width = 8, height = 6, dpi = 300, bg = "white")
  }
  ego
}

ego_up <- run_enrich_go(sig_up_entrez, "higher_in_keto")
ego_down <- run_enrich_go(sig_down_entrez, "lower_in_keto")

message("Running ranked GO GSEA...")
rank_tbl <- keto_vs_control %>%
  filter(!is.na(stat), !is.na(entrez_id)) %>%
  group_by(entrez_id) %>%
  slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(stat))
gene_ranks <- rank_tbl$stat
names(gene_ranks) <- rank_tbl$entrez_id

gsea_go <- gseGO(
  geneList = gene_ranks,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.2,
  pAdjustMethod = "BH",
  verbose = FALSE
)
gsea_go_df <- as.data.frame(gsea_go)
write_csv(gsea_go_df, file.path(tables_dir, "gsea_go_bp_keto_vs_control.csv"))
if (nrow(gsea_go_df) > 0) {
  p <- dotplot(gsea_go, showCategory = min(20, nrow(gsea_go_df)), split = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle("GO BP GSEA: keto vs control")
  ggsave(file.path(figures_dir, "gsea_go_bp_dotplot_keto_vs_control.png"), p, width = 10, height = 6, dpi = 300, bg = "white")
}

message("Running Hallmark GSEA with msigdbr...")
hallmark <- msigdbr(species = "Mus musculus", db_species = "MM", collection = "MH") %>%
  select(gs_name, ncbi_gene) %>%
  filter(!is.na(ncbi_gene)) %>%
  mutate(ncbi_gene = as.character(ncbi_gene)) %>%
  distinct()

gsea_hallmark <- GSEA(
  geneList = gene_ranks,
  TERM2GENE = hallmark %>% transmute(term = gs_name, gene = ncbi_gene),
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE
)
gsea_hallmark_df <- as.data.frame(gsea_hallmark)
write_csv(gsea_hallmark_df, file.path(tables_dir, "gsea_hallmark_keto_vs_control.csv"))
if (nrow(gsea_hallmark_df) > 0) {
  top_hallmark <- gsea_hallmark_df %>%
    arrange(p.adjust) %>%
    slice_head(n = 20)
  p <- ggplot(top_hallmark, aes(NES, reorder(Description, NES), fill = p.adjust)) +
    geom_col() +
    scale_fill_viridis_c(option = "mako", direction = -1) +
    labs(
      title = "Hallmark GSEA: keto vs control",
      x = "Normalized enrichment score",
      y = NULL,
      fill = "FDR"
    ) +
    theme_minimal(base_size = 11)
  ggsave(file.path(figures_dir, "gsea_hallmark_barplot_keto_vs_control.png"), p, width = 10, height = 6, dpi = 300, bg = "white")
}

summary_tbl <- tibble(
  metric = c(
    "samples_total_downloaded",
    "samples_analyzed_after_excluding_contaminated",
    "genes_in_raw_count_matrix",
    "genes_after_low_count_filter",
    "keto_vs_control_FDR_0.05",
    "keto_vs_control_FDR_0.05_absLFC_0.58",
    "keto_higher_FDR_0.05_absLFC_0.58",
    "keto_lower_FDR_0.05_absLFC_0.58",
    "hallmark_terms_FDR_0.25",
    "go_gsea_terms_FDR_0.25"
  ),
  value = c(
    ncol(counts_mat),
    nrow(analysis_samples),
    nrow(counts_mat),
    nrow(counts_use),
    sum(!is.na(keto_vs_control$padj) & keto_vs_control$padj < 0.05),
    sum(!is.na(keto_vs_control$padj) & keto_vs_control$padj < 0.05 & abs(keto_vs_control$log2FoldChange) >= 0.58),
    length(sig_up_entrez),
    length(sig_down_entrez),
    sum(!is.na(gsea_hallmark_df$p.adjust) & gsea_hallmark_df$p.adjust < 0.25),
    sum(!is.na(gsea_go_df$p.adjust) & gsea_go_df$p.adjust < 0.25)
  )
)
write_csv(summary_tbl, file.path(tables_dir, "analysis_summary.csv"))

message("Done. Results written to:")
message("  ", tables_dir)
message("  ", figures_dir)
