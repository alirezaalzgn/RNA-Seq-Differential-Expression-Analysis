
# ============================================================
# 03_gsea_hallmark.R
# Hallmark Gene Set Enrichment Analysis
# Comparison: Liver vs Lung metastases
# ============================================================

# 1. Load required packages
library(DESeq2)
library(clusterProfiler)
library(msigdbr)

# 2. Define project paths
project_dir <- getwd()

dds_file <- file.path(
  project_dir, "Results", "DESeq2", "dds.rds"
)

output_dir <- file.path(
  project_dir, "Results", "Enrichment"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Load the fitted DESeq2 object
if (!file.exists(dds_file)) {
  stop("Cannot find dds.rds. Check the project path.")
}

dds <- readRDS(dds_file)

# 4. Calculate Liver vs Lung differential expression
res <- results(
  dds,
  contrast = c("Tissue", "Liver", "Lung")
)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Keep genes with valid Wald statistics
res_df <- res_df[
  !is.na(res_df$stat) & !is.na(res_df$gene_id),
]

# 5. Prepare the ranked gene list
# Remove Ensembl version suffixes if present
res_df$gene_id <- sub("\\..*$", "", res_df$gene_id)

# If IDs become duplicated, retain the entry with the
# largest absolute Wald statistic
res_df <- res_df[
  order(abs(res_df$stat), decreasing = TRUE),
]

res_df <- res_df[!duplicated(res_df$gene_id), ]

gene_list <- res_df$stat
names(gene_list) <- res_df$gene_id

# GSEA requires decreasing order
gene_list <- sort(gene_list, decreasing = TRUE)

cat("Genes ranked for GSEA:", length(gene_list), "\n")
cat("Positive statistics:", sum(gene_list > 0), "\n")
cat("Negative statistics:", sum(gene_list < 0), "\n")

# 6. Retrieve Hallmark gene sets
msig <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

term2gene <- unique(
  data.frame(
    pathway = msig$gs_name,
    gene = msig$ensembl_gene
  )
)

term2gene <- term2gene[
  !is.na(term2gene$gene) & !is.na(term2gene$pathway),
]

cat("Hallmark pathways:", length(unique(term2gene$pathway)), "\n")

# 7. Run preranked GSEA
set.seed(42)

gsea_result <- clusterProfiler::GSEA(
  geneList = gene_list,
  TERM2GENE = term2gene,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 15,
  maxGSSize = 500,
  eps = 0,
  verbose = FALSE,
  seed = TRUE
)

gsea_df <- as.data.frame(gsea_result)

if (nrow(gsea_df) == 0) {
  stop(
    "GSEA returned no pathways. Check gene ID matching and package versions."
  )
}

# 8. Save all pathway results
write.csv(
  gsea_df,
  file.path(output_dir, "GSEA_Hallmark_Liver_vs_Lung_recomputed.csv"),
  row.names = FALSE
)

# 9. Separate pathways by enrichment direction
# Positive NES = enriched toward Liver
# Negative NES = enriched toward Lung

liver_pathways <- gsea_df[
  gsea_df$NES > 0,
]

lung_pathways <- gsea_df[
  gsea_df$NES < 0,
]

write.csv(
  liver_pathways,
  file.path(output_dir, "GSEA_Hallmark_Liver_enriched_recomputed.csv"),
  row.names = FALSE
)

write.csv(
  lung_pathways,
  file.path(output_dir, "GSEA_Hallmark_Lung_enriched_recomputed.csv"),
  row.names = FALSE
)

# 10. Save significant pathways
significant_pathways <- gsea_df[
  !is.na(gsea_df$p.adjust) & gsea_df$p.adjust < 0.05,
]

write.csv(
  significant_pathways,
  file.path(output_dir, "GSEA_Hallmark_significant_recomputed.csv"),
  row.names = FALSE
)

# 11. Record package versions and session information
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "GSEA_sessionInfo.txt")
)

# 12. Summary
cat("\nGSEA completed successfully.\n")
cat("Total pathways tested:", nrow(gsea_df), "\n")
cat("Significant pathways (adjusted p < 0.05):",
    nrow(significant_pathways), "\n")
cat("Liver-enriched pathways:", nrow(liver_pathways), "\n")
cat("Lung-enriched pathways:", nrow(lung_pathways), "\n")
cat("Results saved in:", output_dir, "\n")

# 13. Save the GSEA object for reproducibility
saveRDS(
  gsea_result,
  file.path(output_dir, "GSEA_Hallmark_recomputed.rds")
)

# 14. Create enrichment plots
library(enrichplot)
library(ggplot2)

p_oxphos <- gseaplot2(
  gsea_result,
  geneSetID = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  title = "Hallmark Oxidative Phosphorylation"
)

ggsave(
  file.path(output_dir, "GSEA_Oxidative_Phosphorylation_recomputed.png"),
  plot = p_oxphos,
  width = 9,
  height = 6,
  dpi = 300
)

p_emt <- gseaplot2(
  gsea_result,
  geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  title = "Hallmark Epithelial–Mesenchymal Transition"
)

ggsave(
  file.path(output_dir, "GSEA_EMT_recomputed.png"),
  plot = p_emt,
  width = 9,
  height = 6,
  dpi = 300
)

cat("GSEA plots saved successfully.\n")
