# deanalysis.R
# Input:  results/psme1_vs_scram_normalized.csv
# Output: results/DEA_results.csv, results/psme1_kd_volcano_plot.pdf, results/psme1_kd_heatmap.pdf

library(limma)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)
library(ggrepel)

# 1. Load normalized matrix
cat("[1/6] Loading normalized data...\n")
mat <- read.csv("results/psme1_vs_scram_normalized.csv", row.names = 1)
mat <- as.matrix(mat)
cat(sprintf("      %d proteins x %d samples\n", nrow(mat), ncol(mat)))

# 2. Define groups and design matrix
cat("[2/6] Building design matrix...\n")
group <- factor(c("PSME1_KD", "PSME1_KD", "PSME1_KD", "Scram", "Scram", "Scram"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Contrast: KD vs Scram (positive logFC = higher in KD)
contrast_matrix <- makeContrasts(PSME1_KD - Scram, levels = design)

# 3. Fit limma model
# eBayes borrows variance information across all proteins —
# critical for statistical stability with only 3 replicates per group
cat("[3/6] Fitting limma model with eBayes...\n")
fit  <- lmFit(mat, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 4. Extract and save full results
cat("[4/6] Extracting results...\n")
results <- topTable(fit2, coef = 1, number = Inf, sort.by = "P") %>%
  rownames_to_column("ProteinID")

results <- results %>%
  mutate(
    significant = adj.P.Val < 0.05 & abs(logFC) > 1,
    direction   = case_when(
      adj.P.Val < 0.05 & logFC >  1 ~ "Up in KD",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down in KD",
      TRUE                           ~ "NS"
    )
  )

write.csv(results, "results/DEA_results.csv", row.names = FALSE)
n_sig <- sum(results$significant)
cat(sprintf("      Significant hits (FDR < 0.05, |logFC| > 1): %d proteins\n", n_sig))
cat(sprintf("      Up in KD: %d  |  Down in KD: %d\n",
    sum(results$direction == "Up in KD"),
    sum(results$direction == "Down in KD")))

# 5. Volcano plot
cat("[5/6] Generating volcano plot...\n")

# All significant proteins excluding PSME1/PSME2 — kept separate to avoid duplicates
sig_labels <- results %>%
  filter(significant & !ProteinID %in% c("PSME1", "PSME2"))

# PSME1/PSME2 as their own exclusive set for bold formatting
psme_labels <- results %>%
  filter(ProteinID %in% c("PSME1", "PSME2"))

volcano <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Up in KD"   = "#d62728",
    "Down in KD" = "#1f77b4",
    "NS"         = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.4) +
  # All significant hits labeled with small black pointer lines
  ggrepel::geom_text_repel(
    data              = sig_labels,
    aes(label         = ProteinID),
    size              = 2.5,
    color             = "black",
    segment.color     = "black",
    segment.size      = 0.3,
    segment.alpha     = 0.6,
    box.padding       = 0.3,
    min.segment.length = 0,
    max.overlaps      = Inf
  ) +
  # PSME1/PSME2 bold and slightly larger to stand out as KD targets
  ggrepel::geom_text_repel(
    data          = psme_labels,
    aes(label     = ProteinID),
    size          = 3.2,
    color         = "black",
    segment.color = "black",
    segment.size  = 0.4,
    box.padding   = 0.5,
    max.overlaps  = Inf
  ) +
  labs(
    title    = "PSME1 KD vs Scramble",
    subtitle = sprintf("%d significant proteins (FDR < 0.05, |logFC| > 1)", n_sig),
    x        = "log2 Fold Change",
    y        = "-log10 adjusted p-value",
    color    = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
  legend.position = "top",
  plot.title    = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5)
  )

ggsave("results/psme1_kd_volcano_plot.png", volcano, width = 7, height = 6)
cat("      Saved: results/psme1_kd_volcano_plot.png\n")

# 6. Heatmap of top 50 significant proteins
cat("[6/6] Generating heatmap...\n")

sig_proteins <- results %>%
  filter(significant) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  pull(ProteinID)

if (length(sig_proteins) > 1) {
  heatmap_mat <- mat[sig_proteins, ]

  # Z-score per protein so colour reflects relative change, not absolute intensity
  heatmap_mat_z <- t(scale(t(heatmap_mat)))

  annotation_col <- data.frame(
    Group = c("KD", "KD", "KD", "Scram", "Scram", "Scram"),
    row.names = colnames(heatmap_mat_z)
  )
  ann_colors <- list(Group = c(KD = "#d62728", Scram = "#1f77b4"))

  png("results/psme1_kd_heatmap.png", width = 8, height = 12, units = "in", res = 150)
  pheatmap(
    heatmap_mat_z,
    annotation_col    = annotation_col,
    annotation_colors = ann_colors,
    color             = colorRampPalette(c("#1f77b4", "white", "#d62728"))(100),
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    show_rownames     = TRUE,
    fontsize_row      = 7,
    main              = sprintf("Top %d DEPs: PSME1 KD vs Scram (z-score)", length(sig_proteins))
  )
  dev.off()
  cat("      Saved: results/psme1_kd_heatmap.png\n")
} else {
  cat("      Not enough significant proteins for heatmap.\n")
}

cat("\nDEA complete. Run 02_annotate.py next to add gene descriptions to DEA_results.csv\n")