# =============================================================================
# 01_CSF_DAP_PPMI.R
# Differential abundance and GO/KEGG enrichment of CSF proteins
# (PPMI Project 151, SomaScan 5K) for all comparisons in the manuscript
#
# Manuscript: "Integrated Biofluid Proteomics Identified Dynamic Functional
#              Biomarkers of LRRK2-Linked Parkinson's Disease Progression"
#
# For each comparison in COMPARISONS:
#   1. PCA of all samples (quality control)
#   2. limma: abundance ~ group + age + sex; empirical Bayes moderated t test;
#      Benjamini-Hochberg FDR
#   3. Differentially abundant proteins (DAPs):
#      FDR < 0.05 and |log2FC| > 1.5 x SD of all log2FC in that comparison
#   4. Volcano plot
#   5. GO (BP, CC, MF) and KEGG over-representation for down- and
#      up-regulated DAPs
#
# Input : data/Working_file_ppmi_project_151 with phenotypes7.12.25.csv
#           (first 31 columns = phenotypes; remaining columns = log2 SomaScan
#            abundances, already normalized by the SomaScan pipeline)
#         data/excluded_ids.csv
#           (one column "PATNO"; participants excluded before analysis.
#            Not distributed; see README and Methods for the exclusion rules)
# Output: results/CSF/<comparison>/
# =============================================================================

library(ggplot2)
library(limma)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

## ---- Settings ---------------------------------------------------------------

PPMI_FILE     <- "data/Working_file_ppmi_project_151 with phenotypes7.12.25.csv"
EXCLUDED_FILE <- "data/excluded_ids.csv"
OUT_ROOT      <- "results/CSF"
N_META        <- 31          # number of phenotype columns before the proteins

# Column and values for alpha-synuclein seed amplification assay (SAA) status.
# Edit to match the phenotype columns of the input file.
SAA_COL <- "SAA"
SAA_POS <- "Positive"
SAA_NEG <- "Negative"

# One row per comparison. Positive log2FC = higher in the case group.
#   col    : phenotype column that defines the two groups
#   within : if set, restrict to this GROUP first (used for SAA subgroups)
cmp <- function(table, case, control, col = "GROUP", within = NA,
                case_lab = case, control_lab = control) {
  data.frame(table, col, case, control, within, case_lab, control_lab)
}
COMPARISONS <- rbind(
  cmp("Table S1",  "PD_LRRK2",  "HC"),
  cmp("Table S2",  "sPD",       "HC"),
  cmp("Table S3",  "PD_LRRK2",  "sPD"),
  cmp("Table S4",  "PD_GBA",    "HC"),
  cmp("Table S5",  "PD_GBA",    "sPD"),
  cmp("Table S6",  "PD_LRRK2",  "Pro_LRRK2", control_lab = "LRRK2_NMC"),
  cmp("Table S7",  "Pro_LRRK2", "HC",        case_lab    = "LRRK2_NMC"),
  cmp("Table S8",  "PD_GBA",    "Pro_GBA",   control_lab = "GBA_NMC"),
  cmp("Table S9",  "Pro_GBA",   "HC",        case_lab    = "GBA_NMC"),
  cmp("Table S10", SAA_POS, SAA_NEG, col = SAA_COL, within = "PD_LRRK2",
      case_lab = "LRRK2_PD_SAApos", control_lab = "LRRK2_PD_SAAneg"),
  cmp("Fig. S2h",  SAA_POS, SAA_NEG, col = SAA_COL, within = "PD_GBA",
      case_lab = "GBA_PD_SAApos",   control_lab = "GBA_PD_SAAneg")
)

# Volcano axis limits per comparison; comparisons not listed use automatic limits.
VOLCANO_LIMITS <- list(
  "PD_LRRK2_vs_HC" = list(xlim = c(-0.8, 0.9), ylim = c(0, 28))
)

# Gene symbols containing a hyphen (e.g. HLA-DQA2) appear with "_" in the
# protein names. TRUE restores the official symbol so these proteins are
# labelled correctly and map in GO/KEGG. FALSE = original rule (text before
# the first "_", which turns HLA-DQA2 into "HLA").
FIX_HYPHEN_SYMBOLS <- TRUE

## ---- Helpers -----------------------------------------------------------------

# Protein name -> gene symbol. Names are "GENE" or "GENE_SeqId" (e.g. GPNMB_8240_207_3).
valid_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
to_symbol <- function(x) {
  if (!FIX_HYPHEN_SYMBOLS) return(sub("_.*$", "", x))
  core  <- sub("_[0-9]+_[0-9]+_[0-9]+$", "", x)        # drop SeqId
  hyph  <- gsub("_", "-", core)                        # HLA_DQA2 -> HLA-DQA2
  first <- sub("_.*$", "", core)                       # GSK3A_GSK3B -> GSK3A
  ifelse(hyph %in% valid_symbols, hyph, first)
}

# Vertical-strip theme for GO dot plots
strip_theme <- theme(
  panel.border = element_blank(), panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
  legend.title = element_text(size = 14), legend.text = element_text(size = 12),
  axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()
)

## ---- Load data and remove excluded participants ----------------------------
data <- read.csv(PPMI_FILE, header = TRUE, row.names = 1)   # row names = PATNO
to_remove <- as.character(read.csv(EXCLUDED_FILE, colClasses = "character")$PATNO)
cat("Excluded participants found in file:", sum(rownames(data) %in% to_remove), "\n")
data <- data[!(rownames(data) %in% to_remove), ]
print(table(data$GROUP))

## ---- Run each comparison -------------------------------------------------------
for (i in seq_len(nrow(COMPARISONS))) {

  cm   <- COMPARISONS[i, ]
  name <- paste0(cm$case_lab, "_vs_", cm$control_lab)
  cat("\n=====", cm$table, ":", name, "=====\n")

  # -- Select samples ------------------------------------------------------------
  if (!cm$col %in% colnames(data)) {
    warning(name, ": column '", cm$col, "' not found; comparison skipped."); next
  }
  d <- if (is.na(cm$within)) data else data[data$GROUP == cm$within, ]
  d <- d[!is.na(d[[cm$col]]) & d[[cm$col]] %in% c(cm$case, cm$control), ]
  n_case <- sum(d[[cm$col]] == cm$case); n_ctrl <- sum(d[[cm$col]] == cm$control)
  if (n_case < 2 || n_ctrl < 2) {
    warning(name, ": fewer than 2 samples in a group; comparison skipped."); next
  }

  outdir <- file.path(OUT_ROOT, name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  lab <- sprintf("%s (%d) vs %s (%d)", cm$case_lab, n_case, cm$control_lab, n_ctrl)
  cat(lab, "\n")

  sample_data <- d[, 1:N_META]
  sample_data$CONTRAST <- factor(d[[cm$col]], levels = c(cm$control, cm$case))
  filtered_new_data <- d[, -(1:N_META)]
  colnames(filtered_new_data) <- sapply(strsplit(colnames(filtered_new_data), "_"),
                                        function(x) paste(x[3:length(x)], collapse = "_"))
  filtered_new_data <- t(filtered_new_data)                 # proteins x samples

  # -- PCA (quality control) ---------------------------------------------------
  pca_result <- prcomp(t(filtered_new_data), scale. = TRUE)
  pca_df <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2],
                       condition = sample_data$CONTRAST, sample_names = rownames(sample_data))
  p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    geom_text(aes(label = sample_names), vjust = -1, size = 3) +
    labs(title = paste("PCA of CSF proteome:", lab))
  ggsave(file.path(outdir, "PCA.pdf"), p_pca, width = 8, height = 7)

  # -- limma: abundance ~ group + age + sex ------------------------------------
  design  <- model.matrix(~ CONTRAST + AGE + SexM1.F0, data = sample_data)
  fit     <- eBayes(lmFit(filtered_new_data, design))
  results <- topTable(fit, coef = 2, number = Inf, adjust = "fdr")   # coef 2 = case vs control
  results <- results[order(results$adj.P.Val), ]
  results$genesid <- to_symbol(rownames(results))
  cat("Proteins with FDR < 0.05:", sum(results$adj.P.Val < 0.05, na.rm = TRUE), "\n")
  write.csv(results, file.path(outdir, paste0("Result of all proteins of ", lab, ".csv")))

  # -- DAPs: FDR < 0.05 and |log2FC| > 1.5 x SD --------------------------------
  SD   <- sd(results$logFC)
  DEPs <- subset(results, adj.P.Val < 0.05 & (logFC > 1.5 * SD | logFC < -1.5 * SD))
  cat(sprintf("log2FC cutoff (1.5 x SD) = %.4f; DAPs: %d proteins (%d up, %d down), %d genes\n",
              1.5 * SD, nrow(DEPs), sum(DEPs$logFC > 0), sum(DEPs$logFC < 0),
              length(unique(DEPs$genesid))))
  write.csv(DEPs, file.path(outdir, paste0("Result of sigdiff proteins of ", lab, ".csv")))

  # -- Volcano plot --------------------------------------------------------------
  lims <- VOLCANO_LIMITS[[name]]
  if (is.null(lims)) {
    lims <- list(xlim = range(results$logFC) * 1.05,
                 ylim = c(0, max(-log10(results$adj.P.Val)) * 1.05))
  }
  pdf(file.path(outdir, "Volcano plot.pdf"), 7.5, 4.5)
  print(EnhancedVolcano(results, lab = rownames(results), x = "logFC", y = "adj.P.Val",
                        xlim = lims$xlim, ylim = lims$ylim,
                        pCutoff = 0.05, FCcutoff = 1.5 * SD,
                        gridlines.major = FALSE, gridlines.minor = FALSE,
                        legendLabels = c("NS", "Log2FC", "FDR", "Log2FC & FDR"),
                        xlab = bquote(~Log[2] ~ "fold change"),
                        ylab = bquote(~-Log[10] ~ FDR),
                        pointSize = 0.2,
                        col = c("#B3B3B3", "#FF9900", "#0072FF", "#E60000"),  # NS, FC, FDR, both
                        colAlpha = 1,
                        legendPosition = "right", legendLabSize = 14, legendIconSize = 2,
                        legendDropLevels = FALSE,
                        drawConnectors = TRUE, min.segment.length = 0.0001,
                        selectLab = c(" "),
                        title = NULL, subtitle = lab, max.overlaps = 10))
  dev.off()

  # -- GO and KEGG enrichment: down- and up-regulated DAPs --------------------
  gene_lists <- list(DOWN = unique(DEPs$genesid[DEPs$logFC < 0]),
                     UP   = unique(DEPs$genesid[DEPs$logFC > 0]))

  for (direction in names(gene_lists)) {
    if (length(gene_lists[[direction]]) == 0) next
    gene_entrez <- suppressWarnings(bitr(gene_lists[[direction]], fromType = "SYMBOL",
                                         toType = "ENTREZID", OrgDb = org.Hs.eg.db))
    cat(direction, "- failed to map:",
        paste(setdiff(gene_lists[[direction]], gene_entrez$SYMBOL), collapse = ", "), "\n")
    if (nrow(gene_entrez) == 0) next
    stem <- paste0(direction, " ", lab)

    for (ont in c("CC", "BP", "MF")) {
      ego_obj <- enrichGO(gene = gene_entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID", ont = ont,
                          pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
      if (is.null(ego_obj) || nrow(as.data.frame(ego_obj)) == 0) next
      write.csv(as.data.frame(ego_obj),
                file.path(outdir, sprintf("GO_%s %s.csv", ont, stem)), row.names = FALSE)
      p <- dotplot(ego_obj, showCategory = 10, x = "GeneRatio", font.size = 12) +
        aes(x = 1) + scale_x_continuous(NULL, breaks = NULL) + strip_theme
      ggsave(file.path(outdir, sprintf("GO_%s vertical dotplot %s.pdf", ont, stem)),
             p, width = 5, height = 6)
    }

    # KEGG (queries the KEGG website; record the analysis date)
    kegg_enrichment <- enrichKEGG(gene = gene_entrez$ENTREZID, organism = "hsa",
                                  pvalueCutoff = 0.05, keyType = "ncbi-geneid")
    if (is.null(kegg_enrichment) || nrow(as.data.frame(kegg_enrichment)) == 0) next
    kegg_res <- setReadable(kegg_enrichment, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    write.csv(as.data.frame(kegg_res),
              file.path(outdir, sprintf("KEGG %s.csv", stem)), row.names = FALSE)
    ggsave(file.path(outdir, sprintf("KEGG dotplot %s.pdf", stem)),
           dotplot(kegg_res, x = "GeneRatio", showCategory = 15), width = 7, height = 7)
  }
}

## ---- Record software versions ----------------------------------------------
writeLines(c(paste("Run date:", Sys.Date()), capture.output(sessionInfo())),
           file.path(OUT_ROOT, "sessionInfo.txt"))
