# =============================================================================
# 02_CSF_WGCNA_PPMI.R
# Weighted gene co-expression network analysis (WGCNA) of CSF proteins
# (PPMI Project 151, SomaScan 5K)
#
# Manuscript: "Integrated Biofluid Proteomics Identified Dynamic Functional
#              Biomarkers of LRRK2-Linked Parkinson's Disease Progression"
#
# Steps
#   1. Load data, define groups (HC, LRRK2, GBA, sPD)
#   2. Quality control: remove poor proteins/samples and outlier samples
#      (hierarchical clustering of samples, tree cut at height 40)
#   3. Choose the soft-thresholding power (scale-free topology)
#   4. Build a signed network and detect modules (blockwiseModules)
#   5. Module-trait correlations (module eigengene vs. group indicator)
#   6. Protein-level output: module, gene significance (GS), module membership (MM)
#   7. Plots: dendrogram, module-trait heatmap, MM vs. GS scatter plots
#   8. Hub-protein networks for Cytoscape
#   9. GO/KEGG enrichment of selected module proteins
#
# Input : data/Working_file_ppmi_project_151 with phenotypes7.12.25.csv
#           (first 31 columns = phenotypes; remaining columns = log2 SomaScan
#            abundances, already normalized by the SomaScan pipeline)
# Output: results/WGCNA/
# =============================================================================

library(WGCNA)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

## ---- Settings ---------------------------------------------------------------

PPMI_FILE <- "data/Working_file_ppmi_project_151 with phenotypes7.12.25.csv"
N_META    <- 31                    # phenotype columns before the proteins
OUTDIR    <- "results/WGCNA"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(OUTDIR, f)

# Group definitions. In the published analysis the LRRK2 and GBA groups
# include both PD patients and non-manifesting (prodromal) carriers.
# Set INCLUDE_PRODROMAL <- FALSE to restrict them to PD patients.
INCLUDE_PRODROMAL <- TRUE

TREE_CUT_HEIGHT <- 40              # sample-clustering height for outlier removal

# Network parameters
SOFT_POWER     <- 10               # chosen from the scale-free topology plot
MIN_MODULE     <- 30
DEEP_SPLIT     <- 3
MERGE_HEIGHT   <- 0.2
NETWORK_TYPE   <- "signed"

# MM vs. GS scatter plots to draw: module, trait, proteins to label
SCATTERS <- list(
  list(module = "greenyellow", trait = "sPD",   label = NULL),
  list(module = "blue",        trait = "sPD",   label = NULL),
  list(module = "grey",        trait = "LRRK2",
       label = c("HLA_DQA2", "GPNMB_8240_207_3", "GPNMB_5080_131_3", "GAA", "GAPDH",
                 "PPT1", "CTSA", "PARK7_9845_33_3", "CTSH_8644_46_3", "CTSD"))
)

# Hub networks for Cytoscape: module = edge-weight threshold
HUB_MODULES <- c(turquoise = 0.21)
N_HUB       <- 10

# Enrichment: proteins in a module with GS for a trait above a cutoff
ENRICH_SETS <- list(
  list(module = "turquoise",   trait = "LRRK2", gs_min = 0.1)
  # list(module = "greenyellow", trait = "sPD",   gs_min = 0.1)
)

## ---- Helper: protein name -> gene symbol ------------------------------------
# Protein names are "GENE" or "GENE_SeqId" (e.g. GPNMB_8240_207_3);
# hyphenated symbols appear with "_" (e.g. HLA_DQA2 = HLA-DQA2).
valid_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
to_symbol <- function(x) {
  core <- sub("_[0-9]+_[0-9]+_[0-9]+$", "", x)
  hyph <- gsub("_", "-", core)
  ifelse(hyph %in% valid_symbols, hyph, sub("_.*$", "", core))
}

## ---- 1. Load data and define groups -----------------------------------------

protein_matrix <- read.csv(PPMI_FILE, header = TRUE, row.names = 1)
print(table(protein_matrix$GROUP))

lrrk2 <- if (INCLUDE_PRODROMAL) c("PD_LRRK2", "Pro_LRRK2") else "PD_LRRK2"
gba   <- if (INCLUDE_PRODROMAL) c("PD_GBA",   "Pro_GBA")   else "PD_GBA"
protein_matrix$GROUP <- ifelse(protein_matrix$GROUP %in% lrrk2, "LRRK2",
                        ifelse(protein_matrix$GROUP %in% gba,   "GBA",
                        ifelse(protein_matrix$GROUP == "HC",    "HC",
                        ifelse(protein_matrix$GROUP == "sPD",   "sPD", NA))))
protein_matrix <- protein_matrix[!is.na(protein_matrix$GROUP), ]
print(table(protein_matrix$GROUP))

# Expression matrix: samples in rows, proteins in columns
ProData <- protein_matrix[, -(1:N_META)]
colnames(ProData) <- sapply(strsplit(colnames(ProData), "_"),
                            function(x) paste(x[3:length(x)], collapse = "_"))

## ---- 2. Quality control -------------------------------------------------------

# 2a. Proteins/samples with too many missing values or zero variance
gsg <- goodSamplesGenes(ProData, verbose = 3)
ProData <- ProData[gsg$goodSamples, gsg$goodGenes]

# 2b. Outlier samples: cluster samples, cut the tree, keep the largest cluster
sampleTree <- hclust(dist(ProData), method = "average")
clusters   <- cutree(sampleTree, h = TREE_CUT_HEIGHT)
keep       <- clusters == as.integer(names(which.max(table(clusters))))
cat("Outlier samples removed:", sum(!keep), "\n")

pdf(out("QC sample clustering.pdf"), 12, 5)
plot(sampleTree, main = "Sample clustering", sub = "", xlab = "", cex = 0.3)
abline(h = TREE_CUT_HEIGHT, col = "red")
dev.off()

ProData        <- ProData[keep, ]
protein_matrix <- protein_matrix[rownames(ProData), ]
print(table(protein_matrix$GROUP))

# PCA after outlier removal (visual check)
pca <- prcomp(ProData)
pvar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)
p <- ggplot(data.frame(pca$x[, 1:2], group = protein_matrix$GROUP),
            aes(PC1, PC2, colour = group)) +
  geom_point(size = 1) +
  labs(x = paste0("PC1: ", pvar[1], "%"), y = paste0("PC2: ", pvar[2], "%")) +
  theme_classic()
ggsave(out("QC PCA of samples.pdf"), p, width = 5, height = 4)

# Traits: one 0/1 indicator per group (group vs. all other participants)
traitMat <- model.matrix(~ GROUP - 1, data = protein_matrix)
colnames(traitMat) <- sub("^GROUP", "", colnames(traitMat))
stopifnot(identical(rownames(traitMat), rownames(ProData)))

nSamples <- nrow(ProData)
cat("Samples:", nSamples, " Proteins:", ncol(ProData), "\n")

## ---- 3. Soft-thresholding power ----------------------------------------------

powers <- c(1:10, seq(12, 30, by = 2))
sft <- pickSoftThreshold(ProData, powerVector = powers, networkType = NETWORK_TYPE, verbose = 5)
write.csv(sft$fitIndices, out("Soft threshold fit indices.csv"), row.names = FALSE)

pdf(out("Soft threshold.pdf"), 10, 5)
par(mfrow = c(1, 2))
fit_r2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
plot(sft$fitIndices[, 1], fit_r2, type = "n", xlab = "Soft threshold (power)",
     ylab = "Scale-free topology fit, signed R^2", main = "Scale independence")
text(sft$fitIndices[, 1], fit_r2, labels = powers, cex = 0.9, col = "red")
abline(h = 0.85, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft threshold (power)",
     ylab = "Mean connectivity", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = "red")
dev.off()

## ---- 4. Network construction and module detection ----------------------------

temp_cor <- cor
cor <- WGCNA::cor                   # avoid clash with other packages' cor()
net <- blockwiseModules(ProData, power = SOFT_POWER, networkType = NETWORK_TYPE,
                        minModuleSize = MIN_MODULE, deepSplit = DEEP_SPLIT,
                        reassignThreshold = 1e-3, mergeCutHeight = MERGE_HEIGHT,
                        numericLabels = FALSE, pamRespectsDendro = FALSE,
                        maxBlockSize = 20000, randomSeed = 1234, verbose = 3)
cor <- temp_cor

moduleColors <- net$colors
print(table(moduleColors))
write.csv(as.data.frame(table(module = moduleColors)), out("Module sizes.csv"), row.names = FALSE)

pdf(out("Protein dendrogram and modules.pdf"), 12, 9)
plotDendroAndColors(net$dendrograms[[1]], cbind(net$unmergedColors, moduleColors),
                    c("Unmerged", "Merged"), dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

## ---- 5. Module-trait relationships -------------------------------------------

MEs <- orderMEs(moduleEigengenes(ProData, moduleColors)$eigengenes)
write.csv(MEs, out("Module eigengenes.csv"))

moduleTraitCor    <- cor(MEs, traitMat, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

pdf(out("Module-trait relationships.pdf"), 5, 8)
par(mar = c(6, 8.5, 3, 3))
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(traitMat),
               yLabels = colnames(MEs), ySymbols = colnames(MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix,
               setStdMargins = FALSE, cex.text = 0.8, zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()

## ---- 6. Protein-level results: GS and MM -------------------------------------

modNames <- substring(names(MEs), 3)

# Module membership (MM): correlation of each protein with each module eigengene
MM  <- as.data.frame(cor(ProData, MEs, use = "p"));  names(MM) <- paste0("MM.", modNames)
pMM <- as.data.frame(corPvalueStudent(as.matrix(MM), nSamples)); names(pMM) <- paste0("p.MM.", modNames)

# Gene significance (GS): correlation of each protein with each group indicator
GS  <- as.data.frame(cor(ProData, traitMat, use = "p")); names(GS) <- paste0("GS.", colnames(traitMat))
pGS <- as.data.frame(corPvalueStudent(as.matrix(GS), nSamples)); names(pGS) <- paste0("p.GS.", colnames(traitMat))

mm_cols <- as.vector(rbind(names(MM), names(pMM)))           # MM.x, p.MM.x, MM.y, ...
geneInfo <- data.frame(Protein = colnames(ProData), Gene = to_symbol(colnames(ProData)),
                       moduleColor = moduleColors, GS, pGS, cbind(MM, pMM)[, mm_cols],
                       check.names = FALSE)
geneInfo <- geneInfo[order(geneInfo$moduleColor, -abs(geneInfo$GS.HC)), ]
write.csv(geneInfo, out("WGCNA results all proteins.csv"), row.names = FALSE)

## ---- 7. MM vs. GS scatter plots ----------------------------------------------

for (s in SCATTERS) {
  inMod <- moduleColors == s$module
  x <- MM[inMod, paste0("MM.", s$module)]
  y <- GS[inMod, paste0("GS.", s$trait)]
  pdf(out(sprintf("MM vs GS %s module, %s.pdf", s$module, s$trait)), 4.5, 4.5)
  verboseScatterplot(x, y, col = s$module, xlim = c(-1, 1), ylim = c(-1, 1),
                     xlab = paste("Module membership in", s$module, "module"),
                     ylab = paste("Protein significance for", s$trait),
                     main = paste("Module membership vs.", s$trait, "\n"),
                     cex.main = 1.1, cex.lab = 1.1, cex.axis = 1.1)
  lab <- intersect(s$label, rownames(MM)[inMod])
  if (length(lab) > 0) {
    text(MM[lab, paste0("MM.", s$module)], GS[lab, paste0("GS.", s$trait)],
         labels = lab, pos = 2, cex = 0.7, col = "red")
  }
  dev.off()
}

## ---- 8. Hub-protein networks for Cytoscape ----------------------------------

# Topological overlap of the whole network (same settings as module detection)
TOM <- TOMsimilarityFromExpr(ProData, power = SOFT_POWER, networkType = NETWORK_TYPE)
dimnames(TOM) <- list(colnames(ProData), colnames(ProData))

for (m in names(HUB_MODULES)) {
  inMod  <- moduleColors == m
  modTOM <- TOM[inMod, inMod]
  # hubs = the N_HUB proteins with the highest intramodular connectivity
  # (softConnectivity default power = 6, as in the original analysis)
  top <- rank(-softConnectivity(ProData[, inMod])) <= N_HUB
  exportNetworkToCytoscape(modTOM[top, top], weighted = TRUE, threshold = HUB_MODULES[[m]],
    edgeFile = out(paste0("Cytoscape edges ", m, " top", N_HUB, ".txt")),
    nodeFile = out(paste0("Cytoscape nodes ", m, " top", N_HUB, ".txt")))
}

## ---- 9. GO and KEGG enrichment of selected module proteins ------------------

strip_theme <- theme(
  panel.border = element_blank(), panel.background = element_blank(),
  panel.grid = element_blank(), axis.text.y = element_text(size = 14),
  legend.title = element_text(size = 14), legend.text = element_text(size = 12),
  axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()
)

for (e in ENRICH_SETS) {
  sel   <- geneInfo$moduleColor == e$module & geneInfo[[paste0("GS.", e$trait)]] > e$gs_min
  genes <- unique(geneInfo$Gene[sel])
  stem  <- sprintf("%s module, GS %s > %s", e$module, e$trait, e$gs_min)
  gene_entrez <- suppressWarnings(bitr(genes, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db))
  cat(stem, ":", length(genes), "genes; failed to map:",
      paste(setdiff(genes, gene_entrez$SYMBOL), collapse = ", "), "\n")

  for (ont in c("CC", "BP", "MF")) {
    ego <- enrichGO(gene_entrez$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                    ont = ont, pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) next
    write.csv(as.data.frame(ego), out(sprintf("GO_%s %s.csv", ont, stem)), row.names = FALSE)
    p <- dotplot(ego, showCategory = 10, x = "GeneRatio", font.size = 12) +
      aes(x = 1) + scale_x_continuous(NULL, breaks = NULL) + strip_theme
    ggsave(out(sprintf("GO_%s vertical dotplot %s.pdf", ont, stem)), p, width = 5, height = 6)
  }

  kegg <- enrichKEGG(gene_entrez$ENTREZID, organism = "hsa", keyType = "ncbi-geneid",
                     pvalueCutoff = 0.05)
  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    write.csv(as.data.frame(kegg), out(sprintf("KEGG %s.csv", stem)), row.names = FALSE)
    ggsave(out(sprintf("KEGG dotplot %s.pdf", stem)),
           dotplot(kegg, x = "GeneRatio", showCategory = 15), width = 4.5, height = 4)
  }
}

writeLines(c(paste("Run date:", Sys.Date()), capture.output(sessionInfo())),
           out("sessionInfo.txt"))
