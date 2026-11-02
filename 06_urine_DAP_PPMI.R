# =============================================================================
# 09_urine_DAP_PPMI.R
# Differential abundance and GO/KEGG enrichment of urine proteins
# (PPMI Project 190, mass spectrometry, baseline visit) for all comparisons
# in the manuscript
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
# Input : ./PPMI_Project_190_log2_translated_imputated_data.csv
#           (PATNO, 30 phenotype columns, then protein-group abundances that were
#            log2-transformed, filtered (detected in >70% of samples) and
#            imputed from a normal distribution in Perseus; see Methods)
# Output: results/urine/<comparison>/
# =============================================================================

library(ggplot2)
library(limma)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

data <- read.csv("./PPMI_Project_190_log2_translated_imputated_data.csv", header = T, row.names = 1)

new_data <- data[data$EVENT_ID == "BL" & (data$GROUP == "HC" |data$GROUP == "PD_LRRK2"),] # Please change the group comparisons and below
table(new_data$GROUP)
row.names(new_data) <- new_data$PATNO
new_data$PATNO <- NULL 

sample_data <- new_data[, c(1:30)]

filtered_new_data <- new_data[, -c(1:30)]


head(filtered_new_data)


filtered_new_data <- t(filtered_new_data)

library(ggplot2)
library(stats)

pca_result <- prcomp(t(filtered_new_data), scale. = TRUE)
pca_df <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], condition = sample_data$GROUP)
pca_df$sample_names <- rownames(pca_df) 

ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names), vjust = -1, size = 3) +
  labs(title = "PCA of Proteomic Data", x = "PC1", y = "PC2")



# Load limma package
library(limma)

# Create design matrix for the condition of interest   Control vs Treatment
sample_data$GROUP <- factor(sample_data$GROUP, levels = c("HC", "PD_LRRK2"))
# Convert the 0/1 code into a two‐level factor with meaningful labels
sample_data$SEX <- factor(sample_data$SEX, levels = c(0, 1),labels = c("Female", "Male"))
# Verify
table(sample_data$SEX)

#sample_data$GROUP <- factor(sample_data$GROUP)
design <- model.matrix(~ GROUP + AGE + SEX, data=sample_data)

# Fit the linear model to the log2-transformed data
fit <- lmFit(filtered_new_data, design)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Extract differentially expressed proteins for condition (Treated vs. Control)
results <- topTable(fit, coef = "GROUPPD_LRRK2", number = Inf, adjust = "fdr")


head(results)
summary(results)
sum(results$adj.P.Val < 0.05, na.rm=TRUE)
results <- results[order(results$adj.P.Val),]
head(results)

results$genesid <- sub("_.*$", "", rownames(results))
write.csv(results, file = "Result of all proteins of PD_LRRK2 (151) vs HC (140).csv")


SD <- sd(results$logFC)
DEPs <-subset(results, adj.P.Val < 0.05 & (logFC > 1.5*SD | logFC < - 1.5*SD))

dim(DEPs)
head(DEPs)
write.csv(DEPs,file= "Result of sigdiff proteins of PD_LRRK2 (151) vs HC (140).csv")


pdf('Volcanoplot for PB.pdf', 7.5,4.5)
library(EnhancedVolcano)
EnhancedVolcano(results, lab = row.names(results), x = 'logFC', y = 'adj.P.Val', xlim = c(-2.5, 2.5), ylim=c(0,37),
                pCutoff = 0.05, FCcutoff = 1.5*SD, gridlines.major = F, gridlines.minor = F,
                legendLabels=c('NS','Log2FC','FDR','Log2FC & FDR'),
                xlab = bquote(~Log[2] ~ "fold change"),
                ylab = bquote(~-Log[10] ~ FDR),
                pointSize = 0.2,# NS = gray; Log2FC only = orange; FDR only = blue; both = red
                col = c("#B3B3B3", "#FF9900", "#0072FF", "#E60000"),
                colAlpha = 1,
                legendPosition = "right",
                legendLabSize = 14,
                legendIconSize = 2,
                legendDropLevels = F,
                drawConnectors = T, min.segment.length = 0.0001, 
                selectLab = c(" "),
                title = NULL,subtitle = "LRRK2_DP vs. HC",
                max.overlaps = 30)
dev.off() 


Total_genes <- DEPs
Total_genes$genesid <- sub("_.*$", "", rownames(DEPs))

Up_regulated_genes <- Total_genes[Total_genes$logFC > 1.5*SD & Total_genes$adj.P.Val < 0.05,]
Up_genes <- Up_regulated_genes$genesid

Down_regulated_genes <- Total_genes[Total_genes$logFC < -1.5*SD & Total_genes$adj.P.Val < 0.05, ]
Down_genes <- Down_regulated_genes$genesid

dim(DEPs)
dim(Up_regulated_genes)
dim(Down_regulated_genes)

library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)


# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(Down_genes, fromType = "SYMBOL",  #change to Up_genes for enriched pathways
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

# Find the genes that failed to map
failed_genes <- setdiff(Down_genes, gene_entrez$SYMBOL)

# Print the failed genes
print("Failed to map genes:")
print(failed_genes)

#Go enrichment analysis
library(stringr)
ego_cc<-enrichGO(gene       = gene_entrez$ENTREZID,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'ENTREZID',
                 ont        = "CC",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable      = TRUE)
head(ego_cc, 4)

ego_bp<-enrichGO(gene       = gene_entrez$ENTREZID,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'ENTREZID',
                 ont        = "BP",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable      = TRUE)
head(ego_bp, 4)

ego_mf<-enrichGO(gene       = gene_entrez$ENTREZID,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'ENTREZID',
                 ont        = "MF",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable      = TRUE)
head(ego_mf, 4)

# 1. collect your GO objects
ego_list <- list(
  CC = ego_cc,
  BP = ego_bp,
  MF = ego_mf
)


strip_theme <- theme(
  panel.border     = element_blank(),
  panel.background = element_blank(),
  panel.grid       = element_blank(),
  

  axis.text.y      = element_text(size = 14), 
  axis.title.y     = element_text(size = 16), 
  

  legend.title     = element_text(size = 14),
  legend.text      = element_text(size = 12),
  

  axis.title.x     = element_blank(),
  axis.text.x      = element_blank(),
  axis.ticks.x     = element_blank()
)

for (ont in names(ego_list)) {
  ego_obj <- ego_list[[ont]]
  

  write.csv(
    as.data.frame(ego_obj),
    file = sprintf("The_GO_%s_of_enriched_proteins_in_DOWN_PD_LRRK2 (151) vs HC (140).csv", ont),
    row.names = FALSE
  )
 
  p <- dotplot(
    ego_obj,
    showCategory = 10,
    x            = "GeneRatio",
    font.size    = 12       
  ) +
    aes(x = 1) +          
    scale_x_continuous(NULL, breaks = NULL) +
    strip_theme
  

  pdf(
    sprintf("GO_%s_vertical_dotplot_DOWN_PD_LRRK2 (151) vs HC (140).pdf", ont),
    width  = 5,
    height = 6
  )
  print(p)
  dev.off()
}


# Perform KEGG pathway enrichment
kegg_enrichment <- enrichKEGG(gene      = gene_entrez$ENTREZID,
                              organism  = 'hsa',
                              pvalueCutoff = 0.05, 
                              keyType = 'ncbi-geneid')

# Make results human-readable
kegg_res <- setReadable(kegg_enrichment, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

head(kegg_res)


pdf('The dotplot of KEGG enrichment analysis of DOWN_PD_LRRK2 (151) vs HC (140).pdf')
dotplot(kegg_res, x = "GeneRatio", showCategory = 15)
dev.off()
write.csv(kegg_enrichment, file = "The KEGG of DOWN_PD_LRRK2 (151) vs HC (140).csv")


# Gene Set Enrichment Analysis（GSEA）
sig.gen <- results
sig.gen$genesid <- sub("_.*$", "", rownames(results))

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(sig.gen$genesid, fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

# Find the genes that failed to map
failed_genes <- setdiff(sig.gen$genesid, gene_entrez$SYMBOL)

print("Failed to map genes:")
print(failed_genes)

genelist <- sig.gen$logFC
names(genelist) <- sig.gen$genesid
genelist <- sort(genelist, decreasing = TRUE)

library(dplyr)

df <- tibble(gene = names(genelist), score = genelist)

# collapse by taking the value with largest |score|
df2 <- df %>%
  group_by(gene) %>%
  slice_max(order_by = abs(score), n = 1) %>%
  ungroup()

genelist_unique <- df2$score
names(genelist_unique) <- df2$gene
genelist_unique <- sort(genelist_unique, decreasing = TRUE)

gsebp <- gseGO(
  geneList = genelist_unique,
  OrgDb    = org.Hs.eg.db,
  keyType  = "SYMBOL",
  ont      = "BP",
  pvalueCutoff = 0.05
)

pdf('The GSEA_BP dotplot.pdf', 5,6)
dotplot(gsebp, x = 'NES', showCategory = 20, title = "Top 20 GSEBP Terms")
dev.off()
write.csv(gsebp, file = 'The GSEA_BP genesets.csv')

gsemf <- gseGO(
  geneList = genelist_unique,
  OrgDb    = org.Hs.eg.db,
  keyType  = "SYMBOL",
  ont      = "MF",
  pvalueCutoff = 0.05
)

pdf('The GSEA_MF dotplot.pdf', 5,6)
dotplot(gsemf, x = 'NES', showCategory = 20, title = "Top 20 GSEMF Terms")
dev.off()
write.csv(gsemf, file = 'The GSEA_MF genesets.csv')

gsecc <- gseGO(
  geneList = genelist_unique,
  OrgDb    = org.Hs.eg.db,
  keyType  = "SYMBOL",
  ont      = "CC",
  pvalueCutoff = 0.05
)

pdf('The GSEA_CC dotplot.pdf', 5,6)
dotplot(gsecc, x = 'NES', showCategory = 20, title = "Top 20 GSECC Terms")
dev.off()
write.csv(gsecc, file = 'The GSEA_CC genesets.csv')
