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
# Input : data/ppmi_project_151.csv
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

data <- read.csv("./ppmi_project_151.csv", header = T, row.names = 1)

new_data <- data[data$GROUP == "HC" |data$GROUP == "PD_LRRK2",]   # for different group comparison, please change the group name here and below

to_remove <- c( ) # Please see the sample id in the excluded_ids.csv

new_data <- new_data[!(row.names(new_data) %in% to_remove), ]

table(new_data$GROUP)

table(new_data$GROUP)

sample_data <- new_data[, c(1:31)]

filtered_new_data <- new_data[, -c(1:31)]


colnames(filtered_new_data) <- sapply(strsplit(colnames(filtered_new_data), "_"), function(x) paste(x[3:length(x)], collapse = "_"))

# View the updated dataframe
head(filtered_new_data)


filtered_new_data <- t(filtered_new_data)

library(ggplot2)
library(stats)

pca_result <- prcomp(t(filtered_new_data), scale. = TRUE)
pca_df <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], condition = sample_data$GROUP)
pca_df$sample_names <- rownames(pca_df)  # Add sample names to the PCA dataframe

ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text(aes(label = sample_names), vjust = -1, size = 3) +  # Add sample names as labels
  labs(title = "PCA of Proteomic Data", x = "PC1", y = "PC2")

# Load limma package
library(limma)

# Create design matrix for the condition of interest   Control vs Treatment
sample_data$GROUP <- factor(sample_data$GROUP, levels = c("HC", "PD_LRRK2"))
#sample_data$GROUP <- factor(sample_data$GROUP)
design <- model.matrix(~ GROUP + AGE + SexM1.F0, data=sample_data)

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
write.csv(results, file = "Result of all proteins of PD_LRRK2 (149) vs HC(184).csv")


#get different expressed genes, condition is padj < 0.05, increased by 50% or decreased by 30%
SD <- sd(results$logFC)
DEPs <-subset(results, adj.P.Val < 0.05 & (logFC > 1.5*SD | logFC < - 1.5*SD))

dim(DEPs)
head(DEPs)
write.csv(DEPs,file= "Result of sigdiff proteins of PD_LRRK2 (149) vs HC(184).csv")


pdf('Volcanoplot.pdf', 7.5,4.5)
library(EnhancedVolcano)
EnhancedVolcano(results, lab = row.names(results), x = 'logFC', y = 'adj.P.Val', xlim = c(-0.8, 0.9), ylim=c(0,28),
                pCutoff = 0.05, FCcutoff = 1.5*SD, gridlines.major = F, gridlines.minor = F,
                legendLabels=c('NS','Log2FC','FDR','Log2FC & FDR'),
                xlab = bquote(~Log[2] ~ "fold change"),
                ylab = bquote(~-Log[10] ~ FDR),
                pointSize = 0.2,
                # NS = gray; Log2FC only = orange; FDR only = blue; both = red
                col = c("#B3B3B3", "#FF9900", "#0072FF", "#E60000"),
                colAlpha = 1,
                legendPosition = "right",
                legendLabSize = 14,
                legendIconSize = 2,
                legendDropLevels = F,
                drawConnectors = T, min.segment.length = 0.0001, 
                selectLab = c(" "),
                #selectLab = c(" HLA_DQA2", "ITGB2", "C1QTNF1", "GRN", "GPNMB_8240_207_3", "ENTPD1_3182_38_2", "TMEM106A", "GREM2", "LGALS9", "GAA","LCT"),
                title = NULL,subtitle = "LRRK2_PD vs. HC",
                max.overlaps = 10)
dev.off() 


# File paths for the three lists
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
gene_entrez <- bitr(Down_genes, fromType = "SYMBOL", #change to Up_genes for enriched pathways
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

# 2. define your frameless, vertical‐strip theme with big text
strip_theme <- theme(
  panel.border     = element_blank(),
  panel.background = element_blank(),
  panel.grid       = element_blank(),
  
  # y‐axis labels (the GO terms)
  axis.text.y      = element_text(size = 14),  
  axis.title.y     = element_text(size = 16),   
  
  # legend (color/size scales, if shown)
  legend.title     = element_text(size = 14),
  legend.text      = element_text(size = 12),
  
  # remove x‐axis entirely
  axis.title.x     = element_blank(),
  axis.text.x      = element_blank(),
  axis.ticks.x     = element_blank()
)

for (ont in names(ego_list)) {
  ego_obj <- ego_list[[ont]]
  
  # 3. write out the table
  write.csv(
    as.data.frame(ego_obj),
    file = sprintf("The_GO_%s_of_enriched_proteins_in_DOWN_LRRK2 PD_(155)_vs_HC(184).csv", ont),
    row.names = FALSE
  )
  
  # 4. make a dotplot with larger text  
  p <- dotplot(
    ego_obj,
    showCategory = 10,
    x            = "GeneRatio",
    font.size    = 12         # increases the label text inside the plot
  ) +
    aes(x = 1) +              # collapse to a vertical strip
    scale_x_continuous(NULL, breaks = NULL) +
    strip_theme
  
  # 5. save as PDF
  pdf(
    sprintf("GO_%s_vertical_dotplot_DOWN_LRRK2 PD_(155)_vs_HC(184).pdf", ont),
    width  = 5,
    height = 6
  )
  print(p)
  dev.off()
}


# Perform KEGG pathway enrichment
kegg_enrichment <- enrichKEGG(gene      = gene_entrez$ENTREZID,
                              organism  = 'hsa',  # KEGG organism code for humans
                              pvalueCutoff = 0.05, 
                              keyType = 'ncbi-geneid')

# 2. Make results human-readable
kegg_res <- setReadable(kegg_enrichment, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

head(kegg_res)


pdf('The dotplot of KEGG enrichment analysis of DOWN_LRRK2 PD_(155)_vs_HC(184).pdf')
dotplot(kegg_res, x = "GeneRatio", showCategory = 15)
dev.off()
write.csv(kegg_enrichment, file = "The KEGG of DOWN_LRRK2 PD_(155)_vs_HC(184).csv")

