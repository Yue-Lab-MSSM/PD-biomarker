# =============================================================================
# 10_urine_WGCNA_PPMI.R
# Weighted gene co-expression network analysis (WGCNA) of urine proteins
# (PPMI Project 190, MS)
#
# Manuscript: "Integrated Biofluid Proteomics Identified Dynamic Functional
#              Biomarkers of LRRK2-Linked Parkinson's Disease Progression"
#
# Steps
#   1. Load data, define groups (HC, LRRK2, GBA, sPD)
#   2. Load traits data and Quality control: remove poor proteins/samples and outlier samples
#      (hierarchical clustering of samples, tree cut at height 40)
#   3. Choose the soft-thresholding power (scale-free topology)
#   4. WGCNA construction and detect modules (blockwiseModules)
#   5. Clustering dendrograms for peptides
#   6. Modules analysis
#   7. Peptide significance caculation
#   8. Module-trait relationships
#   9. Correlated eigenpeptide in modules associated with traits
#   10. WGCNA results Output
#   11. Network heatmap plot visulization
#   12. Exporting Network to external visualization software
#   13. Exporting top hub Network to external visualization software (cytoscope)
#   14. Module membership vs. peptide significance
#   15. GO TERM and KEGG analysis
#
# Input : ./ppmi_project_190.csv
#           (first 31 columns = phenotypes; remaining columns = log2 MS
#            abundances)
# Output: results/WGCNA/
# =============================================================================

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(gridExtra)

setwd(".")


#=====================================================================================
#
#  Code chunk 1	Load raw data
#
#=====================================================================================

protein_matrix <- read.csv("./ppmi_project_190.csv.csv", header = T, row.names = 1)
head(protein_matrix, n=20)
protein_matrix <- protein_matrix[protein_matrix$EVENT_ID == "BL", ]
row.names(protein_matrix) <- protein_matrix$PATNO
table(protein_matrix$GROUP)

protein_matrix$GROUP <- with(protein_matrix,
                             ifelse(GROUP %in% c("PD_LRRK2", "ProLRRK2"), "PD_LRRK2",
                                    ifelse(GROUP %in% c("PD_GBA",   "ProGBA"),   "PD_GBA",
                                           ifelse(GROUP %in% c("HC"),                    "HC",
                                                  ifelse(GROUP %in% c("sPD"),                   "sPD",
                                                         NA)))))

table(protein_matrix$GROUP)
protein_matrix <- protein_matrix[protein_matrix$GROUP %in% c('HC','PD_LRRK2', 'PD_GBA', 'sPD'), ]
table(protein_matrix$GROUP)

filtered_expression_data <- protein_matrix[, -c(1:31)]

# View the updated dataframe
data <- filtered_expression_data
dim(data)


#=====================================================================================
#
#  Code chunk 2	Load traits data
#
#=====================================================================================

traits <- protein_matrix[, c(3,16)]

# Generate dummy variables for each group
traitMat <- model.matrix(~ GROUP - 1, data = traits)
# This yields columns:
#   GROUPHC    (1 if HC, else 0)   GROUPPD_GBA   GROUPPD_LRRK2


# Ensure rownames match samples
rownames(traitMat) <- rownames(traits)


# 1. QC - outlier detection ------------------------------------------------
# detect outlier genes the column should be gene names, the rows should be samples.

htree <- hclust(dist(data), method = "average")
plot(htree)
# Cut tree at a certain height (adjust height as needed, e.g., 60)
clusters <- cutree(htree, h = 40)
table(clusters)
# Identify samples in very small clusters (potential outliers)
for(i in c(2:13)){
  # Extract sample names belonging to cluster i
  outlier_samples <- names(clusters[clusters == i])
  
  # Create a file name dynamically for each cluster
  file_name <- paste0("outlier_samples_cluster_", i, ".csv")
  
  # Save the outlier sample names to a CSV file
  write.csv(outlier_samples, file = file_name, row.names = FALSE)
  
}


# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes, should row samples, columne gene/protein names
pdf("Cluster Dendrogram.pdf", 5,4)

gsg <- goodSamplesGenes(data)
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data_z <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(data), method = "average")
plot(htree)

dev.off()

# pca - method 2 should row samples, columne gene/protein names
pdf("PCA of the samples.pdf", 5,4)
pca <- prcomp(data)
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
dev.off()



# Assuming outlier_samples contains the row names of the outliers
outlier_samples <- c("41438", "41664", "41854", "41989", "42293", "42308", "42379",
                     "42415", "42447", "43046", "51544", "51551")


# Remove outliers from the dataset
data <- filtered_expression_data[!rownames(filtered_expression_data) %in% outlier_samples, ]

dim(data)

# Get the row names of the expression dataset
expression_samples <- rownames(data)

# Subset phenotypeData to include only those rows that are in the expression dataset
traitMat <- traitMat[rownames(traitMat) %in% expression_samples, ]

# Check if the row names now match
all(rownames(traitMat) %in% rownames(data))  # Should return TRUE



#=====================================================================================
#
#  Code chunk 3	Load WGCNA package 
#
#=====================================================================================

# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Read in the proteomic data set		
ProData = as.data.frame(data)
# Take a quick look at what is in the data set:
dim(ProData) #should row samples, columne gene/protein names

names(ProData)



#=====================================================================================
#
#  Code chunk 4 SoftThreshold chosen
#
#=====================================================================================

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(ProData, powerVector = powers, verbose = 5, networkType = "signed")
# Plot the results:
pdf("Topology model and mean connectivity.pdf",10,5)

cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#=====================================================================================
#
#  Code chunk 5	WGCNA construction
#
#=====================================================================================

softPower = 12

temp_cor <- cor
cor <- WGCNA::cor

net = blockwiseModules(ProData, power = softPower, 
                       networkType = "signed",
                       minModuleSize = 30,
                       deepSplit        = 3,
                       reassignThreshold = 1e-3, 
                       mergeCutHeight = 0.2,
                       numericLabels = F, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       randomSeed = 1234,
                       maxBlockSize = 20000,
                       verbose = 3)


cor <- temp_cor

table(net$colors)

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- net$MEs


# Print out a preview
head(module_eigengenes)

write.csv(module_eigengenes, file = "The module eigengenes.csv", row.names = TRUE)
# get number of genes for each module
table(net$colors)
modulegenes <- as.data.frame(net$colors)

write.csv(modulegenes, file = "The genes in each module.csv", row.names = TRUE)
write.csv(table(net$colors), file = "The gene numbers in each module.csv", row.names = TRUE)


#=====================================================================================
#
#  Code chunk 6	Clustering dendrograms for peptides
#
#=====================================================================================


pdf("Cluster Dendrogram.pdf",12,9)

mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], cbind(net$unmergedColors, net$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()					


#=====================================================================================
#
#  Code chunk 7	Modules analysis
#
#=====================================================================================

moduleLabels = net$colors
moduleColors = net$colors
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

#=====================================================================================
#
#  Code chunk 8	Peptide significance caculation
#
#=====================================================================================


# Define numbers of genes and samples
nGenes = ncol(ProData);
nSamples = nrow(ProData);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(ProData, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traitMat, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#=====================================================================================
#
#  Code chunk 9	Module-trait relationships
#
#=====================================================================================

# open a graphics window
pdf("Module-trait relationship.pdf",5,8)

textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traitMat),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


dev.off()

#=====================================================================================
#
#  Code chunk 10 Correlated eigenpeptide in modules associated with traits
#
#=====================================================================================

##Compute Gene–Trait Significance (GS) for all traits, not only only trait

modNames = substring(names(MEs), 3)

# Calculate gene module membership: correlation between gene expression and module eigengenes
geneModuleMembership <- as.data.frame(cor(ProData, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

# Rename columns for clarity
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

# Calculate gene-trait significance for all traits
geneTraitSignificance <- as.data.frame(cor(ProData, traitMat, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

# Rename columns so that each trait’s GS is clearly labeled
names(geneTraitSignificance) <- paste("GS.", colnames(traitMat), sep = "")
names(GSPvalue) <- paste("p.GS.", colnames(traitMat), sep = "")


#=====================================================================================
#
#  Code chunk 11	WGCNA results Output 
#
#=====================================================================================

# Create a starting data frame with peptide names, module colors, and gene-trait significance info
geneInfo <- data.frame(
  Peptides = names(ProData),
  moduleColor = moduleColors,
  geneTraitSignificance,
  GSPvalue
)

# Optionally, add module membership information for each module
for (mod in 1:ncol(geneModuleMembership)) {
  geneInfo <- cbind(geneInfo,
                    geneModuleMembership[, mod],
                    MMPvalue[, mod])
  # Rename new columns accordingly:
  names(geneInfo)[(ncol(geneInfo)-1):ncol(geneInfo)] <- c(
    paste("MM.", modNames[mod], sep = ""),
    paste("p.MM.", modNames[mod], sep = "")
  )
}

# Order the genes by module and by significance for one of the traits (if needed)
geneOrder <- order(geneInfo$moduleColor, -abs(geneInfo[, grep("GS.", names(geneInfo))[1]]))
geneInfo <- geneInfo[geneOrder, ]

# Save the combined results if desired
write.csv(geneInfo, file = "All traits of WGCNA result.csv", row.names = FALSE)



#=====================================================================================
#
#  Code chunk 12	Network heatmap plot visulization
#
#=====================================================================================

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(ProData, power = softPower, networkType = "signed" );
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function

pdf("Network heatmap plot.pdf",9,9)
# sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot")
dev.off()

#=====================================================================================
#
#  Code chunk 13 Exporting Network to external visualization software
#
#=====================================================================================

# Cytoscape or VisANT plot and software
probes = names(ProData)
TOM=1-dissTOM
modules = c(moduleColors[!duplicated(moduleColors)])

for (i in seq(1, nrow(table(modules)))){
  inModule=is.finite(match(moduleColors,modules[i]))
  modProbes=probes[inModule]
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile=paste("Proteome network CytoEdge ",paste(modules[i],collapse="-"),".txt",sep=""),
                                 nodeFile=paste("Proteome network CytoNode ",paste(modules[i],collapse="-"),".txt",sep=""),
                                 weighted = TRUE, threshold = 0.02)
  #  vis = exportNetworkToVisANT(modTOM,
  #                              file=paste("Example-1 Proteome network VisANTInput ", modules[i],  sep=""),
  #                              weighted=TRUE,threshold = 0.02)
}

#=====================================================================================
#
#  Code chunk 14 Exporting top hub Network to external visualization software
#
#=====================================================================================

# choose red module
module = c("magenta");
inModule = (moduleColors==module);
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# choose the top 10 node with highest conncetion
nTop = 10;
IMConn = softConnectivity(ProData[, probes[is.finite(match(moduleColors,module))]]);
top = (rank(-IMConn) <= nTop)

# network file for cystoscape software 
cyt = exportNetworkToCytoscape(modTOM[top, top],
                               edgeFile=paste("Proteome network CytoEdge ",paste(module,collapse="-"),"-top10.txt",sep=""),
                               nodeFile=paste("Proteome network CytoNode ",paste(module,collapse="-"),"-top10.txt",sep=""),
                               weighted = TRUE, threshold = 0.21)


# choose pink module
module = c("pink");
inModule = (moduleColors==module);
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# choose the top 10 node with highest conncetion
nTop = 10;
IMConn = softConnectivity(ProData[, probes[is.finite(match(moduleColors,module))]]);
top = (rank(-IMConn) <= nTop)
cyt = exportNetworkToCytoscape(modTOM[top, top],
                               edgeFile=paste("Proteome network CytoEdge ",paste(module,collapse="-"),"-top10.txt",sep=""),
                               nodeFile=paste("Proteome network CytoNode ",paste(module,collapse="-"),"-top10.txt",sep=""),
                               weighted = TRUE, threshold = 0.11)


#=====================================================================================
#
#  Code chunk 15 Module membership vs. peptide significance
#
#=====================================================================================

# choose red module vs. Genotype_Treatment

pdf("The Module membership vs peptide significance for sPD_green.pdf",4,4)
module = "green"								
column = match(module, modNames);
moduleGenes = moduleColors==module;
par(mfrow = c(1,1));
verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance[moduleGenes, 3],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Peptide significance for sPD",
                   main = paste("Module membership vs. sPD\n"),
                   cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, col = module,xlim=c(-1,1),ylim=c(-1,1))
# If plotting into a file, close it
dev.off(); 

# choose pink module vs. Genotype_Treatment

pdf("The Module membership vs peptide significance for LRRK2 PD_black.pdf")
module = "black"								
column = match(module, modNames);
moduleGenes = moduleColors==module;
par(mfrow = c(1,1));
verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance[moduleGenes, 4],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Peptide significance for LRRK2 PD",
                   main = paste("Module membership vs. LRRK2 PD\n"),
                   cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, col = module,xlim=c(-1,1),ylim=c(-1,1))
# If plotting into a file, close it
dev.off(); 

# choose tan module vs. Genotype_Treatment 

pdf("The Module membership vs peptide significance for sPD_blue.pdf",4,4)
module = "blue"								
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
# pdf(file = "Plots/MMvsGS-1.pdf", width=7, height=7);
par(mfrow = c(1,1));

# Assume 'selectedGenes' contains the indices or names of the points to label.
selectedGenes <- c("ENPEP", "MAN2B1", "SMPD1", "CTSC", "TMEM106B","CTSA","PPT1", "CTSD", "PLD3", "CTSO",
                   "MAN2B2", "TPP1_002", "LRP2") # Replace with your labels of interest.


selectedGenes <- c("LAMB2", "COL12A1", "NGFR_003", "ADGRL1_002", "AGER", "LFNG", "HAMP", "PRNP") 


# Get the row names for the current module
moduleGeneNames <- rownames(geneModuleMembership)[moduleGenes]

# Identify selected genes that belong to the current module
selectedInModule <- moduleGeneNames[moduleGeneNames %in% selectedGenes]


verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance[moduleGenes, 4],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Peptide significance for sPD",
                   main = paste("Module membership vs. sPD  response\n"),
                   cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, col = module,xlim=c(-1,1),ylim=c(-1,1))

# Add labels to the selected genes
text(
  geneModuleMembership[selectedInModule, column],
  geneTraitSignificance[selectedInModule, 4],
  labels = selectedInModule, pos = 2, cex = 0.8, col = "red"
)
# If plotting into a file, close it
dev.off();



library(dplyr)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)

#############################################
WGCNA_data <- read.csv("All traits of WGCNA result.csv", row.names = 1)
id_data <- WGCNA_data[WGCNA_data$moduleColor == "green" & WGCNA_data$GS.GROUPsPD > 0.1, ]
genes <- id_data$X

#############################################


# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(genes, fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

# Find the genes that failed to map
failed_genes <- setdiff(genes, gene_entrez$SYMBOL)

# Print the failed genes
print("Failed to map genes:")
print(failed_genes)

library(stringr)
ego_cc<-enrichGO(gene       = gene_entrez$ENTREZID,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'ENTREZID',
                 ont        = "CC",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = T)
head(ego_cc, 4)
ego_bp<-enrichGO(gene       = gene_entrez$ENTREZID,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'ENTREZID',
                 ont        = "BP",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = T)
head(ego_bp, 4)

ego_mf<-enrichGO(gene       = gene_entrez$ENTREZID,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'ENTREZID',
                 ont        = "MF",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = T)
head(ego_mf, 4)


# 1. collect your GO objects
ego_list <- list(
  CC = ego_cc,
  #BP = ego_bp,
  MF = ego_mf
)

# 2. define your frameless, vertical‐strip theme with big text
strip_theme <- theme(
  panel.border     = element_blank(),
  panel.background = element_blank(),
  panel.grid       = element_blank(),
  
  # y‐axis labels (the GO terms)
  axis.text.y      = element_text(size = 14),      # bump category names
  axis.title.y     = element_text(size = 16),      # if you ever add a y‐title
  
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
    file = sprintf("The_GO_%s_of_cluster_proteins_in_green sPD.csv", ont),
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
    sprintf("GO_%s_vertical_dotplot_cluster_proteins_in_black_green sPD.pdf", ont),
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


pdf('The dotplot of KEGG enrichment analysis of cluster_proteins_in_green sPD.pdf', 4.5,4)
dotplot(kegg_res, x = "GeneRatio", showCategory = 15)
dev.off()
write.csv(kegg_enrichment, file = "The KEGG of cluster_proteins_in_green sPD.csv")



