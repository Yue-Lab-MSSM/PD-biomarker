# =============================================================================
# 09_mouse_urine_DAP.R
# Differential abundance and GO/KEGG enrichment of urine proteins
#
# Manuscript: "Integrated Biofluid Proteomics Identified Dynamic Functional
#              Biomarkers of LRRK2-Linked Parkinson's Disease Progression"
#
# For each comparison in COMPARISONS:
#   1. PCA of all samples (quality control)
#   2. limma: abundance ~ group; empirical Bayes moderated t test;
#      Benjamini-Hochberg FDR
#   3. Differentially abundant proteins (DAPs):
#      FDR/p-value < 0.05 and |log2FC| > 1.5 x SD of all log2FC in that comparison
#   4. Volcano plot
#   5. GO (BP, CC, MF) and KEGG over-representation for down- and
#      up-regulated DAPs
#
# Input : ./The mouse urine proteomic data 70% filtered, log2 transferred, normal dis imputation.csv
#           (protein-group abundances that were log2-transformed, filtered (detected in >70% of samples) and
#           imputed from a normal distribution in Perseus; see Methods)
#          ./Phenotype.txt
# Output: results/urine/<comparison>/
# =============================================================================






library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(gridExtra)

#=====================================================================================
#
#  Code chunk 1	Load Proteome data for missing imputation
#
#=====================================================================================
orignal_data <- read.csv("./The mouse urine proteomic data 70% filtered, log2 transferred, normal dis imputation.csv", header = T, row.names = 1)

head(orignal_data)


medians <- apply(orignal_data, 2, median, na.rm = TRUE)
max_med <- max(medians)
data <- sweep(orignal_data, 2, medians, FUN = "/") * max_med
head(data)


phenotypeData <- read.delim('Phenotype.txt', header = T, row.names = 1)

head(phenotypeData)

# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes
pdf("Before batch correction Cluster Dendrogram.pdf", 5,5)

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data_z <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

dev.off()

# pca - method 2
pdf("Before batch correction PCA of the samples.pdf", 5,5)
pca <- prcomp(t(data))
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



library(sva)

# Assuming `data` is your expression matrix and `phenoData$batch` contains batch information
batch <- phenotypeData$Batch

# Apply ComBat for batch effect correction
#data_corrected <- ComBat(dat = data, batch = batch, mod = NULL)
mod <- model.matrix(~ Genotype, data = phenotypeData)
data_corrected <- ComBat(dat = data, batch = phenotypeData$Batch, mod = mod)


data <- data_corrected

# Scale rows (NEW STEP)
data <- t(scale(t(data)))  # Standardize each protein across all samples
head(data)

pdf("After batch correction Cluster Dendrogram.pdf", 5,5)
# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

dev.off()

# pca - method 2
pdf("After batch correction PCA of the samples.pdf", 5,5)
pca <- prcomp(t(data))
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

#=====================================================================================
#
#  Code chunk 2	Load WGCNA package 
#
#=====================================================================================

# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Read in the proteomic data set		
ProData = as.data.frame(t(data))
# Take a quick look at what is in the data set:
dim(ProData);names(ProData);

#=====================================================================================
#
#  Code chunk 3	Load traits data
#
#=====================================================================================

traits <- phenotypeData
traitMat <- model.matrix(~ Genotype - 1, data = traits)
#=====================================================================================
#
#  Code chunk 4 SoftThreshold chosen
#
#=====================================================================================

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(ProData, powerVector = powers, verbose = 5, networkType = "signed")

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

softPower = 14

temp_cor <- cor
cor <- WGCNA::cor


#Test new settings
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
pdf("Module-trait relationship.pdf",5.5,8)

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
dissTOM = 1-TOMsimilarityFromExpr(ProData, power = softPower, networkType = "signed" );
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;

pdf("Network heatmap plot.pdf",9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot")
dev.off()

#=====================================================================================
#
#  Code chunk 14 Exporting Network to external visualization software
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
#  Code chunk 15 Exporting top hub Network to external visualization software
#
#=====================================================================================

module = c("brown");
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


module = c("grey");
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
#  Code chunk 16 Module membership vs. peptide significance
#
#=====================================================================================

pdf("GenotypeLRRK2_GS The Module membership vs peptide significance for mouse urine_grey.pdf",4,4)
module = "grey"								
column = match(module, modNames);
moduleGenes = moduleColors==module;
par(mfrow = c(1,1));
verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance[moduleGenes, 3],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Peptide significance for sPD",
                   main = paste("Module membership in sPD\n"),
                   cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, col = module,xlim=c(-1,1),ylim=c(-1,1))
dev.off(); 


pdf("The Module membership vs peptide significance for mouse urine_pink.pdf",4,4)
module = "pink"								
column = match(module, modNames);
moduleGenes = moduleColors==module;
par(mfrow = c(1,1));
verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance[moduleGenes, 4],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Peptide significance for sPD",
                   main = paste("Module membership in sPD\n"),
                   cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, col = module,xlim=c(-1,1),ylim=c(-1,1))
dev.off(); 


pdf("brown The Module membership vs peptide significance for mouse urine_brown.pdf.pdf")
module = "brown"								
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));

# Assume 'selectedGenes' contains the indices or names of the points to label.
selectedGenes <- c("C1QTNF1","GREM2","VIP","PSAP","ENTPD1_3182_38_2", "TMEM106A", "CTSB_3061_61_2", "CTSH_8465_52_3", 
                   "CTSO","CTSB_8007_19_3","CTSH_3737_6_3", 
                   "CTSE", "CTSZ", "SMPD1","MANBA", "GNS_8314_71_3", "IDS",
                   "ARSA_3583_54_4", "ACP2", "TMEM106B", "CLN5_8874_53_3", "LRP2", "TMEM9_9249_17_3", "PLD3",
                   "NAGPA", "GNPTG", "PLD3", "PLA2G7", "SMPD1", "LRP2" ) # Replace with your labels of interest.


selectedGenes <- c("Ctsz", "Ctsd", "Ctsa", "Ctsh", "Lamp2", "Psap", "Gba1") 


# Get the row names for the current module
moduleGeneNames <- rownames(geneModuleMembership)[moduleGenes]

# Identify selected genes that belong to the current module
selectedInModule <- moduleGeneNames[moduleGeneNames %in% selectedGenes]


verboseScatterplot(geneModuleMembership[moduleGenes, column],
                   geneTraitSignificance[moduleGenes, 2],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Peptide significance for LRRK2 PD",
                   main = paste("Module membership in LRRK2 PD  response\n"),
                   cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, col = module,xlim=c(-1,1),ylim=c(-1,1))

# Add labels to the selected genes
text(
  geneModuleMembership[selectedInModule, column],
  geneTraitSignificance[selectedInModule, 2],
  labels = selectedInModule, pos = 2, cex = 0.8, col = "red"
)
dev.off();


library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
library(topGO)

#############################################
WGCNA_data <- read.csv("All traits of WGCNA result.csv", row.names = 1)
id_data <- WGCNA_data[WGCNA_data$moduleColor == "pink" & WGCNA_data$GS.GenotypehLRRK2GS_LRRK2_KO > 0.1, ]
genes <- row.names(id_data)

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(genes, fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Mm.eg.db)

# Find the genes that failed to map
failed_genes <- setdiff(genes, gene_entrez$SYMBOL)


print("Failed to map genes:")
print(failed_genes)

library(stringr)
ego_cc<-enrichGO(gene       = gene_entrez$ENTREZID,
                 OrgDb      = org.Mm.eg.db,
                 keyType    = 'ENTREZID',
                 ont        = "CC",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = T)
head(ego_cc, 4)
ego_bp<-enrichGO(gene       = gene_entrez$ENTREZID,
                 OrgDb      = org.Mm.eg.db,
                 keyType    = 'ENTREZID',
                 ont        = "BP",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = T)
head(ego_bp, 4)

ego_mf<-enrichGO(gene       = gene_entrez$ENTREZID,
                 OrgDb      = org.Mm.eg.db,
                 keyType    = 'ENTREZID',
                 ont        = "MF",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = T)
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
    file = sprintf("The_GO_%s_of_cluster_proteins_in_pink hLRRK2GS_LRRK2_KO.csv", ont),
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
    sprintf("GO_%s_vertical_dotplot_cluster_proteins_in_pink hLRRK2GS_LRRK2_KO.pdf", ont),
    width  = 5,
    height = 6
  )
  print(p)
  dev.off()
}


# Perform KEGG pathway enrichment
kegg_enrichment <- enrichKEGG(gene      = gene_entrez$ENTREZID,
                              organism  =  'mmu', 
                              pvalueCutoff = 0.05, 
                              keyType = 'ncbi-geneid')

# 2. Make results human-readable
kegg_res <- setReadable(kegg_enrichment, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

head(kegg_res)


pdf('The dotplot of KEGG enrichment analysis of cluster_proteins_in_pink hLRRK2GS_LRRK2_KO .pdf', 4.5,4)
dotplot(kegg_res, x = "GeneRatio", showCategory = 15)
dev.off()
write.csv(kegg_res, file = "The KEGG of cluster_proteins_in_pink hLRRK2GS_LRRK2_KO.csv")

