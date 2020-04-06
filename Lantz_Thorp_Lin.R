#Load Necessary Packages and set working directory ##########################
library(pkgbuild)
devtools::install_github(repo = "satijalab/seurat", ref = "develop")
library(ggplot2)
library(dplyr)
library(Matrix)
library(sctransform)
library(cowplot)
library(monocle)
library(Seurat)
library(umap)
library(RCurl)
library(scales)
library(tidyverse)
library(MAST)
library(RColorBrewer)
library(destiny)
library(SingleCellExperiment)
library(ggbeeswarm)
library(ggthemes)
#BiocManager::install("slingshot")
library(slingshot)
#BiocManager::install("clusterExperiment")
library(clusterExperiment)
library(gam)
library(scater)

#Set working directory and retreive files
setwd("C:/Users/conno/Documents/Driskill Graduate Program/Projects/Efferocytosis Containaminating RNA/scRNAseq/MANUSCRIPT_scRNAseq")

#Creating Seurat Object for each sample############
for (file in c("WT_Ctrl", "WT_2_hr", "WT_6_hr", "MerKD_Ctrl", "MerKD_2hr", "MerKD_6_hr")){
  seurat_data <- Read10X(data.dir = paste0("Data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, min.features = 100, project = file)
  assign(file, seurat_obj)
}

#Merge Seurat Object for QC
merged_seurat <- merge(x = WT_Ctrl, y = c(WT_2_hr, WT_6_hr, MerKD_Ctrl, MerKD_2hr, MerKD_6_hr), add.cell.ids = c("WT_Ctrl", "WT_2_hr", "WT_6_hr", "MerKD_Ctrl", "MerKD_2hr", "MerKD_6_hr"))


#Added metadata for QC analyses
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
merged_seurat$mitoRatio <- PercentageFeatureSet(merged_seurat, pattern = "^mt-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
merged_seurat$riboRatio <- PercentageFeatureSet(merged_seurat, pattern = "^Rp[ls]")
merged_seurat$riboRatio <- merged_seurat@meta.data$riboRatio / 100
metadata <- merged_seurat@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^WT_Ctrl"))] <- "WT_Ctrl"
metadata$sample[which(str_detect(metadata$cells, "^WT_2_hr"))] <- "WT_2_hr"
metadata$sample[which(str_detect(metadata$cells, "^WT_6_hr"))] <- "WT_6_hr"
metadata$sample[which(str_detect(metadata$cells, "^MerKD_Ctrl"))] <- "MerKD_Ctrl"
metadata$sample[which(str_detect(metadata$cells, "^MerKD_2hr"))] <- "MerKD_2hr"
metadata$sample[which(str_detect(metadata$cells, "^MerKD_6_hr"))] <- "MerKD_6_hr"
merged_seurat@meta.data <- metadata
save(merged_seurat, file="Data/merged_filtered_seurat.RData")

#Quality Metrics###############
#NCells
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

#nUMIs
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1000)

#nGenes
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 500)

#nGenes Boxplot
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

#Correlation of genes v. UMI
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 750) +
  geom_hline(yintercept = 500) +
  facet_wrap(~sample)

#mito counts ratio
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.1) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.1)

#ribo counts ratio
metadata %>% 
  ggplot(aes(color=sample, x=riboRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

#novelty gene
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

#Filtering Cells ##############################################################
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 750) & 
                            (nGene >= 500) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.10))

#Keep genes expressed in 10 or more cells
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 0
filtered_counts <- counts[keep_genes, ]
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
save(filtered_seurat, file="Data/seurat_filtered.RData")


#Cell Cycle Scoring ##############################################################
cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv") 
cell_cycle_genes <- read.csv(text = cc_file)


Seurat_CellCycle <- NormalizeData(filtered_seurat)
Seurat_CellCycle <- CellCycleScoring(Seurat_CellCycle,
                                     g2m.features = cc.genes$g2m.genes,
                                     s.features = cc.genes$s.genes
)
#View(Seurat_CellCycle@meta.data)
Seurat_CellCycle <- FindVariableFeatures(Seurat_CellCycle, verbose = F)
Seurat_CellCycle <- ScaleData(Seurat_CellCycle)
Seurat_CellCycle <- RunPCA(Seurat_CellCycle)
DimPlot(Seurat_CellCycle, reduction = "pca", group.by = "Phase", split.by = "Phase")

#SCTransform #######################################################################
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")
split_seurat <- split_seurat[c("WT_Ctrl", "WT_2_hr", "WT_6_hr", "MerKD_Ctrl", "MerKD_2hr", "MerKD_6_hr")]
options(future.globals.maxSize= 993718400)
#scTransform data without regressing out mitochondrial genes
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], verbose = F)
}

#Integration of DataSets ###########################################################
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)
split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = integ_features, verbose = F)
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, normalization.method = "SCT", anchor.features = integ_features, verbose = F)
Integrated_Seurat <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT", verbose = F)

#Run PCA on integrated dataset
Integrated_Seurat <- RunPCA(object = Integrated_Seurat)
PCAPlot(Integrated_Seurat, split.by = "sample", verbose = F)
DimHeatmap(Integrated_Seurat, dims = 1:9, cells = 500, balanced = T)
print(x = Integrated_Seurat[["pca"]], dims = 1:10, nfeatures = 5)
ElbowPlot(object = Integrated_Seurat, ndims = 40)

#Run UMAP on integrated dataset
Integrated_Seurat <- RunUMAP(Integrated_Seurat, dims = 1:40)
DimPlot(Integrated_Seurat)
DimPlot(Integrated_Seurat, split.by = "sample")

#Clustering################################
Integrated_Seurat <- FindNeighbors(object = Integrated_Seurat, dims = 1:40)
Integrated_Seurat <- FindClusters(object = Integrated_Seurat, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.4))

#Resolution 0.2
Idents(object = Integrated_Seurat) <- "integrated_snn_res.0.2"
Integrated_Seurat <- RunUMAP(Integrated_Seurat, reduction = "pca", dims = 1:40)
DimPlot(Integrated_Seurat, reduction = "umap", label = F)
DimPlot(Integrated_Seurat, reduction = "umap", label = F, cols = c("firebrick", "goldenrod", "darkolivegreen3", "deepskyblue3", "darkslategray4", "forestgreen", "orchid", "darkorchid4", "navy"))

#Resolution 0.4
Idents(object = Integrated_Seurat) <- "integrated_snn_res.0.4"
Integrated_Seurat <- RunUMAP(Integrated_Seurat, reduction = "pca", dims = 1:40)
DimPlot(Integrated_Seurat, reduction = "umap", label = T, label.size = 6)
DimPlot(Integrated_Seurat, label = T, split.by = "sample") + NoLegend()

#Resolution 0.6
Idents(object = Integrated_Seurat) <- "integrated_snn_res.0.6"
Integrated_Seurat <- RunUMAP(Integrated_Seurat, reduction = "pca", dims = 1:40)
DimPlot(Integrated_Seurat, reduction = "umap", label = T, label.size = 6)
DimPlot(Integrated_Seurat, label = T, split.by = "sample") + NoLegend()

#Revisit QC metrics
metrics <- c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")
FeaturePlot(Integrated_Seurat, reduction = "umap", features = metrics, pt.size = 0.4, sort.cell = T, min.cutoff = "q10", label = T)

#SingleR#####################################
library(SingleR)
immgen.se <- ImmGenData()

#All cell main labels
counts <- GetAssayData(Integrated_Seurat)
pred.integrated <- SingleR(test = counts, ref = immgen.se, labels = immgen.se$label.main)
Main.labels <- table(pred.integrated$labels)
all.markers <- metadata(pred.integrated)$de.genes

#SingleR QC
to.remove <- pruneScores(pred.integrated)
summary(to.remove)
plotScoreDistribution(pred.integrated, show = "delta.med", ncol = 3, show.nmads = 3)

#Label Seurat object with SingleR annotations
Integrated_Seurat[["SingleR.labels"]] <- pred.integrated$labels
DimPlot(Integrated_Seurat, group.by = "SingleR.labels", cols = c("royalblue", "forestgreen", "darkolivegreen3", "deepskyblue3", "gray", "navy", "orchid", "darkorchid4", "firebrick", "goldenrod", "azure4", "black", "darkgoldenrod4", "olivedrab1",  "tan2", "azure3", "forestgreen", "salmon"))
SingleR_markers <- metadata(pred.integrated)$de.genes

#Create Heatmap based on SingleR classifications
plotScoreHeatmap(pred.integrated)


saveRDS(Integrated_Seurat, "results/Integrated_Seurat.rds")


#Identification of Macrophage Populations at steady state ##########################
Integrated_Seurat <- readRDS("results/Integrated_Seurat.rds")

#Low resolution
Idents(object = Integrated_Seurat) <- "integrated_snn_res.0.2"
DimPlot(Integrated_Seurat, label = F)

DotPlot(Integrated_Seurat, features = c("Itgam", "Adgre1", "Icam2", "Itgax", "H2-Ab1", "Cd19", "Cd79a", "Cr2", "Cd3e", "Cd4", "Cd8a"))

FeaturePlot(Integrated_Seurat, features = c("Itgam", "Adgre1", "Icam2", "Itgax", "H2-Ab1", "Mertk"), order = T)

#Subset macrophages/DCs and then only WT at steady state
Integrated_Seurat_MacsDCs <- SubsetData(Integrated_Seurat, ident.use = c("0", "3"))

Idents(Integrated_Seurat_MacsDCs) <- "sample"
WT_Control <- SubsetData(Integrated_Seurat_MacsDCs, ident.use = "WT_Ctrl")

Idents(Integrated_Seurat_MacsDCs) <- "integrated_snn_res.0.4"
DimPlot(Integrated_Seurat_MacsDCs)
Clus5.markers <- FindConservedMarkers(Integrated_Seurat_MacsDCs, ident.1 = "5", grouping.var = "sample")
head(Clus5.markers)
write.csv(Clus5.markers, file = "./Results/Clus5.markers.csv", row.names = TRUE, quote = FALSE)

#FeaturePlot(Integrated_Seurat_MacsDCs, features = c("F5", "Serpinb2", "Prg4", "Cd93", "Selp", "Anxa5", "Prdx1", "Saa3", "Fabp7"), order = T, pt.size = 2, ncol = 3)

#Find Markers for Macrophage Clusters (Remove Cluster 12 due to enrichment of cell cycle genes)
Idents(WT_Control) <- "integrated_snn_res.0.4"
DimPlot(WT_Control, label = T)
WT_Control_Macs <- SubsetData(WT_Control, ident.use = c("1", "3", "5", "6", "11"))

DefaultAssay(WT_Control_Macs) <- "RNA"
WT_Control_Macs  <- NormalizeData(WT_Control_Macs)
WT_Control_Macs_markers <- FindAllMarkers(WT_Control_Macs , only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
#write.csv(WT_Control_Macs_markers, file = "results/WT_Control_Macs_markers.csv")

WT_Control_Macs_markers_top50 <- WT_Control_Macs_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
#write.csv(WT_Control_Macs_markers_top50, file = "results/WT_Control_Macs_markers_top50.csv")

#Heatmap in Figure 1C
DoHeatmap(WT_Control_Macs, features = WT_Control_Macs_markers_top50$gene, slot = "scale.data", assay = "SCT", group.colors = c("royalblue3", "firebrick", "deepskyblue2", "goldenrod", "darkorchid"),  disp.min = -2) + NoLegend()

#Remove outlier cells for visualization
Macsplot <- DimPlot(WT_Control)
Macs_selection <- CellSelector(Macsplot)
Idents(WT_Control, cells = Macs_selection) <- "Macrophages"
WT_Control_vis <- SubsetData(WT_Control, ident.use = "Macrophages")
DimPlot(WT_Control_vis)

saveRDS(WT_Control_vis, "results/WT_Control_vis.rds")

#Dimplot in Figure 1B
Idents(object = WT_Control_vis) <- "integrated_snn_res.0.4"
DimPlot(WT_Control_vis, cols = c("royalblue3", "firebrick", "deepskyblue2", "goldenrod", "darkorchid", "gray"), pt.size = 2)

#Peritoneal Mac Markers - Feature Plot in Figure 1D
DefaultAssay(WT_Control_vis) <- "RNA"

FeaturePlot(WT_Control_vis, features = c("Itgam", "Adgre1", "Icam2", "Timd4", "Klf2", "Mertk", "Itgax", "H2-Ab1", "Ccr2"), order = T, pt.size = 2, ncol = 3)

#Subset Macrophages further in WT Steady State Condition
Idents(WT_Control) <- "integrated_snn_res.0.4"
WT_Control_Macs <- SubsetData(WT_Control, ident.use = c("1", "3", "5", "6"))

#Cytokine and Marker Gene Expression - VlnPlot in Figure 1E
DefaultAssay(WT_Control_Macs) <- "SCT"
VlnPlot.Macs <- VlnPlot(WT_Control_Macs, features = c("Klf2", "Gata6", "Cd93", "Csf1r", "Anxa5", "Cd36", "Actb", "Arg1", "Irf4", "Retnla", "Cd74", "Tgfbi",  "Il1b", "Il6", "Tnf", "Cxcl2"),  cols = c("royalblue3", "firebrick", "deepskyblue2", "goldenrod"), pt.size = 0, combine = F)
for(i in 1:length(VlnPlot.Macs)) {
  VlnPlot.Macs[[i]] <- VlnPlot.Macs[[i]] + NoLegend() + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
}
cowplot::plot_grid(plotlist = VlnPlot.Macs, ncol = 4)


#Cluster Comparisons for gProfiler###############################
Mac1_3 <- FindMarkers(WT_Control_Macs, ident.1 = c("1"), ident.2 = "3", test.use = "MAST", min.pct = 0.25, logfc.threshold = 0.25)
Mac1_3.Markers <- subset(Mac1_3, Mac1_3$p_val_adj < 0.05)
write.csv(Mac1_3.Markers , file = "./Results/Mac1_3Markers.csv", row.names = TRUE, quote = FALSE)

Mac1_5 <- FindMarkers(WT_Control_Macs, ident.1 = c("1"), ident.2 = "5", test.use = "MAST", min.pct = 0.25, logfc.threshold = 0.25)
Mac1_5.Markers <- subset(Mac1_5, Mac1_5$p_val_adj < 0.05)
write.csv(Mac1_5.Markers, file = "./Results/Mac1_5Markers.csv", row.names = TRUE, quote = FALSE)

Mac1_6 <- FindMarkers(WT_Control_Macs, ident.1 = "1", ident.2 = "6", test.use = "MAST", min.pct = 0.25, logfc.threshold = 0.25)
Mac1_6.Markers <- subset(Mac1_6, Mac1_6$p_val_adj < 0.05)
write.csv(Mac1_6.Markers, file = "./Results/Mac1_6Markers .csv", row.names = TRUE, quote = FALSE)

Mac5_6 <- FindMarkers(WT_Control_Macs, ident.1 = "5", ident.2 = "6", test.use = "MAST", min.pct = 0.25, logfc.threshold = 0.25)
Mac5_6.Markers <- subset(Mac5_6, Mac5_6$p_val_adj < 0.05)
write.csv(Mac5_6.Markers, file = "./Results/Mac5_6Markers .csv", row.names = TRUE, quote = FALSE)


#Wild Type Efferocytosis ################################################

#Isolate Wild Type Macrophages from all time points
Idents(Integrated_Seurat) <- "integrated_snn_res.0.4"
Integrated_Macs <- SubsetData(Integrated_Seurat, ident.use = c("1", "3", "5", "6"))

Integrated_Macsplot <- DimPlot(Integrated_Macs)
Integrated_Macs_selection <- CellSelector(Integrated_Macsplot)
Idents(Integrated_Macs, cells = Integrated_Macs_selection) <- "Macrophages"
Integrated_Macs <- SubsetData(Integrated_Macs, ident.use = "Macrophages")
DimPlot(Integrated_Macs)

saveRDS(Integrated_Macs, "results/Integrated_Macs_vis")

#Isolate Wild Type Resident Peritoneal Macrophages
Idents(Integrated_Macs) <- "sample"
Integrated_Macs_WT <- SubsetData(Integrated_Macs, ident.use = c("WT_Ctrl", "WT_2_hr", "WT_6_hr"))

DimPlot(Integrated_Macs_WT, group.by = "sample", cols = c("grey", "grey", "blue"), pt.size = 3, order = T)
DimPlot(Integrated_Macs_WT, group.by = "sample", cols = c("grey", "blue", "grey"), pt.size = 3, order = F)
DimPlot(Integrated_Macs_WT, group.by = "sample", cols = c("blue", "grey", "grey"), pt.size = 3, order = F)

#Isolate Macrophage Subpopulations
Idents(Integrated_Macs_WT) <- "integrated_snn_res.0.4"
Integrated_Macs_WT_MHCIILO <- SubsetData(Integrated_Macs_WT, ident.use = c("1", "3", "5"))
Integrated_Macs_WT_MHCIIHI <- SubsetData(Integrated_Macs_WT, ident.use = "6")

#Create single cell experiment objects
Integrated_Macs_WT_MHCIILO_SCE <- as.SingleCellExperiment(Integrated_Macs_WT_MHCIILO)
table(Integrated_Macs_WT_MHCIILO_SCE$sample)

Integrated_Macs_WT_MHCIIHI_SCE <- as.SingleCellExperiment(Integrated_Macs_WT_MHCIIHI)
table(Integrated_Macs_WT_MHCIIHI_SCE$sample)

#Re-order sample levels
Integrated_Macs_WT_MHCIILO_SCE$sample <- factor(Integrated_Macs_WT_MHCIILO_SCE$sample, levels = c("WT_2_hr", "WT_6_hr", "WT_Ctrl"))

Integrated_Macs_WT_MHCIIHI_SCE$sample <- factor(Integrated_Macs_WT_MHCIIHI_SCE$sample, levels = c("WT_2_hr", "WT_6_hr", "WT_Ctrl"))

#Assign Cluster IDs in SCE object
Idents(Integrated_Macs_WT_MHCIILO) <- "integrated_snn_res.0.4"
colData(Integrated_Macs_WT_MHCIILO_SCE)$Seurat_clusters <- as.character(Integrated_Macs_WT_MHCIILO@active.ident)

Idents(Integrated_Macs_WT_MHCIIHI) <- "integrated_snn_res.0.4"
colData(Integrated_Macs_WT_MHCIIHI_SCE)$Seurat_clusters <- as.character(Integrated_Macs_WT_MHCIIHI@active.ident)

###Run Slingshot Analysis on Different Populations of WT Macs
#UMAP
Integrated_Macs_WT_MHCIILO_SCE <- slingshot(Integrated_Macs_WT_MHCIILO_SCE, clusterLabels = "Seurat_clusters", reducedDim = "UMAP")

Integrated_Macs_WT_MHCIIHI_SCE <- slingshot(Integrated_Macs_WT_MHCIIHI_SCE, clusterLabels = "Seurat_clusters", reducedDim = "UMAP")

#PCA
Integrated_Macs_WT_MHCIILO_SCE <- slingshot(Integrated_Macs_WT_MHCIILO_SCE, clusterLabels = "Seurat_clusters", reducedDim = "PCA")

Integrated_Macs_WT_MHCIIHI_SCE <- slingshot(Integrated_Macs_WT_MHCIIHI_SCE, clusterLabels = "Seurat_clusters", reducedDim = "PCA")

###Pseudotime Colors
#MHCII LO
colors <- rainbow(100, alpha = 1)
plot(reducedDims(Integrated_Macs_WT_MHCIILO_SCE)$UMAP, col = colors[cut(Integrated_Macs_WT_MHCIILO_SCE$slingPseudotime_1, breaks=100)], pch=16, asp = 1, cex = 0.5, main = "MHCIILO Macrophages")
lines(SlingshotDataSet(Integrated_Macs_WT_MHCIILO_SCE), lwd=4, col ="black")

#MHCII HI
colors <- rainbow(100, alpha = 1)
plot(reducedDims(Integrated_Macs_WT_MHCIIHI_SCE)$UMAP, col = colors[cut(Integrated_Macs_WT_MHCIIHI_SCE$slingPseudotime_1, breaks=100)], pch=16, asp = 1, cex = 0.5, main = "MHCIIHI Macrophages")
lines(SlingshotDataSet(Integrated_Macs_WT_MHCIIHI_SCE), lwd=4, col ="black")

###Pseudotime by Sample
#MHCII LO
ggplot(as.data.frame(colData(Integrated_Macs_WT_MHCIILO_SCE)), aes(x = Integrated_Macs_WT_MHCIILO_SCE$slingPseudotime_1, y = sample, colour = Seurat_clusters)) +  geom_quasirandom(groupOnX = FALSE) +  scale_color_tableau() + theme_classic() +  xlab("Slingshot pseudotime") + ylab("Timepoint") + ggtitle("MHCIILO Macs by Condition")

#MHCII HI
ggplot(as.data.frame(colData(Integrated_Macs_WT_MHCIIHI_SCE)), aes(x = Integrated_Macs_WT_MHCIIHI_SCE$slingPseudotime_1, y = sample, colour = sample)) +  geom_quasirandom(groupOnX = FALSE) +  scale_color_tableau() + theme_classic() +  xlab("Slingshot pseudotime") + ylab("Timepoint") + ggtitle("MHCIIHI Macs by Condition") + xlim(15,22)

###DE Genes
#MHCII LO
DefaultAssay(Integrated_Macs_WT_MHCIILO) <- "RNA"
Integrated_Macs_WT_MHCIILO <- NormalizeData(Integrated_Macs_WT_MHCIILO)
Idents(Integrated_Macs_WT_MHCIILO) <- "sample"
WT_MHCIILO_Markers <- FindMarkers(Integrated_Macs_WT_MHCIILO, ident.1 =  "WT_6_hr", ident.2 = "WT_Ctrl", test.use = "MAST")
WT_MHCIILO_Markers <- subset(WT_MHCIILO_Markers, WT_MHCIILO_Markers$p_val < 0.05 & abs(WT_MHCIILO_Markers$avg_logFC) > 0.4)

write.csv(WT_MHCIILO_Markers, file = "./results/DEgenes_Integrated_Macs_WT_MHCIILO.csv")

#MHCII HI
DefaultAssay(Integrated_Macs_WT_MHCIIHI) <- "RNA"
Integrated_Macs_WT_MHCIIHI <- NormalizeData(Integrated_Macs_WT_MHCIIHI)
Idents(Integrated_Macs_WT_MHCIIHI) <- "sample"
WT_MHCIIHI_Markers <- FindMarkers(Integrated_Macs_WT_MHCIIHI, ident.1 =  "WT_6_hr", ident.2 = "WT_Ctrl", test.use = "MAST")
WT_MHCIIHI_Markers <- subset(WT_MHCIIHI_Markers, WT_MHCIIHI_Markers$p_val < 0.05 & abs(WT_MHCIIHI_Markers$avg_logFC) > 0.4)

write.csv(WT_MHCIIHI_Markers, file = "./results/DEgenes_Integrated_Macs_WT_MHCIIHI.csv")


#Plot Expression of DEgenes as a function of pseudotime
plotExpression(Integrated_Macs_WT_MHCIILO_SCE, features = c("Cd36", "Il1b", "Il6", "Tgfbi", "H2-Ab1", "Il1rn"), x = "slingPseudotime_1", show_smooth = T, point_alpha = 0.4, colour_by = "sample", scales = "free", show_se = F) 

plotExpression(Integrated_Macs_WT_MHCIIHI_SCE, features = c("Cd36", "Il1b", "Il6", "Tgfbi", "H2-Ab1", "Il1rn"), x = "slingPseudotime_1", show_smooth = T, colour_by = "sample", scales = "free", show_se = F) 

#Cluster averages for heatmap using Phantasus
Idents(Integrated_Macs_WT_MHCIILO) <- "sample"
cluster.averages <- AverageExpression(Integrated_Macs_WT_MHCIILO)
write.csv(cluster.averages[["SCT"]], file = "./results/Integrated_Macs_WT_MHCIILO.csv", row.names = TRUE, quote = FALSE)

Idents(Integrated_Macs_WT_MHCIIHI) <- "sample"
cluster.averages <- AverageExpression(Integrated_Macs_WT_MHCIIHI)
write.csv(cluster.averages[["SCT"]], file = "./results/Integrated_Macs_WT_MHCIIHI.csv", row.names = TRUE, quote = FALSE)

#GO Resident v Recruited #########################################
library(GOplot)

WT_MHCII_LO_GO_pathways <- read.csv(file = "./results/gProfiler Results/WT_MHCIILO.csv")
WT_MHCII_LO_GO_pathways <- as.data.frame(WT_MHCII_LO_GO_pathways)

WT_MHCIILO_Markers_genes <- as.data.frame(WT_MHCIILO_Markers)
WT_MHCIILO_Markers_genes$ID <- row.names(WT_MHCIILO_Markers_genes)
names(WT_MHCIILO_Markers_genes)[1] <- "P.Value"
names(WT_MHCIILO_Markers_genes)[2] <- "logFC"
names(WT_MHCIILO_Markers_genes)[5] <- "adj.P.Val"
head(WT_MHCIILO_Markers_genes)

WT_MHCII_LO_circ <- circle_dat(WT_MHCII_LO_GO_pathways, WT_MHCIILO_Markers_genes)

GOBubble(WT_MHCII_LO_circ, ID = F, table.legend = F, colour = c("gray25", "gray25"), labels = 20)

WT_MHCII_HI_GO_pathways <- read.csv(file = "./results/gProfiler Results/WT_MHCIIHI.csv")
WT_MHCII_HI_GO_pathways <- as.data.frame(WT_MHCII_HI_GO_pathways)

WT_MHCIIHI_Markers_genes <- as.data.frame(WT_MHCIIHI_Markers)
WT_MHCIIHI_Markers_genes$ID <- row.names(WT_MHCIIHI_Markers_genes)
names(WT_MHCIIHI_Markers_genes)[1] <- "P.Value"
names(WT_MHCIIHI_Markers_genes)[2] <- "logFC"
names(WT_MHCIIHI_Markers_genes)[5] <- "adj.P.Val"
head(WT_MHCIIHI_Markers_genes)

WT_MHCII_HI_circ <- circle_dat(WT_MHCII_HI_GO_pathways, WT_MHCIIHI_Markers_genes)

GOBubble(WT_MHCII_HI_circ, ID = F, table.legend = F, colour = c("goldenrod", "goldenrod"), labels = 10)


#MerTK KO Steady State ########################################################
#Isolate MerKD Resident Peritoneal Macrophages
Idents(Integrated_Seurat_MacsDCs) <- "integrated_snn_res.0.4"
Integrated_Seurat_Macs <- SubsetData(Integrated_Seurat_MacsDCs, ident.use = c("1", "3", "5", "6", "11"))
Idents(Integrated_Seurat_Macs) <- "sample"
Integrated_Macrophages_MerKD <- SubsetData(Integrated_Seurat_Macs, ident.use = c("MerKD_Ctrl", "MerKD_2hr", "MerKD_6_hr"))
Integrated_Macrophages_Control <- SubsetData(Integrated_Seurat_Macs, ident.use = c("WT_Ctrl", "MerKD_Ctrl"))

#Determine heterogeneity between control conditions of the 2 geneotypes
Idents(Integrated_Macrophages_Control) <- "integrated_snn_res.0.4"
DimPlot(Integrated_Macrophages_Control, group.by = "sample", cols = c("coral", "firebrick3", "grey"))
DimPlot(Integrated_Macrophages_Control, split.by = "sample", cols = c("royalblue3", "firebrick", "deepskyblue2", "goldenrod", "darkorchid"))

DefaultAssay(Integrated_Macrophages_Control) <- "RNA"
Integrated_Macrophages_Control  <- NormalizeData(Integrated_Macrophages_Control)
Idents(Integrated_Macrophages_Control) <- "sample"
Control_Markers <- FindMarkers(Integrated_Macrophages_Control, ident.1 = "WT_Ctrl", ident.2 = "MerKD_Ctrl", test.use = "MAST")
Control_Markers_pos <- subset(Control_Markers, Control_Markers$p_val_adj < 0.05 & Control_Markers$avg_logFC > 0)
Control_Markers_neg <- subset(Control_Markers, Control_Markers$p_val_adj < 0.05 & Control_Markers$avg_logFC < 0)

write.csv(Control_Markers_pos , file = "./results/Control_Markers_pos.csv")
write.csv(Control_Markers_neg , file = "./results/Control_Markers_neg.csv")

#Rename Resident and Recruited Mac Populations
Idents(Integrated_Macrophages_Control) <- "integrated_snn_res.0.4"
Integrated_Macrophages_Control_ResvRec <- SubsetData(Integrated_Macrophages_Control, ident.use = c("1", "3", "5", "6"))
Integrated_Macrophages_Control_ResvRec <- RenameIdents(Integrated_Macrophages_Control_ResvRec, "1" =  "Resident", "3" = "Resident", "5" = "Resident", "6" = "Recruited")

Integrated_Macrophages_Control_ResvRec$sample <- factor(Integrated_Macrophages_Control_ResvRec$sample, levels = c("WT_Ctrl", "MerKD_Ctrl"))
VlnPlot.Macs <- VlnPlot(Integrated_Macrophages_Control_ResvRec, features = c("Il1b", "Tgfbi", "Axl", "H2-Ab1", "Cd74", "Cd36"), pt.size = 0, cols = c("firebrick", "grey"), combine = F, split.by = "sample", sort = "decreasing")
for(i in 1:length(VlnPlot.Macs)) {
  VlnPlot.Macs[[i]] <- VlnPlot.Macs[[i]] + NoLegend()+ theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
}
cowplot::plot_grid(plotlist = VlnPlot.Macs, ncol = 3)

#Compare to WT Control Macrophages
FeaturePlot(Integrated_Macrophages_Control, features = c("Itgam", "Adgre1", "Icam2", "Itgax", "H2-Ab1", "Ccr2", "Timd4", "Klf2", "Mertk"), order = T, pt.size = 1.5, ncol = 3)

Integrated_Macrophages_Control$sample <- factor(Integrated_Macrophages_Control$sample, levels = c("WT_Ctrl", "MerKD_Ctrl"))
VlnPlot.Macs <- VlnPlot(Integrated_Macrophages_Control, features = c("Klf2", "Gata6", "Cd93", "Csf1r", "Anxa5", "Cd36", "Actb", "Arg1", "Irf4", "Retnla", "Cd74", "Tgfbi",  "Il1b", "Il6", "Tnf", "Cxcl2"),  cols = c("firebrick", "grey"), pt.size = 0, combine = F, split.by = "sample")
for(i in 1:length(VlnPlot.Macs)) {
  VlnPlot.Macs[[i]] <- VlnPlot.Macs[[i]] + NoLegend() + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
}
cowplot::plot_grid(plotlist = VlnPlot.Macs, ncol = 4)


#Isolate "Resident" F4/80hi Macs and Recruited MHCIIhi Macs
Idents(Integrated_Macrophages_MerKD) <- "sample"
Integrated_Macrophages_MerKD$sample <- factor(Integrated_Macrophages_MerKD$sample, levels = c("MerKD_Ctrl", "MerKD_2hr", "MerKD_6_hr"))
DimPlot(Integrated_Macrophages_MerKD, group.by = "sample", cols = c("grey", "grey", "firebrick3"), pt.size = 3, order = T)
Integrated_Macrophages_MerKD$sample <- factor(Integrated_Macrophages_MerKD$sample, levels = c("MerKD_Ctrl", "MerKD_6_hr", "MerKD_2hr"))
DimPlot(Integrated_Macrophages_MerKD, group.by = "sample", cols = c("grey", "grey", "firebrick3"), pt.size = 3, order = T)
Integrated_Macrophages_MerKD$sample <- factor(Integrated_Macrophages_MerKD$sample, levels = c("MerKD_2hr", "MerKD_6_hr", "MerKD_Ctrl"))
DimPlot(Integrated_Macrophages_MerKD, group.by = "sample", cols = c("grey", "grey", "firebrick3"), pt.size = 3, order = T)

#Split Resident versus Recruited
Idents(Integrated_Macrophages_MerKD) <- "integrated_snn_res.0.4"
Integrated_Macrophages_MerKD_Resident <- SubsetData(Integrated_Macrophages_MerKD, ident.use = c("1", "3", "5"))
Integrated_Macrophages_MerKD_Recruited <- SubsetData(Integrated_Macrophages_MerKD, ident.use = "6")

#Create single cell experiment objects
Integrated_Macrophages_MerKD_Resident_SCE <- as.SingleCellExperiment(Integrated_Macrophages_MerKD_Resident)
table(Integrated_Macrophages_MerKD_Resident_SCE$sample)

Integrated_Macrophages_MerKD_Recruited_SCE <- as.SingleCellExperiment(Integrated_Macrophages_MerKD_Recruited)
table(Integrated_Macrophages_MerKD_Recruited_SCE$sample)


#Re-order sample levels
Integrated_Macrophages_MerKD_Resident_SCE$sample <- factor(Integrated_Macrophages_MerKD_Resident_SCE$sample, levels = c("MerKD_6_hr", "MerKD_2hr", "MerKD_Ctrl"))

Integrated_Macrophages_MerKD_Recruited_SCE$sample <- factor(Integrated_Macrophages_MerKD_Recruited_SCE$sample, levels = c("MerKD_6_hr", "MerKD_2hr", "MerKD_Ctrl"))


#Assign Cluster IDs in SCE object
Idents(Integrated_Macrophages_MerKD_Resident) <- "integrated_snn_res.0.4"
colData(Integrated_Macrophages_MerKD_Resident_SCE)$Seurat_clusters <- as.character(Integrated_Macrophages_MerKD_Resident@active.ident)

Idents(Integrated_Macrophages_MerKD_Recruited) <- "integrated_snn_res.0.6"
colData(Integrated_Macrophages_MerKD_Recruited_SCE)$Seurat_clusters <- as.character(Integrated_Macrophages_MerKD_Recruited@active.ident)


###Run Slingshot Analysis on Different Populations of MerKD Macs
#UMAP
Integrated_Macrophages_MerKD_Resident_SCE <- slingshot(Integrated_Macrophages_MerKD_Resident_SCE, clusterLabels = "Seurat_clusters", reducedDim = "UMAP")

Integrated_Macrophages_MerKD_Recruited_SCE <- slingshot(Integrated_Macrophages_MerKD_Recruited_SCE, clusterLabels = "Seurat_clusters", reducedDim = "UMAP")

###Pseudotime Colors
#Resident
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(50)
plot(reducedDims(Integrated_Macrophages_MerKD_Resident_SCE)$UMAP, col = colors[cut(Integrated_Macrophages_MerKD_Resident_SCE$slingPseudotime_1, breaks=100)], pch=16, asp = 1, cex = 0.5, main = "Resident Macrophages")
lines(SlingshotDataSet(Integrated_Macrophages_MerKD_Resident_SCE), lwd=4, col ="black", type = "lineages")
legend_image <- as.raster(matrix(colors[cut(Integrated_Macrophages_MerKD_Resident_SCE$slingPseudotime_1,breaks=100)], ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
text(x=1.5, y = seq(1,0,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 0, 0, 1,1)


#Recruited
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(reducedDims(Integrated_Macrophages_MerKD_Recruited_SCE)$UMAP, col = colors[cut(Integrated_Macrophages_MerKD_Recruited_SCE$slingPseudotime_1,breaks=50)], pch=16, asp = 1, cex = 0.5, main = "Recruited Macrophages")
lines(SlingshotDataSet(Integrated_Macrophages_MerKD_Recruited_SCE), lwd=4, col ="black")
legend_image <- as.raster(matrix(colors[cut(Integrated_Macrophages_MerKD_Recruited_SCE$slingPseudotime_1,breaks=50)], ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 0, 0, 1,1)

###PCA
#Resident
Integrated_Macrophages_MerKD_Resident_SCE <- slingshot(Integrated_Macrophages_MerKD_Resident_SCE, clusterLabels = "Seurat_clusters", reducedDim = "PCA", start.clus = "1")
colors <- rainbow(50, alpha = 1)
plot(reducedDims(Integrated_Macrophages_MerKD_Resident_SCE)$PCA, col = colors[cut(Integrated_Macrophages_MerKD_Resident_SCE$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(Integrated_Macrophages_MerKD_Resident_SCE), lwd=4, col ="black", type ="lineages")


#Recruited
Integrated_Macrophages_MerKD_Recruited_SCE <- slingshot(Integrated_Macrophages_MerKD_Recruited_SCE, clusterLabels = "Seurat_clusters", reducedDim = "PCA")
colors <- rainbow(50, alpha = 1)
plot(reducedDims(Integrated_Macrophages_MerKD_Recruited_SCE)$PCA, col = colors[cut(Integrated_Macrophages_MerKD_Recruited_SCE$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(Integrated_Macrophages_MerKD_Recruited_SCE), lwd=4, col ="black")

###Pseudotime by Sample
#Resident
ggplot(as.data.frame(colData(Integrated_Macrophages_MerKD_Resident_SCE)), aes(x = Integrated_Macrophages_MerKD_Resident_SCE$slingPseudotime_1, y = sample, colour = Seurat_clusters)) +  geom_quasirandom(groupOnX = FALSE) +  scale_color_tableau() + theme_classic() +  xlab("Slingshot pseudotime") + ylab("Timepoint") + ggtitle("Resident Macs by Condition") 

#Recruited
ggplot(as.data.frame(colData(Integrated_Macrophages_MerKD_Recruited_SCE)), aes(x = Integrated_Macrophages_MerKD_Recruited_SCE$slingPseudotime_1, y = sample, colour = integrated_snn_res.0.6)) +  geom_quasirandom(groupOnX = F)  + scale_color_tableau()+ theme_classic() +  xlab("Slingshot pseudotime") + ylab("Timepoint") + ggtitle("Recruited Macs by Condition")

#Gene Ontology #################################
DefaultAssay(Integrated_Macrophages_MerKD_Resident) <- "RNA"
Integrated_Macrophages_MerKD_Resident <- NormalizeData(Integrated_Macrophages_MerKD_Resident)
Idents(Integrated_Macrophages_MerKD_Resident) <- "sample"
MerKD_Resident_Markers <- FindMarkers(Integrated_Macrophages_MerKD_Resident, ident.2 =  "MerKD_Ctrl", ident.1 =  "MerKD_6_hr", test.use = "MAST", only.pos = T)
MerKD_Resident_Markers <- subset(MerKD_Resident_Markers, MerKD_Resident_Markers$p_val_adj < 0.05)
write.csv(MerKD_Resident_Markers, file = "./results/DEgenes_MerKD_Resident_Effero.csv")

DefaultAssay(Integrated_Macrophages_MerKD_Recruited) <- "RNA"
Integrated_Macrophages_MerKD_Recruited <- NormalizeData(Integrated_Macrophages_MerKD_Recruited)
Idents(Integrated_Macrophages_MerKD_Recruited) <- "sample"
MerKD_Recruited_Markers <- FindMarkers (Integrated_Macrophages_MerKD_Recruited, ident.2 =  "MerKD_Ctrl", ident.1 =  "MerKD_6_hr", test.use = "MAST", only.pos = T)
MerKD_Recruited_Markers <- subset(MerKD_Recruited_Markers, MerKD_Recruited_Markers$p_val_adj < 0.05)

write.csv(MerKD_Recruited_Markers, file = "./results/DEgenes_MerKD_Recruited_Effero.csv")

#MHCII LO
DefaultAssay(Integrated_Macs_WT_MHCIILO) <- "RNA"
Integrated_Macs_WT_MHCIILO <- NormalizeData(Integrated_Macs_WT_MHCIILO)
Idents(Integrated_Macs_WT_MHCIILO) <- "sample"
WT_MHCIILO_Markers <- FindMarkers(Integrated_Macs_WT_MHCIILO, ident.2 =  "WT_Ctrl", ident.1 = "WT_6_hr", test.use = "MAST", only.pos = T)
WT_MHCIILO_Markers <- subset(WT_MHCIILO_Markers, WT_MHCIILO_Markers$p_val < 0.05)


#MHCII HI
DefaultAssay(Integrated_Macs_WT_MHCIIHI) <- "RNA"
Integrated_Macs_WT_MHCIIHI <- NormalizeData(Integrated_Macs_WT_MHCIIHI)
Idents(Integrated_Macs_WT_MHCIIHI) <- "sample"
WT_MHCIIHI_Markers <- FindMarkers(Integrated_Macs_WT_MHCIIHI, ident.2 =  "WT_Ctrl", ident.1 = "WT_6_hr", test.use = "MAST")
WT_MHCIIHI_Markers <- subset(WT_MHCIIHI_Markers, WT_MHCIIHI_Markers$p_val < 0.05)

write.csv(WT_MHCIIHI_Markers, file = "./results/WT_MHCIIHI_Markers.csv")

#Resident Mac Comparison between genotypes
Resident_Same_DEgenes <- intersect(row.names(MerKD_Resident_Markers), row.names(WT_MHCIILO_Markers))
Resident_Control_unique <- setdiff(row.names(WT_MHCIILO_Markers), row.names(MerKD_Resident_Markers))
Resident_Knockout_unique <- setdiff(row.names(MerKD_Resident_Markers), row.names(WT_MHCIILO_Markers))

write.csv(Resident_Knockout_unique, file = "./results/Resident_Knockout_unique.csv")
write.csv(Resident_Control_unique, file = "./results/Resident_Control_unique.csv")
write.csv(Resident_Same_DEgenes, file = "./results/Resident_Same_DEgenes.csv")

#Recruited Mac Comparison between genotypes
Recruited_Same_DEgenes <- intersect(row.names(MerKD_Recruited_Markers), row.names(WT_MHCIIHI_Markers))
Recruited_Control_unique <- setdiff(row.names(WT_MHCIIHI_Markers), row.names(MerKD_Recruited_Markers))
Recruited_Knockout_unique <- setdiff(row.names(MerKD_Recruited_Markers), row.names(WT_MHCIIHI_Markers))

write.csv(Recruited_Knockout_unique, file = "./results/Recruited_Knockout_unique.csv")
write.csv(Recruited_Control_unique, file = "./results/Recruited_Control_unique.csv")
write.csv(Recruited_Same_DEgenes, file = "./results/Recruited_Same_DEgenes.csv")

#GO Plots########################################################
library(GOplot)

MHCIILO_Common_GO <- read.csv(file = "./results/gProfiler Results/MHCIILO_Common.csv")
MHCIILO_Common_GO <- as.data.frame(MHCIILO_Common_GO)
MerKD_MHCIILO_Unique_GO <- read.csv(file = "./results/gProfiler Results/MerKD_MHCIILO_Unique.csv")
MerKD_MHCIILO_Unique_GO <- as.data.frame(MerKD_MHCIILO_Unique_GO)
WT_MHCIILO_Unique_GO <- read.csv(file = "./results/gProfiler Results/WT_MHCIILO_Unique.csv")
WT_MHCIILO_Unique_GO <- as.data.frame(WT_MHCIILO_Unique_GO)

Resident_Effero_6hr_genes <- read.csv(file = "./results/DEgenes_MerKD_Resident_Effero.csv")
Resident_Effero_6hr_genes <- as.data.frame(Resident_Effero_6hr_genes)
head(Resident_Effero_6hr_genes)

WT_MHCIILO_Markers_genes <- as.data.frame(WT_MHCIILO_Markers)
WT_MHCIILO_Markers_genes$ID <- row.names(WT_MHCIILO_Markers_genes)
names(WT_MHCIILO_Markers_genes)[1] <- "P.Value"
names(WT_MHCIILO_Markers_genes)[2] <- "logFC"
names(WT_MHCIILO_Markers_genes)[5] <- "adj.P.Val"
head(WT_MHCIILO_Markers_genes)

#Common
MHCIILO_Common_GO_circle <- circle_dat(MHCIILO_Common_GO, WT_MHCIILO_Markers_genes)
GOBubble(MHCIILO_Common_GO_circle , ID = F, table.legend = F, colour = c("gray40", "gray40"), labels = 0)

#MerKD Unique
MerKD_MHCIILO_Unique_GO_circle <- circle_dat(MerKD_MHCIILO_Unique_GO, Resident_Effero_6hr_genes)
GOBubble(MerKD_MHCIILO_Unique_GO_circle , ID = F, table.legend = F, colour = c("firebrick", "firebrick"), labels = 10)

#WT Unique
WT_MHCIILO_Unique_GO_circle <- circle_dat(WT_MHCIILO_Unique_GO, WT_MHCIILO_Markers_genes)
GOBubble(WT_MHCIILO_Unique_GO_circle , ID = F, table.legend = F, colour = c("gray25", "gray25"), labels = 0)


Integrated_Macs_WT
DefaultAssay(Integrated_Macs_WT) <- "RNA"
Integrated_Macs_WT <- NormalizeData(Integrated_Macs_WT)
Idents(Integrated_Macs_WT) <- "sample"
WT_Markers <- FindMarkers(Integrated_Macs_WT, ident.1 =  "WT_6_hr", ident.2 = "WT_Ctrl", test.use = "MAST")
WT_Markers_r <- subset(WT_Markers, WT_Markers$p_val < 0.05 & abs(WT_Markers$avg_logFC) > 0.5)

write.csv(WT_Markers_r, file = "./results/EfferocytosisGeneSet.csv")

Integrated_Macs_WT_MHCIILO
DefaultAssay(Integrated_Macs_WT_MHCIILO) <- "RNA"
Integrated_Macs_WT_MHCIILO <- NormalizeData(Integrated_Macs_WT_MHCIILO)
Idents(Integrated_Macs_WT_MHCIILO) <- "sample"
WT_Markers <- FindMarkers(Integrated_Macs_WT_MHCIILO, ident.1 =  "WT_6_hr", ident.2 = "WT_Ctrl", test.use = "MAST")
WT_Markers_r <- subset(WT_Markers, WT_Markers$p_val < 0.05 & abs(WT_Markers$avg_logFC) > 0.5)

write.csv(WT_Markers_r, file = "./results/F480Macs_EfferocytosisGeneSet.csv")

