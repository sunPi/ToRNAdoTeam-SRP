# UMI cell counts ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
# ------------- Check versions so they match the authors -------------------

A_Versions <- c('CellRangerPipeline' = 3.0, 'Seurat3' = 4.0, 'monocle3' = 4.0)
CellRanger <- 3.0
Seurat3 <- 4.0
monocle3<- 4.0

A_Versions <- data.frame(CellRanger,Seurat3,monocle3)
View(A_Versions)
format.data.frame(A_Versions, nsmall = 1)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocCheck")


if (packageVersion("Seurat")>="4.0.0") {
    print("Seurat version is ok!")
}

# ------------- Dependencies -------------------

install.packages("magrittr")
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
install.packages('Seurat')

library("magrittr")
library("monocle3")
library("Seurat")

# ------------- Import Cell Counts -------------------
library(data.table)
library(Seurat)
c_counts <- fread('IPDCO_hg_midbrain_UMI.tsv')

# ------------- Create the Seurat object from cell counts (UMI matrix) -------------------

ipdco <- CreateSeuratObject(counts = c_counts, project = "ipdco_manual", min.cells = 3, min.features = 200)
ipdco

print("A few cells inside the UMI matrix:")
ipdco.data[1:3]

# ------------- Normalizing the data -------------------

ipdco <- NormalizeData(ipdco, normalization.method = "LogNormalize", scale.factor = 10000)

# ------------- Identification of highly variable features -------------------

ipdco <- FindVariableFeatures(ipdco, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ipdco), 10)

# plot variable features with and without lables
plot1 <- VariableFeaturePlot(ipdco)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

# ------------- Scaling the data -------------------

all.genes <- rownames(ipdco)
ipdco     <- ScaleData(ipdco, features = all.genes)

# ------------- Linear dimensional reduction -------------------

ipdco <- RunPCA(ipdco, features = VariableFeatures(object = ipdco))

#Vizualize
VizDimLoadings(ipdco, dims = 1:2, reduction = "pca")

#Save as .jpeg
jpeg("IPDCO_PCA_Dim_Reduction.jpg", width = 1200, height = "800")
VizDimLoadings(ipdco, dims = 1:2, reduction = "pca")
dev.off()

jpeg("IPDCO_PCA_Dim_Plot.jpg", width = 1200, height = "800")
DimPlot(ipdco, reduction = "pca")
dev.off()

jpeg("IPDCO_PCA_Dim_Plot.jpg", width = 1200, height = "800")
DimHeatmap(ipdco, dims = 1:10, cells = 500, balanced = TRUE)
dev.off()

# ------------- Linear dimensional reduction ------------------
### NOTE: This process can take a long time for big datasets, comment out for expediency. More
### approximate techniques such as those implemented in ElbowPlot() can be used to reduce
### computation time
#ipdco <- JackStraw(ipdco, num.replicate = 100)
#ipdco <- ScoreJackStraw(ipdco, dims = 1:20)

#JackStrawPlot(ipdco, dims = 1:15)

#ElbowPlot(ipdco)

# ------------- Cluster the cells ------------------

print("This step is performed using the FindNeighbors function, and takes as input the previously defined dimensionality of the dataset (first 10 PCs)")
pbmc <- FindNeighbors(ipdco, dims = 1:10)

print("Finding clusters using the FindClusters function. Optimal resolution often increases for larger datasets. The clusters can be found using the Idents function.")
pbmc <- FindClusters(ipdco, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(ipdco), 5)

# ------------- Run non-linear dimensional reduction (UMAP/tSNE) ------------------

print("Cells within the graph-based clusters determined above should co-localize on these dimension reduction plots. As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis.")
pbmc <- RunUMAP(ipdco, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

jpeg("IPDCO_PCA_Dim_Plot.jpg", width = 1200, height = "800")
DimPlot(ipdco, reduction = "umap")
dev.off()

saveRDS(ipdco, file = "/scratch/spectre/j/jr429/SRP")
print("Saved .rds file sucessefuly!")

# ------------- Find differentially expressed features (cluster biomarkers) ------------------

# find all markers of cluster 2
cluster2.markers <- FindMarkers(ipdco, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(ipdco, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
ipdco.markers <- FindAllMarkers(ipdco, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ipdco.markers %>%
    group_by(cluster) %>%
    top_n(n = 2, wt = avg_log2FC)

# ROC test for differential expression
cluster0.markers <- FindMarkers(ipdco, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)x

# Vizualizing marker expression; shows expression probability distributions across clusters
VlnPlot(ipdco, features = c("gene1", "gene2"))

# plotting raw counts as well
VlnPlot(ipdco, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(ipdco, features = c("1", "2", "3", "4", "5", "6", "7", "8",
    "9", "10", "11", "12", "13", "14", "15"))
# DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.

ipdco.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(ipdco, features = top20$gene) + NoLegend()

# ------------- Assigning cell type identity to clusters ------------------

new.cluster.ids <- c("1", "2", "3", "4", "5", "6", "7", "8",
    "9", "10", "11", "12", "13", "14", "15")
names(new.cluster.ids) <- levels(ipdco)
ipdco <- RenameIdents(ipdco, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# ------------- JUNK -------------------
#downsampled = SampleUMI(data = c_counts, max.umi = 1000, upsample = FALSE, verbose = TRUE)
#head(x = downsampled)
#c_counts_final <- data.frame(c_counts[,-1])
#c_counts_final <- data.frame(c_counts_final[,-1], row.names = c_counts_final[,-1])
#order_c_counts <- c_counts_final[order(row.names(c_counts_final))]
#c_counts <- as.matrix(x = GetAssayData(object = c_counts, assay = "RNA", slot = "counts"))
