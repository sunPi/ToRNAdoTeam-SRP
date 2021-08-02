# Re-analysis of single cell data using UMI cell counts --------------------------------------------------------------------------------------------------------------------------------------------------------------
#  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ------------- Check versions so they match the authors -------------------

CellRanger <- 3.0
Seurat3 <- 4.0
monocle3<- 4.0

Required_Versions <- data.frame(CellRanger,Seurat3,monocle3)
format.data.frame(Required_Versions, nsmall = 1)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocCheck")


if (packageVersion("Seurat")>="4.0.0") {
    print("Seurat version is ok!")
}

# ------------- Dependencies -------------------

#install.packages("magrittr")
#install.packages("devtools")
#devtools::install_github('cole-trapnell-lab/leidenbase')
#devtools::install_github('cole-trapnell-lab/monocle3')
#install.packages('Seurat')

library("magrittr")
library("dplyr")
library("monocle3")
library("Seurat")
library(data.table)


# ------------- Import Cell Counts -------------------

c_counts <- read.delim('IPDCO_hg_midbrain_UMI.tsv', sep = "\t")
#head(c_counts)
dimz <- dim(c_counts)
paste0("Successefully read the data. Its length is ", dimz, ".",sep="")

# ------------- Create the Seurat object from cell counts (UMI matrix) -------------------

library(data.table)
library(Seurat)
ipdco <- CreateSeuratObject(counts = c_counts, project = "ipdco_janrogel_2021", min.cells = 3, min.features = 200)
ipdco

# ------------- Normalizing the data -------------------

ipdco <- NormalizeData(ipdco, normalization.method = "LogNormalize", scale.factor = 10000)

# ------------- Identification of highly variable features -------------------

ipdco <- FindVariableFeatures(ipdco, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes
top20 <- head(VariableFeatures(ipdco), 20)

# plot variable features with and without lables
plot1 <- VariableFeaturePlot(ipdco)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

# ------------- Scaling the data -------------------

all.genes <- rownames(ipdco)
ipdco     <- ScaleData(ipdco, features = all.genes)

# ------------- Linear dimensional reduction -------------------

ipdco <- RunPCA(ipdco, features = VariableFeatures(object = ipdco))

#Vizualize
VizDimLoadings(ipdco, dims = 1:20, reduction = "pca")
DimPlot(ipdco, reduction = "pca")
DimHeatmap(ipdco, dims = 1:20, cells = 500, balanced = TRUE)

# ------------- Linear dimensional reduction ------------------
### NOTE: This process can take a long time for big datasets, comment out for expediency. More
### approximate techniques such as those implemented in ElbowPlot() can be used to reduce
### computation time
ipdco <- JackStraw(ipdco, num.replicate = 100)
ipdco <- ScoreJackStraw(ipdco, dims = 1:20)

JackStrawPlot(ipdco, dims = 1:20)

ElbowPlot(ipdco)

# ------------- Cluster the cells ------------------

print("This step is performed using the FindNeighbors function, and takes as input the previously defined dimensionality of the dataset (first 10 PCs)")
ipdco <- FindNeighbors(ipdco, dims = 1:20)

print("Finding clusters using the FindClusters function. Optimal resolution often increases for larger datasets. The clusters can be found using the Idents function.")
ipdco <- FindClusters(ipdco, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(ipdco), 5)

# ------------- Run non-linear dimensional reduction (UMAP/tSNE) ------------------

print("Cells within the graph-based clusters determined above should co-localize on these dimension reduction plots. As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis.")
ipdco <- RunUMAP(ipdco, dims = 1:20)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

DimPlot(ipdco, reduction = "umap")

saveRDS(ipdco, file = "/scratch/spectre/j/jr429/SRP/r-analysis-manual(UMI)")
print("Saved .rds file successfuly!")

# ------------- Find differentially expressed features (cluster biomarkers) ------------------

# find all markers of cluster 2
cluster2.markers <- FindMarkers(ipdco, ident.1 = 2, min.pct = 0.25)
print("This is cluster2.markers:")
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
print("This is cluster5.markers:")
cluster5.markers <- FindMarkers(ipdco, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
ipdco.markers <- FindAllMarkers(ipdco, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
print("This is the ipdco.markers table:")
ipdco.markers %>%
    group_by(cluster) %>%
    top_n(n = 2, wt = avg_log2FC)
fwrite(ipdco.markers,file="/scratch/spectre/j/jr429/SRP/jobs/output/ipdco.gene.markers.tsv",sep = "\t")

# ROC test for differential expression
cluster0.markers <- FindMarkers(ipdco, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# View the whole selection of gene markers
print("This is the cluster0.markers table:")
cluster0.markers

# Write out the selection in a .tsv file.

fwrite(cluster0.markers,file="/scratch/spectre/j/jr429/SRP/jobs/output/ipdco..markers.tsv",sep = "\t")
# Vizualizing marker expression; shows expression probability distributions across clusters
VlnPlot(ipdco, features = c("10616", "11322"))

# plotting raw counts as well
VlnPlot(ipdco, features = c("17280", "15549"), slot = "counts", log = TRUE)


FeaturePlot(ipdco, features = c("10616", "11322", "17280", "15549", "16693", "15549", "5456", "10842",
    "1672", "6360"))
# DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.

ipdco.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top20

DoHeatmap(ipdco, features = top20$gene) + NoLegend()

# ------------- Assigning cell type identity to clusters ------------------

new.cluster.ids <- c("10616", "11322", "17280", "15549", "16693", "15549", "5456", "10842",
    "1672", "6360")
names(new.cluster.ids) <- levels(ipdco)
ipdco <- RenameIdents(ipdco, new.cluster.ids)
DimPlot(ipdco, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

print("Script has been run successfully!")
print("An output.txt has been generated.")
saveRDS(ipdco, file = "/scratch/spectre/j/jr429/SRP/jobs/Final")
