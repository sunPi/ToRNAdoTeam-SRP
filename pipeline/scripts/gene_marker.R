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

# ------------- Open the RDS -------------------

ipdco <- readRDS("r-analysis-manual(UMI).rds")

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
# ones (Using the default Wilcoxon rank sum test)
ipdco.markers <- FindAllMarkers(ipdco, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
print("This is the ipdco.markers table:")
ipdco.markers %>%
    group_by(cluster) %>%
    top_n(n = 2, wt = avg_log2FC)
write.table(ipdco.markers,file="/scratch/spectre/j/jr429/SRP/jobs/output/ipdco.gene.markers.tsv",sep = "\t")

# ROC test for differential expression
cluster0.markers <- FindMarkers(ipdco, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# View the whole selection of gene markers
print("This is the cluster0.markers table:")
cluster0.markers

# Write out the selection in a .tsv file.
library(data.table)
#write.table(ipdco,file="/scratch/spectre/j/jr429/SRP/jobs/output/ipdco.gene.markers.tsv",sep = "\t")
write.table(cluster0.markers,file="/scratch/spectre/j/jr429/SRP/jobs/output/ipdco..markers.tsv",sep = "\t")
# Vizualizing marker expression; shows expression probability distributions across clusters
VlnPlot(ipdco, features = c("10616")
VlnPlot(ipdco, features = c("11322")

# plotting raw counts as well
VlnPlot(ipdco, features = c("17280"), slot = "counts", log = TRUE)
VlnPlot(ipdco, features = c("15549"), slot = "counts", log = TRUE)

FeaturePlot(ipdco, features = c("10616", "11322", "17280", "15549", "16693", "15549", "5456", "10842",
    "1672", "6360"))

# DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting the top 10 markers (or all markers if less than 20) for each cluster.

ipdco.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(ipdco, features = top10$gene) + NoLegend()

# ------------- Assigning cell type identity to clusters ------------------

new.cluster.ids <- c("10616", "11322", "17280", "15549", "16693", "15549", "5456", "10842",
    "1672", "6360")
names(new.cluster.ids) <- levels(ipdco)
ipdco <- RenameIdents(ipdco, new.cluster.ids)
DimPlot(ipdco, reduction = "umap", label = TRUE, pt.size = 0.3) + NoLegend()

print("Script has been run successfully!")
print("An output.txt has been generated.")
saveRDS(ipdco, file = "/scratch/spectre/j/jr429/SRP/jobs/Final.rds")
