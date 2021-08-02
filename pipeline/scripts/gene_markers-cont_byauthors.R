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

# ------------- Read RDS --------------

ipdco <- readRDS("./r-analysis-manual(UMI).rds")

# ------------- Find differentially expressed features (cluster biomarkers) ------------------
# ------------- Identifying the marker genes ------------------------------------------------
# find all markers of cluster 2
cluster2.markers <- FindMarkers(ipdco, ident.1 = 2, min.pct = 0.25)
print("This is cluster2.markers:")
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
print("This is cluster5.markers:")
cluster5.markers <- FindMarkers(ipdco, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# ------------------------- ROC test for differential expression -----------------------------
cluster0.markers <- FindAllMarkers(ipdco,assay = `assay`,
				  only.pos = TRUE,
				  min.pct = 0.25,
				  test.use = "roc",
				  logfc.threshold = 0.25)
fwrite(cluster0.markers,file="/scratch/spectre/j/jr429/SRP/jobs/output/ipdco.gene.markers.ROC.tsv",sep = "\t")
# ------------------------- Wilcox method -----------------------------
ipdco.markers <- FindAllMarkers(ipdco, only.pos = TRUE, 
				     logfc.threshold = 0.25,
				     assay = `assay`,
				     min.pct = 0.25,
				     test.use = "wilcox")
print(head(ipdco.markers 2))

# Do some merging
ipdco.markers <- merge(cluster0.markers, ipdco.markers, by = "gene", all = TRUE)

print("This is the ipdco.markers table FOR ALL CELLS:")
ipdco.markers %>%
    group_by(cluster) %>%
    top_n(n = 2, wt = avg_log2FC)
fwrite(ipdco.markers,file="/scratch/spectre/j/jr429/SRP/jobs/output/ipdco.gene.markers.Wilcox.tsv",sep = "\t")

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

# DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting the top 10 markers (or all markers if less than 20) for each cluster.
ipdco.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(ipdco, features = top10$gene) + NoLegend()
# ------------------------------------------------- AUC Heatmap

pdf(paste0("output/janr_tornado_auc_marker_genes_heatmap.pdf"), 
     width = 18, height = 15)
plot(
DoHeatmap(
  ipdco,
  features = top10$gene,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = "SCT",
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE) + scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b")))

dev.off()

jpeg(paste0("output/janr_tornado_wilcoxon_marker_genes_heatmap.jpeg"), 
    width = 1980, height = 980, pointsize = 74, quality = 100)
plot(
DoHeatmap(
  ipdco,
  features = top10$gene,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = "SCT",
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE) + scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b")))

dev.off()

# ------------------------------------------------- Wilcoxon Heatmap

pdf(paste0("output/janr_tornado_auc_marker_genes_heatmap.pdf"), 
     width = 18, height = 15)
plot(
DoHeatmap(
  ipdco.markers,
  features = top10$gene,
  cells = cells,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = "SCT",
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
) + scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b"))
)
dev.off()

jpeg(paste0("output/janr_tornado_wilcoxon_marker_genes_heatmap.jpeg"), 
 	    width = 1980, height = 980, pointsize = 74, quality = 100)
plot(
DoHeatmap(
  ipdco.markers,
  features = top10$gene,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = "SCT",
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE) + scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b")))
dev.off()
# ------------- Assigning cell type identity to clusters ------------------

new.cluster.ids <- c("10616", "11322", "17280", "15549", "16693", "15549", "5456", "10842",
    "1672", "6360")
names(new.cluster.ids) <- levels(ipdco)
ipdco <- RenameIdents(ipdco, new.cluster.ids)
DimPlot(ipdco, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

print("Script has been run successfully!")
print("An output.txt has been generated.")
saveRDS(ipdco, file = "/scratch/spectre/j/jr429/SRP/jobs/Final")

pdf(paste0("output/janr_tornado_cell_type_identities.pdf"), 
     width = 18, height = 15)
plot(
DimPlot(ipdco, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
)
dev.off()

jpeg(paste0("output/janr_tornado_wilcoxon_marker_genes_heatmap.jpeg"), 
 	    width = 1980, height = 980, pointsize = 74, quality = 100)
plot(
DimPlot(ipdco, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
)
dev.off()
