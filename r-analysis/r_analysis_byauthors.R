# Re-analysis of single cell data using UMI cell counts --------------------------------------------------------------------------------------------------------------------------------------------------------------
#  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

"  Run the pipeline for re-analysis of Smajic et. al. data 
Usage: r_analyis_byauthors.R --jobname=<value> [--UMImatrix=<file>] --outfolder=<folder>
Options:
  -h --help            Show this screen.
  --jobname=<file>     Descriptive name for analysis recreation.
  --UMImatrix=<file>   Mandatory file. This is the UMI cell count matrix that the authors shared publicly (on GEO).
  --infolder=<file>    Optional. Path to the single_cell_data.  
  --outfolder=<file>   Path to results folder.
"-> doc

library(docopt)
arguments <- docopt(doc, quoted_args=TRUE)
print(arguments)



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
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
install.packages('Seurat')

pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
          'limma', 'S4Vectors', 'SingleCellExperiment',
          'SummarizedExperiment', 'batchelor', 'devtools', 'ggplot2',
          'cowplot', "dplyr", "Seurat", "future", "sctransform",
          'Rtsne', "monocle3")

if (!require("BiocManager", character.only = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install()
} else {
  ipkgs <- sapply(pkgs, function(...) require(..., character.only = TRUE))
  if (any(!ipkgs)) {
    BiocManager::install(pkgs[!ipkkgs])
  } else {
    message("\n\nCool! your machine has everything is needed.\n\n")
  }
}


library("magrittr")
library("dplyr")
library("monocle3")
library("Seurat")
library(data.table)

# ------------- Import Cell Counts -------------------

c_counts <- read.delim('./geo/IPDCO_hg_midbrain_UMI.tsv', sep = "\t")
#head(c_counts)
dimz <- dim(c_counts)
paste0("Successefully read the data. Its length is ", dimz, ".",sep="")

# ------------- Create the Seurat object from cell counts (UMI matrix) -------------------

library(data.table)
library(Seurat)
ipdco <- CreateSeuratObject(counts = c_counts, project = "ipdco_janrogel_reprod_2021", min.cells = 3, min.features = 200)
ipdco

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

# The only data I found about the authors resolution used was 0.01
# There were also other resolutions tested (0.01 to 3) for visualization purposes

res <- 0.01

print("This step is performed using the FindNeighbors function, and takes as input the previously defined dimensionality of the dataset (first 10 PCs)")
ipdco <- FindNeighbors(ipdco,  assay = `assay`,
			     dims = 1:20, 
			     verbose = TRUE)

print("Finding clusters using the FindClusters function. Optimal resolution often increases for larger datasets. The clusters can be found using the Idents function.")
ipdco <- FindClusters(ipdco, graph.name = NULL,
			    modularity.fxn = 1, 
			    initial.membership = NULL, 
			    weights = NULL,
			    node.sizes = NULL, 
			    resolution = res, 
			    algorithm = 1, 
			    n.start = 10,
			    n.iter = 10, 
			    random.seed = 0, 
			    group.singletons = TRUE,
			    temp.file.location = NULL, 
			    edge.file.name = NULL, 
			    verbose = TRUE)

# Look at cluster IDs of the first 5 cells
head(Idents(ipdco), 5)

# ------------- Run non-linear dimensional reduction (UMAP/tSNE) ------------------

print("Cells within the graph-based clusters determined above should co-localize on these dimension reduction plots. As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis.")
ipdco <- RunUMAP(ipdco, assay = `assay`,
		       dims = 1:10,
		       reduction = "pca",
		       n.neighbors = 30,
		       metric = "cosine",
		       min.dist = 0.3, spread = 1, 
		       set.op.mix.ratio = 1,
		       local.connectivity = 1L, 
		       repulsion.strength = 1)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

DimPlot(ipdco, 
		     reduction = "umap", 
		     pt.size = 0.01,) + 
	  	ggplot2::theme(legend.position = "none")

	pdf("output/janr_tornado_clustering.pdf", height = 11.69, width = 8.27) 
	    plot(DimPlot(ipdco, 
		     reduction = "umap", 
		     pt.size = 0.01,) + 
	  	ggplot2::theme(legend.position = "none"))
	dev.off()

	jpeg(paste0("output/janr_tornado_clustering.jpeg"), 
			   height = 1980, width = 2980, pointsize = 74, quality = 100)
	    plot(DimPlot(ipdco, 
		     reduction = "umap", 
		     pt.size = 0.01,) + 
	  	ggplot2::theme(legend.position = "none"))
	dev.off()

saveRDS(ipdco, file = "/scratch/spectre/j/jr429/SRP/output/r-analysis-manual(UMI)")
print("Saved .rds file successfuly!")

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
head(ipdco.markers)

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

# Visualizing marker expression; shows expression probability distributions across clusters
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
