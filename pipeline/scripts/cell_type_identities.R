# ------------- Assigning cell type identity to clusters ------------------

new.cluster.ids <- c("1", "2", "3", "4", "5", "6", "7", "8",
    "9", "10", "11", "12", "13", "14", "15")
names(new.cluster.ids) <- levels(ipdco)
ipdco <- RenameIdents(ipdco, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

