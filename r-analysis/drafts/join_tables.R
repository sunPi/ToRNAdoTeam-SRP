library(data.table)

gene_markers <- fread("/home/jr429/Documents/UoL/SteeredResearchProject/Tornado-Group-scRNA-seq-Project/r-analysis/ipdco.gene.markers.tsv")
clusters <- fread("/home/jr429/Documents/UoL/SteeredResearchProject/Tornado-Group-scRNA-seq-Project/r-analysis/ipdco.markers.tsv")

head(clusters)
head(gene_markers)

joined.df <- merge(x=gene_markers, y=clusters, by = "gene", all.x = TRUE)
View(joined.df)

fwrite(joined.df, "/home/jr429/Documents/UoL/SteeredResearchProject/Tornado-Group-scRNA-seq-Project/r-analysis/ipdco.janr.final.tsv", sep="\t", col.names = TRUE, na="", append = TRUE)
