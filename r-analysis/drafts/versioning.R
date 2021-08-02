#UMI cell counts --------------------------------------------------------------
#Check versions so they match the authors
A_Versions <- c('CellRanger' = 3.0, 'Seurat3' = 4.0, 'monocle3' = 4.0)
CellRanger <- 3.0
Seurat3 <- 4.0
monocle3<- 4.0

A_Versions <- data.frame(CellRanger,Seurat3,monocle3)
format.data.frame(A_Versions, nsmall = 1)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocCheck")


if (packageVersion("Seurat")>="4.0.0") {
    print("version is ok")
}

# ------------- Dependencies -------------------

#install.packages("magrittr")
#install.packages("devtools")
#devtools::install_github('cole-trapnell-lab/leidenbase')
#devtools::install_github('cole-trapnell-lab/monocle3')
#install.packages('Seurat')

#library("magrittr")
#library("monocle3")
#library("Seurat")
