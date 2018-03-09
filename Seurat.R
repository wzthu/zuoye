 Args <- commandArgs()
 matrixdata <- Args[6]
 barcodes <- Args[7]
 genes <- Args[8]
 outputpath <- Args[9]
 print(matrixdata)
 print(barcodes)
 print(genes)
 print(outputpath)
 library(Seurat)
 library(dplyr)
 library(Matrix)
 #library(stringr)
 # inputpath <- str_replace_all(inputpath, "\\\\", "/") 
 # outputpath <- str_replace_all(outputpath, "\\\\", "/") 


 
 #install.packages('Seurat')

 # Load the dataset
 #mouse.data <- Read10X(data.dir = "F:/HCA/sample3v3/outs/filtered_gene_bc_matrices/mm10")
 #mouse.data <- Read10X(data.dir = "F:/HCA/pbmc4k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/GRCh38")
 #mouse.data <- Read10X(data.dir = inputpath)
 celldata <- Matrix::readMM(file = matrixdata)
 cellname <- read.table(file = barcodes, header = FALSE, colClasses = "character")[[1]]
 genename <- read.table(file = genes, header = FALSE, colClasses = "character")[[1]]
 rownames(celldata) <- genename
 colnames(celldata) <- cellname
 # Examine the memory savings between regular and sparse matrices
 # dense.size <- object.size(x = as.matrix(x = mouse.data))
 # dense.size 
 # sparse.size <- object.size(x = mouse.data)
 # sparse.size
 # Initialize the Seurat object with the raw (non-normalized data).  Keep all
 # genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
 # least 200 detected genes
 cell <- CreateSeuratObject(raw.data = celldata, min.cells = 3, min.genes = 200, 
                             project = "10X_mouse")
 # The number of genes and UMIs (nGene and nUMI) are automatically calculated
 # for every object by Seurat.  For non-UMI data, nUMI represents the sum of
 # the non-normalized values within a cell.
 jpeg(filename = as.character(paste(as.character(outputpath),"/", "VlnPlot.jpg",sep = "")))
 VlnPlot(object = cell, features.plot = c("nGene", "nUMI"), nCol = 2)
 dev.off()
 jpeg(filename = as.character(paste(as.character(outputpath),"/", "GenePlot.jpg",sep = "")))
 GenePlot(object = cell, gene1 = "nUMI", gene2 = "nGene")
 dev.off()
 