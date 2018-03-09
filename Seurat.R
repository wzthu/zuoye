 Args <- commandArgs()
 matrixdata <- Args[6]
 barcodes <- Args[7]
 genes <- Args[8]
 outputpath <- Args[9]
 #gene_lowthread <- Args[10]
 #gene_highthread <- Args[11]
 print(matrixdata)
 print(barcodes)
 print(genes)
 print(outputpath)
 library(Seurat)
 library(dplyr)
 library(Matrix)
 #install.packages('Seurat')
 celldata <- Matrix::readMM(file = matrixdata)
 cellname <- read.table(file = barcodes, header = FALSE, colClasses = "character")[[1]]
 genename <- read.table(file = genes, header = FALSE, colClasses = "character")[[1]]
 rownames(celldata) <- genename
 colnames(celldata) <- cellname
 #VLplot
 cell <- CreateSeuratObject(raw.data = celldata, min.cells = 3, min.genes = 200, 
                             project = "10X_data")
 jpeg(filename = as.character(paste(as.character(outputpath),"/", "VlnPlot.jpg",sep = "")))
 VlnPlot(object = cell, features.plot = c("nGene", "nUMI"), nCol = 2)
 dev.off()
 jpeg(filename = as.character(paste(as.character(outputpath),"/", "GenePlot.jpg",sep = "")))
 GenePlot(object = cell, gene1 = "nUMI", gene2 = "nGene")
 dev.off()
 #cell <- FilterCells(object = cell, subset.names = "nGene", 
                      #low.thresholds = c(gene_lowthread), high.thresholds = c(gene_highthread))
 #Normalizing the data
 cell <- NormalizeData(object = cell, normalization.method = "LogNormalize", 
                        scale.factor = 10000)
 #Detection of variable genes across the single cells
 jpeg(filename = as.character(paste(as.character(outputpath),"/", "FindVariableGenes.jpg",sep = "")))
 cell <- FindVariableGenes(object = cell, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
 dev.off()
 length(x = cell@var.genes)
 #Scaling the data and removing unwanted sources of variation
 cell <- ScaleData(object = cell, vars.to.regress = "nUMI")
 #Perform linear dimensional reduction
 cell <- RunPCA(object = cell, pc.genes = cell@var.genes, do.print = TRUE, pcs.print = 1:5, 
                 genes.print = 5)
 # Examine and visualize PCA results a few different ways
 # Prints a set of genes that most strongly define a set of principal components
 PrintPCA(object = cell, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
 #Visualize top genes associated with principal components
 #VizPCA(object = cell, pcs.use = 1:2)
 #PCAPlot(object = cell, dim.1 = 1, dim.2 = 2)  
 # ProjectPCA scores each gene in the dataset (including genes not included
 # in the PCA) based on their correlation with the calculated components.
 # Though we don't use this further here, it can be used to identify markers
 # that are strongly correlated with cellular heterogeneity, but may not have
 # passed through variable gene selection.  The results of the projected PCA
 # can be explored by setting use.full=T in the functions above
 cell <- ProjectPCA(object = cell, do.print = FALSE)
 # NOTE: This process can take a long time for big datasets, comment out for
 # expediency.  More approximate techniques such as those implemented in
 # PCElbowPlot() can be used to reduce computation time
 cell <- JackStraw(object = cell, num.replicate = 100, do.print = FALSE)
 #JackStrawPlot(object = cell, PCs = 1)
 #PCElbowPlot(object = cell)
 cell <- FindClusters(object = cell, reduction.type = "pca", dims.use = 1:10, 
                       resolution = 0.5, print.output = 0, save.SNN = TRUE, force.recalc = T) 
 PrintFindClustersParams(object = cell) 
 cell <- RunTSNE(object = cell, dims.use = 1:10, do.fast = TRUE)
 jpeg(filename = as.character(paste(as.character(outputpath),"/", "tSNE_findcluster.jpg",sep = "")))
 TSNEPlot(object = cell)
 dev.off()
 
 
 
 