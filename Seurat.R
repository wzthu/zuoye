Args <- commandArgs()
print(Args)

matrix <- Args[0]
barcodes <- Args[1]
genes <- Args[2]
#projectname <- Args[3]
output <- Args[3]

library(Seurat)
library(dplyr)
library(Matrix)

read10x <- function(matrix,barcode,genes)
{
celldata = Matrix::readMM(file = matrix)
cellname = read.table(file = barcode, header = FALSE, colClasses = "character")[[1]]
genename = read.table(file = genes, header = FALSE, colClasses = "character")[[1]]
rownames(celldata) = genename
colnames(celldata) = cellname
return(celldata)
}

#load data
rawdata = read10X(matrix,barcode,genes)
test = CreateSeuratObject(raw.data = rawdata,project = "Frankie so tired")

#processing
mito.genes <- grep(pattern = "^MT-", x = rownames(x = test@data), value = TRUE)
percent.mito <- Matrix::colSums(test@raw.data[mito.genes, ])/Matrix::colSums(test@raw.data)

test <- AddMetaData(object = test, metadata = percent.mito, col.name = "percent.mito")

##########################
jpeg(file=output+"Violinplot.jpeg")
VlnPlot(object = test, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

########################
jpeg(file=output+"geneplot.jpeg")
par(mfrow = c(1, 2))
GenePlot(object = test, gene1 = "nUMI", gene2 = "nGene")
dev.off()


test <- FilterCells(object = test, subset.names = c("nGene","nUMI"), 
                    low.thresholds = c(-Inf, -Inf), high.thresholds = c(3500,10000))
test <- NormalizeData(object = test, normalization.method = "LogNormalize",scale.factor = 10000)
test <- FindVariableGenes(object = test, mean.function = ExpMean,x.low.cutoff = 0.05, x.high.cutoff = 3, y.cutoff = 1)
test <- ScaleData(object = test,vars.to.regress = "nUMI")
test <- RunPCA(object = test, pc.genes = test@var.genes,pcs.compute = 30)

jpeg(file=output+"Elbowplot.jpeg")
PCElbowPlot(object = test,num.pc = 30)
dev.off()

test <- FindClusters(object = test, reduction.type = "pca", dims.use = 1:5, resolution = 0.3,
                      print.output = 1, save.SNN = TRUE)
test <- RunTSNE(object = test, dims.use = 1:5,do.fast = TRUE)

jpeg(file=output+"TSNEplot.jpeg")
TSNEPlot(object = test)
dev.off()

# FeaturePlot(test,features.plot = "S100B",cols.use = c("white","red"))
# cluster3.markers <- FindMarkers(object = test, ident.1 = 1, logfc.threshold = 0.8)
# print(x = head(x = cluster4.markers, n = 5))

 