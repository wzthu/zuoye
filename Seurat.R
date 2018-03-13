# setwd("C:/Users/Fackie/Desktop/HCA task/µÚ¶þÖÜ/filtered_gene_bc_matrices/hg19")

library(Seurat)
library(dplyr)
library(Matrix)

# celldata = Matrix::readMM(file = 'matrix.mtx')
# cellname = read.table(file = 'barcodes.tsv', header = FALSE, colClasses = "character")[[1]]
# genename = read.table(file = 'genes.tsv', header = FALSE, colClasses = "character")[[1]]
# rownames(celldata) = genename
# colnames(celldata) = cellname

#load data(if datatype is 10X)
celldata <- Read10X("F:\\2018\\Researchs\\chordoma\\hg19")
#load data(if datatype is dense matrix)
#input files should be 
celldata <- read.table("~/Downloads/immune_control_expression_matrix.txt.gz", sep = ",")

test = CreateSeuratObject(raw.data = celldata,project = "Frankie so tired")
rm(celldata)

# processing
# mitogenes portion(unnecessary)
# mito.genes <- grep(pattern = "^MT-", x = rownames(x = test1@data), value = TRUE)
# percent.mito <- Matrix::colSums(test@raw.data[mito.genes, ])/Matrix::colSums(test@raw.data)
# test <- AddMetaData(object = test, metadata = percent.mito, col.name = "percent.mito")

##############################################################################################
#file.create("violinplot.jpeg")
jpeg(file="violinplot.jpeg")
VlnPlot(object = test, features.plot = c("nGene", "nUMI"), nCol = 2)
dev.off()

##############################################################################################
#file.create("geneplot.jpeg")
#jpeg(file="geneplot.jpeg")
GenePlot(object = test, gene1 = "nUMI", gene2 = "nGene")
#GenePlot(object = test, gene1 = "nUMI", gene2 = "percent.mito")
# GenePlot(object = test, gene1 = "nUMI", gene2 = "nGene")
#dev.off()

humi <- quantile(test@meta.data$nUMI,0.95)
hgene <- quantile(test@meta.data$nGene,0.95)
lumi <- quantile(test@meta.data$nUMI,0.05)
lgene <- quantile(test@meta.data$nGene,0.05)

test <- FilterCells(object = test, subset.names = c("nGene","nUMI"), low.thresholds = c(lgene, lumi), high.thresholds = c(hgene,humi))


test <- NormalizeData(object = test, normalization.method = "LogNormalize",scale.factor = 10000)

##############################################################################################
#file.create("variableGenes.jpeg")
#jpeg(filename ="variableGenes.jpeg")
test <- FindVariableGenes(object = test, mean.function = ExpMean, dispersion.function = LogVMR)
#dev.off()


test <- ScaleData(object = test,vars.to.regress = "nUMI")
test <- RunPCA(object = test)

##############################################################################################
#file.create("Elbowplot.jpeg")
jpeg(file="Elbowplot.jpeg")
PCElbowPlot(object = test,num.pc = 30)
dev.off()

test <- FindClusters(object = test, reduction.type = "pca", print.output = 1)
test <- RunTSNE(object = test,dims.use = 1:10, perplexity =10, do.fast = TRUE)

##############################################################################################
#file.create("TSNEplot.jpeg")
jpeg(file="TSNEplot.jpeg")
TSNEPlot(object = test)
dev.off()

# FeaturePlot(test,features.plot = "S100B",cols.use = c("white","red"))
# cluster3.markers <- FindMarkers(object = test, ident.1 = 1, logfc.threshold = 0.8)
# print(x = head(x = cluster4.markers, n = 5))

