# Created by: Fackie
# Created on: 2018/3/15

Args <- commandArgs()
print(Args)

input <- Args[6]
datatype <- Args[7]
output <- Args[8]
if (!dir.exists(output))
{
	dir.create(output)
}
setwd(output)

library(Seurat)
library(dplyr)
library(Matrix)

# celldata = Matrix::readMM(file = matrix)
# cellname = read.table(file = barcodes, header = FALSE, colClasses = "character")[[1]]
# genename = read.table(file = genes, header = FALSE, colClasses = "character")[[1]]
# rownames(celldata) = genename
# colnames(celldata) = cellname

#load data
if (datatype == '10x')
{
celldata = Read10X(data.dir = input)
}
if (datatype == 'densemtx')
{
    celldata = t(t(read.table(input, sep = ',',row.names= 1,header=TRUE,stringsAsFactors =FALSE,fileEncoding='utf-8')))
}
test = CreateSeuratObject(raw.data = celldata,project = "Frankie so tired")
rm(celldata)
#processing
#mito.genes <- grep(pattern = "^MT-", x = rownames(x = test@data), value = TRUE)
# percent.mito <- Matrix::colSums(test@raw.data[mito.genes, ])/Matrix::colSums(test@raw.data)
# test <- AddMetaData(object = test, metadata = percent.mito, col.name = "percent.mito")

##############################################################################################
#file.create("violinplot.jpeg")
jpeg(file="violinplot.jpeg")
VlnPlot(object = test, features.plot = c("nGene", "nUMI"), nCol = 2)
dev.off()

##############################################################################################
#file.create("geneplot.jpeg")
jpeg(file="geneplot.jpeg")
GenePlot(object = test, gene1 = "nUMI", gene2 = "nGene")
dev.off()
file.remove('Preprocessing.Rdata')
save(test,file = 'Preprocessing.Rdata')
