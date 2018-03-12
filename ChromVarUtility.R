# # # # # # # # # # # # # # # # # # # # # # # # # 
#
# script for running chromVAR
#
# author: Weizhang
#
# # # # # # # # # # # # # # # # # # # # # # # # # 

library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(Matrix)
library(pheatmap)
library(BiocParallel)
set.seed(2018)


Args <- commandArgs()
bamDir <- Args[6]  # bamfile dictionary
peakfile <- Args[7]  # peakfile path
thread <- Args[8]  # number of core
genomeidx <- Args[9]  # specify a genome
var.tiff <- Args[10]  # pdf file path
clu.tiff <- Args[11]  # pdf file path

if(genomeidx == "hg19"){
    library(BSgenome.Hsapiens.UCSC.hg19)
    genome <- BSgenome.Hsapiens.UCSC.hg19
    motifs <- getJasparMotifs(species = "Homo sapiens")
}else if(genomeidx == "hg38"){
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome <- BSgenome.Hsapiens.UCSC.hg38
    motifs <- getJasparMotifs(species = "Homo sapiens")
}else if(genomeidx == "mm9"){
    library(BSgenome.Mmusculus.UCSC.mm9)
    genome <- BSgenome.Mmusculus.UCSC.mm9
    motifs <- getJasparMotifs(species = "Mus musculus")
}else if(genomeidx == "mm10"){
    library(BSgenome.Mmusculus.UCSC.mm10)
    genome <- BSgenome.Mmusculus.UCSC.mm10
    motifs <- getJasparMotifs(species = "Mus musculus")
}

# extract all bam files
bamfile <- Sys.glob(file.path(bamDir, "*.bam"))

# peak file
peaks <- getPeaks(peakfile, sort_peaks = TRUE)

# set core
register(SnowParam(workers = 4, type = "SOCK"))

# Counts
colData = S4Vectors::DataFrame(celltype = as.character(seq(length(bamfile))))
example_counts <- getCounts(bamfile, peaks, paired =  TRUE, by_rg = FALSE,
                            format = "bam", colData = colData)

example_counts <- addGCBias(example_counts,
                            genome = BSgenome.Hsapiens.UCSC.hg19)

counts_filtered <- filterSamples(example_counts, min_depth = 1500, min_in_peaks = 0.15, shiny = FALSE)
counts_filtered <- filterPeaks(counts_filtered)

motif_ix <- matchMotifs(motifs, counts_filtered, genome = BSgenome.Hsapiens.UCSC.hg19)

# computing deviations
dev <- computeDeviations(object = counts_filtered, annotations = motif_ix)
variability <- computeVariability(dev)

options(bitmapType='cairo')

# plot variability
tiff(var.tiff)
plotVariability(variability, use_plotly = FALSE)
dev.off()


# clustering
sample_cor <- getSampleCorrelation(dev)
pheatmap(as.dist(sample_cor),
         annotation_row = colData(dev),
         clustering_distance_rows = as.dist(1-sample_cor),
         clustering_distance_cols = as.dist(1-sample_cor),
         filename = clu.tiff)


