#import library
library("SC3")
library("SingleCellExperiment")

#import & preprocess  argv
Args <- commandArgs()
print(Args)
sce_file <- Args[6]
clusters_num <- as.numeric(Args[7])
if(clusters_num==0 )
    {
    estimate_k_flag = TRUE
    }else{
    estimate_k_flag =FALSE
    }
result_folder <- Args[8]

processFlag <- Args[9]


# define processing function
cluster <- function(sce,clusters_num,result_folder)
{
    jpeg(paste(result_folder,"Consensus_Cluster.jpg",sep = "/") )
    sc3_plot_consensus(sce, k = clusters_num )
    dev.off()
    sc3_export_results_xls(sce,filename = paste(result_folder,"Result.xls",sep = "/"))
}

de <- function(sce,clusters_num,result_folder)
{
    jpeg(paste(result_folder,"Expression.jpg",sep = "/") )
    sc3_plot_expression(sce, k = clusters_num)
    dev.off()

    jpeg(paste(result_folder,"DE_Gene.jpg",sep = "/") )
    sc3_plot_de_genes(sce, k = clusters_num)
    dev.off()

    jpeg(paste(result_folder,"Gene_Marker.jpg",sep = "/") )
    sc3_plot_markers(sce, k = clusters_num)
    dev.off()
}


# load sce object
load(sce_file)
    
# if need to auto estimate K 
if(estimate_k_flag ==TRUE)
    {
    show("Auto setting clusters number")
    sce <- sc3_prepare(sce)
    sce <- sc3_estimate_k(sce)
    clusters_num <- metadata(sce)$sc3$k_estimation
    }

# generate sc3 object 
sce <- sc3(sce, ks = clusters_num, biology = TRUE)

# switch which function to be used
switch(processFlag,
    "cluster"=cluster(sce,clusters_num,result_folder),
    "de"=de(sce,clusters_num,result_folder),
    show("Wrong ProcessFlag, No process to be done."))

# save sc3 object
save(sce,file=paste(result_folder,paste("sce",processFlag,sep='_'),sep = "/"))
