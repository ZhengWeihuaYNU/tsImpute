#Seurat, factoextra and umap are prerequisite packages for tsImpute, make sure that you have installed
#these packages before you install tsImpute.
#install.packages('Seurat')
#install.packages('factoextra')
#install.packages('umap')
#install.packages('flexclust') 
install.packages("your path of tsImpute_0.1.0.tar.gz", repos = NULL, type = "source")
library(Seurat)
library(factoextra)
library(umap)
library(flexclust) #This package is not necessary for tsImpute, we use it to calculate adjusted Rand Index
library(tsImpute)
#setwd('Your path of the data')
#read data and corresponding labels:
dat<- readRDS('darmanis.rds') 
cell_types<- readRDS('cell_types.rds')
#run tsimpute:
output<- tsimpute(dat, seed= 1) 


#use Seurat package to cluster the cells and calculate adjusted Rand index (ARI):
x.seurat <- CreateSeuratObject(output) 
x.seurat <- NormalizeData(x.seurat)
x.seurat <- ScaleData(x.seurat)
x.seurat <- FindVariableFeatures(x.seurat, verbose = FALSE)
x.seurat <- RunPCA(x.seurat, features = VariableFeatures(object = x.seurat), seed.use = 1)  
x.seurat <- JackStraw(x.seurat, num.replicate = 100)
x.seurat <- ScoreJackStraw(x.seurat, dims = 1:20)
x.seurat <- FindNeighbors(x.seurat, dims = 1:10) #Use the first 10 components for clustering
x.seurat <- FindClusters(x.seurat, resolution = 0.5, random.seed = 1)
final_cluster<- x.seurat$seurat_clusters
ari<- randIndex(cell_types, final_cluster) #ARI
print(ari) #0.800

#Compare it to the clustering results of raw data:
x.seurat <- CreateSeuratObject(dat) 
x.seurat <- NormalizeData(x.seurat)
x.seurat <- ScaleData(x.seurat)
x.seurat <- FindVariableFeatures(x.seurat, verbose = FALSE)
x.seurat <- RunPCA(x.seurat, features = VariableFeatures(object = x.seurat), seed.use = 1)  
x.seurat <- JackStraw(x.seurat, num.replicate = 100)
x.seurat <- ScoreJackStraw(x.seurat, dims = 1:20)
x.seurat <- FindNeighbors(x.seurat, dims = 1:10) #Use the first 10 components for clustering
x.seurat <- FindClusters(x.seurat, resolution = 0.5, random.seed = 1)
raw_cluster<- x.seurat$seurat_clusters
raw_ari<- randIndex(cell_types, raw_cluster) #ARI
print(raw_ari) #0.663
