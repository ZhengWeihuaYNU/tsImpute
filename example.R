#install.packages('Seurat')
#install.packages('factoextra')
#install.packages('umap')
#install.packages('flexclust') 
library(Seurat)
library(factoextra)
library(umap)
library(flexclust) #This package is not necessary for tsImpute, we use it to calculate adjusted Rand Index
setwd('Your path of the data')
tmp<- read.csv('gse67835.csv', row.names = 1, header= T)
dim(tmp) #466 cells
colnames(tmp)<- substr(colnames(tmp), 1,10)
labels<- read.csv('cell_types.csv', header= F) #load the ground truth labels
colnames(tmp)
colnames(tmp)== labels$V1
cell_types<- labels$V2[match(colnames(tmp), labels$V1)]
cell_types<- as.factor(cell_types) #transform the labels into factors
tmp<- as.matrix(tmp)
unique(cell_types) #9 cell types
rownames(tmp)<- gsub('_','-', rownames(tmp)) 
#filter out genes expressing in less than 5% of the cells:
qc_counts<- filter_genes(tmp, round(0.05*ncol(tmp))) 
#select 2000 hvgs with Seurat:
x.seurat <- CreateSeuratObject(qc_counts)
x.seurat <- NormalizeData(x.seurat)
x.seurat <- ScaleData(x.seurat)
x.seurat <- FindVariableFeatures(x.seurat, verbose = FALSE)
filter_ID <- x.seurat@assays$RNA@var.features
filtered<- qc_counts[filter_ID,]
output<- tsimpute(filtered, seed= 1) #run tsImpute


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
x.seurat <- CreateSeuratObject(filtered) 
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
