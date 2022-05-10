---
title: "Analysis of scRNAseq using Seurat"
author: A. Carissimo
output:
  html_document:
    toc: true
    toc_float: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning=FALSE, message=FALSE)
```

# Single cell RNASeq step-by-step 

## Package loading 
 

```{r message=FALSE}
library(dplyr)
library(monocle3)
library(Seurat)
library(patchwork)
library(mclust)
library(SeuratWrappers)
```

### Reading Tabula Muris data and metadata

```{r message=FALSE}
tabData=readRDS("emptyDroplets_doublet_filtered_tabulamuris_mtx.rds")
info=readRDS("tabulamuris_cell_info.rds")
tabmur <- CreateSeuratObject(counts = tabData, project = "tabulamuris", min.cells = round(dim(tabData)[2]*5/100), min.features = 0, meta.data=info)


tabmur
tabmur[["RNA"]]@counts[1:4,1:4]
summary(Matrix::colSums(tabData))
hist(Matrix::colSums(tabData), xlab="Total Counts", breaks=100,  ylab="Number of cells")       
```



```{r}
#tabmur[["percent.mt"]] <- PercentageFeatureSet(tabmur, pattern = "^mt-")

plot1 <- VlnPlot(tabmur, features = c("nFeature_RNA"), ncol=1)
plot1
plot2 <- VlnPlot(tabmur, features = c("nCount_RNA"), ncol = 1)
plot2

#plot4 <- VlnPlot(tabmur, features = c("percent.mt"), ncol = 1)
#plot4

plot3 <- FeatureScatter(tabmur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "tissue")
plot3

```

```{r}

#tabmur <- subset(tabmur, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#tabmur
```

###  Data normalization: how to make gene expression comparable across individual cells

**Seurat** uses log transformation of the CPM to reduce cell depth variability. We next standardized expression values for each gene across all cells (z-score transformation), as is standard prior to running dimensional reduction tools such as principal component analysis. These steps are implemented in the NormalizeData and ScaleData functions in Seurat.  

```{r}
tabmur <- NormalizeData(tabmur, normalization.method = "LogNormalize", scale.factor = 10000)

tabmur[["RNA"]]@data[1:5,1:5]
```

### Feature selection: how to discard “uninformative” genes

Feature selection aims to detect genes with relevant biological information, while excluding the uninformative ones.
**Seurat** first models the mean-variance relationship using a local polynomial regression function. Then, given the expected variance by the fitted curve and the observed mean, the feature values are standardized, and for each gene, the variance across all cells is computed


```{r}
tabmur <- FindVariableFeatures(object = tabmur, selection.method="vst", nfeatures = 2000)
top10 <- head(VariableFeatures(tabmur), 10)
plot1 <- VariableFeaturePlot(tabmur)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
```

### Scale the data
Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData function:

Shifts the expression of each gene, so that the mean expression across cells is 0 setting *do.center=TRUE*
Scales the expression of each gene, so that the variance across cells is 1 setting *do.scale=TRUE*


```{r}
tabmur <- ScaleData(tabmur, do.center = TRUE, do.scale = TRUE)
tabmur[["RNA"]]@scale.data[1:5,1:5]
```

### Linear dimensional reduction: for the summarization of scRNA-seq data

Dimensionality reduction aims to condense the complexity of the data into a lower-dimensional space by optimally preserving its key properties. Dimensionality reduction methods are essential for clustering, visualization, and summarization of scRNA-seq data. PCA is used to summarise a dataset throughout the top N principal components. The number of PCA to use is usually determined by manually inspecting the elbow plot (*ElbowPlot* function), in which principal components are plotted as a function of the variability they account for, and the number of PCA to use is determined by the point in which an “elbow” is observed. Additional methods can be used, including jackstraw 


```{r}

tabmur <- RunPCA(tabmur, features = VariableFeatures(object = tabmur))

plot1 <- DimPlot(tabmur, reduction = "pca")
plot1

plot2 <- VizDimLoadings(tabmur, dims = 1:2, reduction = "pca")
plot2

tabmur <- JackStraw(tabmur, num.replicate = 100)
tabmur <- ScoreJackStraw(tabmur, dims = 1:20)
JackStrawPlot(tabmur, dims = 1:15)
ElbowPlot(tabmur,ndims=100)
```

### Clustering Analysis: how to identify cellular sub-populations

As transcriptionally distinct populations of cells usually correspond to distinct cell types, a key goal of scRNA-seq consists in the identification of cell subpopulations based on their transcriptional similarity.
**Seurat** constructs a KNN graph based on the euclidean distance in PCA space, and refines the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). This step is performed using the *FindNeighbors* function, and takes as input the first 50 PCs.

To cluster the cells, **Seurat** applies a modularity optimization technique such as the Louvain algorithm (default) or SLM. The *FindClusters* function implements this procedure. When running a graph-based clustering, it is necessary to set the resolution parameter for the community detection algorithm based on modularity optimization. The resolution parameter is correlated to the scale of observing communities. In particular, the higher is the resolution parameter, the larger is the number of smaller communities. 
Here, we set the resolution parameter to 0.5

```{r}
tabmur <- FindNeighbors(tabmur, dims = 1:50)
tabmur <- FindClusters(tabmur, resolution = 0.5)
 
```
```{r}
table(Idents(tabmur))
```

### Nonlinear dimensionality reduction for the visualization of scRNA-seq data

Dimensionality reduction for visualization of scRNA-data uses methods that capture the nonlinearity of the scRNA-seq data, avoiding the overcrowding of the representation. Here we use UMAP and t-SNE


```{r}
tabmur <- RunUMAP(tabmur, dims = 1:50)
DimPlot(tabmur, group.by = "ident", reduction = "umap")
DimPlot(tabmur, group.by = "tissue", reduction = "umap")


tabmur <- RunTSNE(tabmur, dims = 1:50)
DimPlot(tabmur, group.by = "ident", reduction = "tsne")
DimPlot(tabmur, group.by = "tissue", reduction = "tsne")


```

### Results evaluation 


```{r}
adjustedRandIndex(Idents(tabmur), tabmur@meta.data$cell_ontology_class)

 
```


### Finding differentially expressed features

Seurat identifies markers of a single cluster compared to all other cells. FindAllMarkers() automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

```{r}
tabmur.markers <- FindAllMarkers(tabmur, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tabmur.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

 
```


### Find all markers of cluster 2

```{r}
cluster2.markers <- FindMarkers(tabmur, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
```

# Find all markers distinguishing cluster 5 from clusters 0 and 3

```{r}
cluster5.markers <- FindMarkers(tabmur, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
```

### Markers visualization

```{r}
plot1<-VlnPlot(tabmur, features = c("Cd79a", "Cd79b"),idents = c(0,3,5))
plot1

plot2<-FeaturePlot(tabmur, features = c("Cd79a", "Cytl1", "Cd3g", "Krt14", "Fabp4"))
plot2

tabmur.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

```

### Trajectories using Monocle

```{r}
monocle_object <- as.cell_data_set(tabmur)
monocle_object <- cluster_cells(cds = monocle_object, reduction_method = "UMAP")

p1<- plot_cells(monocle_object, show_trajectory_graph = FALSE)
p2 <- plot_cells(monocle_object, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

```

```{r}

#integrated.sub <- subset(as.Seurat(monocle_object), monocle3_partitions == 2)
#cds <- as.cell_data_set(integrated.sub)
```


```{r}
monocle_object <- learn_graph(monocle_object)
plot4<- plot_cells(monocle_object, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
plot4

get_correct_root_state <- function(cds,cell_phenotype="cell_ontology_class",time_bin="mesenchymal stem cell"){
  cell_ids <- which(colData(cds)[, cell_phenotype] == time_bin)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}


node_ids = get_correct_root_state(monocle_object,cell_phenotype = "cell_ontology_class", time_bin ="mesenchymal stem cell") 

#where cell_phenotypes indicates the column where cell types are stored, while root_type is the cell_type we want to select.


monocle_object <- order_cells(monocle_object,root_pr_nodes = node_ids)
```

```{r}
plot_cells(monocle_object,
       color_cells_by = "pseudotime",
       graph_label_size=1,
       show_trajectory_graph = TRUE)
```


