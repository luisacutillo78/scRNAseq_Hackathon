# Monocle3 Example

## 1. First we install all the necessary packages.
Additional information on possible errors during installation process can be found [HERE](https://cole-trapnell-lab.github.io/monocle3/docs/installation/)

```R
if(!require(monocle3))
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(version = "3.10")
  
  BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                         'limma', 'S4Vectors', 'SingleCellExperiment',
                         'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
  
  install.packages("devtools")
  devtools::install_github('cole-trapnell-lab/leidenbase')
  devtools::install_github('cole-trapnell-lab/monocle3')
  
  devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")
  
  library(monocle3)
}
```


## 2. We load into the R environment the downloaded Tabula Muris datasets

```R
# Load data

```

## 3. Gene Filtering and normazlization
Data normalization addresses the unwanted biases arisen by count depth variability while preserving true biological differences.

```R
# Gene filtering and data normalization
cds <- new_cell_data_set(as(mtx, "sparseMatrix"),cell_metadata = cell_metadata,gene_metadata = gene_metadata)
rm(mtx);gc()
```

## 4. Data summarization & Dimensionality Reduction
Dimensionality reduction aims to condense the complexity of the data into a lower-dimensional space by optimally preserving its key properties.

```R
# PCA for data summarization
cds <- preprocess_cds(cds, num_dim = 50)

# Dimensionality reduction with UMAP 
cds <- reduce_dimension(cds) #dimensionality reduction, default value is UMAP
```


## 5. Clustering Analysis: how to identify cellular sub-populations
As transcriptionally distinct populations of cells usually correspond to distinct cell types, a key goal of scRNA-seq consists in the identification of cell subpopulations based on their transcriptional similarity. Thus, organizing cells into groups (i.e. clusters) can allow for de novo detection of cell types or identification of different subpopulations in a single cell state.

```R
# Identify clusters
cds = cluster_cells(cds, cluster_method="louvain") #cell clustering using louvain algorithm
```

## 6. Trajectory analysis: Partition cells into supergroups
Rather than forcing all cells into a single developmental trajectory, Monocle 3 enables you to learn a set of trajectories that describe the biological process you are studying.

```R
# Partition cells
cds <- partitionCells(cds) #separate cells into supergroups to enable trajectory inference for each population
```

## 7. Trajectory analysis: Learn the principal graph
Once cells are paritioned, each supergroup can be organized into a separate trajectory. The default method for doing this in Monocle 3 is SimplePPT, which assumes that each trajectory is a tree (albeit one that may have multiple roots). Finally, each cell will be assigned with a pseudotime value. In order to do so, you need to specify the root nodes of the trajectory graph. If you do not provide them as an argument, it will launch a graphical user interface for selecting one or more root nodes.

```R
# Learn the trajectories
cds <- learnGraph(cds,  RGE_method = 'SimplePPT') #trajectory inference using SimplePPT

# Helper function to identify the root principal points:
get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  
  closest_vertex <-
    cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    V(cds@minSpanningTree)$name[as.numeric(names
                                           (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

# Assign pseudotime value
node_ids = get_correct_root_state(cds,cell_phenotype = 'cell_type', root_type = "cell_type_1") #where cell_phenotypes indicates the column where cell types are stored, while root_type is the cell_type we want to select.
cds <- orderCells(cds, root_pr_nodes = node_ids)
```


## 8. Plots
```R
plot_cells(cds,
           color_cells_by = "pseudotime", #change here to color cells by the metadata of choice
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
```
