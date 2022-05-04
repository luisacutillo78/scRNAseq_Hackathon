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


## 2. We load into the R environment C. Elegans embryo datasets from Packer & Zhu et al.

```R
# Load data
expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_expression.rds"))
cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_colData.rds"))
gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_rowData.rds"))
```

## 3. Gene Filtering and normazlization
Data normalization addresses the unwanted biases arisen by count depth variability while preserving true biological differences.

```R
# Gene filtering and data normalization
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
rm(expression_matrix);gc()
```

## 4. Data summarization & Dimensionality Reduction
Dimensionality reduction aims to condense the complexity of the data into a lower-dimensional space by optimally preserving its key properties.

```R
# PCA for data summarization
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")

# Dimensionality reduction with UMAP 
cds <- reduce_dimension(cds) #dimensionality reduction, default value is UMAP
```


## 5. Clustering Analysis: how to identify cellular sub-populations
As transcriptionally distinct populations of cells usually correspond to distinct cell types, a key goal of scRNA-seq consists in the identification of cell subpopulations based on their transcriptional similarity. Thus, organizing cells into groups (i.e. clusters) can allow for de novo detection of cell types or identification of different subpopulations in a single cell state.
Although cells may continuously transition from one state to the next with no discrete boundary between them, Monocle does not assume that all cells in the dataset descend from a common transcriptional "ancestor". In many experiments, there might in fact be multiple distinct trajectories. Therefore, Monocle assignes each cell not only to a cluster but also to a partition. When you are learning trajectories, each partition will eventually become a separate trajectory.

```R
# Identify clusters
cds = cluster_cells(cds, cluster_method="louvain") #cell clustering using louvain algorithm
```

## 6. Trajectory analysis: Learn the principal graph
Once cells are paritioned, each supergroup can be organized into a separate trajectory. The default method for doing this in Monocle 3 is SimplePPT, which assumes that each trajectory is a tree (albeit one that may have multiple roots). Finally, each cell will be assigned with a pseudotime value. In order to do so, you need to specify the root nodes of the trajectory graph. If you do not provide them as an argument, it will launch a graphical user interface for selecting one or more root nodes.

```R
# Learn the trajectories
cds <- learn_graph(cds) #trajectory inference using SimplePPT

# Helper function to identify the root principal points:
get_correct_root_state <- function(cds,cell_phenotype="embryo.time.bin",time_bin="130-170"){
  cell_ids <- which(colData(cds)[, cell_phenotype] == time_bin)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

# Assign pseudotime value
node_ids = get_correct_root_state(cds,cell_phenotype = "embryo.time.bin", time_bin ="130-170") #where cell_phenotypes indicates the column where cell types are stored, while root_type is the cell_type we want to select.
cds <- order_cells(cds, root_pr_nodes = node_ids)
```


## 7. Plots
```R
plot_cells(cds,
           color_cells_by = "pseudotime", #change here to color cells by the metadata of choice
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
```
