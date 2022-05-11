# Monocle3 Example

## 1. We load into the R environment Nestrowa or Pbmc Data (based on the assignment)

```R
# Load data
data <- "" # Based on the assignment, one of "Nestrowa.rds" or "PbmcData.rds"
expression_matrix <- readRDS(data)
```

## 2. Gene/Cell Filtering and normazlization
Genes and cells are filtered in Seurat, then we move to Monocle3 to proceed with the analysis.
Data normalization addresses the unwanted biases arisen by count depth variability while preserving true biological differences.

```R
# Gene filtering and data normalization
seu_data <- CreateSeuratObject(counts = expression_matrix, project = gsub(".rds","",data), min.cells = round(dim(expression_matrix)[2]*5/100), min.features = 0)
seu_data[["percent.mt"]] <- PercentageFeatureSet(seu_data, pattern = "^mt-") #If data is human, use "^MT-" as pattern

seu_data <- subset(seu_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) #Nestrowa data has high sequencing depth, so remove the "nFeature_RNA < 2500" argument

cds <- as.cell_data_set(seu_data)

rm(expression_matrix);gc()
```

## 3. Data summarization & Dimensionality Reduction
Dimensionality reduction aims to condense the complexity of the data into a lower-dimensional space by optimally preserving its key properties.

```R
# PCA for data summarization
cds <- preprocess_cds(cds, num_dim = 50)

# Dimensionality reduction with UMAP 
cds <- reduce_dimension(cds, reduction_method = 'UMAP') #dimensionality reduction, default value is UMAP
```


## 4. Clustering Analysis: how to identify cellular sub-populations
As transcriptionally distinct populations of cells usually correspond to distinct cell types, a key goal of scRNA-seq consists in the identification of cell subpopulations based on their transcriptional similarity. Thus, organizing cells into groups (i.e. clusters) can allow for de novo detection of cell types or identification of different subpopulations in a single cell state.
Although cells may continuously transition from one state to the next with no discrete boundary between them, Monocle does not assume that all cells in the dataset descend from a common transcriptional "ancestor". In many experiments, there might in fact be multiple distinct trajectories. Therefore, Monocle assignes each cell not only to a cluster but also to a partition. When you are learning trajectories, each partition will eventually become a separate trajectory.

```R
# Identify clusters
cds = cluster_cells(cds, cluster_method="louvain") #cell clustering using louvain algorithm
colData(cds)$clusters_louvain <-(clusters(cds)) #add column to metadata with clusters information

plot_cells(cds,
           color_cells_by = "clusters_louvain", #change here to color cells by the metadata of choice
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
           
plot_cells(cds,
           color_cells_by = "partition", #change here to color cells by the metadata of choice
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
```

## 5. Trajectory analysis: Learn the principal graph
Once cells are paritioned, each supergroup can be organized into a separate trajectory. The default method for doing this in Monocle 3 is SimplePPT, which assumes that each trajectory is a tree (albeit one that may have multiple roots). Finally, each cell will be assigned with a pseudotime value. In order to do so, you need to specify the root nodes of the trajectory graph.

```R
# Learn the trajectories
cds <- learn_graph(cds) #trajectory inference using SimplePPT

# Helper function to identify the root principal points:
get_correct_root_state <- function(cds,cell_phenotype,cell_type){
  cell_ids <- which(colData(cds)[, cell_phenotype] == cell_type)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

# Assign pseudotime value
node_ids <- get_correct_root_state(cds,cell_phenotype = "clusters_louvain", cell_type = 1) #where cell_phenotype indicates the column where cell types are stored, while cell_type is the cell type/state we want to select. Try with the cluster of your choice6
cds <- order_cells(cds, root_pr_nodes = node_ids)
```


## 6. Plot pseudotime
```R
plot_cells(cds,
           color_cells_by = "pseudotime", #change here to color cells by the metadata of choice
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
```
