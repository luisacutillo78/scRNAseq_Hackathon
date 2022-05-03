# TO DO

create a pipeline (from preprocessing to trajectories) in monocole 3 (.Rmd) (Francesco)

create a pipeline (from preprocessing to trajectories) in Seurat (.Rmd) (Annamaria)

# Data Availability
Please download:
1. [emptyDroplets_doublet_filtered_tabulamuris_mtx.rds](https://drive.google.com/file/d/1TLYPcawqbtqApGDTXUZuwE9BJDpeUbPX/view?usp=sharing)
2. [Tabulamuris_cmd.txt](https://drive.google.com/file/d/1ngJ45fOzgY6pm9gNZdCKjXdG7lmJmY03/view?usp=sharing)
3. [Tabulamuris_genes.txt](https://drive.google.com/file/d/1L3IGc59iVLwT2HE7sGBu3R_DhzjxT7oF/view?usp=sharing)


## Sofware perequisites
Required packages and software installation list:


### R software: 
* Latest R version (4.1.3 onwards) https://www.r-project.org/
* Rstudio (optional) https://www.rstudio.com/

### R packages:

Open R\Rstudio and run the code highlighted

* Seurat
 ```
install.pakages("Seurat")
devtools::install_github('satijalab/seurat-data')
remotes::install_github('satijalab/seurat-wrappers')
 ```
* DropletUtils
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DropletUtils")
```
* DoubletFinder
```
install.packages("remotes")
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
```
* Monocle3
```
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(version = "3.10")
# Next, install a few Bioconductor dependencies that aren't automatically installed:

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
# Now, install monocle3 through the cole-trapnell-lab GitHub, execute:

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
# If you wish to install the develop branch of monocle3, execute:

devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")
```
# scRNAseq_Hackathon

https://conferences.leeds.ac.uk/bad-hackathon/programme/

****
## Monocole 3 Pipelines
http://cole-trapnell-lab.github.io/monocle-release/monocle3/#tutorial-1-learning-trajectories-with-monocle-3
https://cole-trapnell-lab.github.io/monocle3/docs/installation/

http://cole-trapnell-lab.github.io/monocle-release/monocle3/#step-3-partition-the-cells-into-supergroups


## Seurat Pipelines
https://satijalab.org/seurat/index.html

https://github.com/gambalab/scRNAseq_chapter/blob/master/pipelines/seurat_pipe.md

http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html

****
## Datasets

https://drive.google.com/drive/folders/1eP16g-thyNZ-wTZwu24IPKB5uRWWYFIn?usp=sharing


PREPROCESSING ONLY

https://www.nature.com/articles/s41586-018-0590-4

https://tabula-muris.ds.czbiohub.org/

PREPROCESSING AND TRAJECTORIES

https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/#pre-process

http://bioconductor.org/books/3.14/OSCA.advanced/trajectory-analysis.html (Main example dataset to be converted into a text file to be loaded in our pipelines)

http://bioconductor.org/books/3.14/OSCA.workflows/bach-mouse-mammary-gland-10x-genomics.html
