# Pipelines
During the day we will follow 2 pipelines:

* [Seurat Pipeline (Rmd format)](https://github.com/luisacutillo78/scRNAseq_Hackathon/blob/main/Analysis_TabulaMuris_AC_2.Rmd)
* [Monocle 3 Pipeline (md format)](https://github.com/luisacutillo78/scRNAseq_Hackathon/blob/main/Pipeline2_monocle3.md)



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

* Devools, Remotes R build tools
 ```
install.packages("devtools")
install.packages("remotes")
install.packages("pkgbuild")
pkgbuild::check_build_tools()
```
* Seurat
```
install.packages("Seurat")
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
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
```
* Monocle3 installation on Windows OS
```
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
# Remove this 
# BiocManager::install(version = "3.10")

# Next, install a few Bioconductor dependencies that aren't automatically installed:

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
# Now, install monocle3 through the cole-trapnell-lab GitHub, execute:

devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
# If you wish to install the develop branch of monocle3, execute:

devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")
```
* Monocle3 installation on macOS Monterey Chip M1

R version 4.2.0 (2022-04-22)

Platform: aarch64-apple-darwin20 (64-bit)

Running under: macOS Monterey 12.3.1

1 - Install Xcode from the Appstore (it requires some time, almost one hour with good connection)
2 - After Xcode is installed, accept the license and the installation of additional components
3 - Open the terminal and run the command 'xcode-select --install' (almost 15/20 minutes)

Ignore "xcode-select: error: command line tools are already installed, use "Software Update" to install updates". It means that it is already installed.

4 - Download gfortran from https://github.com/fxcoudert/gfortran-for-macOS/releases/tag/11-arm-alpha2 (M1 processor) or from https://github.com/fxcoudert/gfortran-for-macOS/releases (Intel processor). The user-friendly installation package is gfortran-ARM-11.0-BigSur.pkg
5 - Click on the .pkg file to install gfortan in /usr/local/bin
6 - In the home directory, create a hidden directory called .R and create inside a file called 'Makevars'. In this file, write:
```
F77 = /usr/local/bin/gfortran
FC = $(F77)
FLIBS = -L/usr/local/gfortran/lib/gcc/aarch64-apple-darwin20.2.0
```
Check for the correspondence with your path, it might change from user to user.

7 - Proceed with Monocle3 installation and select 'no' at the following question:

Do you want to install from sources the packages which need compilation? (Yes/no/cancel)



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
