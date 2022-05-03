#Download Tabula Muris raw data

suppressPackageStartupMessages({
  library(ExperimentHub)
  library(SingleCellExperiment)
  library(TabulaMurisData)
})
M=readRDS("C:/Users/f.panariello/Downloads/Tabula.Muris.Raw.10x.Cells.rds")
eh <- ExperimentHub()
query(eh, "TabulaMurisData")
droplet <- eh[["EH1617"]]
M = assay(droplet)
M = M[,Matrix::colSums(M)>1000]
M = M[,Matrix::colSums(M!=0)>500]

library(DropletUtils)
set.seed(100)
e.out <- emptyDrops(M)
e.out$is.cell <- e.out$FDR <= 0.01

saveRDS(as.data.frame(e.out), "empty_drop_results.rds")


library(ggplot2)
pdf("empty_drops_figure.pdf", height=8, width = 12)
ggplot(as.data.frame(e.out), aes(x=Total, y=-LogProb, fill=is.cell))+
  geom_point(pch=21)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  xlab("Total UMI count")+
  ylab("-Log Probability")+
  scale_fill_manual(values=c("white", "red"))
ggplot(as.data.frame(e.out), aes(x=Total, y=-log10(FDR+10^-5), fill=is.cell))+
  geom_point(pch=21)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  xlab("Total UMI count")+
  ylab("-Log10 (FDR+10^-5)")+
  scale_fill_manual(values=c("white", "red"))+
  coord_cartesian(ylim=c(0,5))
dev.off()


e.out
sum(is.cell, na.rm=TRUE)
M=M[,is.cell]
saveRDS(M,"emptyDroplets_filtered_tabulamuris_mtx.rds")

#Find doublets - Doublet Finder
library(DoubletFinder)
library(Seurat)

## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
seu <- CreateSeuratObject(M)
seu <- NormalizeData(seu)
seu <- ScaleData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:10)

## Pre-process Seurat object (sctransform) -----------------------------------------------------------------------------------
seu <- CreateSeuratObject(kidney.data)
seu <- SCTransform(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:10)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(seu, PCs = 1:10, sct = FALSE)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(seu, PCs = 1:10, sct = FALSE)
gt.calls <- seu@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"]
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu@meta.data$ClusteringResults
nExp_poi <- round(0.075*length(seu@cell.names))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
