#Chargement des librairies
library(Seurat)
library(ggplot2)
library(future)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(SeuratData)
library(dplyr)
library(patchwork)
library(hdf5r)
library(scatterpie)
library(ggcorrplot)


#chargement de l'objet Seurat et des données de Spatial
path<-"~/EP2C/10x_6573G/10x-D24-270-outs/"
spatial.obj<-Load10X_Spatial(data.dir=path)
spatial.obj <- subset(spatial.obj, subset = nCount_Spatial > 0)
p1V<-VlnPlot(spatial.obj, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
spatial.obj$log_nCount_Spatial <- log(spatial.obj$nCount_Spatial)
p2V <- SpatialFeaturePlot(spatial.obj, features = "log_nCount_Spatial") + theme(legend.position = "right")
wrap_plots(p1V, p2V)

#Visualisation des données de spatial
p3V <- DimPlot(spatial.obj, reduction = "umap", label = TRUE)
p4V <- SpatialDimPlot(spatial.obj, stroke = 0)
p4V + p3V

#Visualisation certains clusters sur la coupe
SpatialDimPlot(spatial.obj, cells.highlight = CellsByIdentities(object = spatial.obj, idents = c(2,8,7,4)), facet.highlight = TRUE)

#Intégration de scRNAseq référence
ref <- readRDS("~/Repertoire_travail_deconvolution/scRNAseqSCTvf5.rds")
ref <- UpdateSeuratObject(ref)

#Transfert d'étiquettes pour prédire le type de cellule principal pour chaque spots
transfert <- FindTransferAnchors(reference = ref, query = spatial.obj, normalization.method = "SCT",npcs = 50)
predictions.assay <- TransferData(anchorset = transfert, refdata = ref$Annotationvf5, prediction.assay = TRUE, weight.reduction = spatial.obj[["pca"]], dims = 1:50)
spatial.obj[["predictions"]] <- predictions.assay

#Visulisation1
DefaultAssay(spatial.obj) <- "predictions"
p5V<-SpatialFeaturePlot(spatial.obj, features = c("Basalcells","Luminalcells","Bcells","Macrophages","Fibroblasts","Tcells"), alpha = c(0.1, 1))

#Visualisation2
spatial.obj$predicted.id <- GetTransferPredictions(spatial.obj)
Idents(spatial.obj) <- "predicted.id"
SpatialDimPlot(spatial.obj, cells.highlight = CellsByIdentities(object = spatial.obj, idents = c("Fibroblasts","Luminalcells","Macrophages")), facet.highlight = TRUE)

#Identification d'entités spatialement variables
DefaultAssay(spatial.obj) <- "SCT"
spatial.obj <- FindSpatiallyVariableFeatures(spatial.obj, assay = "SCT", slot = "scale.data", features = VariableFeatures(spatial.obj)[1:1000],selection.method = "moransi", x.cuts = 100, y.cuts = 100)

#Déconvolution spatiale à l'aide d'un RCTD
ref <- readRDS("~/Repertoire_travail_deconvolution/scRNAseqSCTvf5.rds")
ref <- UpdateSeuratObject(ref)
Idents(ref) <- "Annotationvf5"
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$Annotationvf5)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

