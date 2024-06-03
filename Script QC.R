# chargement des librairies
library(ggplot2)
library(ggpubr)
library(Seurat)
library(scater)
library(dplyr)
library(gridExtra)
library(grid)
library(scDblFinder)
library(SingleR)
library(celldex)

#Chargement des data SeqWell de l'ensemble des échantillons
Involution <- Read10X(paste("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/ObjetSeqwellRound3.tier1"))

#Création objet Seurat + Traitement des datas
seuratObjet <- V3
seuratObjet <- CreateSeuratObject(counts = seuratObjet, min.cells = 3, project = "VPLI")
seuratObjet <- NormalizeData(seuratObjet, assay = "RNA")
seuratObjet <- FindVariableFeatures(seuratObjet, nfeatures = 2000)
seuratObjet <- ScaleData(seuratObjet)
seuratObjet <- RunPCA(seuratObjet)
seuratObjet <- RunUMAP(seuratObjet, dims = 1:30)
seuratObjet <- FindNeighbors(seuratObjet, dims = 1:30)
seuratObjet <- FindClusters(seuratObjet, resolution = c(0,seq(0.1, 1, by=0.1)))

#Nombre de reads
seuratObjet$QC_ncount <- seuratObjet$nCount_RNA < 50000
p1 <- (ggplot(as.data.frame(seuratObjet$nCount_RNA), aes(x=seuratObjet$nCount_RNA)) + geom_density(color="darkblue", fill="lightblue") +   labs(x="nCount_RNA", y = "Density") + geom_vline(xintercept = 50000))  + ggtitle("E8-M20221021-INV")  + (FeaturePlot(seuratObjet, features = "nCount_RNA")) + DimPlot(seuratObjet, group.by = "QC_ncount")

#Nombre de gènes
seuratObjet$QC_nfeature <- seuratObjet$nFeature_RNA < 8000
p2 <- (ggplot(as.data.frame(seuratObjet$nFeature_RNA), aes(x=seuratObjet$nFeature_RNA)) + geom_density(color="plum4", fill="thistle1") +   labs(x="nFeature_RNA", y = "Density") + geom_vline(xintercept =8000)) +FeaturePlot(seuratObjet, features = "nFeature_RNA") + DimPlot(seuratObjet, group.by = "QC_nfeature")

#% Ribosomaux
seuratObjet <- PercentageFeatureSet(seuratObjet, pattern = "^Rp[sl]",col.name = "percent.rp")
OUtlierPerrp <- isOutlier(seuratObjet$percent.rp, nmads=3, type="both", log=FALSE)
seuratObjet$QC_percent.rp <- seuratObjet$percent.rp < max(attributes(OUtlierPerrp)$thresholds[2],7)
p6 <- (ggplot(as.data.frame(seuratObjet$percent.rp), aes(x=seuratObjet$percent.rp)) + geom_density(color="grey27", fill="grey50") +   labs(x="percent.rp", y = "Density") + geom_vline(xintercept = max(attributes(OUtlierPerrp)$thresholds[2],7))) + DimPlot(seuratObjet, group.by = "QC_percent.rp") + FeaturePlot(seuratObjet, features = "percent.rp", max.cutoff = 25)

#% Mitochondriales
seuratObjet <- PercentageFeatureSet(seuratObjet, pattern = "mt-", col.name = "percent.mt")
OUtlierPerMT <- isOutlier(seuratObjet$percent.mt, nmads=3, type="both", log=FALSE)
seuratObjet$QC_percent.mt <- seuratObjet$percent.mt < max(attributes(OUtlierPerMT)$thresholds[2],7)
p3 <- (ggplot(as.data.frame(seuratObjet$percent.mt), aes(x=seuratObjet$percent.mt)) + geom_density(color="grey27", fill="grey50") +   labs(x="percent.mt", y = "Density") + geom_vline(xintercept = max(attributes(OUtlierPerMT)$thresholds[2],7))) + DimPlot(seuratObjet, group.by = "QC_percent.mt") + FeaturePlot(seuratObjet, features = "percent.mt", max.cutoff = 25)
p3

#Identification markers
Idents(seuratObjet) <- "seurat_clusters"
markers <- FindAllMarkers(seuratObjet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
p8<-((DimPlot(seuratObjet, group.by = "RNA_snn_res.1") + labs(title="",subtitle="Clustering res.1") + NoLegend()) + (DoHeatmap(seuratObjet, features = top5$gene, size = 3)+ theme(text = element_text(size = 5), axis.text.x.top = element_text(size = 6))))
print(p8)

#Recapitulatif count, feature et mt
seuratObjet$to_discard <- seuratObjet$QC_ncount & seuratObjet$QC_nfeature & seuratObjet$QC_percent.mt 
p7<-DimPlot(seuratObjet, group.by = "to_discard") + ggtitle(paste("Percentage of cells failing QC : ",round(table(seuratObjet$to_discard)[1]/sum(table(seuratObjet$to_discard))*100,2),"%", sep=""))
print(p7)



