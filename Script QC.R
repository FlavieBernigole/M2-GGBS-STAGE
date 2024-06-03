
#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)
#args <- "E8-M2-20221021-"
# test if there is at least one argument: if not, return an error
#if (length(args)==0) {
#  stop("At least one argument must be supplied (id sample)", call.=FALSE)
#}

# chargement des librairies
library(ggplot2)
library(ggpubr)
library(Seurat)
library(scater)
#BiocManager::install("scater")
library(dplyr)
#library(nichenetr)
library(gridExtra)
library(grid)
#BiocManager::install("scDblFinder")
library(scDblFinder)
#"BiocManager::install("SingleR")
library(SingleR)
#BiocManager::install("celldex")
library(celldex)


Involution <- Read10X(paste("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/scRNAseq_15062023/0_STARsolo_Alignment/Involution/Solo.out/GeneFull/filtered"))
Virgin<- Read10X("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/scRNAseq_15062023/0_STARsolo_Alignment/Virgin/Solo.out/GeneFull/filtered")
Lactation<- Read10X("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/scRNAseq_15062023/0_STARsolo_Alignment/Lactation/Solo.out/GeneFull/filtered")
Pregnant<- Read10X("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/scRNAseq_15062023/0_STARsolo_Alignment/Pregnant/Solo.out/GeneFull/filtered")
VirginKI<- Read10X("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/scRNAseq_15062023/0_STARsolo_Alignment/Virgin-Ki/Solo.out/GeneFull/filtered")
P1<- Read10X("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/Seqwell_Jen_29022024/0_STARsolo_Alignment/P1_S3/Solo.out/GeneFull/filtered")
P2<- Read10X("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/Seqwell_Jen_29022024/0_STARsolo_Alignment/P2_S4/Solo.out/GeneFull/filtered")
N1<- Read10X("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/Seqwell_Jen_29022024/0_STARsolo_Alignment/N1_S5/Solo.out/GeneFull/filtered")
N2<- Read10X("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/Seqwell_Jen_29022024/0_STARsolo_Alignment/N2_S6/Solo.out/GeneFull/filtered")
V1<- Read10X("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/Seqwell_Jen_29022024/0_STARsolo_Alignment/V1_S1/Solo.out/GeneFull/filtered")
V2<- Read10X("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/Seqwell_Jen_29022024/0_STARsolo_Alignment/V2_S2/Solo.out/GeneFull/filtered")
I1<- Read10X("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/scRNAseq_26042024/I1_S3/Solo.out/GeneFull/filtered")
I1b<- Read10X("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/scRNAseq_26042024/I1bis_S4/Solo.out/GeneFull/filtered")
I2<- Read10X("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/scRNAseq_26042024/I2_S5/Solo.out/GeneFull/filtered")
N3<- Read10X("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/scRNAseq_26042024/N3_S2/Solo.out/GeneFull/filtered")
V3<- Read10X("/LAB-DATA/BiRD/shares/CRCINA/E08_Projet_SingleCell/Seqwell/scRNAseq_26042024/V3_S1/Solo.out/GeneFull/filtered")


seuratObjet <- readRDS("~/EP2C/ObjetSeqwell/CellChat2/objet_seurat.rds")
seuratObjet <- CreateSeuratObject(counts = seuratObjet, min.cells = 3, project = "VPLI")
seuratObjet <- NormalizeData(seuratObjet, assay = "RNA")
seuratObjet <- FindVariableFeatures(seuratObjet, nfeatures = 2000)
seuratObjet <- ScaleData(seuratObjet)
seuratObjet <- RunPCA(seuratObjet)
seuratObjet <- RunUMAP(seuratObjet, dims = 1:14)
seuratObjet <- FindNeighbors(seuratObjet, dims = 1:14)
seuratObjet <- FindClusters(seuratObjet, resolution = c(0,seq(0.1, 1, by=0.1)))

seuratObjet$QC_ncount <- seuratObjet$nCount_RNA < 50000
p1 <- (ggplot(as.data.frame(seuratObjet$nCount_RNA), aes(x=seuratObjet$nCount_RNA)) + geom_density(color="darkblue", fill="lightblue") +   labs(x="nCount_RNA", y = "Density") + geom_vline(xintercept = 50000))  + ggtitle("E8-M20221021-INV")  + (FeaturePlot(seuratObjet, features = "nCount_RNA")) + DimPlot(seuratObjet, group.by = "QC_ncount")
p1


seuratObjet$QC_nfeature <- seuratObjet$nFeature_RNA < 8000
p2 <- (ggplot(as.data.frame(seuratObjet$nFeature_RNA), aes(x=seuratObjet$nFeature_RNA)) + geom_density(color="plum4", fill="thistle1") +   labs(x="nFeature_RNA", y = "Density") + geom_vline(xintercept =8000)) +FeaturePlot(seuratObjet, features = "nFeature_RNA") + DimPlot(seuratObjet, group.by = "QC_nfeature")
p2

seuratObjet <- PercentageFeatureSet(seuratObjet, pattern = "^Rp[sl]",col.name = "percent.rp")
OUtlierPerrp <- isOutlier(seuratObjet$percent.rp, nmads=3, type="both", log=FALSE)
seuratObjet$QC_percent.rp <- seuratObjet$percent.rp < max(attributes(OUtlierPerrp)$thresholds[2],7)
p6 <- (ggplot(as.data.frame(seuratObjet$percent.rp), aes(x=seuratObjet$percent.rp)) + geom_density(color="grey27", fill="grey50") +   labs(x="percent.rp", y = "Density") + geom_vline(xintercept = max(attributes(OUtlierPerrp)$thresholds[2],7))) + DimPlot(seuratObjet, group.by = "QC_percent.rp") + FeaturePlot(seuratObjet, features = "percent.rp", max.cutoff = 25)
p6

seuratObjet <- PercentageFeatureSet(seuratObjet, pattern = "mt-", col.name = "percent.mt")
OUtlierPerMT <- isOutlier(seuratObjet$percent.mt, nmads=3, type="both", log=FALSE)
seuratObjet$QC_percent.mt <- seuratObjet$percent.mt < max(attributes(OUtlierPerMT)$thresholds[2],7)
p3 <- (ggplot(as.data.frame(seuratObjet$percent.mt), aes(x=seuratObjet$percent.mt)) + geom_density(color="grey27", fill="grey50") +   labs(x="percent.mt", y = "Density") + geom_vline(xintercept = max(attributes(OUtlierPerMT)$thresholds[2],7))) + DimPlot(seuratObjet, group.by = "QC_percent.mt") + FeaturePlot(seuratObjet, features = "percent.mt", max.cutoff = 25)
p3


### doublets detection
# scDblFinder
seuratObjet<- readRDS(file="~/EP2C/ObjetSeqwell/Annotation_v3/VPLI_152926_STARsolo_Round3Tier1_optimized.rds")
sce <- as.SingleCellExperiment(DietSeurat(seuratObjet))
library(scDblFinder)
sce <- scDblFinder(sce, samples = sce$Sequencing.Date)
seuratObjet$scDblFinder_score <- sce$scDblFinder.score
seuratObjet$scDblFinder_class <- sce$scDblFinder.class
d <- scDblFinder(sce, verbose=FALSE, returnType="table")

th <- doubletThresholding(d)


OutlierscDblFinder <- th
seuratObjet$QC_scDblFinder <- seuratObjet$scDblFinder_score < th
p4 <- (ggplot(as.data.frame(seuratObjet$scDblFinder_score), aes(x = seuratObjet$scDblFinder_score))+geom_density(color="#92A9BD", fill="#92A9BD") +   labs(x="scDblFinder score", y = "Density") + geom_vline(xintercept = th))+  DimPlot(seuratObjet, group.by = "QC_scDblFinder")
p4

### annotation avec singleR
library(SingleR)
library(celldex)
ImmGen <- ImmGenData()
sce <- as.SingleCellExperiment(seuratObjet)
annot <- SingleR(sce, ref = ImmGen, labels = ImmGen$label.main)

seuratObjet <- AddMetaData(seuratObjet, metadata = annot$pruned.labels, col.name = "pruned.labels") 
p5 <- DimPlot(seuratObjet, group.by = "pruned.labels") + labs(title="")
p5
#saveRDS(seuratObjet, paste("~/EP2C/Data_Spatial/Vierge1.rds", sep=""))

#dev.off()



#pdf(paste(("~/EP2C/Data10X"), "-with-scDblFinder.pdf", sep=""), width=12, height = 6)
#print(p1)
#print(p2)
#print(p3)
#print(p4)
#print(p5)

Idents(seuratObjet) <- "seurat_clusters"
markers <- FindAllMarkers(seuratObjet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
p8<-((DimPlot(seuratObjet, group.by = "RNA_snn_res.1") + labs(title="",subtitle="Clustering res.1") + NoLegend()) + (DoHeatmap(seuratObjet, features = top5$gene, size = 3)+ theme(text = element_text(size = 5), axis.text.x.top = element_text(size = 6))))
print(p8)


seuratObjet$to_discard <- seuratObjet$QC_ncount & seuratObjet$QC_nfeature & seuratObjet$QC_percent.mt 
#& seuratObjet$QC_scDblFinder
p7<-DimPlot(seuratObjet, group.by = "to_discard") + ggtitle(paste("Percentage of cells failing QC : ",round(table(seuratObjet$to_discard)[1]/sum(table(seuratObjet$to_discard))*100,2),"%", sep=""))
print(p7)


df <- data.frame()
tablediscard <- table(seuratObjet$pruned.labels, seuratObjet$to_discard)
Percentage <- round(tablediscard[,2] / (tablediscard[,1]+tablediscard[,2])*100,1)
finalTable <- cbind(t(t(table(seuratObjet$pruned.labels))), table(seuratObjet$pruned.labels,seuratObjet$QC_ncount)[,1],table(seuratObjet$pruned.labels,seuratObjet$QC_percent.mt)[,1],table(seuratObjet$pruned.labels,seuratObjet$to_discard)[,2], Percentage)
#write.csv(finalTable, file="D:/Seqwell_Jen_29022024/finalTableV1.csv")
finalTable
colnames(finalTable) <- c("Total_cells", "QC_nCount","QC_nFeature","QC_percent.mt", "Percent_keep_cells")
print(finalTable)

p11 <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)+ theme_void() + annotation_custom(tableGrob(cbind(tablediscard,Percentage)), xmin=0, xmax=10, ymin=-0, ymax=100)
p11

##### filter cells
seuratObjet <- subset(seuratObjet, subset = to_discard)
##### reprocess with good cells

seuratObjet <- FindVariableFeatures(seuratObjet, nfeatures = 2000)
seuratObjet <- ScaleData(seuratObjet)
seuratObjet <- RunPCA(seuratObjet)
seuratObjet <- RunUMAP(seuratObjet, dims = 1:30)
seuratObjet <- FindNeighbors(seuratObjet, dims = 1:30)
seuratObjet <- FindClusters(seuratObjet, resolution = c(0,seq(0.1, 1, by=0.1)))

Idents(seuratObjet) <- "round3.tier1"
markers <- FindAllMarkers(seuratObjet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.xlsx(markers, file = "D:/STAGEM2/markers.xlsx")

top5 <- markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
p9<-((DimPlot(seuratObjet, group.by = "RNA_snn_res.1") + labs(title="",subtitle=paste("Clustering after QC - res.1 - n=",sum(table(seuratObjet$to_discard)),sep="")) + NoLegend()) + (DoHeatmap(seuratObjet, features = top5$gene, size = 3) + theme(text = element_text(size = 5), axis.text.x.top = element_text(size = 6))))
p9



p10<-((DimPlot(seuratObjet, group.by = "pruned.labels") + labs(title="",subtitle=paste("Clustering after QC - n=",sum(table(seuratObjet$to_discard)),sep=""))))
p11<-((DimPlot(seuratObjet, group.by = "dvpmt.stage") + labs(title="",subtitle=paste("Clustering after QC - n=",sum(table(seuratObjet$to_discard)),sep=""))))
p12<-((DimPlot(seuratObjet, group.by = "RNA_snn_res.0.1") + labs(title="",subtitle=paste("Clustering after QC - res.0.1 - n=",sum(table(seuratObjet$to_discard)),sep=""))))
p13<-((DimPlot(seuratObjet, group.by = "RNA_snn_res.0.2") + labs(title="",subtitle=paste("Clustering after QC - res.0.2 - n=",sum(table(seuratObjet$to_discard)),sep=""))))
p14<-((DimPlot(seuratObjet, group.by = "RNA_snn_res.0.3", label = TRUE) + labs(title="",subtitle=paste("Clustering after QC - res.0.3 - n=",sum(table(seuratObjet$to_discard)),sep=""))))
p15<-((DimPlot(seuratObjet, group.by = "RNA_snn_res.0.4") + labs(title="",subtitle=paste("Clustering after QC - res.0.4 - n=",sum(table(seuratObjet$to_discard)),sep=""))))
p16<-((DimPlot(seuratObjet, group.by = "RNA_snn_res.0.5") + labs(title="",subtitle=paste("Clustering after QC - res.0.5 - n=",sum(table(seuratObjet$to_discard)),sep=""))))
p17<-((DimPlot(seuratObjet, group.by = "RNA_snn_res.0.6") + labs(title="",subtitle=paste("Clustering after QC - res.0.6 - n=",sum(table(seuratObjet$to_discard)),sep=""))))
p18<-((DimPlot(seuratObjet, group.by = "RNA_snn_res.0.7") + labs(title="",subtitle=paste("Clustering after QC - res.0.7 - n=",sum(table(seuratObjet$to_discard)),sep=""))))

pdf(file = "D:/UMAP_resolution_152926_cutoff_mt.pdf")
p10
p11
p12
p13
p14
p15
p16
p17
p18
dev.off()


#saveRDS(Seurat, paste("D:/merged15_29_26_sanscutoff.rds", sep=""))

#dev.off()


