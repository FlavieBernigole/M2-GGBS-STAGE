#Chargement des librairies
library(STdeconvolve)
library(Seurat)
library(SeuratData)
library(SpotClean)
library(S4vectors)
library(ggplot2)
library(spacexr)

#Chargement des données spatial
Test<-Load10X_Spatial(data.dir = "~/EP2C/10x_6573G/10x-D24-270-outs/")
Test<- SCTransform(Test, assay = "Spatial", verbose = FALSE)
Test <- RunPCA(Test, assay = "SCT", verbose = FALSE)
Test <- FindNeighbors(Test, reduction = "pca", dims = 1:30)
Test <- FindClusters(Test, verbose = FALSE)
Test <- RunUMAP(Test, reduction = "pca", dims = 1:30)
Test <- RunTSNE(Test, dims = 1:10)
Test <- SCTransform(Test, assay = "Spatial", verbose = TRUE) %>% RunPCA(verbose = FALSE)

"Mise en place des variables pour la déconvolution
cd <- Test[["Spatial"]]$counts
pos <- GetTissueCoordinates(Test)
colnames(pos) <- c("x", "y")
pos[is.na(colnames(pos))] <- NULL
annot<-Test@meta.data$first_type
#query <- SpatialRNA(coords, counts, colSums(counts))



counts <- cleanCounts(counts = cd,
                      min.lib.size = 10,
                      min.reads = 1,
                      min.detected = 1,
                      verbose = TRUE)


corpus <- restrictCorpus(counts,
                         removeAbove = 1.0,
                         removeBelow = 0.05,
                         alpha = 0.05,
                         plot = TRUE,
                         verbose = TRUE)

ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, 8, by = 1),
               perc.rare.thresh = 0.05,
               plot=TRUE,
               verbose=TRUE)


optLDA <- optimalModel(models = ldas, opt = "min")
results <- getBetaTheta(optLDA,
                        perc.filt = 0.05,
                        betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta

vizAllTopics(deconProp, pos, 
             groups = annot, 
             group_cols = rainbow(length(levels(annot))),
             r=0.4)

