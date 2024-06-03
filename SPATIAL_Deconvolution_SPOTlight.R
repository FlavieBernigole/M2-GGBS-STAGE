#Chargement des libraries
library(SPOTlight)
library(Seurat)
library(ggplot2)
library(NMF)
library(scran) #pour modelGeneVar
library(tidyverse) #pour pipe %>%

#Chargement des data SeqWEll
ref <- readRDS("~/EP2C/ObjetSeqwell/SeqWellAnnotated_Bcells.rds")
ref <- NormalizeData(ref, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

# Gène hautement variable
dec <- modelGeneVar(ref@assays$SCT$counts)
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)

# Selection de 3000 gènes
hvg <- getTopHVGs(dec, n = 3000)
ref_single_cell = ref@assays$SCT$counts
colnames(ref_single_cell) = ref@meta.data$Annotation
# Gène ribosomaux ou mitochondriales
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(ref_single_cell))

# Calculs des gènes marqueurs
mgs <- scoreMarkers(ref_single_cell, groups = ref@meta.data$Annotation, subset.row = genes)

#Selection des plus pertinents
mgs_fil <- lapply(names(mgs), function(i) {  #applique fonction a chaque element de la liste mgs, en prenant i dans mgs, c est des noms de clusters
  x <- mgs[[i]] #mgs est une liste de listes contenant un tableau en gros... et la on met x = le tableau
  # Diltrer et garder les gènes AUC > 0.8
  x <- x[x$mean.logFC.detected > 2, ]
  # Trier les gènes en fonction de leur poids
  x <- x[order(x$mean.logFC.detected, decreasing = TRUE), ]
  
  if (nrow(x) > 0 ) { #le dataframe n est pas vide, on peut y ajouter des colones
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
    
  } else {
    data.frame()
  }
})
mgs_df <- do.call(rbind, mgs_fil) # combiner tous les dataframes de la liste mgs_fil




#genere une sequence d entiers de 1 a ncol, un nombre represente une colonne
#regroupe les indices des colonnes correspondant a un meme type cellulaire
idx <- split(seq(ncol(ref_single_cell)), ref@meta.data$Annotation) #genere une sequence d entiers de 1 a ncol(sce)


# selection 100 markers par cellules
n_cells <- 100
cs_keep <- lapply(idx, function(i) {
  n <- length(i) #nombre de cellules de ce type cellulaire
  if (n < n_cells) #si moins de 100 cellules de ce type cellulaire
    n_cells <- n #on les prendra toutes
  sample(i, n_cells) #choisi n_cells cellules de i (au final : de chaque type cellulaire)
})
ref_single_cell <- ref_single_cell[, unlist(cs_keep)]  # toutes les lignes (genes), mais que les colonnes correspondant aux cellules precedemment choisies


# Application de la Déconvolution

res <- SPOTlight(
  x = ref_single_cell,
  y = objet_seurat4,
  groups = colnames(ref_single_cell),
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC", 
  group_id = "cluster",
  gene_id = "gene")

head(mat <- res$mat)[, seq_len(3)]
#NMF model fit
mod <- res$NMF




##Visualisation
#Matrice de corrélation spatiale
plot1<- plotCorrelationMatrix(mat, method = "circle")
plot1
plot2<- plotCorrelationMatrix(mat)
plot2

#Co_localisation
plot5<-plotInteractions(mat,"heatmap")
plot5
plot4<-plotInteractions(mat, "heatmap",hc.order = TRUE, type = "lower",outline.col = "white")
plot4

#Gènes sélectionnés 
library(NMF)
sign <- basis(mod)
colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
head(sign)


#Visualisation sur coupe
ct <- colnames(mat)
mat[mat < 0.1] <- 0

# Couleurs
library(RColorBrewer)
paletteMartin <- brewer.pal(12, "Paired")
paletteMartine <- c( "#000000", "#FFFFFF", "#FFFFFF","#FFFFFF","#FFFFFF", "#D3D3D3", "#D3D3D3","#FFFFFF","#FFFFFF","#D3D3D3","#FFFFFF","#FFFFFF")
pal <- colorRampPalette(paletteMartine)(length(ct))
names(pal) <- ct

#Visualisation coupe
plotSpatialScatterpie(
  x = objet_seurat4,
  y = mat,
  cell_types = colnames(mat),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.4) +
  scale_fill_manual(
    values = pal,
    breaks = names(pal))



