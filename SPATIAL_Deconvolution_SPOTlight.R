library(SPOTlight)
library(Seurat)
library(ggplot2)
library(NMF)
library(scran) #pour modelGeneVar
library(tidyverse) #pour pipe %>%



##### loading data #####
ref <- readRDS("~/EP2C/ObjetSeqwell/SeqWellAnnotated_Bcells.rds")
ref <- NormalizeData(ref, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)


# highly variable genes du tuto
dec <- modelGeneVar(ref@assays$SCT$counts)

#plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
#curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
# Ajouter le vecteur transformé à l'objet Seurat

# Get the top 3000 genes.
hvg <- getTopHVGs(dec, n = 3000)

ref_single_cell = ref@assays$SCT$counts
colnames(ref_single_cell) = ref@meta.data$Annotation

# Get vector indicating which genes are neither ribosomal or mitochondrial
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(ref_single_cell))

# Compute marker genes
mgs <- scoreMarkers(ref_single_cell, groups = ref@meta.data$Annotation, subset.row = genes)



#garder que ceux pertinents


mgs_fil <- lapply(names(mgs), function(i) {  #applique fonction a chaque element de la liste mgs, en prenant i dans mgs, c est des noms de clusters
  x <- mgs[[i]] #mgs est une liste de listes contenant un tableau en gros... et la on met x = le tableau
  
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.logFC.detected > 2, ]
  
  # Sort the genes from highest to lowest weight
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



##### Cell Downsampling #####

#genere une sequence d entiers de 1 a ncol, un nombre represente une colonne
#regroupe les indices des colonnes correspondant a un meme type cellulaire
idx <- split(seq(ncol(ref_single_cell)), ref@meta.data$Annotation) #genere une sequence d entiers de 1 a ncol(sce)



# downsample to at most 100 per identity & subset
n_cells <- 50
cs_keep <- lapply(idx, function(i) {
  n <- length(i) #nombre de cellules de ce type cellulaire
  if (n < n_cells) #si moins de 100 cellules de ce type cellulaire
    n_cells <- n #on les prendra toutes
  sample(i, n_cells) #choisi n_cells cellules de i (au final : de chaque type cellulaire)
})

ref_single_cell <- ref_single_cell[, unlist(cs_keep)]  # toutes les lignes (genes), mais que les colonnes correspondant aux cellules precedemment choisies



objet_seurat<- readRDS("~/Repertoire_travail_deconvolution/SpatialI.rds")
objet_seurat2<- readRDS("~/Repertoire_travail_deconvolution/SpatialL.rds")
objet_seurat3<- readRDS("~/Repertoire_travail_deconvolution/SpatialP.rds")
objet_seurat4<-readRDS("~/Repertoire_travail_deconvolution/SpatialV.rds")


##### Deconvolution #####

res <- SPOTlight(
  x = ref_single_cell,
  y = objet_seurat4,
  groups = colnames(ref_single_cell),
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC", 
  group_id = "cluster",
  gene_id = "gene")



#Extract data 
head(mat <- res$mat)[, seq_len(3)]

# Extract NMF model fit
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


#saveRDS(mod, file = "modI.rds")
#saveRDS(mat, file = "matI.rds")
#saveRDS(res, file = "resI.rds")
#saveRDS(sign, file = "signI.rds")

#Visualisation sur coupe
ct <- colnames(mat)
mat[mat < 0.1] <- 0

# Define color palette
# (here we use 'paletteMartin' from the 'colorBlindness' package)
paletteMartin <- c(
  "#FF0000", "#FFFF00", "#00FFFF", "#FF00FF", "#FFA500", 
  "#800080", "#FF1493", "#32CD32", "#6A5ACD", 
  "#000000", "#924900", "#DB6D00", "#24FF24", "#D3D3D3"
)


library(RColorBrewer)

# Obtenez une palette de 14 couleurs distinctes
paletteMartin <- brewer.pal(12, "Paired")
print(paletteMartin)
[1] "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6"
[10] "#6A3D9A" "#FFFF99" "#B15928"

paletteMartine <- c( "#000000", "#FFFFFF", "#FFFFFF","#FFFFFF","#FFFFFF", "#D3D3D3", "#D3D3D3","#FFFFFF","#FFFFFF","#D3D3D3","#FFFFFF","#FFFFFF")



pal <- colorRampPalette(paletteMartine)(length(ct))
names(pal) <- ct

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



