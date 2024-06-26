#Chargement des library
library(SpotClean)
library(S4vectors)
library(Seurat)
library(ggplot2)

vignette("SpotClean")

#Chargement des count matrix et des slide information de 10X Space Ranger output

spatial_image_dir<- "~/EP2C/10X_6573G/10x-D24-267-outs/spatial/"
spatial_data_dir<- "~/EP2C/10X_6573G/10x-D24-267-outs/raw_feature_bc_matrix/"

spatial_raw<-read10xRaw(spatial_data_dir)
spatial_slide_info<-read10xSlide(tissue_csv_file = file.path(spatial_image_dir,"tissue_positions.csv"),
                                 tissue_img_file = file.path(spatial_image_dir, "tissue_lowres_image.png"),
                                 scale_factor_file = file.path(spatial_image_dir, "scalefactors_json.json"))
str(spatial_slide_info)

#Création de slide_obj
slide_obj<-createSlide(spatial_raw,spatial_slide_info)
p1<-slide_obj

#Visualisation de slide_obj
p2<-visualizeSlide(slide_obj)
p3<-visualizeLabel(slide_obj, "tissue")

metadata(slide_obj)$slide$total_counts<-Matrix::colSums(spatial_raw)
p4<-visualizeHeatmap(slide_obj, "total_counts")
p5<-visualizeHeatmap(slide_obj,"Krt14")
p6<-visualizeHeatmap(slide_obj,"Krt8")
p7<-visualizeHeatmap(slide_obj,"Krt18")

#Décontamination des data
decont_obj<-spotclean(slide_obj, maxit = 10, candidate_radius = 20)
decont_obj
names(metadata(decont_obj))
p8<-visualizeHeatmap(decont_obj, "Krt14")
p9<-visualizeHeatmap(decont_obj, "Krt18")
p10<-visualizeHeatmap(decont_obj, "Krt8")

#Estimation du niveau de contamination dans les datas
p11<-summary(metadata(decont_obj)$contamination_rate)

#ARC score
p12<-arcScore(slide_obj)

#Conversion en un objet Seurat pour appliquer la déconvolution
seurat_obj_INV_267<-convertToSeurat(decont_obj,image_dir = spatial_image_dir)

#Enregistrement
pdf(file="~/EP2C/10X_6573G/SpotClean_Involution_267.pdf")
p1
p2
p3
p4
p5
p6
p7
p8
p9
p10
p11
p12
dev.off()


