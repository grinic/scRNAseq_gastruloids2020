if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("batchelor")
# BiocManager::install("Seurat")
# BiocManager::install("scater")
# BiocManager::install("umap")
# install.packages("tidyverse")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(irlba))
suppressPackageStartupMessages(library(umap))

data_location <- "server"
data_location <- "local"
# data_location <- "ext_drive"

if(data_location == "server"){
  folder.RData <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\data\\"
  outFolder <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\results\\integration\\diff_expr\\"
} else if(data_location == "local"){
  folder.RData <- "C:\\Users\\nicol\\OneDrive\\Desktop\\scrnaseq_gastruloids\\data\\"
  outFolder <- "C:\\Users\\nicol\\OneDrive\\Desktop\\scrnaseq_gastruloids\\results\\integration\\diff_expr\\"
} else if(data_location == "ext_drive"){
  folder.RData <- "F:\\scrnaseq_gastruloids\\data\\"
  outFolder <- "F:\\scrnaseq_gastruloids\\results\\integration\\diff_expr\\"
}
load(paste0(folder.RData,"anlas_data.RData"))
load(paste0(folder.RData,"pijuan_data.RData"))

pijuan <- NormalizeData(object = pijuan, normalization.method = "LogNormalize", scale.factor = 10000)
# filter pijuan data to only stages needed
# pijuan <- pijuan[,(pijuan$stage=='E6.5')|(pijuan$stage=='E6.75')|(pijuan$stage=='E7.0')|(pijuan$stage=='E7.25')|(pijuan$stage=='E7.5')]
pijuan <- pijuan[rowSums(pijuan)>10,]
pijuan <- pijuan[grep('mt-',rownames(pijuan),invert=T),]

# filter argelaguet data to only good cells and genes
anlas <- anlas[,anlas$merge.ident!='0h']
anlas <- anlas[rowSums(anlas)>10,]
anlas <- anlas[grep('mt-',rownames(anlas),invert=T),]

# append batch ident metadata
pijuan <- AddMetaData(pijuan, gsub(" ","",paste("pijuan_",pijuan$stage,"_sample",pijuan$sample)), col.name="batch.ident")
anlas <- AddMetaData(anlas, gsub(" ","",paste("anlas_",anlas$merge.ident,"_",anlas$replicate)), col.name="batch.ident")

# filter for intersect features
common_features <- intersect(rownames(pijuan), rownames(anlas))
pijuan <- pijuan[common_features,]
anlas <- anlas[common_features,]

# remove integrated assay from anlas
anlas[['integrated']] <- NULL

# extract metadata
meta_pijuan <- pijuan@meta.data
meta_anlas <- anlas@meta.data

# reassign celltype names for pijuan dataset
c <- c()
for(i in 1:length(meta_pijuan$celltype.pijuan)){
  old <- meta_pijuan$celltype.pijuan[i]
  if(old=='Epiblast'){
    c[i] <- 'Epiblast'
  } else if(old=='Primitive Streak'){
    c[i] <- 'Primitive Streak'
  } else if(old=='ExE ectoderm'){
    c[i] <- 'ExE'
  } else if(old=='Visceral endoderm'){
    c[i] <- 'ExE'
  } else if(old=='ExE endoderm'){
    c[i] <- 'ExE'
  } else if(old=='Nascent mesoderm'){
    c[i] <- 'Mesoderm'
  } else if(old=='Rostral neurectoderm'){
    c[i] <- 'Ectoderm'
  } else if(old=='Blood progenitors 2'){
    c[i] <- 'Mesoderm'
  } else if(old=='Mixed mesoderm'){
    c[i] <- 'Mesoderm'
  } else if(old=='ExE mesoderm'){
    c[i] <- 'ExE'
  } else if(old=='Intermediate mesoderm'){
    c[i] <- 'Mesoderm'
  } else if(old=='Pharyngeal mesoderm'){
    c[i] <- 'Mesoderm'
  } else if(old=='Caudal epiblast'){
    c[i] <- 'Epiblast'
  } else if(old=='PGC'){
    c[i] <- 'PGC'
  } else if(old=='Mesenchyme'){
    c[i] <- 'Mesoderm'
  } else if(old=='Haematoendothelial progenitors'){
    c[i] <- 'Mesoderm'
  } else if(old=='Blood progenitors 1'){
    c[i] <- 'Mesoderm'
  } else if(old=='Surface ectoderm'){
    c[i] <- 'Ectoderm'
  } else if(old=='Gut'){
    c[i] <- 'Endoderm'
  } else if(old=='Paraxial mesoderm'){
    c[i] <- 'Mesoderm'
  } else if(old=="Caudal neurectoderm"){
    c[i] <- 'Ectoderm'
  } else if(old=="Notochord"){
    c[i] <- 'Mesoderm'
  } else if(old=="Somitic mesoderm"){
    c[i] <- 'Mesoderm'
  } else if(old=="Caudal Mesoderm"){
    c[i] <- 'Mesoderm'
  } else if(old=="Erythroid1"){
    c[i] <- 'Mesoderm'
  } else if(old=="Def. endoderm"){
    c[i] <- 'Endoderm'
  } else if(old=="Parietal endoderm"){
    c[i] <- 'ExE'
  } else if(old=="Allantois"){
    c[i] <- 'ExE'
  } else if(old=="Anterior Primitive Streak"){
    c[i] <- 'Primitive Streak'
  } else if(old=="Endothelium"){
    c[i] <- 'Mesoderm'
  } else if(old=="Forebrain/Midbrain/Hindbrain"){
    c[i] <- 'Ectoderm'
  } else if(old=="Spinal cord"){
    c[i] <- 'Ectoderm'
  } else if(old=="Cardiomyocytes"){
    c[i] <- 'Mesoderm'
  } else if(old=="Erythroid2"){
    c[i] <- 'Mesoderm'
  } else if(old=="NMP"){
    c[i] <- 'Mesoderm'
  } else if(old=="Erythroid3"){
    c[i] <- 'Mesoderm'
  } else if(old=="Neural crest"){
    c[i] <- 'Ectoderm'
  }
}
pijuan$celltype.general <- c
# use this to change ident
Idents(pijuan) <- pijuan$celltype.general
levels(pijuan)

# reassign celltype names for anlas dataset
meta_anlas <- anlas@meta.data
c <- c()
for(i in 1:length(meta_anlas$celltype.anlas)){
  old <- meta_anlas$celltype.anlas[i]
  if(grepl("Early diff", old)){
    c[i] <- 'G Early diff'
  } else if(grepl("Early vasc", old)){
    c[i] <- 'G Mesodermal'
  } else if(grepl("Primed ", old)){
    c[i] <- 'G Primed pluripotent'
  } else if(grepl("Pluripotent ", old)){
    c[i] <- 'G Pluripotent'
  } else if(grepl("Mesodermal", old)){
    c[i] <- 'G Mesodermal'
  } else if(grepl("Neural ", old)){
    c[i] <- 'G Neural'
  } else if(grepl("neural&", old)){
    c[i] <- 'G Neural'
  } else{
    c[i] <- paste("G",old)
  }
}
anlas$celltype.general <- c
# use this to change ident
Idents(anlas) <- anlas$celltype.general
levels(anlas)

# merge seurat objects
pijuan@project.name <- "pijuan"
anlas@project.name <- "anlas"
pijuan$dataset <- "pijuan"
anlas$dataset <- "anlas"
seurat_all <- merge( x = pijuan, y = anlas )

# clean up memory
remove(pijuan, anlas, meta_anlas, meta_pijuan)
gc()

###
levels(seurat_all)

# features <- c("T","Nanog","Sox2")
# VlnPlot(seurat_all, features = features)

# filter only the cell types you want to run the diff_exp analysis on
seurat_pluripotent <- seurat_all[,(seurat_all$celltype.general=='Epiblast')|(seurat_all$celltype.general=='G Pluripotent')|(seurat_all$celltype.general=='G Primed pluripotent')|(seurat_all$celltype.general=='G Early diff')]
levels(seurat_pluripotent)

############################################
# markers.epi.pluri <- FindMarkers(seurat_all, ident.1 = "Epiblast", ident.2 = "G Pluripotent",
#                                  test.use = "LR", latent.vars = "batch.ident", min.pct = 0.5)
# head(markers.epi.pluri)
# 
# markers.epi.primed <- FindMarkers(seurat_all, ident.1 = "Epiblast", ident.2 = "G Primed pluripotent",
#                                   test.use = "LR", latent.vars = "batch.ident", min.pct = 0.5)
# head(markers.epi.primed)
# 
# write.csv(seurat_pluripotent[labels(markers.epi.pluri)[[1]],(seurat_all$celltype.general=='Epiblast')|(seurat_all$celltype.general=='G Pluripotent')][['RNA']]@data,paste0(outFolder,"expression_markers_epi_pluri.csv"))
# # write.csv(seurat_pluripotent$celltype.general,paste0(outFolder,"cell_id.csv"))
# write.csv(markers.epi.pluri,paste0(outFolder,"markers_epi_pluri.csv"))
# 
# write.csv(seurat_pluripotent[labels(markers.epi.primed)[[1]],(seurat_all$celltype.general=='Epiblast')|(seurat_all$celltype.general=='G Primed pluripotent')][['RNA']]@data,paste0(outFolder,"expression_markers_epi_primed.csv"))
# # write.csv(seurat_pluripotent$celltype.general,paste0(outFolder,"cell_id.csv"))
# write.csv(markers.epi.primed,paste0(outFolder,"markers_epi_primed.csv"))
############################################

############################################
########## make diff exp analysis
############################################

### epiblast vs pluripotent

c_pluri <- seurat_pluripotent[,(seurat_pluripotent$celltype.general=='Epiblast')|
                               (seurat_pluripotent$celltype.general=='G Pluripotent')]
Idents(c_pluri) <- "dataset"
avg.cells.pluri <- log1p(AverageExpression(c_pluri, verbose = FALSE)$RNA)
avg.cells.pluri$gene <- rownames(avg.cells.pluri)

# compute avg_logFC as in FindMarkers of Seurat
logfc <- c()
for (i in 1: length(avg.cells.pluri$gene)){
  a <- avg.cells.pluri$pijuan[i]
  b <- avg.cells.pluri$anlas[i]
  logfc[i] <- a-b
}
avg.cells.pluri$avg_logFC <- logfc

# extract markers based on logfc
markers_epi_pluri <- avg.cells.pluri[(avg.cells.pluri$avg_logFC>1.0)|(avg.cells.pluri$avg_logFC<(-1.0)),]

p1 <- ggplot(avg.cells.pluri, aes(pijuan, anlas)) + geom_point() + ggtitle("Epiblast vs Pluripotent")
p1 <- LabelPoints(plot = p1, points = markers_epi_pluri$gene, repel = TRUE)
p1

write.csv(c_pluri[labels(markers_epi_pluri)[[1]],][['RNA']]@data,paste0(outFolder,"expression_markersDiff_epi_pluri.csv"))
# write.csv(seurat_pluripotent$celltype.general,paste0(outFolder,"cell_id.csv"))
write.csv(markers_epi_pluri,paste0(outFolder,"markersDiff_epi_pluri.csv"))

### epiblast vs primed pluripotent

c_primed <- seurat_pluripotent[,(seurat_pluripotent$celltype.general=='Epiblast')|
                                (seurat_pluripotent$celltype.general=='G Primed pluripotent')]
Idents(c_primed) <- "dataset"
avg.cells.primed <- log1p(AverageExpression(c_primed, verbose = FALSE)$RNA)
avg.cells.primed$gene <- rownames(avg.cells.primed)

# compute avg_logFC as in FindMarkers of Seurat
logfc <- c()
for (i in 1: length(avg.cells.primed$gene)){
  a <- avg.cells.primed$pijuan[i]
  b <- avg.cells.primed$anlas[i]
  logfc[i] <- a-b
}
avg.cells.primed$avg_logFC <- logfc

# extract markers based on logfc
markers_epi_primed <- avg.cells.primed[(avg.cells.primed$avg_logFC>1.0)|(avg.cells.primed$avg_logFC<(-1.0)),]

p2 <- ggplot(avg.cells.primed, aes(pijuan, anlas)) + geom_point() + ggtitle("Epiblast vs Primed")
p2 <- LabelPoints(plot = p2, points = markers_epi_primed$gene, repel = TRUE)
p2

write.csv(c_primed[labels(markers_epi_primed)[[1]],][['RNA']]@data,paste0(outFolder,"expression_markersDiff_epi_primed.csv"))
# write.csv(seurat_pluripotent$celltype.general,paste0(outFolder,"cell_id.csv"))
write.csv(markers_epi_primed,paste0(outFolder,"markersDiff_epi_primed.csv"))

### epiblast vs early diff

c_ediff <- seurat_pluripotent[,(seurat_pluripotent$celltype.general=='Epiblast')|
                                 (seurat_pluripotent$celltype.general=='G Early diff')]
Idents(c_ediff) <- "dataset"
avg.cells.ediff <- log1p(AverageExpression(c_ediff, verbose = FALSE)$RNA)
avg.cells.ediff$gene <- rownames(avg.cells.ediff)

# compute avg_logFC as in FindMarkers of Seurat
logfc <- c()
for (i in 1: length(avg.cells.ediff$gene)){
  a <- avg.cells.ediff$pijuan[i]
  b <- avg.cells.ediff$anlas[i]
  logfc[i] <- a-b
}
avg.cells.ediff$avg_logFC <- logfc

# extract markers based on logfc
markers_epi_ediff <- avg.cells.ediff[(avg.cells.ediff$avg_logFC>1.0)|(avg.cells.ediff$avg_logFC<(-1.0)),]

p3 <- ggplot(avg.cells.ediff, aes(pijuan, anlas)) + geom_point() + ggtitle("Epiblast vs Ediff")
p3 <- LabelPoints(plot = p3, points = markers_epi_ediff$gene, repel = TRUE)
p3

write.csv(c_ediff[labels(markers_epi_primed)[[1]],][['RNA']]@data,paste0(outFolder,"expression_markersDiff_epi_ediff.csv"))
# write.csv(seurat_pluripotent$celltype.general,paste0(outFolder,"cell_id.csv"))
write.csv(markers_epi_primed,paste0(outFolder,"markersDiff_epi_ediff.csv"))

############################################
########## save data for marker genes
############################################

seurat_all_early <- seurat_all[,((seurat_all$stage!='E7.75')&
                        (seurat_all$stage!='E8.0')&
                        (seurat_all$stage!='E8.25')&
                        (seurat_all$stage!='E8.5'))|(seurat_all$dataset=='anlas')]

genes_pluri <- c('Oct4','Nanog','Sox2','Esrrb','Wnt8a','Nkx1-2','Lefty1','Lefty2','Tdgf1','Evx1')

genes_meso <- c('T','Eomes','Six2','Mesp1','Mesp2','Cer1','Msgn1','Hand1',
           'Myl7','Myocd','Tnnt2','Msgn1','Snail1','Foxc2','Pdgfra','Osr1','Tbx6','Mesp1')

genes_endo <- c('Sox17','Cer1','Gata6','Foxa2','Mixl1','Lhx1')

genes_ecto <- c('Sox2','Oct6','Zeb2','Otx2','Zic2','Zic3','Olig2','Mab21l2','Nkx6-2','Neurog2','Pax6',
                'Epha2','Zic1','Zic5','Hoxa2','Msx3','Utf1','Grik3','Slc7a3','Sox1','Sox3','Olig3','Neurod1')

genes_mech <- c('cdh1','cdh2','Mmps')

genes_all <- c('Zfp42','Nanog','Esrrb','Klf4','T','Eomes','Gsc','Lefty1','Mixl1','Mesp1','Mesp2','Lhx1','Hand1',
               'Tbx6','Msgn1','Foxc2','Pdgfra','Osr1','Snai1','Myl7','Tnnt2','Cer1','Six2','Meox1','Sox1','Sox2','Sox3',
               'Pou3f1','Zeb2','Otx2','Zic1','Zic2','Zic3','Zic5','Olig2','Neurog2','Nes','Gbx2','Pax2','Pax6','Ncam1','Utf1',
               'Epha2','Hoxa2','Msx3','Grik3','Slc7a3','Sox17','Foxa2','Gata4')

### pluripotency
c_pluri <- seurat_all_early[,(seurat_all_early$celltype.general=='Epiblast')|
                                (seurat_all_early$celltype.general=='Primitive Streak')|
                                (seurat_all_early$celltype.general=='G Pluripotent')|
                                (seurat_all_early$celltype.general=='G Early diff')|
                                (seurat_all_early$celltype.general=='G Primed pluripotent')]
write.csv(c_pluri[genes_pluri,][['RNA']]@data,paste0(outFolder,"expression_markersGenes_pluri.csv"))
write.csv(c_pluri@meta.data,paste0(outFolder,"expression_markersGenes_pluri_meta.csv"))

### mesoderm
c_meso <- seurat_all_early[,(seurat_all_early$celltype.general=='Mesoderm')|
                       (seurat_all_early$celltype.general=='Primitive Streak')|
                       (seurat_all_early$celltype.general=='G Mesodermal')|
                       (seurat_all_early$celltype.general=='G Mesendodermal')]
write.csv(c_meso[genes_meso,][['RNA']]@data,paste0(outFolder,"expression_markersGenes_meso.csv"))
write.csv(c_meso@meta.data,paste0(outFolder,"expression_markersGenes_meso_meta.csv"))

### endoderm
c_endo <- seurat_all_early[,(seurat_all_early$celltype.general=='Endoderm')|
                       (seurat_all_early$celltype.general=='G Endodermal')]
write.csv(c_endo[genes_endo,][['RNA']]@data,paste0(outFolder,"expression_markersGenes_endo.csv"))
write.csv(c_endo@meta.data,paste0(outFolder,"expression_markersGenes_endo_meta.csv"))

### ectoderm
c_ecto <- seurat_all_early[,(seurat_all_early$celltype.general=='Ectoderm')|
                       (seurat_all_early$celltype.general=='G Neural')]
write.csv(c_ecto[genes_ecto,][['RNA']]@data,paste0(outFolder,"expression_markersGenes_ecto.csv"))
write.csv(c_ecto@meta.data,paste0(outFolder,"expression_markersGenes_ecto_meta.csv"))

#### all markers
c_all <- seurat_all_early[genes_all,]
write.csv(c_all[['RNA']]@data,paste0(outFolder,"expression_markersGenes_ALL.csv"))
write.csv(c_all@meta.data,paste0(outFolder,"expression_markersGenes_ALL_meta.csv"))

