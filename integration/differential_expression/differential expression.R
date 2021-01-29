if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("batchelor")
# BiocManager::install("Seurat")
# BiocManager::install("scater")
# BiocManager::install("umap")

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

if(data_location == "server"){
  folder.RData <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\data\\"
  outFolder <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\results\\integration\\diff_exp\\"
} else if(data_location == "local"){
  folder.RData <- "F:\\scrnaseq_gastruloids\\data\\"
  outFolder <- "F:\\scrnaseq_gastruloids\\data\\integration\\diff_exp\\"
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

# reassign celltype names for pijuan dataset
meta_pijuan <- pijuan@meta.data
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
    c[i] <- 'Endoderm'
  } else if(old=="Allantois"){
    c[i] <- 'Endoderm'
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
  if(grepl("Early ", old)){
    c[i] <- 'G Early diff'
  } else if(grepl("Primed ", old)){
    c[i] <- 'G Primed pluripotent'
  } else if(grepl("Pluripotent ", old)){
    c[i] <- 'G Pluripotent'
  } else if(grepl("Mesodermal", old)){
    c[i] <- 'G Mesoderm'
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
seurat_all <- merge( x = pijuan, y = anlas )

# clean up memory
remove(pijuan, anlas, meta_anlas, meta_pijuan)
gc()

###
levels(seurat_all)
markers <- FindMarkers(seurat_all, ident.1 = "Epiblast", ident.2 = "G Pluripotent", min.pct = 0.5)
head(markers)

features <- c("T","Nanog","Sox2")
VlnPlot(pijuan, features = features)

