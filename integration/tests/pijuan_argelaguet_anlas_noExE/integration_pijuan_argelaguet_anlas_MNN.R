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

######################################################################################################################

getHVGs <- function(sce, block, min.mean = 1e-3){
  decomp <- modelGeneVar(sce, block=block)
  decomp <- decomp[decomp$mean > min.mean,]
  decomp$FDR <- p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$p.value < 0.05])
}

doBatchCorrect <- function(counts, timepoints, samples, timepoint_order, sample_order, npc = 50, pc_override = NULL, BPPARAM = SerialParam()){
  require(BiocParallel)
  
  if(!is.null(pc_override)){
    pca = pc_override
  } else {
    pca = irlba::prcomp_irlba(t(counts), n = npc)$x
    rownames(pca) = colnames(counts)
  }
  
  if(length(unique(samples)) == 1){
    return(pca)
  }
  
  #create nested list
  pc_list    <- lapply(unique(timepoints), function(tp){
    sub_pc   <- pca[timepoints == tp, , drop = FALSE]
    sub_samp <- samples[timepoints == tp]
    list     <- lapply(unique(sub_samp), function(samp){
      sub_pc[sub_samp == samp, , drop = FALSE]
    })
    names(list) <- unique(sub_samp)
    return(list)
  })
  
  names(pc_list) <- unique(timepoints)
  
  #arrange to match timepoint order
  pc_list <- pc_list[order(match(names(pc_list), timepoint_order))]
  pc_list <- lapply(pc_list, function(x){
    x[order(match(names(x), sample_order))]
  })
  
  #perform corrections within list elements (i.e. within stages)
  correct_list <- lapply(pc_list, function(x){
    if(length(x) > 1){
      #return(do.call(scran::fastMNN, c(x, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected)
      return(do.call(reducedMNN, c(x, BPPARAM = BPPARAM))$corrected) # edited 09.02.2020 because of "Error: 'fastMNN' is not an exported object from 'namespace:scran'", 17.02.2020 changed to reducedMNN because otherwise it thinks PCA space is logcounts which would be utter bullcrap
    } else {
      return(x[[1]])
    }
  })
  
  #perform correction over list
  if(length(correct_list)>1){
    #correct <- do.call(scran::fastMNN, c(correct_list, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected
    correct <- do.call(reducedMNN, c(correct_list, BPPARAM = BPPARAM))$corrected # edited 09.02.2020 because of "Error: 'fastMNN' is not an exported object from 'namespace:scran'", 17.02.2020 changed to reducedMNN because otherwise it thinks PCA space is logcounts which would be utter bullcrap
  } else {
    correct <- correct_list[[1]]
  }
  
  correct <- correct[match(colnames(counts), rownames(correct)),]
  
  return(correct)
  
}

#######################################################################################################################

# load dataset
# load dataset
if(Sys.info()['nodename']=='PCBA-TRIVEDI02'){
  folder.RData <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\R_codes\\integration_MNN\\"
  folder <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\R_codes\\integration_MNN\\pijuan_argelaguet_anlas_noExE\\"
} else if(Sys.info()['nodename']=='PCBA-TRIVEDI03'){
  folder.RData <- "C:\\Users\\nicol\\OneDrive\\Desktop\\integration_MNN\\"
  folder <- "C:\\Users\\nicol\\OneDrive\\Desktop\\integration_MNN\\pijuan_argelaguet_anlas_noExE\\"
}
load(paste0(folder.RData,"anlas_pijuan_data.RData"))
load(paste0(folder.RData,"argelaguet_data.RData"))

# filter pijuan data to only stages needed
pijuan <- pijuan[,(pijuan$stage=='E6.5')|(pijuan$stage=='E6.75')|(pijuan$stage=='E7.0')|(pijuan$stage=='E7.25')|(pijuan$stage=='E7.5')]
pijuan <- pijuan[,grep('ExE', pijuan$celltype.pijuan, invert=T)]
pijuan <- pijuan[,grep('Visceral endoderm', pijuan$celltype.pijuan, invert=T)]
pijuan <- pijuan[rowSums(pijuan)>10,]
pijuan <- pijuan[grep('mt-',rownames(pijuan),invert=T),]

# filter argelaguet data to only good cells and genes
argelaguet <- AddMetaData(argelaguet, argelaguet$lineage, col.name="celltype.argelaguet")
argelaguet <- argelaguet[,as.logical(argelaguet$pass_rnaQC)]
argelaguet <- argelaguet[,argelaguet$lineage!='NOIDEA']
argelaguet <- argelaguet[,argelaguet$lineage!='Primitive_endoderm']
argelaguet <- argelaguet[,argelaguet$lineage!='ExE_ectoderm']
argelaguet <- argelaguet[,argelaguet$lineage!='Visceral_endoderm']
argelaguet <- argelaguet[rowSums(argelaguet)>10,]
argelaguet <- argelaguet[grep('mt-',rownames(argelaguet),invert=T),]

# filter argelaguet data to only good cells and genes
anlas <- anlas[rowSums(anlas)>10,]
anlas <- anlas[grep('mt-',rownames(anlas),invert=T),]

# append batch ident metadata
pijuan <- AddMetaData(pijuan, gsub(" ","",paste("pijuan_",pijuan$stage,"_sample",pijuan$sample)), col.name="batch.ident")
argelaguet <- AddMetaData(argelaguet, gsub(" ","",paste("argelaguet_",argelaguet$stage)), col.name="batch.ident")
anlas <- AddMetaData(anlas, gsub(" ","",paste("anlas_",anlas$merge.ident,"_",anlas$replicate)), col.name="batch.ident")

# filter for intersect features
common_features <- intersect(rownames(pijuan),rownames(argelaguet))
common_features <- intersect(common_features, rownames(anlas))
pijuan <- pijuan[common_features,]
argelaguet <- argelaguet[common_features,]
anlas <- anlas[common_features,]

# remove integrated assay from anlas
anlas[['integrated']] <- NULL

# merge seurat objects
pijuan@project.name <- "pijuan"
argelaguet@project.name <- "argelaguet"
anlas@project.name <- "anlas"
seurat_all <- merge( x = pijuan, y = argelaguet )
seurat_all <- merge( x = seurat_all, y = anlas)

# convert into sce and extract metadata
sce_all <- as.SingleCellExperiment(seurat_all)
sce_pijuan <- as.SingleCellExperiment(pijuan)
sce_argelaguet <- as.SingleCellExperiment(argelaguet)
sce_anlas <- as.SingleCellExperiment(anlas)
meta_pijuan <- pijuan@meta.data
meta_argelaguet <- argelaguet@meta.data
meta_anlas <- anlas@meta.data

# clean up memory
remove(pijuan, argelaguet, anlas)
remove(seurat_all)
gc()

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
  } else if(old=='Surface ectoderm'){
    c[i] <- 'Ectoderm'
  } else if(old=='Blood progenitors 1'){
    c[i] <- 'Mesoderm'
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
  } else if(old=="Erythroid3"){
    c[i] <- 'Mesoderm'
  } else if(old=="Erythroid2"){
    c[i] <- 'Mesoderm'
  }
}
meta_pijuan$celltype.general <- c

# reassign celltype names for argelaguet dataset
c <- c()
for(i in 1:length(meta_argelaguet$celltype.argelaguet)){
  old <- meta_argelaguet$celltype.argelaguet[i]
  if(grepl("ExE", old)){
    c[i] <- 'ExE'
  } else if(grepl("Visceral_endoderm", old)){
    c[i] <- 'ExE'
  } else if(grepl("Primitive_endoderm", old)){
    c[i] <- 'ExE'
  } else if(grepl("Primitive_Streak", old)){
    c[i] <- 'Primitive Streak'
  } else {
    c[i] <- old
  }
}
meta_argelaguet$celltype.general <- c

# reassign celltype names for anlas dataset
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
meta_anlas$celltype.general <- c

### test which genes are not expressed in batches
batch.ident_all <- c(meta_pijuan$batch.ident, meta_argelaguet$batch.ident, meta_anlas$batch.ident)
batch.ident_unique <- unique(batch.ident_all)

for(i in 1:length(batch.ident_unique)){
  submat <- counts(sce_all)[,batch.ident_all==batch.ident_unique[i]]
  print(batch.ident_unique[i])
  print(dim(submat))
  #  print(which(rowMeans(submat)<=0))
}

# normalize sce within each batch
sizeFactors(sce_all) <- NULL
big_sce <- multiBatchNorm(sce_all, batch=batch.ident_all)

# extract hvg
hvgs <- getHVGs(big_sce, block=batch.ident_all)

### perform pca
npcs = 50
big_pca <- multiBatchPCA(big_sce,
                         batch=batch.ident_all,
                         subset.row = hvgs,
                         d = npcs,
                         preserve.single = TRUE,
                         assay.type = "logcounts")[[1]]

n.cell_pijuan <- length(meta_pijuan$batch.ident)
n.cell_argelaguet <- length(meta_argelaguet$batch.ident)
n.cell_anlas <- length(meta_anlas$batch.ident)

rownames(big_pca) <- colnames(big_sce) 

pijuan_rows <- 1:n.cell_pijuan
argelaguet_rows <- (n.cell_pijuan+1):(n.cell_pijuan+n.cell_argelaguet)
anlas_rows <- (n.cell_pijuan+n.cell_argelaguet+1):(n.cell_pijuan+n.cell_anlas+n.cell_argelaguet)

pijuan_pca <- big_pca[pijuan_rows,]
argelaguet_pca <- big_pca[argelaguet_rows,]
anlas_pca <- big_pca[anlas_rows,]
message("Done\n")

### batch correction on the pijuan dataset
message("Batch effect correction for the atlas...")  
order_df        <- meta_pijuan[!duplicated(meta_pijuan$sample), c("stage", "sample")]
order_df$ncells <- sapply(order_df$sample, function(x) sum(meta_pijuan$sample == x))
order_df$stage  <- factor(order_df$stage, 
                          levels = rev(c("E7.5","E7.25","E7.0","E6.75","E6.5")))
order_df       <- order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
order_df$stage <- as.character(order_df$stage)

set.seed(42)
correct <- doBatchCorrect(counts         = logcounts(sce_pijuan[hvgs,]), 
                                  timepoints      = meta_pijuan$stage, 
                                  samples         = meta_pijuan$sample, 
                                  timepoint_order = order_df$stage, 
                                  sample_order    = order_df$sample, 
                                  pc_override     = pijuan_pca,
                                  npc             = npcs)
message("Done\n")

### map argelaguet with mnn
correct <- reducedMNN(rbind(correct, argelaguet_pca, anlas_pca),
                      batch=c(rep("ATLAS", n.cell_pijuan), meta_argelaguet$batch.ident, meta_anlas$batch.ident),
                      merge.order=NULL)$corrected
# ### map anlas with mnn
# correct <- reducedMNN(rbind(correct, anlas_pca),
#                       batch=c(rep("ATLAS", n.cell_pijuan+n.cell_argelaguet), meta_anlas$batch.ident),
#                       merge.order=NULL)$corrected

correct_pijuan <- correct[pijuan_rows,]
correct_argelaguet <- correct[argelaguet_rows,]
correct_anlas <- correct[anlas_rows,]

# mapping <- get_meta(correct_atlas = correct_atlas,
#                     atlas_meta = atlas_meta,
#                     correct_map = correct_map,
#                     map_meta = map_meta,
#                     k_map = k)
message("Done\n")

### compute umap
correct.umap <- umap(correct)

### save data
write.table(correct, file= paste0(folder,"integration.csv"), sep=",")
write.table(correct.umap$layout, file= paste0(folder,"integration_umap.csv"), sep=",")
write.table(meta_pijuan, file= paste0(folder,"meta_pijuan.csv"), sep=",")
write.table(meta_argelaguet, file= paste0(folder,"meta_argelaguet.csv"), sep=",")
write.table(meta_anlas, file= paste0(folder,"meta_anlas.csv"), sep=",")
save(correct, correct.umap, meta_pijuan, meta_argelaguet, meta_anlas, pijuan_rows, argelaguet_rows, anlas_rows, file= paste0(folder,"corrected_data.RData"))

rm(list = ls())
gc()

