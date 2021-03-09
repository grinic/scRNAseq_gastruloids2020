if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("batchelor")
#BiocManager::install("Seurat")
#BiocManager::install("scater")
#BiocManager::install("umap")
#BiocManager::install("scRNAseq")
#BiocManager::install("zellkonverter")

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
suppressPackageStartupMessages(library(scRNAseq))
suppressPackageStartupMessages(library(zellkonverter))


######################################################################################################################

getHVGs <- function(sce, block, min.mean = 1e-3){
  decomp <- modelGeneVar(sce, block=block)
  decomp <- decomp[decomp$mean > min.mean,]
  decomp$FDR <- p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$p.value < 0.05])
}

getmode <- function(v, dist) {
  tab <- table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied <- names(tab)[tab == max(tab)]
    sub  <- dist[v %in% tied]
    names(sub) <- v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

getcelltypes <- function(v, dist) {
  tab <- table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied <- names(tab)[tab == max(tab)]
    sub  <- dist[v %in% tied]
    names(sub) <- v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
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

getMappingScore <- function(mapping){
  out <- list()
  celltypes_accrossK <- matrix(unlist(mapping$celltypes.mapped),
                               nrow=length(mapping$celltypes.mapped[[1]]),
                               ncol=length(mapping$celltypes.mapped))
  cellstages_accrossK <- matrix(unlist(mapping$cellstages.mapped),
                                nrow=length(mapping$cellstages.mapped[[1]]),
                                ncol=length(mapping$cellstages.mapped))
  out$celltype.score <- NULL
  for (i in 1:nrow(celltypes_accrossK)){
    p <- max(table(celltypes_accrossK[i,]))
    index <- which(table(celltypes_accrossK[i,]) == p)
    p <- p/length(mapping$celltypes.mapped)
    out$celltype.score <- c(out$celltype.score,p)
  }
  out$cellstage.score <- NULL
  for (i in 1:nrow(cellstages_accrossK)){
    p <- max(table(cellstages_accrossK[i,]))
    index <- which(table(cellstages_accrossK[i,]) == p)
    p <- p/length(mapping$cellstages.mapped)
    out$cellstage.score <- c(out$cellstage.score,p)
  }
  return(out)  
}

get_meta <- function(correct_atlas, atlas_meta, correct_map, map_meta, k_map = 10){
  knns <- BiocNeighbors::queryKNN(correct_atlas, correct_map, k = k_map, get.index = TRUE,
                                  get.distance = FALSE)
  #get closest k matching cells
  k.mapped  <- t(apply(knns$index, 1, function(x) atlas_meta$cell[x]))
  celltypes <- t(apply(k.mapped, 1, function(x) atlas_meta$celltype.pijuan[match(x, atlas_meta$cell)]))
  stages    <- t(apply(k.mapped, 1, function(x) atlas_meta$stage[match(x, atlas_meta$cell)]))
  celltype.mapped <- apply(celltypes, 1, function(x) getmode(x, 1:length(x)))
  stage.mapped    <- apply(stages, 1, function(x) getmode(x, 1:length(x)))
  out <- lapply(1:length(celltype.mapped), function(x){
    list(cells.mapped     = k.mapped[x,],
         celltype.mapped  = celltype.mapped[x],
         stage.mapped     = stage.mapped[x],
         celltypes.mapped = celltypes[x,],
         stages.mapped    = stages[x,])
  })
  names(out) <- map_meta$cell
  return(out)  
}

#######################################################################################################################

data_location <- "server"
# data_location <- "local"
# data_location <- "ext_drive"

if(data_location == "server"){
  anlas.RData <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\raw_data\\gastruloids\\"
  folder.RData <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\data\\"
  outFolder <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\results\\integration\\integration_withGastrIntegratedData\\"
} else if(data_location == "local"){
  folder.RData <- "C:\\Users\\nicol\\OneDrive\\Desktop\\scrnaseq_gastruloids\\data\\"
  outFolder <- "C:\\Users\\nicol\\OneDrive\\Desktop\\scrnaseq_gastruloids\\results\\integration\\integration_withGastrIntegratedData\\"
} else if(data_location == "ext_drive"){
  folder.RData <- "F:\\scrnaseq_gastruloids\\data\\"
  outFolder <- "F:\\scrnaseq_gastruloids\\results\\integration\\integration_withGastrIntegratedData\\"
}
anlas <- readRDS(paste0(anlas.RData,"all_data.combined_2.rds"))
DefaultAssay(anlas) <- "RNA"
# remove integrated assay from anlas
anlas[['integrated']] <- NULL
load(paste0(folder.RData,"pijuan_data.RData"))

# filter pijuan data to only stages needed
pijuan <- NormalizeData(object = pijuan, normalization.method = "LogNormalize", scale.factor = 10000)
anlas <- NormalizeData(object = anlas, normalization.method = "LogNormalize", scale.factor = 10000)
# pijuan <- pijuan[,(pijuan$stage=='E6.5')|(pijuan$stage=='E6.75')|(pijuan$stage=='E7.0')|(pijuan$stage=='E7.25')|(pijuan$stage=='E7.5')]
pijuan <- pijuan[rowSums(pijuan)>10,]
pijuan <- pijuan[grep('mt-',rownames(pijuan),invert=T),]

# filter argelaguet data to only good cells and genes
# anlas <- anlas[,anlas$merge.ident!='0h'] # there are no 0h cells in this dataset
anlas <- anlas[rowSums(anlas)>10,]
anlas <- anlas[grep('mt-',rownames(anlas),invert=T),]

# append batch ident metadata
pijuan <- AddMetaData(pijuan, gsub(" ","",paste("pijuan_",pijuan$stage,"_sample",pijuan$sample)), col.name="batch_ident")
batch.ident <- c()
new.stage <- c()
new.rep <- c()
m <- anlas$merge.ident
r <- anlas$replicate
for(i in 1:length(m)){
  stage <- m[[i]]
  rep <- substr(r[[i]],1,4)
  if(is.na(rep)){
    rep <- "rep2"
  }
  batch.ident[i] <- paste0("anlas_",stage,"_",rep)
  new.stage[i] <- stage
  new.rep[i] <- rep
}
anlas <- AddMetaData(anlas, batch.ident, col.name="batch_ident")
anlas <- AddMetaData(anlas, new.stage, col.name="anlas_stage")
anlas <- AddMetaData(anlas, new.rep, col.name="anlas_rep")

# filter for intersect features
common_features <- intersect(rownames(pijuan), rownames(anlas))
pijuan <- pijuan[common_features,]
anlas <- anlas[common_features,]

# merge seurat objects
pijuan@project.name <- "pijuan"
anlas@project.name <- "anlas"
seurat_all <- merge( x = pijuan, y = anlas )

# convert into sce and extract metadata
sce_all <- as.SingleCellExperiment(seurat_all)
sce_pijuan <- as.SingleCellExperiment(pijuan)
sce_anlas <- as.SingleCellExperiment(anlas)
meta_pijuan <- pijuan@meta.data
meta_anlas <- anlas@meta.data

# clean up memory
remove(pijuan, anlas)
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
meta_pijuan$celltype_simplified <- c

# reassign celltype names for anlas dataset
cl.name <- c()
cl.name.simplified <- c()
for(i in 1:length(meta_anlas$seurat_clusters)){
  old <- meta_anlas$seurat_clusters[i]
  if(old==0){
    cl.name[i] <- "Mesodermal (Bra+)"
    cl.name.simplified[i] <- "Mesodermal"
  } else if(old==1){
    cl.name[i] <- "Pluripotent (naive-primed)"
    cl.name.simplified[i] <- "Pluripotent"
  } else if(old==2){
    cl.name[i] <- "Mesodermal (paraxial, intermediate, lateral plate)"
    cl.name.simplified[i] <- "Mesodermal"
  } else if(old==3){
    cl.name[i] <- "Early differentiated (neural, mesodermal)"
    cl.name.simplified[i] <- "Early diff"
  } else if(old==4){
    cl.name[i] <- "Mesendodermal (Eomes+)"
    cl.name.simplified[i] <- "Mesendodermal"
  } else if(old==5){
    cl.name[i] <- "Primed pluripotent, Early differentiated"
    cl.name.simplified[i] <- "Primed pluripotent"
  } else if(old==6){
    cl.name[i] <- "Pluripotent (naive)"
    cl.name.simplified[i] <- "Pluripotent"
  } else if(old==7){
    cl.name[i] <- "Neural, mesodermal (posterior)"
    cl.name.simplified[i] <- "Neural"
  } else if(old==8){
    cl.name[i] <- "Endodermal"
    cl.name.simplified[i] <- "Endodermal"
  } else if(old==9){
    cl.name[i] <- "Primitive streak, Early differentiated"
    cl.name.simplified[i] <- "Primitive streak"
  } else if(old==10){
    cl.name[i] <- "Early differentiated"
    cl.name.simplified[i] <- "Early diff"
  }
}
meta_anlas$celltype <- cl.name
meta_anlas$celltype_simplified <- cl.name.simplified

### test which genes are not expressed in batches
batch.ident_all <- c(meta_pijuan$batch_ident, meta_anlas$batch_ident)
batch.ident_unique <- unique(batch.ident_all)

# sanity check that all batches are populated with cells
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

######################################

### perform pca
npcs = 50
big_pca <- multiBatchPCA(big_sce,
                         batch=batch.ident_all,
                         subset.row = hvgs,
                         d = npcs,
                         preserve.single = TRUE,
                         assay.type = "logcounts")[[1]]

n.cell_pijuan <- length(meta_pijuan$batch_ident)
n.cell_anlas <- length(meta_anlas$batch_ident)

rownames(big_pca) <- colnames(big_sce) 

pijuan_rows <- 1:n.cell_pijuan
anlas_rows <- (n.cell_pijuan+1):(n.cell_pijuan+n.cell_anlas)

pijuan_pca <- big_pca[pijuan_rows,]
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

### map anlas with mnn
correct <- reducedMNN(rbind(correct, anlas_pca),
                      batch=c(rep("ATLAS", n.cell_pijuan), meta_anlas$batch_ident),
                      merge.order=NULL)$corrected

correct_pijuan <- correct[pijuan_rows,]
correct_anlas <- correct[anlas_rows,]

mapping <- get_meta(correct_atlas = correct_pijuan,
                   atlas_meta = meta_pijuan,
                   correct_map = correct_anlas,
                   map_meta = meta_anlas,
                   k_map = 30)
message("Done\n")

message("Computing mapping scores...") 
out <- list()
for (i in seq(from = 1, to = 30)) {
  out$closest.cells[[i]]     <- sapply(mapping, function(x) x$cells.mapped[i])
  out$celltypes.mapped[[i]]  <- sapply(mapping, function(x) x$celltypes.mapped[i])
  out$cellstages.mapped[[i]] <- sapply(mapping, function(x) x$stages.mapped[i])
}  
multinomial.prob <- getMappingScore(out)
message("Done\n")


cell_mapped <- c()
celltype_mapped <- c()
stage_mapped <- c()
celltype_score <- c()
stage_score <- c()
for(i in 1:length(mapping)){
  cell_mapped[i] <- mapping[[i]]$cells.mapped[1]
  celltype_mapped[i] <- mapping[[i]]$celltype.mapped
  stage_mapped[i] <- mapping[[i]]$stage.mapped
  celltype_score[i] <- multinomial.prob[[1]][i]
  stage_score[i] <- multinomial.prob[[2]][i]
}

out$cell_mapped       <- cell_mapped
out$celltype_mapped   <- celltype_mapped
out$stage_mapped      <- stage_mapped
out$celltype_score    <- celltype_score
out$stage_score       <- stage_score


# message("Writing output...") 
# out$correct_pijuan <- correct_pijuan
# out$correct_anlas <- correct_anlas
# ct <- sapply(mapping, function(x) x$celltype.mapped); is.na(ct) <- lengths(ct) == 0
# st <- sapply(mapping, function(x) x$stage.mapped); is.na(st) <- lengths(st) == 0
# cm <- sapply(mapping, function(x) x$cells.mapped[1]); is.na(cm) <- lengths(cm) == 0
# out$mapping <- data.frame(
#   cell            = names(mapping), 
#   celltype.mapped = unlist(ct),
#   stage.mapped    = unlist(st),
#   closest.cell    = unlist(cm))
# 
# out$mapping <- cbind(out$mapping,multinomial.prob)
# out$pca <- big_pca
# message("Done\n")


### save data
write.table(out, file= paste0(outFolder,"mapping.csv"), sep=",")
write.table(correct, file= paste0(outFolder,"integration.csv"), sep=",")
write.table(meta_pijuan, file= paste0(outFolder,"meta_pijuan.csv"), sep=",")
write.table(meta_anlas, file= paste0(outFolder,"meta_anlas.csv"), sep=",")
save(correct, meta_pijuan, meta_anlas, pijuan_rows, anlas_rows, file= paste0(outFolder,"corrected_data.RData"))

### save the buig_sce as loom to open in python
b <- as.Seurat(big_sce)[hvgs]
c <- as.loom(b, filename = paste0(outFolder,"data_all.loom"), verbose = FALSE)


rm(list = ls())
gc()

