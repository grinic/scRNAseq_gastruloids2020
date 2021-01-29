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


data_location <- "server"
data_location <- "local"
data_location <- "ext_drive"

if(data_location == "server"){
  folder.RData <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\data\\"
  outFolder <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\results\\integration\\pijuan\\"
} else if(data_location == "local"){
  folder.RData <- "C:\\Users\\nicol\\OneDrive\\Desktop\\scrnaseq_gastruloids\\data\\"
  outFolder <- "C:\\Users\\nicol\\OneDrive\\Desktop\\scrnaseq_gastruloids\\results\\integration\\pijuan\\"
} else if(data_location == "ext_drive"){
  folder.RData <- "F:\\scrnaseq_gastruloids\\data\\"
  outFolder <- "F:\\scrnaseq_gastruloids\\results\\integration\\pijuan\\"
}
load(paste0(folder.RData,"pijuan_data.RData"))

# filter pijuan data to only stages needed
pijuan <- NormalizeData(object = pijuan, normalization.method = "LogNormalize", scale.factor = 10000)
# pijuan <- pijuan[,(pijuan$stage=='E6.5')|(pijuan$stage=='E6.75')|(pijuan$stage=='E7.0')|(pijuan$stage=='E7.25')|(pijuan$stage=='E7.5')]
pijuan <- pijuan[rowSums(pijuan)>10,]
pijuan <- pijuan[grep('mt-',rownames(pijuan),invert=T),]

# append batch ident metadata
pijuan <- AddMetaData(pijuan, gsub(" ","",paste("pijuan_",pijuan$stage,"_sample",pijuan$sample)), col.name="batch.ident")

# merge seurat objects
pijuan@project.name <- "pijuan"

# convert into sce and extract metadata
sce_pijuan <- as.SingleCellExperiment(pijuan)
meta_pijuan <- pijuan@meta.data

# clean up memory
remove(pijuan)
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
meta_pijuan$celltype.general <- c

### test which genes are not expressed in batches
batch.ident_all <- c(meta_pijuan$batch.ident)
batch.ident_unique <- unique(batch.ident_all)

for(i in 1:length(batch.ident_unique)){
  submat <- counts(sce_pijuan)[,batch.ident_all==batch.ident_unique[i]]
  print(batch.ident_unique[i])
  print(dim(submat))
  #  print(which(rowMeans(submat)<=0))
}

# normalize sce within each batch
sizeFactors(sce_pijuan) <- NULL
big_sce <- multiBatchNorm(sce_pijuan, batch=batch.ident_all)

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

rownames(big_pca) <- colnames(big_sce) 

pijuan_rows <- 1:n.cell_pijuan

pijuan_pca <- big_pca[pijuan_rows,]
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

# mapping <- get_meta(correct_atlas = correct_atlas,
#                     atlas_meta = atlas_meta,
#                     correct_map = correct_map,
#                     map_meta = map_meta,
#                     k_map = k)
message("Done\n")

### save data
write.table(correct, file= paste0(outFolder,"integration.csv"), sep=",")
write.table(meta_pijuan, file= paste0(outFolder,"meta_pijuan.csv"), sep=",")
save(correct, meta_pijuan, pijuan_rows, file= paste0(outFolder,"corrected_data.RData"))

rm(list = ls())
gc()

