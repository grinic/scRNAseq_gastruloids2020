if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("batchelor")
BiocManager::install("Seurat")
BiocManager::install("scater")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(scater))

######
# define arguments
######
# load dataset
if(Sys.info()['nodename']=='PCBA-TRIVEDI02'){
  folder.RData <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\data\\"
} else if(Sys.info()['nodename']=='PCBA-TRIVEDI03'){
  folder.RData <- "C:\\Users\\nicol\\OneDrive\\Desktop\\scRNAseq_Gastruloids\\data\\"
}

data.file <- list(
  "g00h"=paste0(folder.RData,"anlas\\g_00h.rds"),
  "g24h"=paste0(folder.RData,"anlas\\g_24h.rds"),
  "g48h"=paste0(folder.RData,"anlas\\g_48h.rds"),
  "g72h"=paste0(folder.RData,"anlas\\g_72h.rds"),
  "g48hWT"=paste0(folder.RData,"anlas\\g_48h_WT.rds"),
  "pijuan"=paste0(folder.RData,"pijuan"),
  "arg"=paste0(folder.RData,"argelaguet\\SingleCellExperiment.rds")
)

###
# Load pijuan dataset into Seurat object
###

counts = readMM(gsub(" ","",paste(data.file$pijuan,"\\raw_counts.mtx")))
genes = read.table(gsub(" ","",paste(data.file$pijuan,"\\genes.tsv")), stringsAsFactors = F)
meta = read.table(gsub(" ","",paste(data.file$pijuan,"\\meta.csv")), header = TRUE, sep = ",", stringsAsFactors = FALSE, comment.char = "$")

rownames(counts) = genes[,2]
colnames(counts) = meta$cells

sce = SingleCellExperiment(assays = list("counts" = counts))

# remove doublet
sce = sce[,!meta$doublet]
meta = meta[!meta$doublet,]

# remove stripped
sce = sce[,!meta$stripped]
meta = meta[!meta$stripped,]

#get order: oldest to youngest; most cells to least cells
order_df = meta[!duplicated(meta$sample), c("stage", "sample")]
order_df$ncells = sapply(order_df$sample, function(x) sum(meta$sample == x))
order_df$stage = factor(order_df$stage, 
                        levels = rev(c("E8.5", 
                                       "E8.25", 
                                       "E8.0", 
                                       "E7.75", 
                                       "E7.5", 
                                       "E7.25", 
                                       "mixed_gastrulation", 
                                       "E7.0", 
                                       "E6.75", 
                                       "E6.5")))
order_df = order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
order_df$stage = as.character(order_df$stage)

# create Seurat object
pijuan <- as.Seurat(sce,assay='RNA',data=NULL)
pijuan <- AddMetaData(pijuan, meta$cell, col.name="cell")
pijuan <- AddMetaData(pijuan, meta$barcode, col.name="barcode")
pijuan <- AddMetaData(pijuan, meta$sample, col.name="sample")
pijuan <- AddMetaData(pijuan, meta$stage, col.name="stage")
pijuan <- AddMetaData(pijuan, meta$cluster, col.name="cluster")
pijuan <- AddMetaData(pijuan, meta$celltype, col.name="celltype.pijuan")

remove(counts, sce)

### save dataset
save(pijuan, file = paste0(folder.RData,"pijuan_data.RData"))

#####
# load anlas data
#####
g00 <- readRDS(data.file$g00h)
g00$merge.ident <- "0h"
g00$replicate <- "rep2"
g00 <- AddMetaData(object = g00, metadata = gsub(" ","",paste(g00$replicate,"_",g00$merge.ident)), col.name = "integrate.ident")

# assign name to clusters
new.cluster.ids <- c("Early diff (neural priming)", 
                     "Pluripotent", "Early diff (cyto upreg)", "Primed pluripotent", 
                     "Early diff (translation reg)", "Primed pluripotent (cyto upreg)", "Primed pluripotent (splicing reg)")
names(new.cluster.ids) <- levels(g00)
#g00 <- RenameIdents(g00, new.cluster.ids)

x_1 <- c()
for (i in 1 : length(g00$seurat_clusters)) {
  idx <- g00$seurat_clusters[[i]]
  name <- new.cluster.ids[[idx]]
  x_1[i] <- name
}
g00 <- AddMetaData(object = g00, metadata = x_1, col.name = "celltype.anlas")

### load and filter reference dataset for 24h

g24 <- readRDS(data.file$g24h)
g24 <- AddMetaData(object = g24, metadata = gsub(" ","",paste(g24$replicate,"_",g24$merge.ident)), col.name = "integrate.ident")

# assign name to clusters
new.cluster.ids <- c("Early diff (neural&meso priming)", 
                     "Pluripotent II", "Primed pluripotent & early diff", "Mesodermal", "Mesendodermal", "Primed pluripotent", 
                     "Mesodermal & epithelial", "Neural", "Pluripotent")
names(new.cluster.ids) <- levels(g24)
#g24 <- RenameIdents(g24, new.cluster.ids)

x_1 <- c()
for (i in 1 : length(g24$seurat_clusters)) {
  idx <- g24$seurat_clusters[[i]]
  name <- new.cluster.ids[[idx]]
  x_1[i] <- name
}
g24 <- AddMetaData(object = g24, metadata = x_1, col.name = "celltype.anlas")

### load and filter reference dataset for 48h

g48 <- readRDS(data.file$g48h)
g48 <- AddMetaData(object = g48, metadata = gsub(" ","",paste(g48$replicate,"_",g48$merge.ident)), col.name = "integrate.ident")

# assign name to clusters
new.cluster.ids <- c("Mesodermal (posterior)", 
                     "Pluripotent II", "Mesodermal (paraxial&intermediate)", 
                     "Mesendodermal", "Primed pluripotent & early diff","Neural & mesodermal (posterior)",
                     "Endodermal", "Pluripotent", "Mixed (neural&mesodermal)","Neural")
names(new.cluster.ids) <- levels(g48)
#g48 <- RenameIdents(g48, new.cluster.ids)

x_1 <- c()
for (i in 1 : length(g48$seurat_clusters)) {
  idx <- g48$seurat_clusters[[i]]
  name <- new.cluster.ids[[idx]]
  x_1[i] <- name
}
g48 <- AddMetaData(object = g48, metadata = x_1, col.name = "celltype.anlas")

### load and filter reference dataset for 72h
g72 <- readRDS(data.file$g72h)
g72$merge.ident <- "72h"
g72$replicate <- "rep1"
g72 <- AddMetaData(object = g72, metadata = gsub(" ","",paste(g72$replicate,"_",g72$merge.ident)), col.name = "integrate.ident")

# assign name to clusters
new.cluster.ids <- c("Neural & mesodermal (posterior)", "Mesodermal trunk and neural",
                     "Mesodermal (paraxial)", "Mesodermal", "Mesodermal (posterior)",
                     "Pluripotent", "Early diff (SMAD, ERK activ)", "Mesodermal (mesench & cardiac)", "Endodermal", "Early vascular & endo")
names(new.cluster.ids) <- levels(g72)
#g72 <- RenameIdents(g72, new.cluster.ids)

x_1 <- c()
for (i in 1 : length(g72$seurat_clusters)) {
  idx <- g72$seurat_clusters[[i]]
  name <- new.cluster.ids[[idx]]
  x_1[i] <- name
}
g72 <- AddMetaData(object = g72, metadata = x_1, col.name = "celltype.anlas")

### merge data and split according to replicate and timepoint

anlas <- merge(g00, y=c(g24,g48,g72), add.cell.ids = c("h00", "h24", "h48", "h72"))
remove(g00,g24,g48,g72)

### save dataset
save(anlas, file = paste0(folder.RData,"anlas_data.RData"))

#####
# load and filter reference dataset
#####

argelaguet <- readRDS(data.file$arg)
argelaguet
# filter out genes with less than 10 total counts
argelaguet <- argelaguet[rowSums(counts(argelaguet))>10,]
argelaguet
# filter out cells that didn't pass rnaQC
argelaguet <- argelaguet[,as.logical(colData(argelaguet)$pass_rnaQC)]
argelaguet
# add ensembl id as rowData entry and make gene symbol the row name
rowData(argelaguet)$ens_id <- rownames(argelaguet)
rownames(argelaguet) <- rowData(argelaguet)$symbol
argelaguet
# remove mito genes
argelaguet <- argelaguet[grep('mt-',rownames(argelaguet),invert=T),]
argelaguet

argelaguet <- as.Seurat(argelaguet)#, counts = "counts", data = "counts")

save(argelaguet, file = paste0(folder.RData,"argelaguet_data.RData"))

### clean memory

remove(list=ls())
gc()

#####
# load anlas WT data
#####
anlas_WT <- readRDS(data.file$g48hWT)
anlas_WT$merge.ident <- "WT_48h"
anlas_WT$replicate <- "rep1"

# assign name to clusters
# new.cluster.ids <- c("Early diff (neural priming)", 
#                     "Pluripotent", "Early diff (cyto upreg)", "Primed pluripotent", 
#                     "Early diff (translation reg)", "Primed pluripotent (cyto upreg)", "Primed pluripotent (splicing reg)")
# names(new.cluster.ids) <- levels(g48wt)
# g48wt <- RenameIdents(g48wt, new.cluster.ids)

# x_1 <- c()
# for (i in 1 : length(g48wt$seurat_clusters)) {
#   idx <- g48wt$seurat_clusters[[i]]
#   name <- new.cluster.ids[[idx]]
#   x_1[i] <- name
# }
# g48wt <- AddMetaData(object = g48wt, metadata = x_1, col.name = "celltype.anlas")

anlas_WT <- AddMetaData(object = anlas_WT, metadata = paste0(anlas_WT$replicate,"_",anlas_WT$merge.ident), col.name = "integrate.ident")

### save dataset
save(anlas_WT, file = paste0(folder.RData,"anlas_WT_data.RData"))
