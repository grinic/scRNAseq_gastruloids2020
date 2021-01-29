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

if(Sys.info()['nodename']=='PCBA-TRIVEDI02'){
  folder <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\data\\anlas\\"
} else if(Sys.info()['nodename']=='PCBA-TRIVEDI03'){
  folder <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\data\\anlas\\"
}


### 24h
anlas24 <- readRDS(paste0(folder,"g_24h.rds"))

outFolder <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\data\\for_scvelo\\24h\\"
seurat.object <- anlas24[,anlas24$replicate=='rep1']

write.csv(Cells(seurat.object), file = paste0(outFolder,"cellID_obs.csv"), row.names = FALSE)
write.csv(Embeddings(seurat.object, reduction = "umap"), file = paste0(outFolder,"cell_embeddings.csv"))
write.csv(seurat.object@meta.data$seurat_clusters, file = paste0(outFolder,"clusters.csv"))

### 48h
anlas48 <- readRDS(paste0(folder,"g_48h.rds"))

outFolder <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\data\\for_scvelo\\48h\\"
seurat.object <- anlas48[,anlas48$replicate=='rep1']

write.csv(Cells(seurat.object), file = paste0(outFolder,"cellID_obs.csv"), row.names = FALSE)
write.csv(Embeddings(seurat.object, reduction = "umap"), file = paste0(outFolder,"cell_embeddings.csv"))
write.csv(seurat.object@meta.data$seurat_clusters, file = paste0(outFolder,"clusters.csv"))

### 72h
anlas72 <- readRDS(paste0(folder,"g_72h.rds"))

outFolder <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\data\\for_scvelo\\72h\\"
seurat.object <- anlas72

write.csv(Cells(seurat.object), file = paste0(outFolder,"cellID_obs.csv"), row.names = FALSE)
write.csv(Embeddings(seurat.object, reduction = "umap"), file = paste0(outFolder,"cell_embeddings.csv"))
write.csv(seurat.object@meta.data$seurat_clusters, file = paste0(outFolder,"clusters.csv"))

