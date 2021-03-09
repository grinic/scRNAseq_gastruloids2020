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

data_location <- "server"
data_location <- "local"
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
load(paste0(outFolder,"corrected_data.RData"))

col2hex <- function(rcolor) {
  cols <- c(rgb(t(col2rgb(colors())), maxColorValue=255))
  colCheck <- colors() %in% rcolor
  cols[colCheck]
}

custom.config = umap.defaults
custom.config$n_neighbors = 20
custom.config$min_dist = 0.7

set.seed(42)
correct.umap <- umap(correct, config=custom.config)

write.table(correct.umap$layout, file= paste0(outFolder,"integration_umap.csv"), sep=",")
save(correct.umap, file= paste0(outFolder,"corrected_data_umap.RData"))

