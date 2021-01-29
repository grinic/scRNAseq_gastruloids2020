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

# load dataset
if(Sys.info()['nodename']=='PCBA-TRIVEDI02'){
  folder.RData <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\data\\"
  outFolder <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\results\\integration\\pijuan_BraGFP_WT\\"
} else if(Sys.info()['nodename']=='PCBA-TRIVEDI03'){
  folder.RData <- "C:\\Users\\nicol\\OneDrive\\Desktop\\scRNAseq_Gastruloids\\data\\"
  outFolder <- "C:\\Users\\nicol\\OneDrive\\Desktop\\scRNAseq_Gastruloids\\results\\integration\\ppijuan_BraGFP_WT\\"
}
load(paste0(outFolder,"corrected_data.RData"))

custom.config = umap.defaults
custom.config$n_neighbors = 20
custom.config$min_dist = 0.7

set.seed(42)
correct.umap <- umap(correct, config=custom.config)

write.table(correct.umap$layout, file= paste0(outFolder,"integration_umap.csv"), sep=",")
save(correct.umap, file= paste0(outFolder,"corrected_data_umap.RData"))

