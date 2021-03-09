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
  folder <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\raw_data\\gastruloids\\"
} else if(Sys.info()['nodename']=='PCBA-TRIVEDI03'){
  folder <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\raw_data\\gastruloids\\"
}


### 24h
anlas <- readRDS(paste0(folder,"all_data.combined_2.rds"))

meta_anlas <- anlas@meta.data

cl.name <- c()
cl.name.simplified <- c()
for(i in 1:length(meta_anlas$seurat_clusters)){
  old <- meta_anlas$seurat_clusters[i]
  if(old==0){
    cl.name[i] <- "Mesodermal (Bra+)"
    cl.name.simplified[i] <- "Mesodermal I"
  } else if(old==1){
    cl.name[i] <- "Pluripotent (naive-primed)"
    cl.name.simplified[i] <- "Pluripotent II"
  } else if(old==2){
    cl.name[i] <- "Mesodermal (paraxial, intermediate, lateral plate)"
    cl.name.simplified[i] <- "Mesodermal II"
  } else if(old==3){
    cl.name[i] <- "Early differentiated (neural, mesodermal)"
    cl.name.simplified[i] <- "Early diff II"
  } else if(old==4){
    cl.name[i] <- "Mesendodermal (Eomes+)"
    cl.name.simplified[i] <- "Mesodermal III"
  } else if(old==5){
    cl.name[i] <- "Primed pluripotent, Early differentiated"
    cl.name.simplified[i] <- "Primed pluripotent"
  } else if(old==6){
    cl.name[i] <- "Pluripotent (naive)"
    cl.name.simplified[i] <- "Pluripotent I"
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
    cl.name.simplified[i] <- "Early diff I"
  }
}
meta_anlas$celltype <- cl.name
anlas$celltype_simplified <- cl.name.simplified

outFolder <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\data\\for_scvelo\\"
seurat.object <- anlas[,(anlas$orig.ident=='scRNA_gast24h')|(anlas$orig.ident=='scRNA_gast48h')|(anlas$orig.ident=='scRNA_gast72h')]

write.csv(Cells(seurat.object), file = paste0(outFolder,"cellID_obs.csv"), row.names = FALSE)
write.csv(Embeddings(seurat.object, reduction = "umap"), file = paste0(outFolder,"cell_embeddings.csv"))
write.csv(seurat.object@meta.data$celltype_simplified, file = paste0(outFolder,"clusters.csv"))
write.csv(seurat.object@meta.data$orig.ident, file = paste0(outFolder, "orig_ident.csv"))

