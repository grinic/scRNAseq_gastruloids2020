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

################################## plot

rm(list = ls())
gc()
if(Sys.info()['nodename']=='PCBA-TRIVEDI02'){
  folder.RData <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\data\\"
  outFolder <- "Y:\\Nicola_Gritti\\analysis_code\\scRNAseq_Gastruloids\\new_codes\\results\\integration\\pijuan6.5-7.5_anlas\\"
} else if(Sys.info()['nodename']=='PCBA-TRIVEDI03'){
  folder.RData <- "C:\\Users\\nicol\\OneDrive\\Desktop\\scRNAseq_Gastruloids\\data\\"
  outFolder <- "C:\\Users\\nicol\\OneDrive\\Desktop\\scRNAseq_Gastruloids\\results\\integration\\pijuan6.5-7.5_anlas\\"
}

load(paste0(outFolder,"corrected_data.RData"))
load(paste0(outFolder,"corrected_data_umap.RData"))

correct_df <- as.data.frame(correct.umap$layout)
correct_df_pijuan <- correct_df[pijuan_rows,]
correct_df_anlas <- correct_df[anlas_rows,]

# define colors

stage_color_Publication <- function(...){
  library(scales)
  discrete_scale("colour", "Publication",
                 manual_pal(values = c(
                   "#000000", "#404040", "#707070",
                   "#330000", "#660000", "#990000", "#cc0000", "#ff0000", "#ff3232", "#ff4c4c", "#ff7f7f", "#ffb2b2", "#ffcccc")))
}

lineage_anlas_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour", "Publication",
                 manual_pal(values = c(
                   "#000000", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                   "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                   "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                   "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                   "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                   "#372101", "#FFB500", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                   "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                   "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                   "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                   "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                   "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                   "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                   "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")), ...)
}

lineage_pijuan_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour", "Publication",
                 manual_pal(values = c(
                   "#000000", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                   "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                   "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                   "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                   "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                   "#372101", "#FFB500", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                   "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                   "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                   "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                   "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                   "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                   "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                   "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")), ...)
}

lineage_color_Publication <- function(...){
  library(scales)
  discrete_scale("colour", "Publication",
                 manual_pal(values = c(
                   "#000000", "#404040", "#707070", 
                   "#008941", "#006FA6", "#A30059", "#7A4900", 
                   "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF")))
}

library(cowplot)

############################################################################################################

p1 <- ggplot(correct_df_pijuan, aes(x=V1, y=V2, col=c(meta_pijuan$celltype.pijuan))) + 
  geom_point(size=0.01) + 
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black')) + 
  theme(legend.title=element_blank()) +
#  xlim(-16,16) + ylim(-16,16) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(y= "UMAP 2", x = "UMAP 1") + 
  lineage_pijuan_colour_Publication()

p2 <- ggplot(correct_df_anlas, aes(x=V1, y=V2, col=c(meta_anlas$celltype.anlas))) + 
  geom_point(size=0.01) + 
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black')) + 
  theme(legend.title=element_blank()) +
#  xlim(-16,16) + ylim(-16,16) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(y= "UMAP 2", x = "UMAP 1") + 
  lineage_anlas_colour_Publication()

p1+p2

p3 <- ggplot(correct_df, aes(x=V1, y=V2, col=c(meta_pijuan$stage, meta_anlas$merge.ident))) + 
  geom_point(size=0.01) + 
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black')) + 
  theme(legend.title=element_blank()) +
#  xlim(-16,16) + ylim(-16,16) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(y= "UMAP 2", x = "UMAP 1") + 
  stage_color_Publication()

p4 <- ggplot(correct_df, aes(x=V1, y=V2, col=c(meta_pijuan$celltype.general, meta_anlas$merge.ident))) + 
  geom_point(size=0.01) + 
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black')) + 
  theme(legend.title=element_blank()) +
#  xlim(-16,16) + ylim(-16,16) +
  guides(colour = guide_legend(override.aes = list(size=5))) + 
  labs(y= "UMAP 2", x = "UMAP 1") + 
  lineage_color_Publication()

p3 + p4

p5 <- ggplot(correct_df, aes(x=V1, y=V2, col=c(meta_pijuan$celltype.general, meta_anlas$celltype.general))) + 
  geom_point(size=0.01) + 
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black')) + 
  theme(legend.title=element_blank()) +
#  xlim(-16,16) + ylim(-16,16) +
  guides(colour = guide_legend(override.aes = list(size=5))) + 
#  labs(y= "UMAP 2", x = "UMAP 1") + 
  lineage_pijuan_colour_Publication()

p5

#########################################################################################################

ggsave(
  paste0(outFolder, "p1_pijuanCelltype.pdf"),
  plot = p1,
  device = NULL,
  path = NULL,
#  scale = 1,
  width = 10,
  height = 5,
#  units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = FALSE,
)

ggsave(
  paste0(outFolder, "p2_anlasCelltype.pdf"),
  plot = p2,
  device = NULL,
  path = NULL,
  #  scale = 1,
  width = 10,
  height = 5,
  #  units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = FALSE,
)

ggsave(
  paste0(outFolder, "p3_pijuanStage_anlasStage.pdf"),
  plot = p3,
  device = NULL,
  path = NULL,
  #  scale = 1,
  width = 6,
  height = 5,
  #  units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = FALSE,
)

ggsave(
  paste0(outFolder, "p4_pijuanCelltype_anlasStage.pdf"),
  plot = p4,
  device = NULL,
  path = NULL,
  #  scale = 1,
  width = 6,
  height = 5,
  #  units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = FALSE,
)

ggsave(
  paste0(outFolder, "p5_pijuanCelltype_anlasCelltype.pdf"),
  plot = p5,
  device = NULL,
  path = NULL,
  #  scale = 1,
  width = 6,
  height = 5,
  #  units = c("in", "cm", "mm"),
  dpi = 300,
  limitsize = FALSE,
)
