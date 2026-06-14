rm(list = ls())
Sys.setenv(R_MAX_NUM_DLLS = 999)
options(stringsAsFactors = FALSE)

set.seed(1234)

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggpubr)
library(reshape2)

sample_dirs <- list(
  age12 = "./data/age12/",
  age32 = "./data/age32/",
  age34 = "./data/age34/",
  age42 = "./data/age42/"
)

sce_list <- lapply(names(sample_dirs), function(s) {
  counts <- Read10X(sample_dirs[[s]])
  CreateSeuratObject(
    counts = counts,
    project = s,
    min.cells = 3,
    min.features = 200
  )
})

names(sce_list) <- names(sample_dirs)

sce <- merge(
  sce_list[[1]],
  y = sce_list[-1],
  project = "Ovary_merged"
)

sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")

sce <- subset(
  sce,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)

sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, nfeatures = 2000)
sce <- ScaleData(sce)
sce <- RunPCA(sce)

ElbowPlot(sce, ndims = 30)

sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.5)
sce <- RunUMAP(sce, dims = 1:15)

DotPlot(
  sce,
  features = c(
    "COL1A1","DCN","LUM",
    "PECAM1","VWF",
    "ACTA2","RGS5",
    "PTPRC","LST1","TYROBP",
    "FOXL2","AMH","CYP19A1",
    "ZP1","ZP2","ZP3","GDF9","BMP15"
  )
) + RotatedAxis()

new.cluster.ids <- c(
  "Theca & stroma","Theca & stroma","Theca & stroma",
  "Theca & stroma","Endothelial cell","Endothelial cell",
  "Endothelial cell","Smooth muscle","Theca & stroma",
  "Endothelial cell","Granulosa cell","Theca & stroma",
  "Macrophage","Endothelial cell","Theca & stroma"
)

names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)

sce$celltype <- Idents(sce)

macrophage <- subset(sce, subset = celltype == "Macrophage")

pyro_genes <- c(
  "CASP1","GSDMD","IL1B","NLRP3","PYCARD",
  "IL18","NLRP1"
)

macrophage <- AddModuleScore(
  macrophage,
  features = list(pyro_genes),
  name = "pyroptosis"
)

macrophage$group <- macrophage$orig.ident

ggplot(macrophage@meta.data,
       aes(x = group, y = pyroptosis1, fill = group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_classic() +
  xlab("Group") +
  ylab("Pyroptosis score")
