.libPaths("/lustre/user/liclab/liocean/maosl/Rpackages")
setwd("/lustre/user/liclab/liocean/maosl/CR/")

library(Seurat)
library(ggplot2)
library(Matrix)
library(ggsignif)
library(ggsci)
library(dplyr)
library(DoubletFinder)
library(stringr)

# preprocessing
#read in data
for(directory in dir("./data/", pattern = "GSM")){
  print(directory)
  tissue <- unlist(str_split(directory, pattern = "_"))[2]
  gender <- unlist(str_split(tissue, pattern = "-"))[2]
  age <- unlist(str_split(tissue, pattern = "-"))[3]
  sample.counts <- Read10X(data.dir = paste0("./data/", directory))
  sample <- CreateSeuratObject(counts = sample.counts, min.cells = 5)
  sample@meta.data$gender<- gender
  sample@meta.data$age<- age
  sample@meta.data$sample <- paste0(gender, ".", age)
  sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^Mt-")
  #VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  # 这行要测试一下
  if (startsWith(tissue, "Kidney")){
    sample <- subset(sample, subset = nFeature_RNA > 500 & percent.mt < 50)
  }else{
    sample <- subset(sample, subset = nFeature_RNA > 500 & percent.mt < 10)
  }
  
  sample <- NormalizeData(object = sample)
  sample <- FindVariableFeatures(object = sample)
  sample <- ScaleData(object = sample)
  sample <- RunPCA(object = sample)
  sample <- FindNeighbors(object = sample)
  sample <- FindClusters(object = sample)
  sample <- RunTSNE(object = sample)
  p <- DimPlot(object = sample, reduction = "tsne")
  ggsave(paste0("./figure/tsne/", tissue, ".pdf"), p,  height = 5, width = 6)
  #doublet find
  sweep.res.list <- paramSweep_v3(sample, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_value <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  annotations <- sample@meta.data$res.0.6
  homotypic.prop <- modelHomotypic(annotations)  
  nExp_poi <- round(0.031*length(sample@meta.data$orig.ident)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  pN_value <- 0.25
  pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
  sample <- doubletFinder_v3(sample, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
  sample <- doubletFinder_v3(sample, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
  sample@meta.data$Doublet <- sample[[paste0("DF.classifications_",pN_value,"_",pK_value,'_',nExp_poi)]]
  saveRDS(sample, paste0("./CR_rds/", tissue, ".rds"))
}


# integration
# BAT
BAT_F_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/BAT-F-CR.rds")
BAT_F_CR$gender <- "F"
BAT_F_CR$age <- "CR"
BAT_F_CR$sample <- "F.CR"
BAT_F_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/BAT-F-O.rds")
BAT_F_O$gender <- "F"
BAT_F_O$age <- "O"
BAT_F_O$sample <- "F.O"
BAT_F_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/BAT-F-Y.rds")
BAT_F_Y$gender <- "F"
BAT_F_Y$age <- "Y"
BAT_F_Y$sample <- "F.Y"
BAT_M_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/BAT-M-CR.rds")
BAT_M_CR$gender <- "M"
BAT_M_CR$age <- "CR"
BAT_M_CR$sample <- "M.CR"
BAT_M_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/BAT-M-O.rds")
BAT_M_O$gender <- "M"
BAT_M_O$age <- "O"
BAT_M_O$sample <- "M.O"
BAT_M_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/BAT-M-Y.rds")
BAT_M_Y$gender <- "M"
BAT_M_Y$age <- "Y"
BAT_M_Y$sample <- "M.Y"

BAT_list <- list(BAT_F_CR, BAT_F_O, BAT_F_Y,
                 BAT_M_CR, BAT_M_O, BAT_M_Y)

BAT_list <- lapply(X = BAT_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = BAT_list)
BAT_anchors <- FindIntegrationAnchors(object.list = BAT_list, anchor.features = features)
BAT <- IntegrateData(anchorset = BAT_anchors)


DefaultAssay(BAT) <- "integrated"

BAT <- ScaleData(BAT, verbose = FALSE)
BAT <- RunPCA(BAT, npcs = 30, verbose = FALSE)
BAT <- RunTSNE(BAT)
BAT <- RunUMAP(BAT, reduction = "pca", dims = 1:30)
BAT <- FindNeighbors(BAT, reduction = "pca", dims = 1:30)
BAT <- FindClusters(BAT, resolution = 0.5)
p1 <- DimPlot(BAT, reduction = "umap", group.by = "sample")
p2 <- DimPlot(BAT, reduction = "tsne", group.by = "sample")
p3 <- DimPlot(BAT, reduction = "umap", group.by = "age")
p4 <- DimPlot(BAT, reduction = "tsne", group.by = "age")


ggsave(plot = p1+p2, filename = "./figure/BAT/BAT_integral_1.pdf", width = 12, height = 5)
ggsave(plot = p4+p3, filename = "./figure/BAT/BAT_integral_2.pdf", width = 12, height = 5)

saveRDS(BAT, "./CR_integral_rds/BAT.rds")

# Aorta
Aorta_F_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Aorta-F-CR.rds")
Aorta_F_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Aorta-F-O.rds")
Aorta_F_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Aorta-F-Y.rds")
Aorta_M_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Aorta-M-CR.rds")
Aorta_M_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Aorta-M-O.rds")
Aorta_M_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Aorta-M-Y.rds")

Aorta_list <- list(Aorta_F_CR, Aorta_F_O, Aorta_F_Y,
                   Aorta_M_CR, Aorta_M_O, Aorta_M_Y)

Aorta_list <- lapply(X = Aorta_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = Aorta_list)
Aorta_anchors <- FindIntegrationAnchors(object.list = Aorta_list, anchor.features = features)
Aorta <- IntegrateData(anchorset = Aorta_anchors)


DefaultAssay(Aorta) <- "integrated"

Aorta <- ScaleData(Aorta, verbose = FALSE)
Aorta <- RunPCA(Aorta, npcs = 30, verbose = FALSE)
Aorta <- RunTSNE(Aorta)
Aorta <- RunUMAP(Aorta, reduction = "pca", dims = 1:30)
Aorta <- FindNeighbors(Aorta, reduction = "pca", dims = 1:30)
Aorta <- FindClusters(Aorta, resolution = 0.5)
p1 <- DimPlot(Aorta, reduction = "umap", group.by = "sample")
p2 <- DimPlot(Aorta, reduction = "tsne", group.by = "sample")
p3 <- DimPlot(Aorta, reduction = "umap", group.by = "age")
p4 <- DimPlot(Aorta, reduction = "tsne", group.by = "age")

ggsave(plot = p1+p2, filename = "./figure/Aorta/Aorta_integral_1.pdf", width = 12, height = 5)
ggsave(plot = p4+p3, filename = "./figure/Aorta/Aorta_integral_2.pdf", width = 12, height = 5)

saveRDS(Aorta, "./CR_integral_rds/Aorta.rds")
rm(list = ls())

# BM
BM_F_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/BM-F-CR.rds")
BM_F_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/BM-F-O.rds")
BM_F_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/BM-F-Y.rds")
BM_M_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/BM-M-CR.rds")
BM_M_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/BM-M-O.rds")
BM_M_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/BM-M-Y.rds")

BM_list <- list(BM_F_CR, BM_F_O, BM_F_Y,
                BM_M_CR, BM_M_O, BM_M_Y)

BM_list <- lapply(X = BM_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = BM_list)
BM_anchors <- FindIntegrationAnchors(object.list = BM_list, anchor.features = features)
BM <- IntegrateData(anchorset = BM_anchors)


DefaultAssay(BM) <- "integrated"

BM <- ScaleData(BM, verbose = FALSE)
BM <- RunPCA(BM, npcs = 30, verbose = FALSE)
BM <- RunTSNE(BM)
BM <- RunUMAP(BM, reduction = "pca", dims = 1:30)
BM <- FindNeighbors(BM, reduction = "pca", dims = 1:30)
BM <- FindClusters(BM, resolution = 0.5)
p1 <- DimPlot(BM, reduction = "umap", group.by = "sample")
p2 <- DimPlot(BM, reduction = "tsne", group.by = "sample")
p3 <- DimPlot(BM, reduction = "umap", group.by = "age")
p4 <- DimPlot(BM, reduction = "tsne", group.by = "age")

ggsave(plot = p1+p2, filename = "./figure/BM/BM_integral_1.pdf", width = 12, height = 5)
ggsave(plot = p4+p3, filename = "./figure/BM/BM_integral_2.pdf", width = 12, height = 5)

saveRDS(BM, "./CR_integral_rds/BM.rds")
rm(list = ls())


# Kidney
Kidney_F_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Kidney-F-CR.rds")
Kidney_F_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Kidney-F-O.rds")
Kidney_F_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Kidney-F-Y.rds")
Kidney_M_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Kidney-M-CR.rds")
Kidney_M_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Kidney-M-O.rds")
Kidney_M_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Kidney-M-Y.rds")

Kidney_list <- list(Kidney_F_CR, Kidney_F_O, Kidney_F_Y,
                    Kidney_M_CR, Kidney_M_O, Kidney_M_Y)

Kidney_list <- lapply(X = Kidney_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = Kidney_list)
Kidney_anchors <- FindIntegrationAnchors(object.list = Kidney_list, anchor.features = features)
Kidney <- IntegrateData(anchorset = Kidney_anchors)


DefaultAssay(Kidney) <- "integrated"

Kidney <- ScaleData(Kidney, verbose = FALSE)
Kidney <- RunPCA(Kidney, npcs = 30, verbose = FALSE)
Kidney <- RunTSNE(Kidney)
Kidney <- RunUMAP(Kidney, reduction = "pca", dims = 1:30)
Kidney <- FindNeighbors(Kidney, reduction = "pca", dims = 1:30)
Kidney <- FindClusters(Kidney, resolution = 0.5)
p1 <- DimPlot(Kidney, reduction = "umap", group.by = "sample")
p2 <- DimPlot(Kidney, reduction = "tsne", group.by = "sample")
p3 <- DimPlot(Kidney, reduction = "umap", group.by = "age")
p4 <- DimPlot(Kidney, reduction = "tsne", group.by = "age")


ggsave(plot = p1+p2, filename = "./figure/Kidney/Kidney_integral_1.pdf", width = 12, height = 5)
ggsave(plot = p4+p3, filename = "./figure/Kidney/Kidney_integral_2.pdf", width = 12, height = 5)

saveRDS(Kidney, "./CR_integral_rds/Kidney.rds")
rm(list = ls())


# Liver
Liver_F_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Liver-F-CR.rds")
Liver_F_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Liver-F-O.rds")
Liver_F_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Liver-F-Y.rds")
Liver_M_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Liver-M-CR.rds")
Liver_M_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Liver-M-O.rds")
Liver_M_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Liver-M-Y.rds")

Liver_list <- list(Liver_F_CR, Liver_F_O, Liver_F_Y,
                   Liver_M_CR, Liver_M_O, Liver_M_Y)

Liver_list <- lapply(X = Liver_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = Liver_list)
Liver_anchors <- FindIntegrationAnchors(object.list = Liver_list, anchor.features = features)
Liver <- IntegrateData(anchorset = Liver_anchors)


DefaultAssay(Liver) <- "integrated"

Liver <- ScaleData(Liver, verbose = FALSE)
Liver <- RunPCA(Liver, npcs = 30, verbose = FALSE)
Liver <- RunTSNE(Liver)
Liver <- RunUMAP(Liver, reduction = "pca", dims = 1:30)
Liver <- FindNeighbors(Liver, reduction = "pca", dims = 1:30)
Liver <- FindClusters(Liver, resolution = 0.5)
p1 <- DimPlot(Liver, reduction = "umap", group.by = "sample")
p2 <- DimPlot(Liver, reduction = "tsne", group.by = "sample")
p3 <- DimPlot(Liver, reduction = "umap", group.by = "age")
p4 <- DimPlot(Liver, reduction = "tsne", group.by = "age")

ggsave(plot = p1+p2, filename = "./figure/Liver/Liver_integral_1.pdf", width = 12, height = 5)
ggsave(plot = p4+p3, filename = "./figure/Liver/Liver_integral_2.pdf", width = 12, height = 5)

saveRDS(Liver, "./CR_integral_rds/Liver.rds")
rm(list = ls())


# Skin
Skin_F_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Skin-F-CR.rds")
Skin_F_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Skin-F-O.rds")
Skin_F_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Skin-F-Y.rds")
Skin_M_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Skin-M-CR.rds")
Skin_M_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Skin-M-O.rds")
Skin_M_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Skin-M-Y.rds")

Skin_list <- list(Skin_F_CR, Skin_F_O, Skin_F_Y,
                  Skin_M_CR, Skin_M_O, Skin_M_Y)

Skin_list <- lapply(X = Skin_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = Skin_list)
Skin_anchors <- FindIntegrationAnchors(object.list = Skin_list, anchor.features = features)
Skin <- IntegrateData(anchorset = Skin_anchors)


DefaultAssay(Skin) <- "integrated"

Skin <- ScaleData(Skin, verbose = FALSE)
Skin <- RunPCA(Skin, npcs = 30, verbose = FALSE)
Skin <- RunTSNE(Skin)
Skin <- RunUMAP(Skin, reduction = "pca", dims = 1:30)
Skin <- FindNeighbors(Skin, reduction = "pca", dims = 1:30)
Skin <- FindClusters(Skin, resolution = 0.5)
p1 <- DimPlot(Skin, reduction = "umap", group.by = "sample")
p2 <- DimPlot(Skin, reduction = "tsne", group.by = "sample")
p3 <- DimPlot(Skin, reduction = "umap", group.by = "age")
p4 <- DimPlot(Skin, reduction = "tsne", group.by = "age")

ggsave(plot = p1+p2, filename = "./figure/Skin/Skin_integral_1.pdf", width = 12, height = 5)
ggsave(plot = p4+p3, filename = "./figure/Skin/Skin_integral_2.pdf", width = 12, height = 5)

saveRDS(Skin, "./CR_integral_rds/Skin.rds")
rm(list = ls())


# WAT
WAT_F_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/WAT-F-CR.rds")
WAT_F_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/WAT-F-O.rds")
WAT_F_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/WAT-F-Y.rds")
WAT_M_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/WAT-M-CR.rds")
WAT_M_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/WAT-M-O.rds")
WAT_M_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/WAT-M-Y.rds")

WAT_list <- list(WAT_F_CR, WAT_F_O, WAT_F_Y,
                 WAT_M_CR, WAT_M_O, WAT_M_Y)

WAT_list <- lapply(X = WAT_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = WAT_list)
WAT_anchors <- FindIntegrationAnchors(object.list = WAT_list, anchor.features = features)
WAT <- IntegrateData(anchorset = WAT_anchors)


DefaultAssay(WAT) <- "integrated"

WAT <- ScaleData(WAT, verbose = FALSE)
WAT <- RunPCA(WAT, npcs = 30, verbose = FALSE)
WAT <- RunTSNE(WAT)
WAT <- RunUMAP(WAT, reduction = "pca", dims = 1:30)
WAT <- FindNeighbors(WAT, reduction = "pca", dims = 1:30)
WAT <- FindClusters(WAT, resolution = 0.5)
p1 <- DimPlot(WAT, reduction = "umap", group.by = "sample")
p2 <- DimPlot(WAT, reduction = "tsne", group.by = "sample")
p3 <- DimPlot(WAT, reduction = "umap", group.by = "age")
p4 <- DimPlot(WAT, reduction = "tsne", group.by = "age")

ggsave(plot = p1+p2, filename = "./figure/WAT/WAT_integral_1.pdf", width = 12, height = 5)
ggsave(plot = p4+p3, filename = "./figure/WAT/WAT_integral_2.pdf", width = 12, height = 5)

saveRDS(WAT, "./CR_integral_rds/WAT.rds")
rm(list = ls())


# Brain
Brain_F_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Brain-F-CR.rds")
Brain_F_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Brain-F-O.rds")
Brain_F_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Brain-F-Y.rds")
Brain_M_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Brain-M-CR.rds")
Brain_M_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Brain-M-O.rds")
Brain_M_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Brain-M-Y.rds")

Brain_list <- list(Brain_F_CR, Brain_F_O, Brain_F_Y,
                   Brain_M_CR, Brain_M_O, Brain_M_Y)

Brain_list <- lapply(X = Brain_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = Brain_list)
Brain_anchors <- FindIntegrationAnchors(object.list = Brain_list, anchor.features = features)
Brain <- IntegrateData(anchorset = Brain_anchors)


DefaultAssay(Brain) <- "integrated"

Brain <- ScaleData(Brain, verbose = FALSE)
Brain <- RunPCA(Brain, npcs = 30, verbose = FALSE)
Brain <- RunTSNE(Brain)
Brain <- RunUMAP(Brain, reduction = "pca", dims = 1:30)
Brain <- FindNeighbors(Brain, reduction = "pca", dims = 1:30)
Brain <- FindClusters(Brain, resolution = 0.5)
p1 <- DimPlot(Brain, reduction = "umap", group.by = "sample")
p2 <- DimPlot(Brain, reduction = "tsne", group.by = "sample")
p3 <- DimPlot(Brain, reduction = "umap", group.by = "age")
p4 <- DimPlot(Brain, reduction = "tsne", group.by = "age")

ggsave(plot = p1+p2, filename = "./figure/Brain/Brain_integral_1.pdf", width = 12, height = 5)
ggsave(plot = p4+p3, filename = "./figure/Brain/Brain_integral_2.pdf", width = 12, height = 5)

saveRDS(Brain, "./CR_integral_rds/Brain.rds")
rm(list = ls())


# Muscle
Muscle_F_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Muscle-F-CR.rds")
Muscle_F_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Muscle-F-O.rds")
Muscle_F_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Muscle-F-Y.rds")
Muscle_M_CR <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Muscle-M-CR.rds")
Muscle_M_O <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Muscle-M-O.rds")
Muscle_M_Y <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_rds/Muscle-M-Y.rds")

Muscle_list <- list(Muscle_F_CR, Muscle_F_O, Muscle_F_Y,
                    Muscle_M_CR, Muscle_M_O, Muscle_M_Y)

Muscle_list <- lapply(X = Muscle_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = Muscle_list)
Muscle_anchors <- FindIntegrationAnchors(object.list = Muscle_list, anchor.features = features)
Muscle <- IntegrateData(anchorset = Muscle_anchors)


DefaultAssay(Muscle) <- "integrated"

Muscle <- ScaleData(Muscle, verbose = FALSE)
Muscle <- RunPCA(Muscle, npcs = 30, verbose = FALSE)
Muscle <- RunTSNE(Muscle)
Muscle <- RunUMAP(Muscle, reduction = "pca", dims = 1:30)
Muscle <- FindNeighbors(Muscle, reduction = "pca", dims = 1:30)
Muscle <- FindClusters(Muscle, resolution = 0.5)
p1 <- DimPlot(Muscle, reduction = "umap", group.by = "sample")
p2 <- DimPlot(Muscle, reduction = "tsne", group.by = "sample")
p3 <- DimPlot(Muscle, reduction = "umap", group.by = "age")
p4 <- DimPlot(Muscle, reduction = "tsne", group.by = "age")

ggsave(plot = p1+p2, filename = "./figure/Muscle/Muscle_integral_1.pdf", width = 12, height = 5)
ggsave(plot = p4+p3, filename = "./figure/Muscle/Muscle_integral_2.pdf", width = 12, height = 5)

saveRDS(Muscle, "./CR_integral_rds/Muscle.rds")
rm(list = ls())

# calculate aging score
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE){
  library(plyr)
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  datac <- plyr::rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

aging_gene_glmnet_all_tissue <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/aging_gene_glmnet_all_tissue_add_50_corrected.rds")
rat_mouse_gene_transform <- readRDS("/lustre/user/liclab/liocean/maosl/CR/rat_mouse_gene_transform.rds")

# BAT
BAT <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_integral_rds/BAT.rds")

aging_genes_mouse_down <- c(as.character(aging_gene_glmnet_all_tissue[["facs-BAT"]][1:50, "down"]))
aging_genes_rat_down <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_down]
aging_genes_rat_down <- aging_genes_rat_down[aging_genes_rat_down %in% rownames(BAT@assays$RNA@counts)]
aging_genes_mouse_up <- c(as.character(aging_gene_glmnet_all_tissue[["facs-BAT"]][1:50, "up"]))
down_counts <- length(aging_genes_rat_down)
aging_genes_rat_up <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_up]
aging_genes_rat_up <- aging_genes_rat_up[aging_genes_rat_up %in% rownames(BAT@assays$RNA@counts)]
up_counts <- length(aging_genes_rat_up)
aging_genes_rat <- c(as.character(aging_genes_rat_down), as.character(aging_genes_rat_up))


weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(BAT@assays$RNA@counts[aging_genes_rat,]>0)/dim(BAT@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(BAT@assays$RNA@data[aging_genes_rat,]),
    center = T, scale = T))
aging.score = t(z_matrix[aging_genes_rat,]) %*% weights 

df_aging_score <- as.data.frame(cbind(BAT$age, aging.score))
colnames(df_aging_score) <- c("age", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$age <- factor(df_aging_score$age, 
                             levels = c("Y", "O", "CR"),
                             labels = c("Young", "Old", "CR"))

a <- summarySE(df_aging_score, measurevar="aging_score", groupvars=c("age"))

ggplot(a, aes(x=age, y=aging_score)) + 
  geom_errorbar(aes(ymin=aging_score-se, ymax=aging_score+se), width=.1) +
  geom_line() +
  geom_point()

ggplot(df_aging_score, aes(x = age, y = aging_score)) +
  #geom_line() +
  #geom_violin(aes(fill = age)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Young", "Old"), c("Old", "CR")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  scale_fill_npg() + 
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  labs(x="Age", y="Aging score", title="BAT")
p
ggsave("./figure/BAT/BAT_aging_score.pdf", width = 6, height = 4)
rm("BAT")


# Aorta 
Aorta <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_integral_rds/Aorta.rds")

# facs
aging_genes_mouse_down <- c(as.character(aging_gene_glmnet_all_tissue[["facs-Aorta"]][1:50, "down"]))
aging_genes_rat_down <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_down]
aging_genes_rat_down <- aging_genes_rat_down[aging_genes_rat_down %in% rownames(Aorta@assays$RNA@counts)]
down_counts <- length(aging_genes_rat_down)
aging_genes_mouse_up <- c(as.character(aging_gene_glmnet_all_tissue[["facs-Aorta"]][1:50, "up"]))
aging_genes_rat_up <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_up]
aging_genes_rat_up <- aging_genes_rat_up[aging_genes_rat_up %in% rownames(Aorta@assays$RNA@counts)]
up_counts <- length(aging_genes_rat_up)

aging_genes_rat <- as.character(c(as.character(aging_genes_rat_down), as.character(aging_genes_rat_up)))

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(Aorta@assays$RNA@counts[aging_genes_rat,]>0)/dim(Aorta@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(Aorta@assays$RNA@data[aging_genes_rat,]),
    center = T, scale = T))
aging.score = t(z_matrix[aging_genes_rat,]) %*% weights 

df_aging_score <- as.data.frame(cbind(Aorta$age, aging.score))
colnames(df_aging_score) <- c("age", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$age <- factor(df_aging_score$age, 
                             levels = c("Y", "O", "CR"),
                             labels = c("Young", "Old", "CR"))

p <- ggplot(df_aging_score, aes(x = age, y = aging_score)) +
  geom_violin(aes(fill = age)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Young", "Old"), c("Old", "CR")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  scale_fill_npg() + 
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  labs(x="Age",y="Aging score", title="Aorta")
ggsave("./figure/Aorta/Aorta_aging_score_facs.pdf", width = 6, height = 4)

# droplet
aging_genes_mouse_down <- c(as.character(aging_gene_glmnet_all_tissue[["droplet-Heart_and_Aorta"]][1:50, "down"]))
aging_genes_rat_down <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_down]
aging_genes_rat_down <- aging_genes_rat_down[aging_genes_rat_down %in% rownames(Aorta@assays$RNA@counts)]
down_counts <- length(aging_genes_rat_down)
aging_genes_mouse_up <- c(as.character(aging_gene_glmnet_all_tissue[["droplet-Heart_and_Aorta"]][1:50, "up"]))
aging_genes_rat_up <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_up]
aging_genes_rat_up <- aging_genes_rat_up[aging_genes_rat_up %in% rownames(Aorta@assays$RNA@counts)]
up_counts <- length(aging_genes_rat_up)

aging_genes_rat <- as.character(c(as.character(aging_genes_rat_down), as.character(aging_genes_rat_up)))

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(Aorta@assays$RNA@counts[aging_genes_rat,]>0)/dim(Aorta@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(Aorta@assays$RNA@data[aging_genes_rat,]),
    center = T, scale = T))
aging.score = t(z_matrix[aging_genes_rat,]) %*% weights 

df_aging_score <- as.data.frame(cbind(Aorta$age, aging.score))
colnames(df_aging_score) <- c("age", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$age <- factor(df_aging_score$age, 
                             levels = c("Y", "O", "CR"),
                             labels = c("Young", "Old", "CR"))

p <- ggplot(df_aging_score, aes(x = age, y = aging_score)) +
  geom_violin(aes(fill = age)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Young", "Old"), c("Old", "CR")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  scale_fill_npg() + 
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  labs(x="Age",y="Aging score", title = "Aorta")
ggsave("./figure/Aorta/Aorta_aging_score_droplet.pdf", width = 6, height = 4)

rm("Aorta")

# BM-bone marrow
BM <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_integral_rds/BM.rds")

# facs
aging_genes_mouse_down <- c(as.character(aging_gene_glmnet_all_tissue[["facs-Marrow"]][1:50, "down"]))
aging_genes_rat_down <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_down]
aging_genes_rat_down <- aging_genes_rat_down[aging_genes_rat_down %in% rownames(BM@assays$RNA@counts)]
down_counts <- length(aging_genes_rat_down)
aging_genes_mouse_up <- c(as.character(aging_gene_glmnet_all_tissue[["facs-Marrow"]][1:50, "up"]))
aging_genes_rat_up <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_up]
aging_genes_rat_up <- aging_genes_rat_up[aging_genes_rat_up %in% rownames(BM@assays$RNA@counts)]
up_counts <- length(aging_genes_rat_up)

aging_genes_rat <- as.character(c(as.character(aging_genes_rat_down), as.character(aging_genes_rat_up)))

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(BM@assays$RNA@counts[aging_genes_rat,]>0)/dim(BM@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(BM@assays$RNA@data[aging_genes_rat,]),
    center = T, scale = T))
aging.score = t(z_matrix[aging_genes_rat,]) %*% weights 

df_aging_score <- as.data.frame(cbind(BM$age, aging.score))
colnames(df_aging_score) <- c("age", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$age <- factor(df_aging_score$age, 
                             levels = c("Y", "O", "CR"),
                             labels = c("Young", "Old", "CR"))

p <- ggplot(df_aging_score, aes(x = age, y = aging_score)) +
  geom_violin(aes(fill = age)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Young", "Old"), c("Old", "CR")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  scale_fill_npg() + 
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  labs(x="Age",y="Aging score", title="Bone marrow")
p
ggsave("./figure/BM/BM_aging_score_facs.pdf", width = 6, height = 4)

# droplet
aging_genes_mouse_down <- c(as.character(aging_gene_glmnet_all_tissue[["droplet-Marrow"]][1:50, "down"]))
aging_genes_rat_down <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_down]
aging_genes_rat_down <- aging_genes_rat_down[aging_genes_rat_down %in% rownames(BM@assays$RNA@counts)]
down_counts <- length(aging_genes_rat_down)
aging_genes_mouse_up <- c(as.character(aging_gene_glmnet_all_tissue[["droplet-Marrow"]][1:50, "up"]))
aging_genes_rat_up <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_up]
aging_genes_rat_up <- aging_genes_rat_up[aging_genes_rat_up %in% rownames(BM@assays$RNA@counts)]
up_counts <- length(aging_genes_rat_up)

aging_genes_rat <- as.character(c(as.character(aging_genes_rat_down), as.character(aging_genes_rat_up)))

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(BM@assays$RNA@counts[aging_genes_rat,]>0)/dim(BM@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(BM@assays$RNA@data[aging_genes_rat,]),
    center = T, scale = T))
aging.score = t(z_matrix[aging_genes_rat,]) %*% weights 

df_aging_score <- as.data.frame(cbind(BM$age, aging.score))
colnames(df_aging_score) <- c("age", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$age <- factor(df_aging_score$age, 
                             levels = c("Y", "O", "CR"),
                             labels = c("Young", "Old", "CR"))

p <- ggplot(df_aging_score, aes(x = age, y = aging_score)) +
  geom_violin(aes(fill = age)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Young", "Old"), c("Old", "CR")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  scale_fill_npg() + 
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) + 
  labs(x="Age",y="Aging score", title="Bone marrow")
p
ggsave("./figure/BM/BM_aging_score_droplet.pdf", width = 6, height = 4)

rm("BM")

# Kidney
Kidney <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_integral_rds/Kidney.rds")

# facs
aging_genes_mouse_down <- c(as.character(aging_gene_glmnet_all_tissue[["facs-Kidney"]][1:50, "down"]))
aging_genes_rat_down <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_down]
aging_genes_rat_down <- aging_genes_rat_down[aging_genes_rat_down %in% rownames(Kidney@assays$RNA@counts)]
down_counts <- length(aging_genes_rat_down)
aging_genes_mouse_up <- c(as.character(aging_gene_glmnet_all_tissue[["facs-Kidney"]][1:50, "up"]))
aging_genes_rat_up <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_up]
aging_genes_rat_up <- aging_genes_rat_up[aging_genes_rat_up %in% rownames(Kidney@assays$RNA@counts)]
up_counts <- length(aging_genes_rat_up)

aging_genes_rat <- as.character(c(as.character(aging_genes_rat_down), as.character(aging_genes_rat_up)))

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(Kidney@assays$RNA@counts[aging_genes_rat,]>0)/dim(Kidney@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(Kidney@assays$RNA@data[aging_genes_rat,]),
    center = T, scale = T))
aging.score = t(z_matrix[aging_genes_rat,]) %*% weights 

df_aging_score <- as.data.frame(cbind(Kidney$age, aging.score))
colnames(df_aging_score) <- c("age", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$age <- factor(df_aging_score$age, 
                             levels = c("Y", "O", "CR"),
                             labels = c("Young", "Old", "CR"))

p <- ggplot(df_aging_score, aes(x = age, y = aging_score)) +
  geom_violin(aes(fill = age)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Young", "Old"), c("Old", "CR")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  scale_fill_npg() + 
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) + 
  labs(x="Age",y="Aging score", title="Kidney")
p
ggsave("./figure/Kidney/Kidney_aging_score_facs.pdf", width = 6, height = 4)

# droplet
aging_genes_mouse_down <- c(as.character(aging_gene_glmnet_all_tissue[["droplet-Kidney"]][1:50, "down"]))
aging_genes_rat_down <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_down]
aging_genes_rat_down <- aging_genes_rat_down[aging_genes_rat_down %in% rownames(Kidney@assays$RNA@counts)]
down_counts <- length(aging_genes_rat_down)
aging_genes_mouse_up <- c(as.character(aging_gene_glmnet_all_tissue[["droplet-Kidney"]][1:50, "up"]))
aging_genes_rat_up <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_up]
aging_genes_rat_up <- aging_genes_rat_up[aging_genes_rat_up %in% rownames(Kidney@assays$RNA@counts)]
up_counts <- length(aging_genes_rat_up)

aging_genes_rat <- as.character(c(as.character(aging_genes_rat_down), as.character(aging_genes_rat_up)))

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(Kidney@assays$RNA@counts[aging_genes_rat,]>0)/dim(Kidney@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(Kidney@assays$RNA@data[aging_genes_rat,]),
    center = T, scale = T))
aging.score = t(z_matrix[aging_genes_rat,]) %*% weights 

df_aging_score <- as.data.frame(cbind(Kidney$age, aging.score))
colnames(df_aging_score) <- c("age", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$age <- factor(df_aging_score$age, 
                             levels = c("Y", "O", "CR"),
                             labels = c("Young", "Old", "CR"))

p <- ggplot(df_aging_score, aes(x = age, y = aging_score)) +
  geom_violin(aes(fill = age)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Young", "Old"), c("Old", "CR")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  scale_fill_npg() + 
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  labs(x="Age",y="Aging score", title="Kideny")
p
ggsave("./figure/Kidney/Kidney_aging_score_droplet.pdf", width = 6, height = 4)

rm("Kidney")


# Liver
Liver <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_integral_rds/Liver.rds")

# facs
aging_genes_mouse_down <- c(as.character(aging_gene_glmnet_all_tissue[["facs-Liver"]][1:50, "down"]))
aging_genes_rat_down <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_down]
aging_genes_rat_down <- aging_genes_rat_down[aging_genes_rat_down %in% rownames(Liver@assays$RNA@counts)]
down_counts <- length(aging_genes_rat_down)
aging_genes_mouse_up <- c(as.character(aging_gene_glmnet_all_tissue[["facs-Liver"]][1:50, "up"]))
aging_genes_rat_up <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_up]
aging_genes_rat_up <- aging_genes_rat_up[aging_genes_rat_up %in% rownames(Liver@assays$RNA@counts)]
up_counts <- length(aging_genes_rat_up)

aging_genes_rat <- as.character(c(as.character(aging_genes_rat_down), as.character(aging_genes_rat_up)))

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(Liver@assays$RNA@counts[aging_genes_rat,]>0)/dim(Liver@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(Liver@assays$RNA@data[aging_genes_rat,]),
    center = T, scale = T))
aging.score = t(z_matrix[aging_genes_rat,]) %*% weights 

df_aging_score <- as.data.frame(cbind(Liver$age, aging.score))
colnames(df_aging_score) <- c("age", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$age <- factor(df_aging_score$age, 
                             levels = c("Y", "O", "CR"),
                             labels = c("Young", "Old", "CR"))

p <- ggplot(df_aging_score, aes(x = age, y = aging_score)) +
  geom_violin(aes(fill = age)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Young", "Old"), c("Old", "CR")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  scale_fill_npg() + 
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  labs(x="Age",y="Aging score", title="Liver")
p
ggsave("./figure/Liver/Liver_aging_score_facs.pdf", width = 6, height = 4)

# droplet
aging_genes_mouse_down <- c(as.character(aging_gene_glmnet_all_tissue[["droplet-Liver"]][1:50, "down"]))
aging_genes_rat_down <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_down]
aging_genes_rat_down <- aging_genes_rat_down[aging_genes_rat_down %in% rownames(Liver@assays$RNA@counts)]
down_counts <- length(aging_genes_rat_down)
aging_genes_mouse_up <- c(as.character(aging_gene_glmnet_all_tissue[["droplet-Liver"]][1:50, "up"]))
aging_genes_rat_up <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_up]
aging_genes_rat_up <- aging_genes_rat_up[aging_genes_rat_up %in% rownames(Liver@assays$RNA@counts)]
up_counts <- length(aging_genes_rat_up)

aging_genes_rat <- as.character(c(as.character(aging_genes_rat_down), as.character(aging_genes_rat_up)))

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(Liver@assays$RNA@counts[aging_genes_rat,]>0)/dim(Liver@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(Liver@assays$RNA@data[aging_genes_rat,]),
    center = T, scale = T))
aging.score = t(z_matrix[aging_genes_rat,]) %*% weights 

df_aging_score <- as.data.frame(cbind(Liver$age, aging.score))
colnames(df_aging_score) <- c("age", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$age <- factor(df_aging_score$age, 
                             levels = c("Y", "O", "CR"),
                             labels = c("Young", "Old", "CR"))

p <- ggplot(df_aging_score, aes(x = age, y = aging_score)) +
  geom_violin(aes(fill = age)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Young", "Old"), c("Old", "CR")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  scale_fill_npg() + 
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  labs(x="Age",y="Aging score", title="Liver")
p
ggsave("./figure/Liver/Liver_aging_score_droplet.pdf", width = 6, height = 4)

rm("Liver")


# Skin
Skin <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_integral_rds/Skin.rds")

# facs
aging_genes_mouse_down <- c(as.character(aging_gene_glmnet_all_tissue[["facs-Skin"]][1:50, "down"]))
aging_genes_rat_down <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_down]
aging_genes_rat_down <- aging_genes_rat_down[aging_genes_rat_down %in% rownames(Skin@assays$RNA@counts)]
down_counts <- length(aging_genes_rat_down)
aging_genes_mouse_up <- c(as.character(aging_gene_glmnet_all_tissue[["facs-Skin"]][1:50, "up"]))
aging_genes_rat_up <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_up]
aging_genes_rat_up <- aging_genes_rat_up[aging_genes_rat_up %in% rownames(Skin@assays$RNA@counts)]
up_counts <- length(aging_genes_rat_up)

aging_genes_rat <- as.character(c(as.character(aging_genes_rat_down), as.character(aging_genes_rat_up)))

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(Skin@assays$RNA@counts[aging_genes_rat,]>0)/dim(Skin@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(Skin@assays$RNA@data[aging_genes_rat,]),
    center = T, scale = T))
aging.score = t(z_matrix[aging_genes_rat,]) %*% weights 

df_aging_score <- as.data.frame(cbind(Skin$age, aging.score))
colnames(df_aging_score) <- c("age", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$age <- factor(df_aging_score$age, 
                             levels = c("Y", "O", "CR"),
                             labels = c("Young", "Old", "CR"))

p <- ggplot(df_aging_score, aes(x = age, y = aging_score)) +
  geom_violin(aes(fill = age)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Young", "Old"), c("Old", "CR")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  scale_fill_npg() + 
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  labs(x="Age",y="Aging score", title="Skin")
p
ggsave("./figure/Skin/Skin_aging_score_facs.pdf", width = 6, height = 4)

# droplet
aging_genes_mouse_down <- c(as.character(aging_gene_glmnet_all_tissue[["droplet-Skin"]][1:50, "down"]))
aging_genes_rat_down <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_down]
aging_genes_rat_down <- aging_genes_rat_down[aging_genes_rat_down %in% rownames(Skin@assays$RNA@counts)]
down_counts <- length(aging_genes_rat_down)
aging_genes_mouse_up <- c(as.character(aging_gene_glmnet_all_tissue[["droplet-Skin"]][1:50, "up"]))
aging_genes_rat_up <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_up]
aging_genes_rat_up <- aging_genes_rat_up[aging_genes_rat_up %in% rownames(Skin@assays$RNA@counts)]
up_counts <- length(aging_genes_rat_up)

aging_genes_rat <- as.character(c(as.character(aging_genes_rat_down), as.character(aging_genes_rat_up)))

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(Skin@assays$RNA@counts[aging_genes_rat,]>0)/dim(Skin@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(Skin@assays$RNA@data[aging_genes_rat,]),
    center = T, scale = T))
aging.score = t(z_matrix[aging_genes_rat,]) %*% weights 

df_aging_score <- as.data.frame(cbind(Skin$age, aging.score))
colnames(df_aging_score) <- c("age", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$age <- factor(df_aging_score$age, 
                             levels = c("Y", "O", "CR"),
                             labels = c("Young", "Old", "CR"))

p <- ggplot(df_aging_score, aes(x = age, y = aging_score)) +
  geom_violin(aes(fill = age)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Young", "Old"), c("Old", "CR")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  scale_fill_npg() + 
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  labs(x="Age",y="Aging score", title="Skin")
p
ggsave("./figure/Skin/Skin_aging_score_droplet.pdf", width = 6, height = 4)

rm("Skin")


# Brain
Brain <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_integral_rds/Brain.rds")

# facs
aging_genes_mouse_down <- c(as.character(aging_gene_glmnet_all_tissue[["facs-Brain_Myeloid"]][1:50, "down"]))
aging_genes_rat_down <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_down]
aging_genes_rat_down <- aging_genes_rat_down[aging_genes_rat_down %in% rownames(Brain@assays$RNA@counts)]
down_counts <- length(aging_genes_rat_down)
aging_genes_mouse_up <- c(as.character(aging_gene_glmnet_all_tissue[["facs-Brain_Myeloid"]][1:50, "up"]))
aging_genes_rat_up <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_up]
aging_genes_rat_up <- aging_genes_rat_up[aging_genes_rat_up %in% rownames(Brain@assays$RNA@counts)]
up_counts <- length(aging_genes_rat_up)

aging_genes_rat <- as.character(c(as.character(aging_genes_rat_down), as.character(aging_genes_rat_up)))

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(Brain@assays$RNA@counts[aging_genes_rat,]>0)/dim(Brain@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(Brain@assays$RNA@data[aging_genes_rat,]),
    center = T, scale = T))
aging.score = t(z_matrix[aging_genes_rat,]) %*% weights 

df_aging_score <- as.data.frame(cbind(Brain$age, aging.score))
colnames(df_aging_score) <- c("age", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$age <- factor(df_aging_score$age, 
                             levels = c("Y", "O", "CR"),
                             labels = c("Young", "Old", "CR"))

p <- ggplot(df_aging_score, aes(x = age, y = aging_score)) +
  geom_violin(aes(fill = age)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Young", "Old"), c("Old", "CR")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  scale_fill_npg() + 
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  labs(x="Age",y="Aging score", title="Brain")
p
ggsave("./figure/Brain/Brain_aging_score_facs.pdf", width = 6, height = 4)

# droplet
aging_genes_mouse_down <- c(as.character(aging_gene_glmnet_all_tissue[["facs-Brain_Non"]][1:50, "down"]))
aging_genes_rat_down <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_down]
aging_genes_rat_down <- aging_genes_rat_down[aging_genes_rat_down %in% rownames(Brain@assays$RNA@counts)]
down_counts <- length(aging_genes_rat_down)
aging_genes_mouse_up <- c(as.character(aging_gene_glmnet_all_tissue[["facs-Brain_Non"]][1:50, "up"]))
aging_genes_rat_up <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_up]
aging_genes_rat_up <- aging_genes_rat_up[aging_genes_rat_up %in% rownames(Brain@assays$RNA@counts)]
up_counts <- length(aging_genes_rat_up)

aging_genes_rat <- as.character(c(as.character(aging_genes_rat_down), as.character(aging_genes_rat_up)))

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(Brain@assays$RNA@counts[aging_genes_rat,]>0)/dim(Brain@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(Brain@assays$RNA@data[aging_genes_rat,]),
    center = T, scale = T))
aging.score = t(z_matrix[aging_genes_rat,]) %*% weights 

df_aging_score <- as.data.frame(cbind(Brain$age, aging.score))
colnames(df_aging_score) <- c("age", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$age <- factor(df_aging_score$age, 
                             levels = c("Y", "O", "CR"),
                             labels = c("Young", "Old", "CR"))

p <- ggplot(df_aging_score, aes(x = age, y = aging_score)) +
  geom_violin(aes(fill = age)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Young", "Old"), c("Old", "CR")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  scale_fill_npg() + 
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  labs(x="Age",y="Aging score", title="Brain")
p
ggsave("./figure/Brain/Brain_aging_score_droplet.pdf", width = 6, height = 4)

rm("Brain")


# Muscle
Muscle <- readRDS("/lustre/user/liclab/liocean/maosl/CR/CR_integral_rds/Muscle.rds")

# facs
aging_genes_mouse_down <- c(as.character(aging_gene_glmnet_all_tissue[["facs-Limb_Muscle"]][1:50, "down"]))
aging_genes_rat_down <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_down]
aging_genes_rat_down <- aging_genes_rat_down[aging_genes_rat_down %in% rownames(Muscle@assays$RNA@counts)]
down_counts <- length(aging_genes_rat_down)
aging_genes_mouse_up <- c(as.character(aging_gene_glmnet_all_tissue[["facs-Limb_Muscle"]][1:50, "up"]))
aging_genes_rat_up <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_up]
aging_genes_rat_up <- aging_genes_rat_up[aging_genes_rat_up %in% rownames(Muscle@assays$RNA@counts)]
up_counts <- length(aging_genes_rat_up)

aging_genes_rat <- as.character(c(as.character(aging_genes_rat_down), as.character(aging_genes_rat_up)))

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(Muscle@assays$RNA@counts[aging_genes_rat,]>0)/dim(Muscle@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(Muscle@assays$RNA@data[aging_genes_rat,]),
    center = T, scale = T))
aging.score = t(z_matrix[aging_genes_rat,]) %*% weights 

df_aging_score <- as.data.frame(cbind(Muscle$age, aging.score))
colnames(df_aging_score) <- c("age", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$age <- factor(df_aging_score$age, 
                             levels = c("Y", "O", "CR"),
                             labels = c("Young", "Old", "CR"))

p <- ggplot(df_aging_score, aes(x = age, y = aging_score)) +
  geom_violin(aes(fill = age)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Young", "Old"), c("Old", "CR")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  scale_fill_npg() + 
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  labs(x="Age",y="Aging score", title="Muscle")
p
ggsave("./figure/Muscle/Muscle_aging_score_facs.pdf", width = 6, height = 4)

# droplet
aging_genes_mouse_down <- c(as.character(aging_gene_glmnet_all_tissue[["droplet-Limb_Muscle"]][1:50, "down"]))
aging_genes_rat_down <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_down]
aging_genes_rat_down <- aging_genes_rat_down[aging_genes_rat_down %in% rownames(Muscle@assays$RNA@counts)]
down_counts <- length(aging_genes_rat_down)
aging_genes_mouse_up <- c(as.character(aging_gene_glmnet_all_tissue[["droplet-Limb_Muscle"]][1:50, "up"]))
aging_genes_rat_up <- rat_mouse_gene_transform$rat[rat_mouse_gene_transform$mouse %in% aging_genes_mouse_up]
aging_genes_rat_up <- aging_genes_rat_up[aging_genes_rat_up %in% rownames(Muscle@assays$RNA@counts)]
up_counts <- length(aging_genes_rat_up)

aging_genes_rat <- as.character(c(as.character(aging_genes_rat_down), as.character(aging_genes_rat_up)))

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(Muscle@assays$RNA@counts[aging_genes_rat,]>0)/dim(Muscle@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(Muscle@assays$RNA@data[aging_genes_rat,]),
    center = T, scale = T))
aging.score = t(z_matrix[aging_genes_rat,]) %*% weights 

df_aging_score <- as.data.frame(cbind(Muscle$age, aging.score))
colnames(df_aging_score) <- c("age", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$age <- factor(df_aging_score$age, 
                             levels = c("Y", "O", "CR"),
                             labels = c("Young", "Old", "CR"))

p <- ggplot(df_aging_score, aes(x = age, y = aging_score)) +
  geom_violin(aes(fill = age)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Young", "Old"), c("Old", "CR")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  scale_fill_npg() + 
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  labs(x="Age",y="Aging score", title="Muscle")
p
ggsave("./figure/Muscle/Muscle_aging_score_droplet.pdf", width = 6, height = 4)

rm("Muscle")

