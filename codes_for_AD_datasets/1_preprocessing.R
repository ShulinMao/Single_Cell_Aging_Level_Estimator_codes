.libPaths("/lustre/user/liclab/liocean/maosl/Rpackages")
setwd("/lustre/user/liclab/liocean/maosl/AD_mouse_single_nucleus/")

library(Seurat)
library(ggplot2)
library(Matrix)
library(ggsignif)
library(ggsci)

for(directory in dir("./data/")){
  print(directory)
  sample.counts <- Read10X(data.dir = paste0("./data/", directory))
  sample <- CreateSeuratObject(counts = sample.counts, min.cells = 5)
  sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^mt-")
  p1 <- VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  sample <- subset(sample, subset = nFeature_RNA > 300 & nFeature_RNA < 5600 &
                     nCount_RNA > 300 & nCount_RNA < 9000 &
                     percent.mt < 5)
  p2 <- VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  sample <- NormalizeData(object = sample, scale.factor = 10000)
  sample <- FindVariableFeatures(object = sample, nfeatures = 3000)
  sample <- ScaleData(object = sample)
  sample <- RunPCA(object = sample)
  sample <- FindNeighbors(object = sample)
  sample <- FindClusters(object = sample, resolution = 0.6)
  sample <- RunTSNE(object = sample, dims = 1:20)
  p3 <- DimPlot(object = sample, reduction = "tsne")
  ggsave(paste("./figure/filter/", directory, "_1.pdf", sep = ""), plot = p1, width = 7, height = 4)
  ggsave(paste("./figure/filter/", directory, "_2.pdf", sep = ""), plot = p2, width = 7, height = 4)
  ggsave(paste("./figure/tsne/", directory, ".pdf", sep = ""), plot = p3, width = 7, height = 5)
  saveRDS(sample, paste("./data_rds/", directory, ".rds", sep = ""))
}

# integration
`5xfad_1` <- readRDS("/lustre/user/liclab/liocean/maosl/AD_mouse_single_nucleus/data_rds/5xfad_1.rds")
`5xfad_1`$sample <- "5xfad_1"
`5xfad_1`$type <- "AD"
`5xfad_2` <- readRDS("/lustre/user/liclab/liocean/maosl/AD_mouse_single_nucleus/data_rds/5xfad_2.rds")
`5xfad_2`$sample <- "5xfad_2"
`5xfad_2`$type <- "AD"
`5xfad_3` <- readRDS("/lustre/user/liclab/liocean/maosl/AD_mouse_single_nucleus/data_rds/5xfad_3.rds")
`5xfad_3`$sample <- "5xfad_3"
`5xfad_3`$type <- "AD"
wt_1 <- readRDS("/lustre/user/liclab/liocean/maosl/AD_mouse_single_nucleus/data_rds/wt_1.rds")
wt_1$sample <- "wt_1"
wt_1$type <- "Normal"
wt_2 <- readRDS("/lustre/user/liclab/liocean/maosl/AD_mouse_single_nucleus/data_rds/wt_2.rds")
wt_2$sample <- "wt_2"
wt_2$type <- "Normal"
wt_3 <- readRDS("/lustre/user/liclab/liocean/maosl/AD_mouse_single_nucleus/data_rds/wt_3.rds")
wt_3$sample <- "wt_3"
wt_3$type <- "Normal"


AD_list <- list(`5xfad_1`, `5xfad_2`, `5xfad_3`,
                wt_1, wt_2, wt_3)

AD_list <- lapply(X = AD_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

features <- SelectIntegrationFeatures(object.list = AD_list)
AD_anchors <- FindIntegrationAnchors(object.list = AD_list, anchor.features = features)
saveRDS(AD_anchors, "./AD_anchor_temp.rds")
AD <- IntegrateData(anchorset = AD_anchors)
saveRDS(AD, "./AD_temp.rds")

DefaultAssay(AD) <- "integrated"

AD <- ScaleData(AD, verbose = FALSE)
AD <- RunPCA(AD)
ElbowPlot(AD)
AD <- FindNeighbors(AD, reduction = "pca", dims = 1:10)
AD <- FindClusters(AD, resolution = 0.1)
AD <- RunTSNE(AD, reduction = "pca")
p1 <- DimPlot(AD, reduction = "tsne", group.by = "sample")
p2 <- DimPlot(AD, reduction = "tsne", group.by = "type")
ggsave(plot = p1+p2, "./figure/integrated/AD_tsne.pdf", width = 14, height = 5)
DimPlot(AD, reduction = "tsne")
ggsave("./figure/integrated/AD_cluster.pdf", width = 7, height = 5)
saveRDS(AD, "AD_integral.rds")

# Microglia 单独取出来用Brain_Myeloid
# 其他用Brain_Non_Myeloid
# Microglia marker "C1qa" & "Csf1r"
p3 <- VlnPlot(AD, features = c("C1qa", "Csf1r"), pt.size = 0)
p4 <- FeaturePlot(AD, features = c("C1qa", "Csf1r"))
ggsave(plot = p3+p4, "./figure/integrated/Microglia_marker_gene.pdf", width = 10, height = 5)

AD_Microglia <- subset(AD, idents = c(8, 14))
AD_Non_Myeloid <- subset(AD, idents = c(0:7, 9:13))
saveRDS(AD_Non_Myeloid, "AD_Non_Myeloid.rds")
saveRDS(AD_Microglia, "AD_Microglia.rds")