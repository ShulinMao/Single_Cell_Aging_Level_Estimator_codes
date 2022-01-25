.libPaths("/lustre/user/liclab/liocean/maosl/Rpackages")
setwd("/lustre/user/liclab/liocean/maosl/AD_mouse_single_nucleus/")

library(Seurat)
library(ggplot2)
library(Matrix)
library(ggsignif)
library(ggsci)

# preprocessing
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

# Microglia marker "C1qa" & "Csf1r"
p3 <- VlnPlot(AD, features = c("C1qa", "Csf1r"), pt.size = 0)
p4 <- FeaturePlot(AD, features = c("C1qa", "Csf1r"))
ggsave(plot = p3+p4, "./figure/integrated/Microglia_marker_gene.pdf", width = 10, height = 5)

AD_Microglia <- subset(AD, idents = c(8, 14))
AD_Non_Myeloid <- subset(AD, idents = c(0:7, 9:13))
saveRDS(AD_Non_Myeloid, "AD_Non_Myeloid.rds")
saveRDS(AD_Microglia, "AD_Microglia.rds")


# caluclate Aging score
aging_gene_glmnet_all_tissue <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/aging_gene_glmnet_all_tissue_add_50_corrected.rds")

# Microglia
AD_Microglia <- readRDS("/lustre/user/liclab/liocean/maosl/AD_mouse_single_nucleus/data_rds/AD_Microglia.rds")
aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Brain_Myeloid`[1:50, "down"])
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(AD_Microglia@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Brain_Myeloid`[1:50, "up"])
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(AD_Microglia@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)


weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(AD_Microglia@assays$RNA@counts[aging_gene,])>0)/dim(AD_Microglia@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(AD_Microglia@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

df_aging_score <- as.data.frame(cbind(AD_Microglia$type, aging.score))
colnames(df_aging_score) <- c("Type", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$Type <- factor(df_aging_score$Type, levels = c("Normal", "AD"))

AD_aging_score <- df_aging_score[df_aging_score$Type=="AD", "aging_score"]
Normal_aging_score <- df_aging_score[df_aging_score$Type=="Normal", "aging_score"]
AD_normal_ttest <- t.test(AD_aging_score, Normal_aging_score, alternative = "greater")

AD_normal_summary <- summarySE(df_aging_score, measurevar="aging_score", groupvars=c("Type"))

ggplot(AD_normal_summary, aes(x=Type, y=aging_score, group = 1)) + 
  geom_errorbar(aes(ymin=aging_score-se, ymax=aging_score+se), width=.1) +
  geom_point() +
  geom_line() +
  geom_signif(annotations = "***", y_position = 0.2, xmin=1, xmax=2, textsize = 7) +
  coord_cartesian(ylim = c(-0.25, 0.25)) +
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  labs(x="Sample",y="Aging score",title="Brain (Myeloid)")
ggsave("./figure/aging_score/Microglia_line.pdf", width = 5, height = 3.5)

ggplot(df_aging_score, aes(x = Type, y = aging_score)) +
  geom_violin(aes(fill = Type)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("AD", "Normal")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  theme_classic() + 
  scale_fill_npg() + 
  labs(x="Sample",y="Aging score",title="Brain_Myeloid")
ggsave("./figure/aging_score/Microglia.pdf", width = 5, height = 3.5)


# Non_Myeloid
AD_Non_Myeloid <- readRDS("/lustre/user/liclab/liocean/maosl/AD_mouse_single_nucleus/data_rds/AD_Non_Myeloid.rds")
aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Brain_Non`[1:50, "down"])
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(AD_Non_Myeloid@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Brain_Non`[1:50, "up"])
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(AD_Non_Myeloid@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)


weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(AD_Non_Myeloid@assays$RNA@counts[aging_gene,]>0)/dim(AD_Non_Myeloid@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(AD_Non_Myeloid@assays$RNA@data[aging_gene,]),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

df_aging_score <- as.data.frame(cbind(AD_Non_Myeloid$type, aging.score))
colnames(df_aging_score) <- c("Type", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))
df_aging_score$Type <- factor(df_aging_score$Type, levels = c("Normal", "AD"))

AD_aging_score <- df_aging_score[df_aging_score$Type=="AD", "aging_score"]
Normal_aging_score <- df_aging_score[df_aging_score$Type=="Normal", "aging_score"]
AD_normal_ttest <- t.test(AD_aging_score, Normal_aging_score, alternative = "greater")

AD_normal_summary <- summarySE(df_aging_score, measurevar="aging_score", groupvars=c("Type"))

ggplot(AD_normal_summary, aes(x=Type, y=aging_score, group = 1)) + 
  geom_errorbar(aes(ymin=aging_score-se, ymax=aging_score+se), width=.1) +
  geom_point() +
  geom_line() +
  geom_signif(annotations = "***", y_position = 0.15, xmin=1, xmax=2, textsize = 7) +
  coord_cartesian(ylim = c(-0.1, 0.17)) +
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  labs(x="Sample",y="Aging score",title="Brain (Non-Myeloid)")
ggsave("./figure/aging_score/Non_Microglia_line.pdf", width = 5, height = 3.5)

ggplot(df_aging_score, aes(x = Type, y = aging_score)) +
  geom_violin(aes(fill = Type)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("AD", "Normal")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  geom_line(data=AD_normal_summary, aes(x=Type, y=aging_score, group = 1)) +
  theme_classic() + 
  scale_fill_npg() + 
  labs(x="Sample",y="Aging score",title="Brain_Non-Myeloid")
ggsave("./figure/aging_score/Non_Myeloid.pdf", width = 5, height = 3.5)

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