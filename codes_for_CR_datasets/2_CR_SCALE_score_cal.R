.libPaths("/lustre/user/liclab/liocean/maosl/Rpackages")
setwd("/lustre/user/liclab/liocean/maosl/CR/")

library(Seurat)
library(ggplot2)
library(Matrix)
library(ggsignif)
library(ggsci)
library(dplyr)


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

ggplot(df_aging_score, aes(x = age, y = aging_score)) +
  geom_violin(aes(fill = age)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Young", "Old"), c("Old", "CR")),
              step_increase = 0.1, map_signif_level = T, test = t.test, textsize = 7) +
  scale_fill_npg() + 
  theme_classic() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  coord_cartesian(ylim = c(-8, 11)) +
  labs(x="Age", y="SCALE score", title="BAT")
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
  labs(x="Age",y="SCALE score", title="Aorta")
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
  labs(x="Age",y="SCALE score", title = "Aorta")
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
  labs(x="Age",y="SCALE score", title="Bone marrow")
p
ggsave("./figure/BM/BM_aging_score_facs.pdf", width = 6, height = 4)

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
  labs(x="Age",y="SCALE score", title="Kidney")
p
ggsave("./figure/Kidney/Kidney_aging_score_facs.pdf", width = 6, height = 4)

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
  labs(x="Age",y="SCALE score", title="Liver")
p
ggsave("./figure/Liver/Liver_aging_score_facs.pdf", width = 6, height = 4)

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
  labs(x="Age",y="SCALE score", title="Skin")
p
ggsave("./figure/Skin/Skin_aging_score_facs.pdf", width = 6, height = 4)


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
  labs(x="Age",y="SCALE score", title="Brain")
p
ggsave("./figure/Brain/Brain_aging_score_facs.pdf", width = 6, height = 4)

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
  labs(x="Age",y="SCALE score", title="Muscle")
p
ggsave("./figure/Muscle/Muscle_aging_score_facs.pdf", width = 6, height = 4)



rm("Muscle")

