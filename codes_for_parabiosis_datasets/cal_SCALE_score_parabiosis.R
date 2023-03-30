.libPaths("/lustre/user/liclab/liocean/maosl/Rpackages")
setwd("/lustre/user/liclab/liocean/maosl/parabiosis/")

library(Seurat)
library(ggplot2)
library(stringr)
library(ggsignif)
library(ggsci)

aging_gene_glmnet_all_tissue <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/aging_gene_glmnet_all_tissue_add_50_corrected.rds")

# limb_muscle
limb_muscle <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/limb_muscle.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Limb_Muscle`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(limb_muscle@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Limb_Muscle`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(limb_muscle@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(limb_muscle@assays$RNA@counts[aging_gene,])>0)/dim(limb_muscle@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(limb_muscle@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

limb_muscle[["condition_character"]] = as.character(limb_muscle$condition)

df_aging_score <- as.data.frame(cbind(limb_muscle$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "IA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))


ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*2) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(6, 5, 10), textsize = 5, tip_length = 0, vjust = 0.5) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() + 
  labs(x="Conditions",y="SCALE score",title="Limb-Muscle")
ggsave("./figures/limb_muscle.pdf", width = 12, height = 12, unit = "cm")

# liver
liver <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/liver.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Liver`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(liver@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Liver`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(liver@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(liver@assays$RNA@counts[aging_gene,])>0)/dim(liver@assays$RNA@counts)[2]

z_matrix = t(
  scale(
        t(as.matrix(liver@assays$RNA@data[aging_gene,])),
        center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

liver[["condition_character"]] = as.character(liver$condition)

df_aging_score <- as.data.frame(cbind(liver$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*2) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(19, 15, 21.5), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() + 
  labs(x="Conditions",y="SCALE score",title="Liver")
ggsave("./figures/liver.pdf", width = 6, height = 6)

# kidney
kidney <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/kidney.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Kidney`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(kidney@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Kidney`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(kidney@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(kidney@assays$RNA@counts[aging_gene,])>0)/dim(kidney@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(kidney@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

kidney[["condition_character"]] = as.character(kidney$condition)

df_aging_score <- as.data.frame(cbind(kidney$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "A", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*2) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(8, 7, 9.5), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="Kidney")
ggsave("./figures/kidney.pdf", width = 6, height = 6)

# lung
lung <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/lung.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Lung`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(lung@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Lung`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(lung@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(lung@assays$RNA@counts[aging_gene,])>0)/dim(lung@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(lung@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

lung[["condition_character"]] = as.character(lung$condition)

df_aging_score <- as.data.frame(cbind(lung$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*2) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(10, 13, 15), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="Lung")
ggsave("./figures/lung.pdf", width = 6, height = 6)


# large_intestine
large_intestine <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/large_intestine.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Large_Intestine`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(large_intestine@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Large_Intestine`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(large_intestine@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(large_intestine@assays$RNA@counts[aging_gene,])>0)/dim(large_intestine@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(large_intestine@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

large_intestine[["condition_character"]] = as.character(large_intestine$condition)

df_aging_score <- as.data.frame(cbind(large_intestine$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*2) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(15, 13, 20), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="Large Intestine")
ggsave("./figures/large_intestine.pdf", width = 6, height = 6)


# marrow
marrow <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/marrow.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Marrow`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(marrow@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Marrow`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(marrow@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(marrow@assays$RNA@counts[aging_gene,])>0)/dim(marrow@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(marrow@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

marrow[["condition_character"]] = as.character(marrow$condition)

df_aging_score <- as.data.frame(cbind(marrow$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])


df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*2) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(11, 10, 15), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="Marrow")
ggsave("./figures/marrow.pdf", width = 6, height = 6)

# bladder
bladder <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/bladder.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Bladder`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(bladder@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Bladder`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(bladder@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(bladder@assays$RNA@counts[aging_gene,])>0)/dim(bladder@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(bladder@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

bladder[["condition_character"]] = as.character(bladder$condition)

df_aging_score <- as.data.frame(cbind(bladder$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])


df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*2) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(28, 20, 35), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="Bladder")
ggsave("./figures/bladder.pdf", width = 6, height = 6)


# heart
heart <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/heart.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Heart`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(heart@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Heart`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(heart@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(heart@assays$RNA@counts[aging_gene,])>0)/dim(heart@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(heart@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

heart[["condition_character"]] = as.character(heart$condition)

df_aging_score <- as.data.frame(cbind(heart$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])


df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*2) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(10, 7, 12), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="Heart")
ggsave("./figures/heart.pdf", width = 6, height = 6)

# tongue
tongue <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/tongue.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Tongue`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(tongue@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Tongue`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(tongue@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(tongue@assays$RNA@counts[aging_gene,])>0)/dim(tongue@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(tongue@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

tongue[["condition_character"]] = as.character(tongue$condition)

df_aging_score <- as.data.frame(cbind(tongue$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])


df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*2) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(20, 15, 25), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="Tongue")
ggsave("./figures/tongue.pdf", width = 6, height = 6)

# brain
brain <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/brain.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Brain_Myeloid`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(brain@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Brain_Myeloid`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(brain@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(brain@assays$RNA@counts[aging_gene,])>0)/dim(brain@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(brain@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

brain[["condition_character"]] = as.character(brain$condition)

df_aging_score <- as.data.frame(cbind(brain$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])


df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*2) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(5, 3, 6), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="Brain")
ggsave("./figures/brain.pdf", width = 6, height = 6)

# brain <- subset(brain, cells = colnames(brain)[brain$condition_character %in% c("A", "HA", "HY", "Y")])
# brain <- NormalizeData(object = brain)
# brain <- FindVariableFeatures(object = brain)
# brain <- ScaleData(object = brain)
# brain <- RunPCA(object = brain)
# brain <- FindNeighbors(object = brain)
# brain <- FindClusters(object = brain)
# brain <- RunUMAP(object = brain, dims = 1:15)
# brain[["parabiosis_condition"]] <- factor(brain$condition_character, levels = c("A", "HA", "HY", "Y"), 
#                                           labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))
# DimPlot(brain, reduction = "umap", group.by = "parabiosis_condition")
# saveRDS(brain, "../GR_figure_scripts/data/paabiosis_brain.rds")

# trachea
trachea <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/trachea.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Trachea`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(trachea@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Trachea`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(trachea@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(trachea@assays$RNA@counts[aging_gene,])>0)/dim(trachea@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(trachea@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

trachea[["condition_character"]] = as.character(trachea$condition)

df_aging_score <- as.data.frame(cbind(trachea$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])


df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*2) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(20, 15, 25), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="Trachea")
ggsave("./figures/trachea.pdf", width = 6, height = 6)

# diaphragm
diaphragm <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/diaphragm.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Diaphragm`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(diaphragm@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Diaphragm`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(diaphragm@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(diaphragm@assays$RNA@counts[aging_gene,])>0)/dim(diaphragm@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(diaphragm@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

diaphragm[["condition_character"]] = as.character(diaphragm$condition)

df_aging_score <- as.data.frame(cbind(diaphragm$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])


df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*10) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(20, 15, 25), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="Diaphragm")
ggsave("./figures/diaphragm.pdf", width = 6, height = 6)

# skin
skin <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/skin.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Skin`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(skin@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Skin`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(skin@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(skin@assays$RNA@counts[aging_gene,])>0)/dim(skin@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(skin@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

skin[["condition_character"]] = as.character(skin$condition)

df_aging_score <- as.data.frame(cbind(skin$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])


df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*4) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(35, 25, 45), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="Skin")
ggsave("./figures/skin.pdf", width = 6, height = 6)

# spleen
spleen <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/spleen.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Spleen`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(spleen@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Spleen`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(spleen@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(spleen@assays$RNA@counts[aging_gene,])>0)/dim(spleen@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(spleen@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

spleen[["condition_character"]] = as.character(spleen$condition)

df_aging_score <- as.data.frame(cbind(spleen$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])


df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*2) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(12, 10, 15), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="Spleen")
ggsave("./figures/spleen.pdf", width = 6, height = 6)

# pancreas
pancreas <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/pancreas.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Pancreas`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(pancreas@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Pancreas`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(pancreas@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(pancreas@assays$RNA@counts[aging_gene,])>0)/dim(pancreas@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(pancreas@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

pancreas[["condition_character"]] = as.character(pancreas$condition)

df_aging_score <- as.data.frame(cbind(pancreas$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])


df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*2) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(10, 10, 12), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="Pancreas")
ggsave("./figures/pancreas.pdf", width = 6, height = 6)

# BAT
BAT <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/BAT.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-BAT`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(BAT@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-BAT`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(BAT@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(BAT@assays$RNA@counts[aging_gene,])>0)/dim(BAT@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(BAT@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

BAT[["condition_character"]] = as.character(BAT$condition)

df_aging_score <- as.data.frame(cbind(BAT$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])


df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*5) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(6, 2, 8), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="BAT")
ggsave("./figures/BAT.pdf", width = 6, height = 6)

# MAT
MAT <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/MAT.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-MAT`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(MAT@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-MAT`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(MAT@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(MAT@assays$RNA@counts[aging_gene,])>0)/dim(MAT@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(MAT@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

MAT[["condition_character"]] = as.character(MAT$condition)

df_aging_score <- as.data.frame(cbind(MAT$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])


df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*3) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(15, 8, 20), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="MAT")
ggsave("./figures/MAT.pdf", width = 6, height = 6)

# GAT
GAT <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/GAT.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-GAT`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(GAT@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-GAT`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(GAT@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(GAT@assays$RNA@counts[aging_gene,])>0)/dim(GAT@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(GAT@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

GAT[["condition_character"]] = as.character(GAT$condition)

df_aging_score <- as.data.frame(cbind(GAT$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])


df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*2) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(30, 15, 35), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="GAT")
ggsave("./figures/GAT.pdf", width = 6, height = 6)

# SCAT
SCAT <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/SCAT.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-SCAT`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(SCAT@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-SCAT`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(SCAT@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(SCAT@assays$RNA@counts[aging_gene,])>0)/dim(SCAT@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(SCAT@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

SCAT[["condition_character"]] = as.character(SCAT$condition)

df_aging_score <- as.data.frame(cbind(SCAT$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "IA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "HA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])


df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))

ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*6) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(24, 15, 28), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() +
  labs(x="Conditions",y="SCALE score",title="SCAT")
ggsave("./figures/SCAT.pdf", width = 6, height = 6)

# thymus
thymus <- readRDS("/lustre/user/liclab/liocean/maosl/parabiosis/data/thymus.rds")

aging_gene_down <- as.character(aging_gene_glmnet_all_tissue$`facs-Thymus`[1:50, "down"])
aging_gene_down <- str_to_lower(aging_gene_down)
aging_gene_down <- aging_gene_down[aging_gene_down %in% row.names(thymus@assays$RNA@counts)]
down_counts <- length(aging_gene_down)
aging_gene_up <- as.character(aging_gene_glmnet_all_tissue$`facs-Thymus`[1:50, "up"])
aging_gene_up <- str_to_lower(aging_gene_up)
aging_gene_up <- aging_gene_up[aging_gene_up %in% row.names(thymus@assays$RNA@counts)]
up_counts <- length(aging_gene_up)
aging_gene <- c(aging_gene_down, aging_gene_up)

weights = c(rep(-1, down_counts), rep(1, up_counts))*
  Matrix::rowSums(as.matrix(thymus@assays$RNA@counts[aging_gene,])>0)/dim(thymus@assays$RNA@counts)[2]

z_matrix = t(
  scale(
    t(as.matrix(thymus@assays$RNA@data[aging_gene,])),
    center = T, scale = T))
aging.score = t(z_matrix[aging_gene,]) %*% weights 

thymus[["condition_character"]] = as.character(thymus$condition)

df_aging_score <- as.data.frame(cbind(thymus$condition_character, aging.score))
colnames(df_aging_score) <- c("Conditions", "aging_score")
df_aging_score$aging_score <- as.numeric(as.character(df_aging_score$aging_score))

t.test(df_aging_score[df_aging_score$Conditions == "IA", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "IY", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HY", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "A", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "HA", "aging_score"])

t.test(df_aging_score[df_aging_score$Conditions == "Y", "aging_score"], 
       df_aging_score[df_aging_score$Conditions == "A", "aging_score"])

df_aging_score = df_aging_score[(df_aging_score$Conditions != "IA") & (df_aging_score$Conditions != "IY"),]
df_aging_score$Conditions = factor(df_aging_score$Conditions, levels = c("A", "HA", "HY", "Y"), 
                                   labels = c("Old", "Old Heterochronic", "Young Heterochronic", "Young"))


ggplot(df_aging_score, aes(x = Conditions, y = aging_score)) +
  geom_boxplot(aes(fill = Conditions), outlier.shape = NA) +
  coord_cartesian(ylim = boxplot.stats(df_aging_score$aging_score)$stats[c(1, 5)]*3) + 
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("Old", "Old Heterochronic"), c("Young Heterochronic", "Young"), c("Old", "Young")), 
              map_signif_level = T, test = t.test, y_position = c(28, 20, 36), textsize = 5, tip_length = 0) +
  theme_classic() +
  theme_set(theme_classic(base_size = 15)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  guides(fill = F) +
  scale_fill_npg() + 
  labs(x="Conditions",y="SCALE score",title="Thymus")
ggsave("./figures/thymus.pdf", width = 6, height = 6)
