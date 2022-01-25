.libPaths("/lustre/user/liclab/liocean/maosl/Rpackages")
setwd("/lustre/user/liclab/liocean/maosl/AD_mouse_single_nucleus/")

library(Matrix)
library(Seurat)
library(ggplot2)
library(ggsignif)
library(ggsci)

# Microglia单独取出来用Brain_Myeloid
# 其他用Brain_Non_Myeloid
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
  labs(x="Sample",y="SCALE score",title="Brain (Myeloid)")
ggsave("./figure/aging_score/Microglia_line.pdf", width = 5, height = 3.5)

ggplot(df_aging_score, aes(x = Type, y = aging_score)) +
  geom_violin(aes(fill = Type)) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("AD", "Normal")),
              step_increase = 0.1, map_signif_level = F, test = t.test) +
  theme_classic() + 
  scale_fill_npg() + 
  labs(x="Sample",y="SCALE score",title="Brain_Myeloid")
ggsave("./figure/aging_score/Microglia.pdf", width = 5, height = 3.5)


# other tissue
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
  
  # 计算长度
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # 以 groupvars 为组,计算每组的长度,均值,以及标准差
  # ddply 就是 dplyr 中的 group_by + summarise
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # 重命名  
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  # 计算标准偏差
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  # 计算置信区间
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
