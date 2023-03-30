.libPaths("/lustre/user/liclab/liocean/maosl/Rpackages")
setwd("/lustre/user/liclab/liocean/maosl/methy_analysis/")

library(Seurat)
library(stringr)
library(dplyr)
library(zoo)
library(ggplot2)
library(ggsignif)
library(ggsci)

# read counts matrix
raw_counts <- read.table("./dataset/muscle_cell_stem_cell/RNAseq_data/GSE121364_Table_raw_counts_RNA_QSC.txt",
                         header = T, stringsAsFactors = F)
cell_name <- colnames(raw_counts)[-1:-3]
cell_name <- str_extract(cell_name, "[OY][0-9]_[A-Z]{2}_[A-Za-z]{2}_[0-9]{1,2}")
colnames(raw_counts) <- c(colnames(raw_counts)[1:3], cell_name)
raw_counts <- raw_counts[-2:-3]

unique_gene_counts <- raw_counts[!duplicated(raw_counts$geneName),]
row.names(unique_gene_counts) <- unique_gene_counts$geneName
unique_gene_counts <- unique_gene_counts[-1]

non_unique_gene <- unique(raw_counts$geneName[duplicated(raw_counts$geneName)])
for (gene in non_unique_gene){
  unique_gene_counts[gene,] <- colSums(raw_counts[raw_counts$geneName == gene, -1])
}
raw_counts <- unique_gene_counts
rm(unique_gene_counts)

# load metadata
meta_data <- read.csv("./dataset/muscle_cell_stem_cell/meta_data/Muscle_stem_cells_meta_sc_rna.txt",
                        stringsAsFactors = F, header = T)
cell_name_to_meta_data <- read.csv("./dataset/muscle_cell_stem_cell/meta_data/cell_barcodes_scRNAseq.csv",
                                   stringsAsFactors = F, header = F)
colnames(cell_name_to_meta_data) <- c("Sample.Name", "cell_name", "data_type")
meta_data <- left_join(meta_data, cell_name_to_meta_data[,1:2])
cell_name <- colnames(raw_counts)
meta_data <- meta_data[meta_data$cell_name %in% cell_name,]
rm(cell_name_to_meta_data)
row.names(meta_data) <- meta_data$cell_name

# create Seurat object
data <- CreateSeuratObject(counts = raw_counts, project = "muscle_stem_cell")
data <- AddMetaData(data, meta_data)
data <- NormalizeData(object = data)
data <- FindVariableFeatures(object = data)
data <- ScaleData(object = data)
data <- RunPCA(object = data)
data <- FindNeighbors(object = data)
data <- FindClusters(object = data)
data <- RunTSNE(object = data)
DimPlot(object = data, reduction = "tsne")
data <- RunUMAP(object = data, dims = 1:30)
DimPlot(object = data, reduction = "umap")

DimPlot(object = data, reduction = "tsne", group.by = "Age")
DimPlot(object = data, reduction = "umap", group.by = "Age")

# compute SCALE score
aging_gene_glmnet_all_tissue <- readRDS("/lustre/user/liclab/liocean/maosl/tms/tms_h5ad/aging_gene_glmnet_all_tissue_add_50_corrected.rds")
aging_genes <- c(as.character(aging_gene_glmnet_all_tissue[["droplet-Limb_Muscle"]][1:50, "up"]),
                 as.character(aging_gene_glmnet_all_tissue[["droplet-Limb_Muscle"]][1:50, "down"]))

coefficient <- c(rep(1, 50), rep(-1, 50))
coefficient <- coefficient[aging_genes %in% row.names(raw_counts)]
aging_genes <- aging_genes[aging_genes %in% row.names(raw_counts)]

weights <- coefficient *
  Matrix::rowSums(GetAssayData(data, "counts")[aging_genes,]>0)/dim(GetAssayData(data, "counts"))[2]
aging_genes <- aging_genes[weights != 0]
weights <- weights[weights != 0]

z_matrix <- t(scale(t(as.matrix(GetAssayData(data, "data")[aging_genes,])), center = T, scale = T))

aging.score = t(z_matrix[aging_genes,]) %*% weights  

data$aging_scores <- as.numeric(aging.score)                                  

data[["age_value"]] <- factor(data$Age, levels = c("2 months", "24 months"), labels = c(2, 24))
SCALE_result <- data[[c("cell_name", "Age", "age_value", "aging_scores")]]
SCALE_result$age_value <- as.numeric(as.character(SCALE_result$age_value))
regression <- summary(lm(age_value~aging_scores, data = SCALE_result))
regression
cor(SCALE_result$age_value, SCALE_result$aging_scores)

ggplot(SCALE_result, aes(x = Age, y = aging_scores)) +
  geom_boxplot(aes(fill = Age)) +
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("2 months", "24 months")), 
              step_increase = 0.1, map_signif_level = F, test = t.test, textsize = 5) +
  theme_classic() +
  scale_fill_npg() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  coord_cartesian(ylim = c(-6, 8)) +
  labs(x="Chronological age",y="SCALE score", title = "Muscle stem cells \n(Hernando et al. 2019)")
#ggsave(filename = "./dataset/muscle_cell_stem_cell/figure/SCALE_score_all_cells.pdf", width = 5, height = 5)

SCALE_result["aging_scores_scaled"] <- as.numeric(scale(SCALE_result$aging_scores, center = T, scale = T))

ggplot(SCALE_result, aes(x = Age, y = aging_scores_scaled)) +
  geom_boxplot(aes(fill = Age)) +
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("2 months", "24 months")), 
              step_increase = 0.1, map_signif_level = F, test = t.test, textsize = 5) +
  theme_classic() +
  scale_fill_npg() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  coord_cartesian(ylim = c(-6, 8)) +
  labs(x="Chronological age",y="SCALE score (scaled)", title = "Muscle stem cells \n(Hernando et al. 2019)")

# integrate scBS-seq and scRNA-seq data
meta_data_methy <- read.csv("./dataset/muscle_cell_stem_cell/meta_data/Muscle_stem_cells_meta_sc_methy.txt",
                      stringsAsFactors = F, header = T)
scAge_result <- read.table("./dataset/muscle_cell_stem_cell/scAge_result/Muscle_stem_cells-train(Thompson_Liver_BL6)-mode(percentile)-param(top_1_pct).tsv",
                           sep = "\t", header = T, stringsAsFactors = F)
colnames(scAge_result)[1] <- "Sample.Name"
scAge_result["scAge_scaled"] <- as.numeric(scale(scAge_result$PredictedAge, center = T, scale = T))
meta_data_methy <- left_join(meta_data_methy, scAge_result)

cell_name_to_meta_data <- read.csv("./dataset/muscle_cell_stem_cell/meta_data/cell_barcodes_scBSseq.csv",
                                   stringsAsFactors = F, header = F)
colnames(cell_name_to_meta_data) <- c("Sample.Name", "cell_name", "data_type")
meta_data_methy <- left_join(meta_data_methy, cell_name_to_meta_data[,1:2])
meta_data_methy[startsWith(meta_data_methy$cell_name, "O"), "Age"] <- "24 months"

SCALE_scAge_result <- full_join(SCALE_result, meta_data_methy)

regression <- summary(lm(age_value~aging_scores, data = SCALE_scAge_result))
regression

regression_scAge <- summary(lm(age_value~PredictedAge, data = SCALE_scAge_result))
regression_scAge

ggplot(SCALE_scAge_result, aes(x = Age, y = aging_scores)) +
  geom_boxplot(aes(fill = Age)) +
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("2 months", "24 months")), 
              step_increase = 0.1, map_signif_level = F, test = t.test, textsize = 5) +
  theme_classic() +
  scale_fill_npg() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  #coord_cartesian(ylim = c(-6, 8)) +
  labs(x="Chronological age",y="SCALE score", title = "Muscle stem cells \n(Hernando et al. 2019)")
ggsave(filename = "./dataset/muscle_cell_stem_cell/figure/SCALE_score_multi_omics_cells.pdf", width = 5, height = 5)

ggplot(SCALE_scAge_result, aes(x = Age, y = PredictedAge)) +
  geom_boxplot(aes(fill = Age)) +
  #stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange") +
  geom_signif(comparisons = list(c("2 months", "24 months")), 
              step_increase = 0.1, map_signif_level = F, test = t.test, textsize = 5) +
  theme_classic() +
  scale_fill_npg() + 
  theme_set(theme_classic(base_size = 20)) +
  guides(fill = F) +
  labs(x="Chronological age",y="scAge predicted age", title = "Muscle stem cells \n(Hernando et al. 2019)")
ggsave(filename = "./dataset/muscle_cell_stem_cell/figure/scAge_multi_omics_cell.pdf", width = 5, height = 5.5)

cor(SCALE_scAge_result$age_value, SCALE_scAge_result$PredictedAge)
cor(SCALE_scAge_result$age_value, SCALE_scAge_result$aging_scores)
cor(SCALE_scAge_result$aging_scores, SCALE_scAge_result$PredictedAge)

cor(SCALE_scAge_result$age_value, SCALE_scAge_result$scAge_scaled)
cor(SCALE_scAge_result$age_value, SCALE_scAge_result$aging_scores_scaled)
cor(SCALE_scAge_result$aging_scores_scaled, SCALE_scAge_result$scAge_scaled)

prob_set = 0.5
t.test(SCALE_scAge_result[SCALE_scAge_result$age_value == 2 & 
                            SCALE_scAge_result$aging_scores > quantile(SCALE_scAge_result[SCALE_scAge_result$age_value == 2, "aging_scores"], probs = prob_set), c("PredictedAge")],
       SCALE_scAge_result$PredictedAge)
t.test(SCALE_scAge_result[SCALE_scAge_result$age_value == 24 & 
                            SCALE_scAge_result$aging_scores > quantile(SCALE_scAge_result[SCALE_scAge_result$age_value == 24, "aging_scores"], probs = prob_set), c("PredictedAge")],
       SCALE_scAge_result$PredictedAge)
t.test(SCALE_scAge_result[SCALE_scAge_result$age_value == 2 & 
                            SCALE_scAge_result$PredictedAge > quantile(SCALE_scAge_result[SCALE_scAge_result$age_value == 2, "PredictedAge"], probs = prob_set), c("aging_scores")],
       SCALE_scAge_result$aging_scores)
t.test(SCALE_scAge_result[SCALE_scAge_result$age_value == 24 & 
                            SCALE_scAge_result$PredictedAge > quantile(SCALE_scAge_result[SCALE_scAge_result$age_value == 24, "PredictedAge"], probs = prob_set), c("aging_scores")],
       SCALE_scAge_result$aging_scores)

SCALE_scAge_result_2_month = SCALE_scAge_result[SCALE_scAge_result$age_value == 2,]
SCALE_scAge_result_24_month = SCALE_scAge_result[SCALE_scAge_result$age_value == 24,]

SCALE_scAge_result_2_month$aging_scores_scaled <- as.numeric(scale(SCALE_scAge_result_2_month$aging_scores))
SCALE_scAge_result_2_month$scAge_scaled <- as.numeric(scale(SCALE_scAge_result_2_month$PredictedAge))

SCALE_scAge_result_24_month$aging_scores_scaled <- as.numeric(scale(SCALE_scAge_result_24_month$aging_scores))
SCALE_scAge_result_24_month$scAge_scaled <- as.numeric(scale(SCALE_scAge_result_24_month$PredictedAge))

cor(SCALE_scAge_result_2_month$aging_scores_scaled, SCALE_scAge_result_2_month$scAge_scaled)
cor(SCALE_scAge_result_24_month$aging_scores_scaled, SCALE_scAge_result_24_month$scAge_scaled)

t.test(SCALE_scAge_result_2_month[SCALE_scAge_result_2_month$aging_scores_scaled > quantile(SCALE_scAge_result_2_month$aging_scores_scaled, probs = 0.7), c("scAge_scaled")],
       SCALE_scAge_result$scAge_scaled)
t.test(SCALE_scAge_result_24_month[SCALE_scAge_result_24_month$aging_scores_scaled > quantile(SCALE_scAge_result_24_month$aging_scores_scaled, probs = 0.7), c("scAge_scaled")],
       SCALE_scAge_result$scAge_scaled)

t.test(SCALE_scAge_result_2_month[SCALE_scAge_result_2_month$scAge_scaled > quantile(SCALE_scAge_result_2_month$scAge_scaled, probs = 0.7), c("aging_scores_scaled")],
       SCALE_scAge_result$aging_scores_scaled)
t.test(SCALE_scAge_result_24_month[SCALE_scAge_result_24_month$aging_scores_scaled > quantile(SCALE_scAge_result_24_month$aging_scores_scaled, probs = 0.7), c("aging_scores_scaled")],
       SCALE_scAge_result$aging_scores_scaled)

save.image("./dataset/muscle_cell_stem_cell/data_analysis.Rdata")
