
### Network analysis 
## Preprocessing 

panda_output_pred <- read.table("sisana/network/Predicted/norm/results/panda_output.txt")
raw_edges_pred <- read.csv("raw_edges_pred.csv", sep = ",")
Pred_net_ind <- read.csv("sisana/network/Predicted/norm/results/lioness_transformed_edges_indegree.csv", sep = ",")

panda_output_real <- read.table("sisana/network/Real/norm/results/panda_output.txt")
raw_edges_real <- read.csv("raw_edges_real.csv", sep = ",")
Real_net_ind <- read.csv("sisana/network/Real/norm/results/lioness_transformed_edges_indegree.csv", sep = ",")

## Data extraction from sisana 

indegrees_pred <- read.csv("sisana/network_filt/Predicted/results/lioness_indegree.csv", header = TRUE)
indegrees_real <- read.csv("sisana/network_filt/Real/results/lioness_indegree.csv", header = TRUE)

indegrees_pred <- indegrees_pred %>% rename(Target = target)
indegrees_real <- indegrees_real %>% rename(Target = target)


## Calculation and quantile normalization of indegrees 
## Real Network

raw_edges_real <- raw_edges_real %>% dplyr::select(-X)

colnames <- colnames(Real_net_ind[,2:77]) 

colnames(raw_edges_real) <- colnames 


raw_edges_real$TF <- panda_output_real$V1

raw_edges_real$Target <- panda_output_real$V2

## Indegree calculation 

indegrees_real <- raw_edges_real %>% group_by(Target) %>% summarise(across(!c(TF), sum))
indegrees_real <- indegrees_real[,!colnames(indegrees_real) %in% colnames_filt]
targets <- indegrees_real$Target

## Quantile normalization

indegrees_real <- read.csv("sisana/network_filt/Real/results/lioness_indegree.csv", header = TRUE)

colnames(indegrees_real) <- c("Target", colnames(indegrees_real)[-1])

group_labels <- sapply(strsplit(colnames(indegrees_real), "\\."), function(x) x[1]) 

qs <- qsmooth(object = as.matrix(indegrees_real[,-1]), group_factor = group_labels[-1])

indegrees_real <- as.data.frame(qsmoothData(qs))

indegrees_real$Target <- targets 

indegrees_real <- write.table(indegrees_real, "sisana/network/Real/results/indegrees_real.txt", sep = "\t", row.names = FALSE, col.names = TRUE) 


## Predicted 

raw_edges_pred <- raw_edges_pred %>% dplyr::select(-X)

colnames <- colnames(Pred_net_ind[,2:77]) 

colnames(raw_edges_pred) <- colnames 


raw_edges_pred$TF <- panda_output_pred$V1

raw_edges_pred$Target <- panda_output_pred$V2

## Indegree calculation 

indegrees_pred <- raw_edges_pred %>% group_by(Target) %>% summarise(across(!c(TF), sum))
indegrees_pred <- indegrees_pred[,!colnames(indegrees_pred) %in% colnames_filt]
targets <- indegrees_pred$Target

## Quantile normalization

indegrees_pred <- read.csv("sisana/network_filt/Predicted/results/lioness_indegree.csv", header = TRUE)

colnames(indegrees_pred) <- c("Target", colnames(indegrees_pred)[-1]) 

group_labels <- sapply(strsplit(colnames(indegrees_pred), "\\."), function(x) x[1]) 

qs <- qsmooth(object = as.matrix(indegrees_pred[,-1]), group_factor = group_labels[-1])

indegrees_pred <- as.data.frame(qsmoothData(qs))

indegrees_pred$Target <- targets 

indegrees_pred <- write.table(indegrees_pred, "sisana/network/Predicted/results/indegrees_pred.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

## Long format transformation

com <- intersect(indegrees_pred$Target, indegrees_real$Target)
real_net <- indegrees_real %>% as.data.frame() %>%
pivot_longer(-Target, names_to = "Sample", values_to = "Real")
predict_net <- indegrees_pred %>% as.data.frame() %>% 
pivot_longer(-Target, names_to = "Sample", values_to = "Predicted")
predict_net <- predict_net %>% mutate(type = "Predicted")
real_net <- real_net %>% mutate(type = "Real")
matrix_fin_net <- inner_join(real_net, predict_net, by = c("Target", "Sample"))

matrix_fin_net <- matrix_fin_net %>%
  separate(Sample, into = c("CellType", "Patient"), sep = "\\.", remove = FALSE) 


## Summary 

stats_real <- summary(matrix_fin_net$Real)
stats_pred <- summary(matrix_fin_net$Predicted)

## Regression 

for (n in unique(matrix_fin_net$CellType)) {

  reg <- ggplot(matrix_fin_net[matrix_fin_net$CellType == n,], aes(x = Real, y = Predicted)) +
    geom_point(alpha = 0.6, color = "salmon") + 
    geom_smooth(method = "lm", se = FALSE, color = "black") +  
    stat_poly_eq(
    formula = y ~ x,
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    parse = TRUE,
    label.x = "right", label.y = "bottom"
  ) +
    labs(title = paste("Linear regression: Real vs. Predicted -", n),
         x = "Real", y = "Predicted") +
    theme_minimal()  +
    theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

  ggsave(paste0("Results/Networks/regression/REGRESSION_", n, "_net.pdf"), 
         plot = reg, width = 10, height = 6)
}


## Calculation of general statistical metrics 

metrics_by_patient_cell_net <- matrix_fin_net %>% 
  group_by(Patient, CellType) %>% 
  summarise(
    cor_pear = cor(Real, Predicted, method = "pearson", use = "complete.obs"),
    cor_spear = cor(Real, Predicted, method = "spearman", use = "complete.obs"),
    pearson_pval = cor.test(Real, Predicted, method = "pearson", use = "complete.obs")$p.value,
    spearman_pval = cor.test(Real, Predicted, method = "spearman", use = "complete.obs")$p.value,
    MSE = mean((Real - Predicted)^2, na.rm = TRUE),  
    RMSE = sqrt(MSE), 
    MAE = mean(abs(Real - Predicted), na.rm = TRUE),
    mean_r = mean(Real, na.rm = TRUE),
    mean_p = mean(Predicted, na.rm = TRUE),
    sd_r = sd(Real, na.rm = TRUE),
    sd_p = sd(Predicted, na.rm = TRUE)
  ) %>% 
  ungroup()

write.csv2(metrics_by_patient_cell_net, "metrics_by_patient_cell_net.csv", row.names = FALSE)

## PCA and UMAP construction  
## Predicted network 

indegrees_real <- as.data.frame(indegrees_real)
indegrees_pred <- as.data.frame(indegrees_pred)

rownames(indegrees_real) <- indegrees_real$Target
rownames(indegrees_pred) <- indegrees_pred$Target

gene_variance <- apply(indegrees_real[, -ncol(indegrees_real)], 1, var)

top_genes <- names(sort(gene_variance, decreasing = TRUE))[1:3000]

pred_pca <-  indegrees_real[top_genes,]
pred_pca <- t(pred_pca[, -ncol(indegrees_real)])
pred_pca <- as.matrix(pred_pca)
pred_pca <- scale(pred_pca)

## Add read number variability to metadata

n_reads_pred <- expr_tot[,-1] %>%          
  summarise(across(everything(), sum)) %>% 
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "n_reads")

n_reads_real <- final_matrix[,-1] %>%          
  summarise(across(everything(), sum)) %>% 
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "n_reads") 

metadata <- data.frame(Sample = rownames(pred_pca)) %>%
  mutate(CellType = str_split_fixed(Sample, "\\.", 2)[,1])  %>% 
  mutate(Patient = str_split_fixed(Sample, "\\.", 2)[,2]) 

metadata <- merge(metadata, n_reads_pred, by = "Sample")
metadata <- merge(metadata, n_reads_real, by = "Sample")

pca_result <- prcomp(pred_pca, center = TRUE, scale. = FALSE)


pca_harmony <- harmony::RunHarmony(
  data_mat = pca_result$x,  
  meta_data = metadata,   
  vars_use = "Patient"    
)

pca_df <- data.frame(pca_result$x)

pca_df$Sample <- rownames(pca_df)


pca_df <- merge(pca_df, metadata, by = "Sample")

p7 <- ggplot(pca_df, aes(PC1, PC2, color = CellType)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA: color = CellType") +
  theme_minimal()

ggsave("Results/Networks_filt/regression/pca_pred_cell(1-2).pdf", plot = p7, width = 10, height = 6)

pca_reduced <- pca_df[, 2:8]

umap_result <- umap(pca_reduced, n_components = 4)

df_umap <- as.data.frame(umap_result$layout)

df_umap$Sample <- pca_df$Sample


df_umap <- merge(df_umap, metadata, by = "Sample")

  p5 <- ggplot(df_umap, aes(V1, V2, color = CellType, label = CellType)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "UMAP") +
  theme_minimal()

ggsave("umap_real_cell.pdf", plot = p5, width = 10, height = 6)

## Ranking comparison, consistency between genes with more difference in exp and ind 
## Ranking creation in the indegrees 

matrix_fin_net <- matrix_fin_net %>% group_by(Patient, CellType) %>%
 mutate(ranking_ind_real = rank(-Real, ties.method = "average")) %>% 
 mutate(ranking_ind_pred = rank(-Predicted, ties.method = "average")) %>%
 mutate(diff_rank_ind = ranking_ind_real - ranking_ind_pred) %>% 
 ungroup()
  
## Selection of the top genes

top_genes_ind <- matrix_fin_net %>%
  group_by(Patient, CellType) %>%
  arrange(desc(abs(diff_rank_ind))) %>% slice_head(n = 100) %>%
  ungroup() %>%
  dplyr::select(Target, Patient, CellType)

## Ranking creation in the expression 

expr_tot <- read.table("sisana/data/Pred_exp_mat_norm.tsv", header = TRUE, sep = "\t")
final_matrix <- read.table("sisana/data/Real_exp_mat_norm.tsv", header = TRUE, sep = "\t")

## Select same genes due to the sisana preprocess 

matrix_final <- matrix_final[matrix_final$Target %in% com, ]
matrix_final <- matrix_final %>% group_by(Patient, CellType) %>%
mutate(ranking_exp_real = rank(-Real, ties.method = "average")) %>% 
mutate(ranking_exp_pred = rank(-Predicted, ties.method = "average")) %>% 
mutate(diff_rank_exp = ranking_exp_real - ranking_exp_pred) %>%
mutate(Residues = Real - Predicted) %>%
ungroup()  

## Filter on the top genes 

matrix_fin_net <- inner_join(matrix_fin_net, top_genes_ind, by = c("Target", "Patient", "CellType"))
matrix_final   <- inner_join(matrix_final,   top_genes_ind, by = c("Target", "Patient", "CellType"))

data_comp <- inner_join(matrix_final, matrix_fin_net, by = c("Patient", "CellType", "Target")) %>%
mutate(diff_rank_exp_ind = (diff_rank_exp - diff_rank_ind)) 

data_comp <- data_comp %>% mutate(dif_ranking = ranking_exp - ranking_ind) %>% mutate(group = paste(CellType, Patient, sep = "_")) %>% 
dplyr::select(Target, Patient, CellType, group, diff_rank_exp, diff_rank_ind, diff_rank_exp_ind) 

correl_exp_ind <- data_comp %>% group_by(Patient, CellType) %>% 
summarise(cor_spear = cor(diff_rank_exp, diff_rank_ind, method = "spearman"),
          kendall = cor(diff_rank_exp, diff_rank_ind, method = "kendall"))
        

## Heatmap rankings 

data_heatmap <- data_comp %>% filter(CellType == "Luminal", Patient == "patient12") %>% arrange(diff_rank_ind) %>% 
mutate(Target = factor(Target, levels = pull(., Target)))

data_long <- data_heatmap %>%
  pivot_longer(
    cols = c(diff_rank_exp, diff_rank_ind, diff_rank_exp_ind),
    names_to = "rankingType",
    values_to = "rankingValue"
  ) %>%
  mutate(rankingType = factor(rankingType, levels = c("diff_rank_exp", "diff_rank_ind", "dif_rank_exp_ind")))

p2 <- ggplot(data_long, aes(x = grupo, y = Target, fill = rankingValue)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "Ranking"
  ) +
  facet_wrap(~ rankingType, ncol = 3) +
  labs(
    title = "Ranking heatmaps",
    x = "Group (CellType_Patient)",
    y = "Gene"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6)
  )

ggsave("Results/Networks/ranks/heatmap_ranking_luminal12_res.pdf", plot = p2, width = 10, height = 6, device = cairo_pdf)


## Analysis of consistencies between celltypes

rank_ct <-  matrix_fin_net %>%
  dplyr::select(Target, Patient, CellType, diff_rank_ind) %>%
  group_by(CellType) %>%
  group_split()

 get_spearman_corr <- function(df_celltype) {
  mat <- df_celltype %>% dplyr::select(-CellType) %>% 
    pivot_wider(names_from = Patient, values_from = diff_rank_ind) %>%
    column_to_rownames("Target") %>% 
    as.matrix()
  
  cor(mat, method = "spearman", use = "pairwise.complete.obs")
}

cor_matrices <- map(rank_ct, get_spearman_corr)

names(cor_matrices) <- map_chr(rank_ct, ~ unique(.x$CellType))

my_colors <- colorRampPalette(c("blue", "white", "red"))(200)

## Correlation matrixes

for (n in names(cor_matrices)) {
  
  pdf(file = paste0("Correlation_", n, ".pdf"),
      width = 10, height = 6)
  
 
  par(mar = c(1, 1, 4, 1))
  
  
  corrplot(
    cor_matrices[[n]],
    method = "color",
    type = "full",
    col = my_colors,
    tl.col = "black",
    tl.cex = 0.8,
    addCoef.col = "black"
  )


  dev.off()
}

## Correlation between genes with highest residuals and genes with expression not registered in the signature matrix or the real matrix 

genes_ord <- matrix_fin_net %>% ungroup() %>% filter(CellType == "T_cell", Patient == "patient10") %>%
  arrange(desc(abs(diff_rank_ind))) %>% dplyr::select(Target)

genes_ord <- as.vector(genes_ord$Target)
genes_ord_mat <- genes_only_p_no_sign[["patient10"]][["T_cell"]]

pos <- which(genes_ord %in% genes_ord_mat)


## Same with expression 

rank_ct <- matrix_final %>%
  dplyr::select(Target, Patient, CellType, diff_rank_exp) %>%
  group_by(Patient) %>%
  group_split()

 get_spearman_corr <- function(df_celltype) {
  mat <- df_celltype %>% dplyr::select(-Patient) %>% 
    pivot_wider(names_from = CellType, values_from = diff_rank_exp) %>%
    column_to_rownames("Target") %>% 
    as.matrix()
  
  cor(mat, method = "spearman", use = "pairwise.complete.obs")
}

cor_matrices <- map(rank_ct, get_spearman_corr)

names(cor_matrices) <- map_chr(rank_ct, ~ unique(.x$Patient))

my_colors <- colorRampPalette(c("blue", "white", "red"))(200)

## Correlation matrixes

for (n in names(cor_matrices)) {
  
  pdf(file = paste0("Correlation_exp_", n, ".pdf"),
      width = 10, height = 6)
  
 
  par(mar = c(1, 1, 4, 1))
  
  
  corrplot(
    cor_matrices[[n]],
    method = "color",
    type = "full",
    col = my_colors,
    tl.col = "black",
    tl.cex = 0.8,
    addCoef.col = "black"
  )


  dev.off()
}
