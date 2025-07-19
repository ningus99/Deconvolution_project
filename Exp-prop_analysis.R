
## Cell type proportion analysis
## Real proportion calculation 

table_cheng <- matrix_dec@meta.data %>% group_by(sample, cell_types) %>% summarise(num_cells = n(), .groups = 'drop_last') %>%
mutate(percent = (num_cells / sum(num_cells))) %>%
ungroup()
cells <- c("Luminal", "Basal", "Fibroblast", "Macrophage", "T_cell", "Mastocyte", "Endothelial", "B_cell")
real_prop <- list()
for (n in cells){ 
real_prop[[n]] <- table_cheng %>% filter(cell_types == n) %>% dplyr::select(sample, percent) %>% 
rename(!!n := percent)
}

real_propor <- real_prop[[1]]
for (i in 2:length(real_prop)) {
  real_propor <- full_join(real_propor, real_prop[[i]], by = "sample")
}

real_propor[,-1] <- replace(real_propor[,-1], is.na(real_propor[,-1]), 0)
real_propor <- as.data.frame(real_propor)
real_propor$sample <- c(paste0("patient_", 1:13))
real_propor <- real_propor %>% mutate(type = "Real")

write.table(real_propor, file = "GT_Data/real_propor.txt")

## Data combination and long format transformation 

real_propor <- read.table("GT_Data/real_propor.txt", header = TRUE)

combined_data <- bind_rows(real_propor, predict_prop)
combined_data_long <- combined_data %>%
  pivot_longer(
    cols = -c(sample, type), 
    names_to = "cell_type", 
    values_to = "percentage" 
  )
combined_data_long <- combined_data_long %>%
mutate(sample_type = paste(sample, type, sep = "_")) %>% arrange(sample)

## Proportion comparison 

 p3 <- ggplot(combined_data_long, aes(x = sample, y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") + 
  coord_flip() +
  facet_wrap(~type) +
  labs(
    title = "Comparison of Cellular Proportions (True vs Predicted)",
    x = "Patients",
    y = "Proportion (%)",
    fill = "Cell Types"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8, margin = margin(t = 5)), 
    panel.grid.major = element_blank()
  )

ggsave("Results/Proportions/Proportions_comparison_exp2.pdf", plot = p3, width = 10, height = 6, dpi = 300)

## Correlation test 

predict_cor <- predict_prop[, !names(predict_prop) %in% c("sample", "type")]
real_cor <- real_propor[, !names(real_propor) %in% c("sample", "type")]

res_cor <- list()
for (tipo_celular in names(predict_cor)){
  res_cor[[tipo_celular]] <- cor.test(real_cor[[tipo_celular]], predict_cor[[tipo_celular]], method = "pearson")
}

res_cor <- do.call(rbind, res_cor)

write.table(res_cor, file = "Results/Proportions/res_cor_exp2.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

### Gene expresion analysis
## Gene expression extraction from BayesPrism output

expres <- list()

for (n in names(table(cell.type.labels))) {
  expres[[n]] <- get(paste0(n, "_exp")) %>%
    as.data.frame() %>%
    mutate(sample = paste0("patient", 1:13))
  
  rownames(expres[[n]]) <- expres[[n]]$sample
  expres[[n]] <- expres[[n]] %>% dplyr::select(-sample)
}

expr_tot <- do.call(rbind, expres)
expr_tot <- t(expr_tot)
expr_tot <- as.data.frame(expr_tot)

write.table(expr_tot, file = "sisana/data/Pred_exp_mat_exp.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)


## Calculation of real expression per patient per celltype

matrix_real_exp <- list()

for (n in 1:13) { 
  
  matrix_pat <- list()  
  
  for (m in names(table(cell.type.labels))) { 
    
    matrix_pat[[m]] <- matrix_dec@assays$RNA$counts[, 
      matrix_dec@meta.data$cell_name[
        matrix_dec@meta.data$cell_types == m & 
        matrix_dec@meta.data$sample == n
      ]
    ] %>% 
      as.data.frame() %>% 
      mutate(!!paste0(m, ".", "patient", n) := rowSums(as.matrix(.))) %>% 
      dplyr::select(!!paste0(m, ".", "patient", n))
  }
  
  
  matrix_real_exp[[n]] <- do.call(cbind, matrix_pat)
}

matrix_real_exp <- lapply(matrix_real_exp, function(x) replace(x, is.na(x), 0))

final_matrix <- as.data.frame(do.call(cbind, matrix_real_exp))
nan_cols <- which(colSums(final_matrix[,!names(final_matrix) %in% "Target"]) == 0)
colnames_with_nan <- colnames(final_matrix)[nan_cols]
colnames_with_nan <- colnames_with_nan[!colnames_with_nan %in% c("Target")]

## Remove non existing cell_types and filter common genes and genes with 0 expression 

final_matrix <- read.table('GT_Data/final_matrix_all.txt', header = TRUE, row.names = 1)
final_matrix$Target <- rownames(final_matrix)
expr_tot <- read.table('sisana/data/Pred_exp_mat_exp.tsv', header = TRUE, row.names = 1)
expr_tot$Target <- rownames(expr_tot)
nan_cols <- which(colSums(final_matrix[,!names(final_matrix) %in% "Target"]) == 0)
colnames_with_nan <- colnames(final_matrix)[nan_cols]
colnames_with_nan <- colnames_with_nan[!colnames_with_nan %in% c("Target")]
final_matrix <- final_matrix[,!colnames(final_matrix) %in% colnames_with_nan]
expr_tot <- expr_tot[,!colnames(expr_tot) %in% colnames_with_nan] 

genes_disc <- setdiff(rownames(expr_tot), rownames(final_matrix))
final_matrix <- final_matrix[!rownames(final_matrix) %in% genes_disc,]
genes_0 <- rowSums(expr_tot > 0) == 0 & rowSums(final_matrix > 0) == 0
final_matrix <- final_matrix[!genes_0,]
expr_tot <- expr_tot[!genes_0,]

## Format problem 

rownames(final_matrix) <- gsub("-", ".", rownames(final_matrix))
rownames(final_matrix) <- gsub("/", ".", rownames(final_matrix))
rownames(final_matrix) <- gsub("|", ".", rownames(final_matrix))
rownames(final_matrix) <- gsub("\\\\", ".", rownames(final_matrix))

identical(rownames(final_matrix), rownames(expr_tot))

matches1 <- sapply(genes_disc, function(gene) {
  grep(paste0("^", gene, "$"), rownames(final_matrix), value = TRUE)
})

matches <- unname(matches1)

final_matrix <- final_matrix[!rownames(final_matrix) %in% matches, ]
rownames(expr_tot) <- rownames(final_matrix)

## Depth sequencing normalization 

rownames(final_matrix) <- final_matrix$Target
rownames(expr_tot) <- expr_tot$Target

expr_tot <- expr_tot[, -1]
final_matrix <- final_matrix[, -1]

final_matrix <- sweep(final_matrix, 2, colSums(final_matrix), FUN = "/") * 1e6

expr_tot <- sweep(expr_tot, 2, colSums(expr_tot), FUN = "/") * 1e6

expr_tot$Target <- rownames(expr_tot)
expr_tot <- expr_tot[,c("Target", setdiff(names(expr_tot), "Target"))]

final_matrix$Target <- rownames(final_matrix)
final_matrix <- final_matrix[, c("Target", setdiff(names(final_matrix), "Target"))]

## Long format transformation

final_matrix <- read.table('sisana/data/Real_exp_mat_norm.tsv', header = TRUE)
expr_tot <- read.table('sisana/data/Pred_exp_mat_norm.tsv', header = TRUE)

real_mat <- final_matrix %>% as.data.frame() %>% pivot_longer(-Target, names_to = "Sample", values_to = "Real")
predict_mat <- expr_tot %>% as.data.frame() %>% pivot_longer(-Target, names_to = "Sample", values_to = "Predicted")
predict_mat <- predict_mat %>% mutate(type = "Predicted")
real_mat <- real_mat %>% mutate(type = "Real")
matrix_final <- inner_join(real_mat, predict_mat, by = c("Target", "Sample"))

## Calculation of several statistical metrics  

matrix_final <- matrix_final %>%
  separate(Sample, into = c("CellType", "Patient"), sep = "\\.", remove = FALSE)  
                                                                                   

metrics_by_patient_cell <- matrix_final %>% 
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

write.csv2(metrics_by_patient_cell, file = "metrics_by_patient_cell_exp.csv", row.names = FALSE)

## Cell Type proportion adittion  

real_prop <- read.table('real_propor.txt', header = TRUE)
predict_prop <- read.table('predict_prop.txt', header = TRUE)
real_prop[, sapply(real_prop, is.numeric)] <- real_prop[, sapply(real_prop, is.numeric)] * 100
predict_prop[, sapply(predict_prop, is.numeric)] <- predict_prop[, sapply(predict_prop, is.numeric)] * 100
real_prop <- real_prop %>% pivot_longer(cols = -c(sample, type), values_to = "real_prop", names_to = "CellType")
predict_prop <- predict_prop %>% pivot_longer(cols = -c(sample, type), values_to = "predict_prop", names_to = "CellType")
props <- inner_join(predict_prop, real_prop, by = c("CellType", "sample"))
props <- props[, !colnames(props) %in% c("type.x", "type.y")]
colnames(props)[colnames(props) == "sample"] <- "Patient"
metrics_by_patient_cell <- metrics_by_patient_cell[order(metrics_by_patient_cell$Patient, metrics_by_patient_cell$CellType),]
props <- props[order(props$Patient, props$CellType),]
props <- props %>% mutate(Patient = str_replace(Patient, "_", ""))
props <- props %>% unite("Pat_cell", CellType, Patient, sep = "." )
rownames(props) <- props$Pat_cell
props <- props[!rownames(props) %in% colnames_with_nan, ]
metrics_by_patient_cell$predict_prop <- props$predict_prop
metrics_by_patient_cell$real_prop <- props$real_prop  

## Analysis of not shared gens between predicted, real expression and signature matrix 

# Matrix genes extraction 

sc.dat <- read.table("BayesPrism/sc.dat.txt")
cell.type.labels <- read.table("BayesPrism/cell.type.labels.txt")
cell.type.labels$x <- make.unique(cell.type.labels$x)
rownames(sc.dat) <- cell.type.labels$x
sc.dat <- t(sc.dat)
sc.dat <- as.data.frame(sc.dat)
sc.dat$Gene <- rownames(sc.dat)

matrix_bulk_sign <- list()

for (m in names(table(cell.type.labels))) {  
    matrix_bulk_sign[[m]] <- sc.dat %>% as.data.frame %>% 
      dplyr::select(Gene, starts_with(m)) %>%  
      mutate(!!sym(m) := rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%  
      filter(!!sym(m) != 0) %>% dplyr::select(Gene)
  }

## Shared and unique expressed genes 

n_genes_int <- data.frame(gene_p = numeric(), gene_r = numeric(), int = numeric(), n_genes_p1 = numeric(), n_genes_r1 = numeric())
genes_only_p <- list()
genes_only_r <- list()
gene_p_in_sign = numeric()
gene_r_in_sign = numeric()
n_genes_tot = data.frame()
genes_only_p_no_sign <- list()

for (m in paste0("patient", 1:13)) { 
  for (n in names(table(cell.type.labels))) {  
    
   
  n_genes_p <- matrix_final$Target[matrix_final$Patient == m & matrix_final$CellType == n & matrix_final$Predicted != 0]
  n_genes_r <- matrix_final$Target[matrix_final$Patient == m & matrix_final$CellType == n & matrix_final$Real != 0] 

  genes_only_p[[m]][[n]] <- setdiff(n_genes_p, intersect(n_genes_p, n_genes_r))
  genes_only_r[[m]][[n]] <- setdiff(n_genes_r, intersect(n_genes_p, n_genes_r)) 
  genes_only_p_no_sign[[m]][[n]] <- setdiff(genes_only_p[[m]][[n]], matrix_bulk_sign[[n]]$Gene)

   # n_genes_int <- data.frame(
   #   gene_p = length(n_genes_p),
    #  gene_r = length(n_genes_r),
    #  int = length(intersect(n_genes_p, n_genes_r)),
    #  gene_p1 = length(n_genes_p) - length(intersect(n_genes_p, n_genes_r)),
    #  gene_r1 = length(n_genes_r) - length(intersect(n_genes_p, n_genes_r)),
    #  gene_p_in_sign = length(intersect(genes_only_p[[m]][[n]], matrix_bulk_sign[[n]]$Gene)),
    #  gene_r_in_sign = length(intersect(genes_only_r[[m]][[n]], matrix_bulk_sign[[n]]$Gene)),
    #  plot1 = ((length(n_genes_p) - length(intersect(n_genes_p, n_genes_r))) / (length(n_genes_p))) * 100, 
    #  plot2 = ((length(intersect(genes_only_p[[m]][[n]], matrix_bulk_sign[[n]]$Gene))) / (length(n_genes_p) - length(intersect(n_genes_p, n_genes_r)))) * 100,
    #  propor_predicted = length(n_genes_p) / length(n_genes_r),
    #  patient = m,
    #  celltype = n
    #)
   
    #n_genes_tot <- rbind(n_genes_tot, n_genes_int)
  }
}

## Unification of all the information

n_genes_tot <- n_genes_tot %>% drop_na()
n_genes_tot <- n_genes_tot[order(n_genes_tot$patient, n_genes_tot$celltype),]
n_genes_tot <- n_genes_tot %>% unite("Pat_cell", celltype, patient, sep = "." )
n_genes_tot <- n_genes_tot[!duplicated(n_genes_tot$Pat_cell), ]
rownames(n_genes_tot) <- n_genes_tot$Pat_cell
n_genes_tot <- n_genes_tot[!rownames(n_genes_tot) %in% colnames_with_nan, ]
metrics_by_patient_cell$prop_pred_noreal <- n_genes_tot$plot1
metrics_by_patient_cell$prop_pred_noreal_sign <- n_genes_tot$plot2 
metrics_by_patient_cell$gene_p <- n_genes_tot$gene_p 
metrics_by_patient_cell$gene_only_p <- n_genes_tot$gene_p1
metrics_by_patient_cell$gene_p_in_sign <- n_genes_tot$gene_p_in_sign
metrics_by_patient_cell$gene_proport <- n_genes_tot$propor_predicted

write.table(metrics_by_patient_cell, file = "Results/Exp_results/metrics_by_patient_cell_exp.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

## Plot construction 
## Reggresion plots 

for (n in unique(matrix_final$CellType)) {

  p4 <- ggplot(matrix_final[matrix_final$CellType == n,], aes(x = Real, y = Predicted)) +
    geom_point(alpha = 0.6, color = "salmon") + 
    geom_smooth(method = "lm", se = FALSE, color = "black") +  
    stat_poly_eq(
    formula = y ~ x,
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    parse = TRUE,
    label.x = "right", label.y = "bottom"
  ) 
    labs(title = paste("Linear regression: Real vs. Predicted -", n),
         x = "Real", y = "Predicted") +
    theme_minimal()

  ggsave(paste0("REGRESSION_", n, ".pdf"), 
         plot = p4, width = 10, height = 6)
}


## Barplots of not shared and shared genes between the three matrixes 

metrics_by_patient_cell_long <- metrics_by_patient_cell %>%
  pivot_longer(cols = c(predict_prop, gene_p, gene_only_p),
               names_to = "Metric", values_to = "Value")


metrics_by_patient_cell_long$Metric <- factor(metrics_by_patient_cell_long$Metric, levels = c("gene_only_p", "gene_p", "predict_prop"))

custom_colors <- c("gene_p" = "#1f78b4",         
                   "gene_only_p" = "#33a02c", 
                   "gene_p_in_sign" = "#e31a1c", 
                   "predict_prop" = "#ff7f00")  

metrics_by_patient_cell_long <- metrics_by_patient_cell_long %>%
  mutate(
    Scaled_Value = case_when(
      Metric == "predict_prop" ~ -Value * 300, 
      TRUE ~ Value  
    )
  )


p1 <- ggplot(metrics_by_patient_cell_long, aes(x = CellType, y = Scaled_Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "stack", data = metrics_by_patient_cell_long %>% 
             filter(Metric %in% c("gene_p", "gene_only_p"))) +
  geom_bar(stat = "identity", position = "identity", data = metrics_by_patient_cell_long %>% 
             filter(Metric == "predict_prop")) +
  scale_y_continuous(labels = abs) +
  facet_wrap(~ Patient) + 
  scale_fill_manual(values = custom_colors) +
  labs(x = "Cell Type", y = "Proportion", fill = "Metric") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

  metrics_by_patient_cell_long <- metrics_by_patient_cell %>%
  pivot_longer(cols = c(predict_prop, gene_only_p, gene_p_in_sign),
               names_to = "Metric", values_to = "Value")

metrics_by_patient_cell_long$Metric <- factor(metrics_by_patient_cell_long$Metric, levels = c("gene_p_in_sign", "gene_only_p", "predict_prop"))
  metrics_by_patient_cell_long <- metrics_by_patient_cell_long %>%
  mutate(
    Scaled_Value = case_when(
      Metric == "predict_prop" ~ -Value * 300, 
      TRUE ~ Value  
    )
  )

p2 <- ggplot(metrics_by_patient_cell_long, aes(x = CellType, y = Scaled_Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "stack", data = metrics_by_patient_cell_long %>% 
             filter(Metric %in% c("gene_only_p", "gene_p_in_sign"))) +
  geom_bar(stat = "identity", position = "identity", data = metrics_by_patient_cell_long %>% 
             filter(Metric == "predict_prop")) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(labels = abs) +  
  facet_wrap(~ Patient) + 
  labs(x = "Cell Type", y = "Proportion", fill = "Metric") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

ggsave("bar_plot_by_patient_pred.pdf", plot = p1, width = 12, height = 6)
ggsave("bar_plot_by_patient_sign.pdf", plot = p2, width = 12, height = 6)

## UMAP construction   
## Predicted expression 

sample_names <- colnames(expr_tot)


CellType <- sub("\\..*$", "", sample_names)


Paciente <- sub("^.*\\.", "", sample_names)


metadata <- data.frame(Sample = sample_names, CellType = CellType, Paciente = Paciente)


matrix_pred <- CreateSeuratObject(counts = expr_tot, meta.data = metadata)
matrix_pred <- NormalizeData(matrix_pred) 
matrix_pred <- FindVariableFeatures(matrix_pred, selection.method = "vst", nfeatures = 2000)
matrix_pred <- ScaleData(matrix_pred, features = VariableFeatures(object = matrix_pred))
matrix_pred <- RunPCA(matrix_pred, features = VariableFeatures(object = matrix_pred))
PCA_pred <- DimPlot(matrix_pred, reduction = "harmony", group.by = "CellType")
matrix_pred <- RunHarmony(matrix_pred, group.by.vars = "Paciente") 
matrix_pred <- RunUMAP(matrix_pred, reduction = "harmony", dims = 1:30)
UMAP_pred <- DimPlot(matrix_pred, reduction = "umap", group.by = "CellType")

ggsave("pca_pred.pdf", plot = PCA_pred, width = 12, height = 6)
ggsave("umap_pred.pdf", plot = UMAP_pred, width = 12, height = 6)

## Real expression 

sample_names <- colnames(final_matrix)


CellType <- sub("\\..*$", "", sample_names)


Paciente <- sub("^.*\\.", "", sample_names)


metadata <- data.frame(Sample = sample_names, CellType = CellType, Paciente = Paciente)
matrix_real <- CreateSeuratObject(counts = final_matrix, meta.data = metadata)
matrix_real <- NormalizeData(matrix_real) 
matrix_real <- FindVariableFeatures(matrix_real, selection.method = "vst", nfeatures = 2000)
matrix_real <- ScaleData(matrix_real, features = VariableFeatures(object = matrix_pred))
matrix_real <- RunPCA(matrix_real, features = VariableFeatures(object = matrix_pred))
PCA_real <- DimPlot(matrix_real, reduction = "harmony", group.by = "CellType")
matrix_real <- RunHarmony(matrix_real, group.by.vars = "Paciente") 
matrix_real <- RunUMAP(matrix_real, reduction = "harmony", dims = 1:30)
UMAP_real <- DimPlot(matrix_real, reduction = "umap", group.by = "CellType")

ggsave("pca_real.pdf", plot = PCA, width = 12, height = 6)
ggsave("umap_real.pdf", plot = UMAP, width = 12, height = 6)

## Test counts in genes only predicted 

ggplot(data.frame(expr_tot[rownames(expr_tot) %in% genes_only_p[["patient10"]][["B_cell"]],]), aes(x = Luminal.patient1)) +
  geom_density(fill = "blue", alpha = 0.5) +  
  theme_minimal() + 
  labs(title = "Gráfico de Densidad", x = "Valor", y = "Densidad")

