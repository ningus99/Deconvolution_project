
## Network enrichment analysis 
## Manual rankings 

results <- list()

matrix_fin_net_gsea <- matrix_fin_net %>% pivot_longer(cols = c(Real, Predicted), names_to = "Net_type", values_to = "value_ind") %>% 
  mutate(Net_type = factor(Net_type, levels = c("Real", "Predicted"))) 

for (n in unique(matrix_fin_net$CellType)) {

df_ct <- matrix_fin_net_gsea %>% filter(CellType == n) 
df_all <- matrix_fin_net_gsea %>% filter(CellType != n)

avg_ct <- df_ct %>% group_by(Net_type, Target) %>% summarise(avg_ct = median(value_ind), .groups = "drop") 

avg_all <- df_all %>% group_by(Net_type, Target) %>% summarise(avg_all = median(value_ind), .groups = "drop")

result <- avg_ct %>% inner_join(avg_all, by = c("Net_type", "Target")) %>% mutate(Celltype = n, residue = avg_ct - avg_all)

results[[n]] <- result 
} 

results_res <- bind_rows(results) 

results_res <- results_res %>% group_by(Net_type, Celltype) %>%  mutate(ranking_ind = rank(-residue, ties.method = "average")) %>% ungroup()

## Limma ranking calculation 
## Limma analysis 
## Predicted network 

rownames(indegrees_pred) <- indegrees_pred$Target

indegrees_pred <- indegrees_pred[, -length(colnames(indegrees_pred))]

metadata <- data.frame(Sample = colnames(indegrees_pred)) %>%
  mutate(CellType = str_split_fixed(Sample, "\\.", 2)[,1])  %>% 
  mutate(Patient = str_split_fixed(Sample, "\\.", 2)[,2]) 

metadata <- metadata[match(colnames(indegrees_pred), metadata$Sample), ] 

indegrees_pred <- as.matrix(indegrees_pred)

CellTypes <- factor(metadata$CellType)

design <- model.matrix(~ 0 + CellTypes)

colnames(design) <- levels(CellTypes)

fit <- lmFit(indegrees_pred, design)

types <- colnames(design)

rank_list <- list()

for (n in types) {
  others <- setdiff(types, n)
  
  contrast_expr <- paste0(n, " - (", paste(others, collapse = " + "), ")/", length(others))
  contrast_matrix <- makeContrasts(contr = contrast_expr, levels = design)
  
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  tab <- topTable(fit2, coef = 1, number = Inf, sort.by = "t")
  
  rank_list[[n]] <- tab
}

## Real network 

rownames(indegrees_real) <- indegrees_real$Target

indegrees_real <- indegrees_real[, -length(colnames(indegrees_real))]

metadata <- metadata[match(colnames(indegrees_real), metadata$Sample), ]

indegrees_real <- as.matrix(indegrees_real)

CellTypes <- factor(metadata$CellType)

design <- model.matrix(~ 0 + CellTypes)

colnames(design) <- levels(CellTypes)

fit <- lmFit(indegrees_real, design)

types <- colnames(design)

rank_list_real <- list()

for (n in types) {
  others <- setdiff(types, n)
  
  contrast_expr <- paste0(n, " - (", paste(others, collapse = " + "), ")/", length(others))
  contrast_matrix <- makeContrasts(contr = contrast_expr, levels = design)
  
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  tab <- topTable(fit2, coef = 1, number = Inf, sort.by = "t")
  
  rank_list_real[[n]] <- tab
}

convert_results_list_to_df <- function(results_list, net_type_label) {
  lapply(names(results_list), function(CellType) {
    df <- results_list[[CellType]]
    tibble(
      Target = rownames(df),
      Celltype = CellType,
      t_stat = df$t,
      FDR = df$adj.P.Val
    )
  }) |> bind_rows() |> mutate(Net_type = net_type_label)
}

df_real <- convert_results_list_to_df(rank_list_real, "Real")

df_pred <- convert_results_list_to_df(rank_list, "Predicted")

results_res <- bind_rows(df_real, df_pred)

results_res <- results_res %>% rename(residue = t_stat)

## gsea object creation 

msigdb_BP_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
msigdb_MF_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")

msigdb_BP <- split(msigdb_BP_df$gene_symbol, msigdb_BP_df$gs_name)
msigdb_MF <- split(msigdb_MF_df$gene_symbol, msigdb_MF_df$gs_name)

run_fgsea <- function(celltype, results_df, pathways_bp, pathways_mf, Inf_type) {

     ranks_real <- results_df %>%
    filter(Celltype == celltype, .data[[Inf_type]] == "Real") %>%
    dplyr::select(Target, residue) %>%
    deframe()

    ranks_pred <- results_df %>%
    filter(Celltype == celltype, .data[[Inf_type]] == "Predicted") %>%
    dplyr::select(Target, residue) %>%
    deframe()

  fgsea_real_bp <- fgsea(pathways = pathways_bp, stats = ranks_real, minSize = 10, maxSize = 500, nproc = 1)
  fgsea_pred_bp <- fgsea(pathways = pathways_bp, stats = ranks_pred, minSize = 10, maxSize = 500, nproc = 1)
  fgsea_real_mf <- fgsea(pathways = pathways_mf, stats = ranks_real, minSize = 10, maxSize = 500, nproc = 1)
  fgsea_pred_mf <- fgsea(pathways = pathways_mf, stats = ranks_pred, minSize = 10, maxSize = 500, nproc = 1)
  
  
  fgsea_real_bp$Celltype <- celltype; fgsea_real_bp[[Inf_type]] <- "Real"; fgsea_real_bp$Ontology <- "BP"
  fgsea_pred_bp$Celltype <- celltype; fgsea_pred_bp[[Inf_type]]  <- "Predicted"; fgsea_pred_bp$Ontology <- "BP"
  fgsea_real_mf$Celltype <- celltype; fgsea_real_mf[[Inf_type]]  <- "Real"; fgsea_real_mf$Ontology <- "MF"
  fgsea_pred_mf$Celltype <- celltype; fgsea_pred_mf[[Inf_type]]  <- "Predicted"; fgsea_pred_mf$Ontology <- "MF"
  bind_rows(fgsea_real_bp, fgsea_pred_bp, fgsea_real_mf, fgsea_pred_mf)
} 

## Object creation 
## enrichment in each network 

for (n in c("T_cell", "Fibroblast", "Basal", "Luminal", "Endothelial", "Macrophage")) {

assign(paste0("fgsea_", n, "_net"), run_fgsea(celltype = n, results_res, msigdb_BP, msigdb_MF, Inf_type = "Net_type"))
} 

## Enrichment comparison between Real and Predicted networks

plot_data_num <- list()
plot_data <- list()

for (n in c("T_cell", "Fibroblast", "Basal", "Luminal", "Endothelial", "Macrophage")) {

fgsea_net <- get(paste0("fgsea_", n, "_net")) %>% mutate(Sign = ifelse(NES > 0, "Positive", "Negative"))

plot_real <- fgsea_net %>%
  filter(Net_type == "Real", Ontology == "BP", padj < 0.01) %>%
  mutate(logFDR = -log10(padj)) %>% filter(NES > 1 | NES < -1) %>%
  group_by(Sign) %>%
  arrange(desc(abs(NES)))  %>%
  mutate(NAME = factor(pathway, levels = pathway)) %>% ungroup ()

plot_pred <- fgsea_net %>%
  filter(Net_type == "Predicted", Ontology == "BP", padj < 0.01) %>%
  mutate(logFDR = -log10(padj)) %>% 
  group_by(Sign) %>%
  arrange(desc(abs(NES))) %>%
  mutate(NAME = factor(pathway, levels = pathway)) %>% ungroup()

com_path <- intersect(plot_real$pathway, plot_pred$pathway) 

plot_com_num <- data.frame(path_pred = length(plot_pred$pathway), path_real = length(plot_real$pathway), path_com = length(com_path), cell_type = n)

plot_com <- bind_rows(plot_real, plot_pred) %>% filter(pathway %in% com_path) %>% group_by(Net_type, Sign) %>% 
mutate(Net_type = factor(Net_type, levels = c("Real", "Predicted"))) %>% group_by(Net_type, Sign) %>% arrange(desc(abs(NES))) %>% 
slice_head (n =20) %>% arrange(logFDR) %>% mutate(
          pathway = gsub("_", " ", as.character(pathway)),
          pathway = str_wrap(pathway, width = 50, whitespace_only = FALSE), 
          pathway = gsub(" ", "_", pathway)
            ) %>%  mutate(NAME = factor(pathway, levels = unique(pathway))) %>% ungroup()

plot_data[[n]] <- plot_com

plot_data_num[[n]] <- plot_com_num

    p9 <- ggplot(plot_com, aes(x = logFDR, y = NAME, size = size, fill = NES)) +
   geom_point(shape = 21, color = "grey20") +
   scale_size_area(max_size = 10) +
   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
   theme_minimal() +
   facet_wrap(Sign ~ Net_type, scales = "free_y") +
   labs(x = "−log10(FDR)", y = NULL, title = paste0("Common pathways in ", n)) +
   theme_minimal() +
   theme(
     axis.text.y = element_text(size = 12),        
     axis.text.x = element_text(size = 12),         
     axis.title.x = element_text(size = 12),        
     plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
     strip.text = element_text(size = 12, face = "bold"),               
     legend.text = element_text(size = 12),         
     legend.title = element_text(size = 12)        
   )

ggsave(paste0("GSEA_common_filt", n, ".pdf"), width = 32, height = 16, plot = p9)
}

plot_data_num <- bind_rows(plot_data_num)

## plot total number of pathways 

plot_data_num <- plot_data_num %>%
  pivot_longer(cols = c(path_real, path_pred, path_com),
               names_to = "Group", values_to = "Count") %>%
  mutate(Group = recode(Group,
                        path_real = "Real",
                        path_pred = "Predicted",
                        path_com  = "Common"))

p8 <- ggplot(plot_data_num, aes(x = cell_type, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(
    title = "Number of Significant Pathways per Cell Type",
    x = "Cell Type",
    y = "Number of Pathways",
    fill = "Source"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14)

  )

ggsave("GSEA_common_filt_net.pdf", width = 12, height = 8, plot = p8)

## Enrichment analysis for expression 

results <- list()

matrix_final_gsea <- matrix_final %>% pivot_longer(cols = c(Real, Predicted), names_to = "Exp_type", values_to = "value_exp") %>% 
  mutate(Exp_type = factor(Exp_type, levels = c("Real", "Predicted"))) 
  
for (n in unique(matrix_final$CellType)) {

df_ct <- matrix_final_gsea %>% filter(CellType == n) 
df_all <- matrix_final_gsea %>% filter(CellType != n)

avg_ct <- df_ct %>% group_by(Exp_type, Target) %>% summarise(avg_ct = mean(value_exp), .groups = "drop") 

avg_all <- df_all %>% group_by(Exp_type, Target) %>% summarise(avg_all = mean(value_exp), .groups = "drop")

result <- avg_ct %>% inner_join(avg_all, by = c("Exp_type", "Target")) %>% mutate(Celltype = n, residue = avg_ct - avg_all)

results[[n]] <- result 
} 

results_res <- bind_rows(results) 

results_res <- results_res %>% group_by(Exp_type, Celltype) %>%  mutate(ranking_exp = rank(-residue, ties.method = "average")) %>% ungroup()

## Limma analysis 
## Predicted exp

rownames(expr_tot) <- expr_tot$Target

expr_tot <- expr_tot[, -1]

metadata <- data.frame(Sample = colnames(expr_tot)) %>%
  mutate(CellType = str_split_fixed(Sample, "\\.", 2)[,1])  %>% 
  mutate(Patient = str_split_fixed(Sample, "\\.", 2)[,2]) 

metadata <- metadata[match(colnames(expr_tot), metadata$Sample), ] 

expr_tot <- as.matrix(expr_tot)

CellTypes <- factor(metadata$CellType)

design <- model.matrix(~ 0 + CellTypes)

colnames(design) <- levels(CellTypes)

fit <- lmFit(expr_tot, design)

types <- colnames(design)

rank_list <- list()

for (n in types) {
  others <- setdiff(types, n)
  
  contrast_expr <- paste0(n, " - (", paste(others, collapse = " + "), ")/", length(others))
  contrast_matrix <- makeContrasts(contr = contrast_expr, levels = design)
  
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  tab <- topTable(fit2, coef = 1, number = Inf, sort.by = "t")
  
  rank_list[[n]] <- tab
}

## Real exp 

rownames(final_matrix) <- final_matrix$Target

final_matrix <- final_matrix[, -1]

metadata <- metadata[match(colnames(final_matrix), metadata$Sample), ]

final_matrix <- as.matrix(final_matrix)

CellTypes <- factor(metadata$CellType)

design <- model.matrix(~ 0 + CellTypes)

colnames(design) <- levels(CellTypes)

fit <- lmFit(final_matrix, design)

types <- colnames(design)

rank_list_real <- list()

for (n in types) {
  others <- setdiff(types, n)
  
  contrast_expr <- paste0(n, " - (", paste(others, collapse = " + "), ")/", length(others))
  contrast_matrix <- makeContrasts(contr = contrast_expr, levels = design)
  
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  tab <- topTable(fit2, coef = 1, number = Inf, sort.by = "t")
  
  rank_list_real[[n]] <- tab
}

convert_results_list_to_df <- function(results_list, net_type_label) {
  lapply(names(results_list), function(CellType) {
    df <- results_list[[CellType]]
    tibble(
      Target = rownames(df),
      Celltype = CellType,
      t_stat = df$t,
      FDR = df$adj.P.Val
    )
  }) |> bind_rows() |> mutate(Exp_type = net_type_label)
}

df_real <- convert_results_list_to_df(rank_list_real, "Real")

df_pred <- convert_results_list_to_df(rank_list, "Predicted")

results_res <- bind_rows(df_real, df_pred)

results_res <- results_res %>% rename(residue = t_stat)

## Object creation 

for (n in c("T_cell", "Fibroblast", "Basal", "Luminal", "Endothelial", "Macrophage")) {

assign(paste0("fgsea_", n, "_exp"), run_fgsea(celltype = n, results_res, msigdb_BP, msigdb_MF, Inf_type = "Exp_type"))
} 


for (n in c("T_cell", "Fibroblast", "Basal", "Luminal", "Endothelial", "Macrophage")) {
plot_data <- fgsea_Fibroblast %>%
  filter(Exp_type == "Predicted", Ontology == "BP", padj < 0.05) %>%
  mutate(logFDR = -log10(padj)) %>% 
  arrange(desc(logFDR)) %>%
  mutate(NAME = factor(pathway, levels = pathway))

top_pos <- plot_data %>% filter(NES > 0) %>% arrange(desc(NES)) %>% slice_head(n = 30) 
top_neg <- plot_data %>% filter(NES < 0) %>% arrange(NES) %>% slice_head(n = 30)

plot_data_filt <- bind_rows(top_pos, top_neg) %>% arrange(desc(logFDR)) %>% 
mutate(NAME = factor(NAME, levels = rev(NAME)))

ggplot(plot_data_filt, aes(x = logFDR, y = NAME, size = size, fill = NES)) +
  geom_point(shape = 21, color = "grey20") +
  scale_size_area(max_size = 10) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  labs(x = "−log10(FDR)", y = NULL, title = "GSEA - Predicted network (GO:BP) - Fibroblast") +
  theme(axis.text.y = element_text(size = 10))
}

## Enrichment analysis for expression - comparison between Real and Predicted 

plot_data <- list()
plot_com_num_exp <- list()

for (n in c("T_cell", "Fibroblast", "Basal", "Luminal", "Endothelial", "Macrophage")) {

fgsea_net <- get(paste0("fgsea_", n, "_exp")) %>% mutate(Sign = ifelse(NES > 0, "Positive", "Negative"))

plot_real <- fgsea_net %>%
  filter(Exp_type == "Real", Ontology == "BP", padj < 0.001) %>% filter(NES > 1 | NES < -1) %>%
  mutate(logFDR = -log10(padj)) %>% group_by(Sign) %>%
  mutate(NAME = factor(pathway, levels = pathway)) %>% arrange(desc(abs(NES))) %>% ungroup()

plot_pred <- fgsea_net %>%
  filter(Exp_type == "Predicted", Ontology == "BP", padj < 0.001) %>% filter(NES > 1 | NES < -1) %>% group_by(Sign) %>%
  mutate(logFDR = -log10(padj)) %>% arrange(desc(abs(NES))) %>% ungroup()
 
com_path <- intersect(plot_real$pathway, plot_pred$pathway)

plot_com_num <- data.frame(path_pred = length(plot_pred$pathway), path_real = length(plot_real$pathway), path_com = length(com_path), cell_type = n)

plot_com_num_exp[[n]] <- plot_com_num


plot_com <- bind_rows(plot_real, plot_pred) %>% 
  filter(pathway %in% com_path) %>% 
  mutate(Exp_type = factor(Exp_type, levels = c("Real", "Predicted"))) %>% 
  group_by(Exp_type, Sign) %>% slice_head(n = 20) %>%
  arrange(logFDR) %>% 
   mutate(
          pathway = gsub("_", " ", as.character(pathway)),
          pathway = str_wrap(pathway, width = 50, whitespace_only = FALSE), 
          pathway = gsub(" ", "_", pathway)
            ) %>% 
  mutate(NAME = factor(pathway, levels = unique(pathway))) %>% 
  ungroup()

plot_data[[n]] <- plot_com

p9 <- ggplot(plot_com, aes(x = logFDR, y = NAME, size = size, fill = NES)) +
  geom_point(shape = 21, color = "grey20") +
  scale_size_area(max_size = 10) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  facet_wrap(Sign ~ Exp_type, scales = "free_y") +
  labs(x = "−log10(FDR)", y = NULL, title = paste0("Common pathways in ", n)) +
  theme_minimal() + theme(
    axis.text.y = element_text(size = 12),        
    axis.text.x = element_text(size = 12),         
    axis.title.x = element_text(size = 12),        
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold"),               
    legend.text = element_text(size = 12),         
    legend.title = element_text(size = 12)        
  )

ggsave(paste0("GSEA_common_", n, "_exp.pdf"), width = 32, height = 16, plot = p9)
}

plot_com_num_exp <- bind_rows(plot_com_num_exp) 

plot_com_num_exp <- plot_com_num_exp %>%
  pivot_longer(cols = c(path_real, path_pred, path_com),
               names_to = "Group", values_to = "Count") %>%
  mutate(Group = recode(Group,
                        path_real = "Real",
                        path_pred = "Predicted",
                        path_com  = "Common"))

p6 <- ggplot(plot_com_num_exp, aes(x = cell_type, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(
    title = "Number of Significant Pathways per Cell Type",
    x = "Cell Type",
    y = "Number of Pathways",
    fill = "Source"
  ) +
  theme_minimal() +
     theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
  axis.text.y = element_text(size = 13),
  axis.title = element_text(size = 14),
  plot.title = element_text(size = 16, face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 12)
)

ggsave("GSEA_common_exp.pdf", width = 12, height = 8, plot = p6)

## Correlation difference analysis between expression and indegrees 

metrics_by_patient_cell_net$type <- "filter" 
metrics_by_patient_cell_net_bef$type <- "no_filter"
real_propor <- read.table("GT_Data/real_propor.txt", 
                 header = TRUE, 
                 sep = "",         
                 quote = "\"",     
                 stringsAsFactors = FALSE,
                 check.names = FALSE,            
                ) 

predict_propor <- read.table("Results/Proportions/predict_prop.txt", 
                 header = TRUE, 
                 sep = "",         
                 quote = "\"",    
                 stringsAsFactors = FALSE,
                 check.names = FALSE,
                ) 


real_propor <- real_propor %>% pivot_longer(cols = -c(sample, type), names_to = "CellType", values_to = "real_cell_pror") %>% 
mutate(sample = gsub("_", "", sample))
predict_propor <- predict_propor %>% pivot_longer(cols = -c(sample, type), names_to = "CellType", values_to = "pred_cell_pror") %>% 
mutate(sample = gsub("_", "", sample))

colnames(real_propor)[1] <- "Patient"
colnames(predict_propor)[2] <- "Patient"

combined_data <- inner_join(predict_propor, metrics_by_patient_cell, by = c("Patient", "CellType"))
combined_data_fin <- inner_join(combined_data, metrics_by_patient_cell_net, by = c("Patient", "CellType")) 
combined_data_fin <- combined_data %>% mutate(diff_prop = real_cell_pror - pred_cell_pror) %>% mutate(diff_prop_perc = abs(diff_prop) / real_cell_pror) %>%
  dplyr::select(Patient, CellType, cor_spear, cor_pear, pred_cell_pror, real_cell_pror, diff_prop, diff_prop_perc)
combined_data_fin_sel <- inner_join(combined_data_fin, metrics_by_patient_cell, by = c("Patient", "CellType")) 
combined_data_fin_sel <- combined_data_fin_sel %>% mutate(diff_spear = cor_spear.y - cor_spear.x) %>%
  dplyr::select(Patient, CellType, cor_spear.x, cor_spear.y, pred_cell_pror, real_cell_pror, diff_prop, diff_spear, diff_prop_perc) 

## test bad samples bigger difference on proportions 

bad_samples <- combined_data_fin %>% arrange(desc(diff_prop_perc)) %>% slice_head(n = 30) %>% dplyr::select(diff_prop_perc)
combined_data_fin$highlight <- combined_data_fin$diff_prop_perc %in% bad_samples$diff_prop_perc

p3 <- ggplot(combined_data, aes(x = pred_cell_pror, y = cor_pear, color = Patient, shape = CellType)) +
geom_point(size = 2) + 
geom_hline(yintercept = 0.75, color = "red", linetype = "dashed") +
geom_vline(xintercept = 0.02, color = "red", linetype = "dotted") +
scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7)) +
  labs(title = "Regression: Pearson correlation vs predicted proportions",
       x = "Predicted proportions",
       y = "Pearson correlation")


p3 <- ggsave("pear_correlation_vs_prop_ind_filt.pdf", width = 12, height = 8)

## Filter analysis 

length(which(combined_data_fin$pred_cell_pror < 0.02)) 
length(which(combined_data_fin$cor_spear < 0.8 & combined_data_fin$pred_cell_pror > 0.02))
length(which(combined_data_fin$cor_spear > 0.8 & combined_data_fin$pred_cell_pror < 0.02))

filt_samp <- combined_data_fin %>% filter(pred_cell_pror < 0.02) %>% dplyr::select(Patient, CellType) %>% mutate(sample = paste0(CellType, ".", Patient)) %>%
dplyr::select(sample)

## Correlations 

cor.test(combined_data_fin$cor_spear, combined_data_fin$real_cell_pror, method = "pearson")

modelo_base <- lm(cor_pear ~ real_cell_pror + diff_prop_perc, data = combined_data_fin)
summary(modelo_base)
combined_data_fin$resid_cor_vs_prop <- resid(modelo_base)

outliers <- combined_data_fin %>%
  filter(abs(resid_cor_vs_prop) > 0.12) %>% dplyr::select(Patient, CellType, cor_pear, real_cell_pror, resid_cor_vs_prop, diff_prop)

## Correlation difference analysis between expression and indegrees

dif <- inner_join(combined_data_fin, correl_exp_ind, by = c("Patient", "CellType"))
dif <- dif[1:16,]
dif <- dif %>% arrange(cor_pear)
dif <- dif %>%
mutate(cor_pear_exp = cor_pear, 
      cor_spear_exp = cor_spear,
      ) %>%
dplyr::select(Patient, CellType, cor_pear_exp, pred_cell_pror, real_cell_pror, cor_spear_exp) %>% arrange(desc(cor_pear_exp)) %>%
  mutate(
    id = paste0(Patient, "_", CellType),
    id = factor(id, levels = unique(id))  
  ) %>%
  pivot_longer(cols = c(cor_pear_exp, cor_spear_exp),
               names_to = "cor_type", values_to = "cor") %>%
  mutate(cor_type = factor(cor_type, levels = c("cor_pear_exp", "cor_spear_exp")))


 p3 <- ggplot(dif, aes(x = 1, y = id, fill = cor)) +
  geom_tile(color = "white", size = 0.5, width = 0.8, height = 0.8) +  
  geom_text(aes(label = round(cor, 2)), color = "black", size = 5) +
  scale_fill_gradient2(low = "white", mid = "lightpink", high = "red", midpoint = 0.5) +
  facet_wrap(~ cor_type, ncol = 2, strip.position = "top") +  
  theme_minimal(base_size = 14) +
  labs(
    title = "Correlation Comparison between Real and Predicted Expression",
    x = "Correlation Value",
    y = "Patient - Cell Type",
    fill = "Correlation"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.y = unit(0.5, "lines"),
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text.y = element_blank(),  
    axis.text.y = element_text()    
  )
ggsave("correlation_comparison_expp.pdf", width = 12, height = 8, plot = p3)

dif <- combined_data_fin
dif <- dif[1:16,]
dif <- dif %>% arrange (desc(cor_pear)) %>%
  mutate(
    id = paste0(Patient, "_", CellType),
    id = factor(id, levels = unique(id)), 
  ) %>% pivot_longer(cols = c(cor_pear, cor_spear),
               names_to = "cor_type", values_to = "cor") %>%
  mutate(cor_type = factor(cor_type, levels = c("cor_pear", "cor_spear")))


  p3 <- ggplot(dif, aes(x = 1, y = id, fill = cor)) +
  geom_tile(color = "white", size = 0.5, width = 0.8, height = 0.8) +  
  geom_text(aes(label = round(cor, 2)), color = "black", size = 5) +
  facet_wrap(~ cor_type, ncol = 2, strip.position = "top") +
  scale_fill_gradient2(low = "white", mid = "lightpink", high = "red", midpoint = 0.5) + 
  theme_minimal(base_size = 14) +
  labs(
    title = "Correlation between Real and Predicted Expression",
    x = "Correlation Value",
    y = "Patient - Cell Type",
    fill = "Correlation"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.y = unit(0.5, "lines"),
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text.y = element_blank(),  
    axis.text.y = element_text()    
  )

ggsave("correlation_exp1.pdf", width = 28, height = 16, plot = p3)

## Limma analysis 
## Preanalysis 

rownames(indegrees_pred) <- indegrees_pred$Target
indegrees_pred <- indegrees_pred[, -length(colnames(indegrees_pred))]
indegrees_pred <- as.matrix(indegrees_pred)

## Check distribution
hist(as.vector(indegrees_pred), breaks = 200, main = "Distribución general de in-degree")

## Continued data 
unique_vals <- unique(as.vector(indegrees_pred))
length(unique_vals)

metadata <- data.frame(Sample = setdiff(colnames(pred_indegrees), "Target")) %>%
  mutate(CellType = str_split_fixed(Sample, "\\.", 2)[,1])  %>% 
  mutate(Patient = str_split_fixed(Sample, "\\.", 2)[,2])

indegrees_pred <- indegrees_pred[, metadata$Sample]

## Homocedasticity test
set.seed(42)
sample_genes <- sample(rownames(indegrees_pred), 2000)

var_data <- lapply(sample_genes, function(g) {
  values <- as.numeric(indegrees_pred[g, ])
  data.frame(
    gene = g,
    value = values,
    CellType = metadata$CellType
  ) |>
    group_by(CellType) |>
    summarise(var = var(value), .groups = "drop")
})

var_df <- bind_rows(var_data)

boxplot(var ~ CellType, data = var_df,
        main = "Varianza del in-degree por tipo celular (500 genes)",
        ylab = "Varianza")

## Residues per cell type, comparative analysis between and inside groups 

results <- lapply(rownames(indegrees_pred), function(gene) {
  values <- as.numeric(indegrees_pred[gene, ])
  df <- data.frame(value = values, CellType = metadata$CellType)
  model <- lm(value ~ 0 + CellType, data = df)

  list(
    gene = gene,
    coefs = coef(model),             
    residuals = residuals(model),     
    model = model                    
  )
})

names(results) <- rownames(indegrees_pred)

genes_to_plot <- c("A1BG", "A1BG-AS1", "A2M") 

res_plot_data <- lapply(genes_to_plot, function(g) {
  data.frame(
    Gene = g,
    Residual = results[[g]]$residuals,
    CellType = metadata$CellType
  )
}) 

res_plot_data <- bind_rows(res_plot_data) 

ggplot(res_plot_data, aes(x = CellType, y = Residual)) +
  geom_boxplot(outlier.size = 0.5, fill = "grey90") +
  facet_wrap(~ Gene, scales = "free_y") +
  theme_minimal() +
  labs(title = "Residuos por tipo celular para genes seleccionados",
       y = "Residuo", x = "Tipo celular")


metrics_tot <- read.table("metrics_by_patient_cell_net_total.txt", header = TRUE)

## Significant differences between networks with and without filtering 

p7 <- ggplot(metrics_by_patient_merge, aes(x = filter, y = cor_pear, fill = filter)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Correlación real vs filter",
       y = "Correlación de Pearson",
       x = "") 

ggsave("dif_filters_total.pdf", width = 12, height = 8, plot = p7)

write.csv2(metrics_by_patient_merge, "metrics_by_patient_cell_net_comp.csv", row.names = FALSE)

