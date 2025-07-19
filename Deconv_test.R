
## Pseudobulks and signature matrix without Mastocytes and B cells

matrix_dec_test <- matrix[,!colnames(matrix) %in% signature]
matrix_dec_test <- matrix_dec_test[,!colnames(matrix_dec_test) %in% matrix@meta.data$cell_name[matrix@meta.data$cell_types %in% c("Mastocyte", "B_cell")]]
matrix_cont_dec_test <- list()

for (n in 1:13) {
    matrix_cont_dec_test[[n]] <- matrix_dec_test@assays$RNA$counts[, 
        matrix_dec_test@meta.data$cell_name[matrix_dec_test@meta.data$sample == n]] %>%
        as.data.frame() %>%
        mutate(!!paste0("Bulk", n) := rowSums(.)) %>%
        dplyr::select(!!paste0("Bulk", n))
}

pseudobulks <- do.call(cbind, matrix_cont_dec_test)
pseudobulks$Gene <- rownames(pseudobulks)
pseudobulks <- pseudobulks[, c("Gene", colnames(pseudobulks)[-ncol(pseudobulks)])]
rownames(pseudobulks) <- pseudobulks$Gene 
pseudobulks <- pseudobulks[,-1]
pseudobulks <- t(pseudobulks)
write.table(pseudobulks, file = "pseudobulks.txt", row.names = TRUE)

## SM 
signature_matrix <- matrix[,(colnames(matrix) %in% signature)]
signature_matrix <- signature_matrix[,!colnames(signature_matrix) %in% signature_matrix@meta.data$cell_name[signature_matrix@meta.data$cell_types %in% c("Mastocyte", "B_cell")]]
cell.type.labels <- signature_matrix@meta.data$cell_types
write.table(cell.type.labels, file = "cell.type.labels.txt")
sc.dat <- signature_matrix@assays$RNA$counts 
sc.dat <- t(sc.dat) 
sc.dat <- as.matrix(sc.dat) 
write.table(sc.dat, file = "sc.dat.txt", row.names = TRUE)

## Real propor 
table_cheng <- matrix_dec_test@meta.data %>% group_by(sample, cell_types) %>% summarise(num_cells = n(), .groups = 'drop_last') %>%
mutate(percent = (num_cells / sum(num_cells))) %>%
ungroup()
cells <- c("Luminal", "Basal", "Fibroblast", "Macrophage", "T_cell", "Endothelial")
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
write.table(real_propor, file = "real_propor.txt")

## Gene expresion prediction 
expr_tot <- read.table("expr_tot.txt")
