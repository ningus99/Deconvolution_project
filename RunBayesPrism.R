
# Deconvolution with BayesPrism

## Signature matrix and cell type label creation  

## Random selection 

set.seed(46)

## Shenling dataset 

set.seed(43)
set.seed(20)

signature <- colnames(matrix)

write.table(signature, file = "BayesPrism/signature_exp.txt")

## Signature matrix and cell type label 

signature_matrix <- matrix[, colnames(matrix) %in% signature]
cell.type.labels <- signature_matrix@meta.data$cells
write.table(cell.type.labels, file = "BayesPrism/cell.type.labels_exp.txt")
sc.dat <- signature_matrix@assays$RNA$counts 
sc.dat <- t(sc.dat) 
sc.dat <- as.matrix(sc.dat) 
write.table(sc.dat, file = "BayesPrism/sc.dat_exp.txt", row.names = TRUE)

## Pseudobulk construction 

matrix_dec <- matrix[,!(colnames(matrix) %in% signature)]
matrix_cont_dec <- list()

for (n in 1:13){
     matrix_cont_dec[[n]] <- matrix_dec@assays$RNA$counts[, matrix_dec@meta.data$cell_name[matrix_dec@meta.data$sample == n]] %>% as.data.frame() %>% 
     mutate(!!paste0("Bulk", n) := rowSums(.)) %>% dplyr::select(!!paste0("Bulk", n))
}
pseudobulks <- do.call(cbind, matrix_cont_dec)
pseudobulks$Gene <- rownames(pseudobulks)
pseudobulks <- pseudobulks[, c("Gene", colnames(pseudobulks)[-ncol(pseudobulks)])]

rownames(pseudobulks) <- pseudobulks$Gene 
pseudobulks <- pseudobulks[,-1]
pseudobulks <- t(pseudobulks)
write.table(pseudobulks, file = "BayesPrism/Data/pseudobulks.txt", row.names = TRUE)
pseudobulks <- read.table("BayesPrism/pseudobulks.txt", check.names = FALSE)
cell.type.labels <- cell.type.labels$x

cell.type.labels <- read.table("BayesPrism/cell.type.labels_exp.txt")
sc.dat <- read.table("BayesPrism/sc.dat_exp.txt", check.names = FALSE, row.names = 1)

## Run BayesPrism 

myPrism <- new.prism(
  reference=sc.dat, 
  mixture=pseudobulks,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.type.labels,
  outlier.cut=0.01,
  outlier.fraction=0.1,
  key = NULL
)

## Create BayesPrism object 

bp_res <- run.prism(prism = myPrism, n.cores=50)



predict_prop <- theta %>% mutate(type = "Predicted")
predict_prop$sample <- c(paste0("patient_", 1:13))
theta <- as.data.frame(theta)

## Cell proportion extraction 

theta <- get.fraction (bp=bp_res,
                       which.theta="final",
                       state.or.type="type")

write.table(predict_prop, file = "predicted_prop_exp.txt", row.names = FALSE)

## Predicted gene expresison extraction 

for (n in names(table(cell.type.labels))) {
  assign(paste0(n, "_exp"),
         get.exp(bp = bp_res,
                 state.or.type = "type",
                 cell.name = n))
}


## Check identity
## Cell proportion 

identical_df_prop <- data.frame(
  seed = integer(),
  identical = logical(),
  stringsAsFactors = FALSE
)

for(i in seeds) {
  result <- identical(theta[[i]], theta[[1]])
  
  identical_df_prop <- rbind(identical_df_prop, data.frame(
    seed = i,
    identical = result
  ))
}

## Expression

identical_df <- data.frame(
  cell_type = character(),
  seed = integer(),
  identical = logical(),
  stringsAsFactors = FALSE
)

for(n in names(table(cell.type.labels))) {
    result <- identical(Pred_exp[[i]][[n]], Pred_exp[[1]][[n]])
    
    identical_df <- rbind(identical_df, data.frame(
      cell_type = n,
      seed = i,
      identical = result
    ))
  }
}

write.table(identical_df, file = "BayesPrism/Results/identical_df.txt", row.names = FALSE)
write.table(identical_df_prop, file = "BayesPrism/Results/identical_df_prop.txt", row.names = FALSE)
