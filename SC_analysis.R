#### SC analysis #### 

## CHENG datasets

counts <- readMM("/storage/kuijjerarea/dcosorioh/CancerHallmarks/Data_Chen2021_Prostate/Exp_data_UMIcounts.mtx")
genes <- read.table("/storage/kuijjerarea/dcosorioh/CancerHallmarks/Data_Chen2021_Prostate/Genes.txt")
cells <- read.csv("/storage/kuijjerarea/dcosorioh/CancerHallmarks/Data_Chen2021_Prostate/Cells.csv")
metadata <- read.csv("/storage/kuijjerarea/dcosorioh/CancerHallmarks/Data_Chen2021_Prostate/Meta-data.csv")
colnames(counts) <- cells$cell_name
rownames(counts) <- genes$V1
meta_data <- merge(cells, metadata, by = "sample")
rownames(meta_data) <- meta_data$cell_name
matrix <- CreateSeuratObject(counts = counts, meta.data = meta_data, min.cells = 3, min.features = 200) 

## Shenling dataset 

metadata <- read.csv(gzfile("SC_data/GSE181294_scRNAseq.ano (1).csv.gz"))

files <- list.files("SC_data/", pattern = "count.csv.gz")

combined_data <- read.csv(gzfile(paste0("SC_data/",files[1])))


for (f in files[-1]){

df <- read.csv(gzfile(paste0("SC_data/",(f))))

combined_data <- merge(combined_data, df, by = "X", all = TRUE)
}

rownames(combined_data) <- combined_data$X

combined_data <- combined_data[,-1]

## Transfomation in the same format 

metadata$X <- gsub("-", ".", metadata$X)

## Common genes 

metadata <- metadata[metadata$X %in% colnames(combined_data),]
combined_data <- combined_data[, colnames(combined_data) %in% metadata$X]
rownames(metadata) <- metadata$X
metadata <- metadata[,-1]

## Correction of filters 

metadata$cells <- gsub("Macrophage[1-4]", "Macrophage", metadata$cells)
metadata$cells <- gsub("Epitheial_Luminal|Tumor", "Luminal", metadata$cells)
metadata$cells <- gsub("Epitheial_Basal", "Basal", metadata$cells)
metadata$cells <- gsub("Endothelial cells-1|Endothelial cells-2", "Endothelial", metadata$cells)
metadata$cells <- gsub("Fibroblasts", "Fibroblast", metadata$cells)
metadata$cells <- gsub("B cells", "B_cell", metadata$cells)
metadata$cells <- gsub("Mast cells", "Mastocyte", metadata$cells)

# Agrupar todos los subtipos T en una sola categoría

metadata$cells <- gsub("CTL-1|CTL-2|Cycling T|Naive Th|Th1|Th17|Treg", "T_cell", metadata$cells)

## Seurat object creation 

matrix <- CreateSeuratObject(counts = combined_data, meta.data = metadata, min.cells = 3, min.features = 200)

## Filter cell types of the analysis 

matrix <- subset(matrix, subset = cells %in% c("Basal", "Luminal", "T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Mastocyte"))
matrix <- matrix[rownames(matrix) %in% ey]

## Create Seurat object 

matrix <- CreateSeuratObject(counts = counts, meta.data = meta_data, min.cells = 3, min.features = 200) 

## Quality filters

matrix[['percent_mt']] <- PercentageFeatureSet(matrix, pattern = '^MT-')
qualCheng <- VlnPlot(matrix, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"), ncol = 3)
FeatureScatter(matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
ggsave("results/qualCheng.png", plot = qualCheng, width = 8, height = 6)

## Filters for doublet searching

matrix <- subset(matrix, subset = nFeature_RNA > 250 & nCount_RNA > 500)

## DONG dataset 

counts1 <- readMM("/storage/kuijjerarea/dcosorioh/CancerHallmarks/Data_Dong2020_Prostate/Exp_data_UMIcounts.mtx")
genes1 <- read.table("/storage/kuijjerarea/dcosorioh/CancerHallmarks/Data_Dong2020_Prostate/Genes.txt")
cells1 <- read.csv("/storage/kuijjerarea/dcosorioh/CancerHallmarks/Data_Dong2020_Prostate/Cells.csv")
metadata1 <- read.csv("//storage/kuijjerarea/dcosorioh/CancerHallmarks/Data_Dong2020_Prostate/Meta-data.csv")

colnames(counts1) <- cells1$cell_name
rownames(counts1) <- genes1$V1
meta_data1 <- merge(cells1, metadata1, by = "sample")
rownames(meta_data1) <- meta_data1$cell_name
matrix1 <- CreateSeuratObject(counts = counts1, meta.data = meta_data1, min.cells = 3, min.features = 200)

## Filters for doublet search

matrix1[['percent_mt']] <- PercentageFeatureSet(matrix1, pattern = '^MT-')
qualDong <- VlnPlot(matrix1, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"), ncol = 3)
FeatureScatter(matrix1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
ggsave("results/qualDong.png", plot = qualDong, width = 8, height = 6)
matrix1 <- subset(matrix1, subset = nFeature_RNA > 250 & nCount_RNA > 500 & percent_mt < 10)

## Doublet search 
## CHENG 

sce <- as.SingleCellExperiment(matrix)
sce <- scDblFinder(sce, samples = "sample")
table(sce$scDblFinder.class, sce$sample)
sce@colData@listData %>% as.data.frame() %>% head()
meta_scdblfinder <- sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) 
  rownames(meta_scdblfinder) <- sce@colData@rownames
matrix <- AddMetaData(object = matrix, metadata = meta_scdblfinder %>% dplyr::select('scDblFinder.class'))
doublet_qual <- VlnPlot(matrix, group.by = 'sample', split.by = "scDblFinder.class",
        features = c("nFeature_RNA", "nCount_RNA"), 
        ncol = 2, pt.size = 0) + theme(legend.position = 'right')
ggsave("results/doublet_qual_Cheng.png", plot = doublet_qual, width = 30, height = 8)

## DONG 

sce1 <- as.SingleCellExperiment(matrix1)
sce1 <- scDblFinder(sce1, samples = "sample") 
table(sce1$scDblFinder.class, sce1$sample)
sce1@colData@listData %>% as.data.frame() %>% head()
meta_scdblfinder1 <- sce1@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) 
  rownames(meta_scdblfinder1) <- sce1@colData@rownames
matrix1 <- AddMetaData(object = matrix1, metadata = meta_scdblfinder1 %>% dplyr::select('scDblFinder.class'))
doublet_qual_Dong <-VlnPlot(matrix1, group.by = 'sample', split.by = "scDblFinder.class",
        features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
ggsave("results/doublet_qual_Dong.png", plot = doublet_qual_Dong, width = 30, height = 8)

## Remove doublets 

matrix <- subset(matrix, subset = scDblFinder.class == "singlet")
matrix1 <- subset(matrix1, subset = scDblFinder.class == "singlet")

## NORMALIZE 

matrix <- NormalizeData(matrix) 
matrix1 <- NormalizeData(matrix1)

## VARIABLE GENES

matrix <- FindVariableFeatures(matrix, selection.method = "vst", nfeatures = 2000)
all_genes <- rownames(matrix) 

matrix1 <- FindVariableFeatures(matrix1, selection.method = "vst", nfeatures = 2000)
all_genes1 <- rownames(matrix1)

## SCALE WITH CELL SCORE 

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

matrix <- CellCycleScoring(matrix, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
matrix <- ScaleData(matrix, vars.to.regress = c("S.Score", "G2M.Score"), features = all_genes)
matrix1 <- CellCycleScoring(matrix1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
matrix1 <- ScaleData(matrix1, vars.to.regress = c("S.Score", "G2M.Score"), features = all_genes1)
matrix <- ScaleData(matrix, features = all_genes)

## PCA DONG 

matrix1 <- RunPCA(matrix1, features = VariableFeatures(matrix1))
print(matrix1[['pca']], dims = 1:5, nfeatures = 5)
DimHeatmap(matrix1, dims = 1, cells = 500, balanced = TRUE)
DimPlot(matrix1, reduction = "pca") + NoLegend()

pca_var <- matrix1[["pca"]]@stdev^2
explained_variance <- pca_var / sum(pca_var) * 100
PC_15 <- sum(explained_variance[1:27])

## PCA CHENG

matrix <- RunPCA(matrix, features = VariableFeatures(matrix))
print(matrix[['pca']], dims = 1:5, nfeatures = 5)
DimHeatmap(matrix, dims = 1, cells = 500, balanced = TRUE)
DimPlot(matrix, reduction = "pca") + NoLegend()

pca_var <- matrix[["pca"]]@stdev^2
explained_variance <- pca_var / sum(pca_var) * 100
PC_15 <- sum(explained_variance[1:30])

# Determine number of components 

if(interactive()){
  ElbowPlot(matrix, ndims = 70)
  ElbowPlot(matrix1)
}

## Batch correction 
## CHENG 

matrix <- RunHarmony(
    object = matrix,
    group.by.vars = "patient",                
    dims.use = 1:30                   
)

## DONG 

matrix1 <- RunHarmony(
    object = matrix1,
    group.by.vars = "patient",                
    dims.use = 1:27             
)

## Clustering 
## CHENG 

matrix <- FindNeighbors(matrix, reduction = "pca", dims = 1:30)
matrix <- FindClusters(matrix, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
DimPlot(matrix, group.by = 'RNA_snn_res.0.3', label = TRUE)
Idents(matrix) <- 'RNA_snn_res.0.1'

# DONG 

matrix1 <- FindNeighbors(matrix1, reduction = "harmony", dims = 1:27)
matrix1 <- FindClusters(matrix1, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
DimPlot(matrix1, group.by = 'RNA_snn_res.0.5', label = TRUE)
Idents(matrix1) <- 'RNA_snn_res.0.1'

## UMAP 
## CHENG 

matrix <- RunUMAP(matrix, reduction = "pca", dims = 1:30)
colors <- c("blue", "red", "grey", "yellow", "orange", "pink", "brown", "purple", "black", "limegreen", "darkblue", "gold", "darkred", "skyblue", "coral")
DimPlot(matrix, reduction = 'umap', group.by = "cells")
ggsave("umap_Cheng_shen.pdf", plot = last_plot(), width = 10, height = 6)
## DONG 

matrix1 <- RunUMAP(matrix1, reduction = "harmony", dims = 1:27)
colors1 <- c("blue", "red", "grey", "yellow", "orange", "pink", "brown", "purple", "black", "limegreen", "darkblue", "gold", "darkred")
DimPlot(matrix1, reduction = 'umap', cols = colors, group.by = "cell_types")

## Visualization individual variability

met <- rownames(matrix@meta.data[matrix@meta.data$cell_type == "T_cell" & matrix@meta.data$sample == "1",])
met1 <- rownames(matrix1@meta.data[matrix1@meta.data$cell_types == "Fibroblast" & matrix1@meta.data$sample == "patient #2",])
FeaturePlot(matrix1, features = c("ACTA2"))
DimPlot(matrix1, reduction = 'umap', cells.highlight = met1, cols.highlight = "red", label = FALSE)
DimPlot(matrix, reduction = 'umap', cols = colors1, cells.highlight = met, cols.highlight = "red", label = FALSE)
DimPlot(matrix1, reduction = 'umap', cols = colors1, group.by = "patient", label = FALSE)

## Cluster by cluster finding markers 
## CHENG  

clusters <- levels(Idents(matrix))
results <- list()
for (i in 1:length(clusters)) {
    for (j in 1:length(clusters)) {
        if (i != j) { 
            cluster1 <- clusters[i]
            cluster2 <- clusters[j]
            
            
            markers <- FindMarkers(
                matrix,
                ident.1 = cluster1,
                ident.2 = cluster2,
                only.pos = FALSE, 
                min.pct = 0, 
                logfc.threshold = 0
            )
            
           
            markers$comparison <- paste(cluster1, "vs", cluster2, sep = "_")
            markers$gene <- rownames(markers)
            
           
            results[[paste(cluster1, "vs", cluster2)]] <- markers
        }
    }
 
}
 
 results_fin <- list()
for (n in 0:14){
  cluster <- grep(paste0("^", n, "\\b vs"), names(results), value = TRUE)
  cluster_dfs <- lapply(results[cluster], function(df) {
   for (s in colnames(df)) {
    if (s != "gene") {
  colnames(df)[colnames(df) == s] <- paste0(s, "_cluster", n)
    }
  }
  return(df)
  })
results_fin[[n + 1]] <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), cluster_dfs)}
resultsfin <- list() 
for (n in 0:length(clusters)){
  resultsfin[[n]] <- results %>% {[names(results)[grepl(paste0(n, "^"), names())]]} %>% list() 
  resultsfin[[n]] <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), results[n])
  
}

## Finding markers average of counts 
## CHENG

marcadores_Cheng <- FindAllMarkers(matrix,
               logfc.threshold = 0.25, 
               min.pct = 0.25, 
               only.pos = FALSE,
               test.use = 'wilcox',   
               )

markersfilt <- subset(marcadores_Cheng, (avg_log2FC < -1.5 | avg_log2FC > 1.5) & p_val_adj < 0.05 & pct.1 > 0.2)

markersfilt <- markersfilt %>% group_by(cluster) 

# DONG 

marcadores_Dong <- FindAllMarkers(matrix1,
               logfc.threshold = 0.25, 
               min.pct = 0.25, 
               only.pos = FALSE,
               test.use = 'wilcox',   
               )

markersfilt1 <- subset(marcadores_Dong, (avg_log2FC < -1.5 | avg_log2FC > 1.5) & p_val_adj < 0.05 & pct.1 > 0.2)
markersfilt1 <- markersfilt1 %>% group_by(cluster) 

## Marker catalog  

genes <- c("TP63", "PROM1", "CDH1", "CD44", "S100A6", "KLK3")
cell_type <- c("Basal cell", "Cancer stem cell", "Epithelial cell", "Cancer stem cell", "Basal cell", "Luminal")
marcado <- data.frame(genes, cell_type)
inm_cells <- data.frame(
    genes = c(
        "CD3E", "CD4", "CD8A", "FOXP3", "PDCD1", "CTLA4", "CD19", "CD20", "CD27",
        "CD68", "NOS2", "IL1B", "ARG1", "CD163", "IL10", "NCAM1", "KLRD1", 
        "NKG7", "GZMB", "S100A8", "S100A9", "CD33", "KIT", "TPSAB1"
    ),
    cell_type = c(
        "Lymphocyte T", "Lymphocyte T Helper", "Lymphocyte T Cytotoxic", "Regulatory T Cell", 
        "Exhausted T Cell", "Regulatory T Cell", "Lymphocyte B", "Lymphocyte B Mature", 
        "Lymphocyte B Activated", "Macrophage General", "Macrophage M1", "Macrophage M1", 
        "Macrophage M2", "Macrophage M2", "Macrophage M2", "Natural Killer Cell", 
        "Natural Killer Cell", "Natural Killer Cell", "Natural Killer Cell", 
        "Myeloid Suppressor Cell", "Myeloid Suppressor Cell", "Myeloid Suppressor Cell", 
        "Mast Cell", "Mast Cell"
    ),
    stringsAsFactors = FALSE
)
fibroblast_markers <- data.frame(
    genes = c(
        "VIM", "S100A4", "COL1A1", "PDGFRA", "PDGFRB", 
        "ACTA2", "FN1", "THY1"
    ),
    cell_type = c(
        "Fibroblast", "Fibroblast", "Fibroblast", "Fibroblast", "Fibroblast",
        "Cancer-Associated Fibroblast", "Cancer-Associated Fibroblast", "Cancer-Associated Fibroblast"
    ),
    stringsAsFactors = FALSE
)


new_rows <- data.frame(
  genes = c(
    "PECAM1", "VWF", "ENG", "CMA1", "MS4A2", "TPSB2", 
    "AR", "KRT19", "KRT18", "KRT8", "KRT14", "KRT5", "KRT17", "KRT15", "DKK1","CYR61",
    "LYZ", "FCGR3A", "CSF1R", "UCHL1", "HAVCR2", "SELL",
    "BTLA", "IL2RA", "IL7R", "CCR7", "CD28", "SLAMF1",
    "DPP4", "CD7", "CD2", "CD3G", "CD3D", "ENO2", "ASCL1", "INSM1", "REST", "NKX3-1", "PSCA", "ALDH1A3", "LMO7", "CLDN3", "AMACR", "KLK3"
  ),
  cell_type = c(
    "Endothelial Cell", "Endothelial Cell", "Endothelial Cell", 
    "Mast Cell", "Mast Cell", "Mast Cell", 
    "Prostate Luminal Cell", "Prostate Luminal Cell", "Prostate Luminal Cell", 
    "Prostate Luminal Cell", "Basal Cell", "Basal cell",  "Basal Cell", "Basal cell", "Basal cell", "Basal cell",
    "Monocyte", "Monocyte", "Monocyte", "Neuronal Cell", 
    "Exhausted T Cell", "Naive T Cell", 
    "Regulatory T Cell", "Activated T Cell", "T Cell Memory", 
    "Naive T Cell", "T Cell Activation", "Monocyte",
    "Endothelial Cell", "T Cell", "T Cell", "T Cell", "T Cell", "Neuroendocrine", 
    "Neuroendocrine", "Neuroendocrine", "Neuroendocrine", "Luminal", "Luminal", "Luminal", "Luminal", "Luminal", "Luminal", "Luminal"
  ),
  stringsAsFactors = FALSE
)

new_rowss <- data.frame(
  genes = c(
    "MKI67", "CCND1", "CCNE1", "CDKN1A", "CDKN1B", 
    "BAX", "BCL2", "CASP3", "CASP9", "BBC3", 
    "EGFR", "AKT1", "AKT2", "MTOR", "CTNNB1", 
    "TGFB1", "MMP9", "ITGB1", 
    "SNAI1", "TP53", "HIF1A", "FOXO1", "MYC", 
    "NFKB1", "CDK4", "CCNA2", "PLK1", "AURKA", "AURKB", 
    "E2F1", "PCNA", "FOXM1", "CASP8", "BID", 
    "FAS", "DIABLO", "PIK3CA", "MAPK1", "ERBB2", 
    "SRC", "MMP2", "CDH2", "ITGA6", 
    "RB1", "PTEN", "CDKN2A", "SMAD4", "VEGFA", 
    "SOD2", "GADD45A", "ATF4"
  ),
  cell_type = c(
    "Proliferation", "Cell Cycle", "Cell Cycle", "Cell Cycle", "Cell Cycle", 
    "Apoptosis", "Apoptosis", "Apoptosis", "Apoptosis", "Apoptosis", 
    "Signaling", "Signaling", "Signaling", "Signaling", 
    "Signaling", "Metastasis", "Adhesion", "Adhesion", 
    "Metastasis", "Tumor Suppression", "Hypoxia Response", "Stress Response", "Proliferation", 
    "Inflammation", "Cell Cycle", "Cell Cycle", "Cell Cycle", "Cell Cycle", "Cell Cycle", 
    "Proliferation", "Proliferation", "Proliferation", "Apoptosis", 
    "Apoptosis", "Apoptosis", "Apoptosis", "Signaling", "Signaling", 
    "Signaling", "Signaling", "Metastasis", "Adhesion",
    "Adhesion", "Tumor Suppression", "Tumor Suppression", 
    "Tumor Suppression", "Tumor Suppression", "Hypoxia Response", 
    "Stress Response", "Stress Response", "Stress Response"
  ),
  stringsAsFactors = FALSE # Para evitar problemas con factores
)

marcado <- rbind(marcado, new_rows)
marcado <- rbind(marcado, inm_cells)
marcado <- rbind(marcado, fibroblast_markers)
marcado <- rbind(marcado, new_rowss)
marcado <- marcado %>% rename(gene = genes)

## Characterization of clusters 

markersfilt$cell_type <- ""

for (n in marcado$gene){
  if (n %in% markersfilt$gene){
cell <- marcado %>% filter(gene == n) %>% pull(cell_type) 
markersfilt$cell_type[markersfilt$gene == n] <- cell
  }
}

markersfilt1$cell_type <- ""

for (n in marcado$gene){
  if (n %in% markersfilt1$gene){
cell1 <- marcado %>% filter(gene == n) %>% pull(cell_type) 
markersfilt1$cell_type[markersfilt1$gene == n] <- cell1
  }
}

## Veryfing duplicates

markersfilt1 %>%
  group_by(gene) %>%
  filter(n() > 1) %>%
  arrange(gene)

marcado %>%
  group_by(gene) %>%
  filter(n() > 1) %>%
  arrange(gene)

## Dataset with markers 

marks_fin <- markersfilt[markersfilt$cell_type != "", ]
marks_fin1 <- markersfilt1[markersfilt1$cell_type != "", ]

marks_fin1 <- marks_fin1 %>% group_by(cluster) 
marks_fin1 <- marks_fin1 %>% group_by(cluster) %>% arrange(desc(avg_log2FC))

#GSEA

gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
gmt <- split(gene_sets$gene_symbol, gene_sets$gs_name)

## CHENG

rank_genes <- list()
for (n in 0:11){
rank_genes[[paste0("cluster", n)]] <- markersfilt %>% filter(cluster == n) %>% arrange(desc(avg_log2FC)) %>% 
select(gene, avg_log2FC) %>%  { setNames(.$avg_log2FC, .$gene) }   
}

fgsea_Cheng <- lapply(rank_genes, function(gene_rank) {
  fgsea(
    pathways = gmt,          
    stats = gene_rank,        
    minSize = 10,            
    maxSize = 500         
  )
})

# DONG 

rank_genes1 <- list()
for (n in 0:11){
rank_genes1[[paste0("cluster", n)]] <- markersfilt1 %>% filter(cluster == n) %>% arrange(desc(avg_log2FC)) %>% 
select(gene, avg_log2FC) %>%  { setNames(.$avg_log2FC, .$gene) }   
}

fgsea_Dong <- lapply(rank_genes1, function(gene_rank) {
  fgsea(
    pathways = gmt,          
    stats = gene_rank,        
    minSize = 10,            
    maxSize = 500         
  )
})

## GSEA Filters

for (n in names(rank_genes)){
fgsea_Cheng[[n]] <- fgsea_Cheng[[n]] %>% filter(padj < 0.05) 
} 

for (n in names(rank_genes1)){
fgsea_Dong[[n]] <- fgsea_Dong[[n]] %>% filter(padj < 0.05)
} 

## GSEA running 

gsea <- do.call(rbind, lapply(names(fgsea_Cheng), function(cluster){
  results <- fgsea_Cheng[[cluster]]
  results$cluster <- cluster 
  return(results) 
})) 

gsea1 <- do.call(rbind, lapply(names(fgsea_Dong), function(cluster){
  results <- fgsea_Dong[[cluster]]
  results$cluster <- cluster 
  return(results) 
})) 
gseam <- gsea[gsea$cluster %in% c("cluster1", "cluster7","cluster10", "cluster11", "cluster9"),]
gseam1 <- gsea1[gsea1$cluster %in% c("cluster2", "cluster6"),]
ggplot(gseam, aes(x = NES, y = reorder(pathway, NES), color = cluster)) +
  geom_point(size = 3) +
  facet_wrap(~ cluster, scales = "free_y") +
  labs(
    title = "Enriched_pathways",
    x = "NES (Normalized Enrichment Score)",
    y = "Pathway"
  ) +
  theme_minimal()

ggplot(gseam1, aes(x = NES, y = reorder(pathway, NES), color = cluster)) +
  geom_point(size = 3) +
  facet_wrap(~ cluster, scales = "free_y") +
  labs(
    title = "Enriched_pathways",
    x = "NES (Normalized Enrichment Score)",
    y = "Pathway"
  ) +
  theme_minimal() 
 
 ggsave("enrich_neuroend.png", plot = enrich_neuroend, width = 8, height = 6, bg ="white")

#CNV analysis 
##Cheng
##Annotation data

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
genelist_Cheng <- rownames(matrix@assays$RNA$counts) 
gene_positions_Cheng <- getBM(
    attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
    filters = "hgnc_symbol",
    values = genelist_Cheng,
    mart = ensembl 
)

colnames(gene_positions_Cheng) <- c("Gene", "Chromosome", "Start", "End")
selected_cells <- matrix@meta.data$cell_name[matrix@meta.data$RNA_snn_res.0.1 %in% c(0, 1, 3, 5, 6, 7, 9, 11, 13, 14)]
annotations_Cheng <- data.frame(Cell = selected_cells, Cluster = matrix@meta.data[selected_cells, "RNA_snn_res.0.1"])
cnv_cheng <- matrix@assays$RNA$counts[,selected_cells]

write.table(
    cnv_cheng, 
    file = "cnv_cheng.txt", 
    sep = "\t", 
    row.names = TRUE,  # Nombres de los genes
    col.names = TRUE,  # Nombres de las células
    quote = FALSE
)
write.table(
    annotations_Cheng, 
    file = "annotations_Cheng.txt", 
    sep = "\t", 
    row.names = TRUE,  # Nombres de los genes
    col.names = TRUE,  # Nombres de las células
    quote = FALSE
)
write.table(
    gene_positions_Cheng, 
    file = "gene_positions_Cheng.txt", 
    sep = "\t", 
    row.names = TRUE,  # Nombres de los genes
    col.names = TRUE,  # Nombres de las células
    quote = FALSE
)
 
## Annotations
## CHENG

matrix@meta.data$cell_types <- NA

matrix@meta.data$cell_types[matrix@meta.data$RNA_snn_res.0.1 == 2] <- "T_cell"
matrix@meta.data$cell_types[matrix@meta.data$RNA_snn_res.0.1 == 4] <- "Macrophage"
matrix@meta.data$cell_types[matrix@meta.data$RNA_snn_res.0.1 == 8] <- "B_cell"
matrix@meta.data$cell_types[matrix@meta.data$RNA_snn_res.0.1 == 10] <- "Luminal"
matrix@meta.data$cell_types[matrix@meta.data$RNA_snn_res.0.1 == 1] <- "Luminal"
matrix@meta.data$cell_types[matrix@meta.data$RNA_snn_res.0.1 == 0] <- "Luminal"
matrix@meta.data$cell_types[matrix@meta.data$RNA_snn_res.0.1 == 3] <- "Endothelial"
matrix@meta.data$cell_types[matrix@meta.data$RNA_snn_res.0.1 == 5] <- "Fibroblast"
matrix@meta.data$cell_types[matrix@meta.data$RNA_snn_res.0.1 == 7] <- "Basal"
matrix@meta.data$cell_types[matrix@meta.data$RNA_snn_res.0.1 == 6] <- "Mastocyte"
matrix@meta.data$cell_types[matrix@meta.data$RNA_snn_res.0.1 == 9] <- "Luminal"
matrix@meta.data$cell_types[matrix@meta.data$RNA_snn_res.0.1 == 11] <- "Fibroblast"


## DONG 

matrix1@meta.data$cell_types <- NA

matrix1@meta.data$cell_types[matrix1@meta.data$RNA_snn_res.0.1 == 1] <- "T_cell"
matrix1@meta.data$cell_types[matrix1@meta.data$RNA_snn_res.0.1 == 10] <- "Fibroblast"
matrix1@meta.data$cell_types[matrix1@meta.data$RNA_snn_res.0.1 == 6] <- "Macrophage"
matrix1@meta.data$cell_types[matrix1@meta.data$RNA_snn_res.0.1 == 7] <- "Mastocyte"
matrix1@meta.data$cell_types[matrix1@meta.data$RNA_snn_res.0.1 == 11] <- "B_cell"
matrix1@meta.data$cell_types[matrix1@meta.data$RNA_snn_res.0.1 == 3] <- "Luminal_trans"
matrix1@meta.data$cell_types[matrix1@meta.data$RNA_snn_res.0.1 == 8] <- "Endothelial"
matrix1@meta.data$cell_types[matrix1@meta.data$RNA_snn_res.0.1 == 0] <- "Luminal"
matrix1@meta.data$cell_types[matrix1@meta.data$RNA_snn_res.0.1 == 2] <- "Fibroblast"
matrix1@meta.data$cell_types[matrix1@meta.data$RNA_snn_res.0.1 == 4] <- "Basal"
matrix1@meta.data$cell_types[matrix1@meta.data$RNA_snn_res.0.1 == 5] <- "Neuroendocrine"
matrix1@meta.data$cell_types[matrix1@meta.data$RNA_snn_res.0.1 == 9] <- "Luminal"
matrix1@meta.data$cell_types[matrix1@meta.data$RNA_snn_res.0.1 == 12] <- "B_cell"

## Basal or luminal analysis  
# CHENG

basal_luminal <- marcado[grep("Basal|Luminal", marcado$cell_type, ignore.case = TRUE), ]
clusters_of_interest <- c(0, 7, 10)
basorlum <- marks_fin[marks_fin$cluster %in% clusters_of_interest &
                                marks_fin$gene %in% basal_luminal$gene, ]

basorlum <- basorlum[order(basorlum$cell_type), ]
custom_order <- c()
basorlum$gene <- factor(basorlum$gene, levels = unique(basorlum$gene))
basal_rows <- grep("Basal", basorlum$cell_type, ignore.case = TRUE)
last_basal_row <- max(basal_rows)
basorlum <- basorlum[!basorlum$gene == "KRT19",]

# Plots 

ggplot(basorlum, aes(x = as.factor(cluster), y = gene)) +
  geom_point(aes(size = pct.1, color = avg_log2FC)) +
  scale_color_gradient(low = "blue", high = "red", name = "Fold Change") +
  scale_size_continuous(range = c(3, 10), name = "% Cells Expressing") +
  labs(
    x = "Cluster",
    y = "Gene",
    title = "Cluster Marker Expression (Basal & Luminal Markers)"
  ) + 
  theme_minimal() + theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    axis.text.y = element_text(size = 10, face = ifelse(levels(as.factor(basorlum$gene)) == "KRT19", "bold", "plain"), 
                               color = ifelse(levels(as.factor(basorlum$gene)) == "KRT19", "red", "black"))
  ) +
  annotate("segment", x = 0.5, xend = length(unique(basorlum$cluster)) + 0.5, 
           y = 5 + 0.5, yend = 5 + 0.5,
           linetype = "dashed", color = "gray") +
  annotate("segment", x = 0.5, xend = length(unique(basorlum$cluster)) + 0.5, 
           y = 11 + 0.5, yend = 11 + 0.5,
           linetype = "dashed", color = "gray")
 
 
 ## DONG 

basal_luminal <- marcado[grep("Basal|Luminal", marcado$cell_type, ignore.case = TRUE), ]
clusters_of_interest1 <- c(0, 3, 11)
basorlum1 <- marks_fin1[marks_fin1$cluster %in% clusters_of_interest1 &
                                marks_fin1$gene %in% basal_luminal$gene, ]

basorlum1 <- basorlum1[order(basorlum1$cell_type), ]
custom_order <- c("S100A6", "KRT15", "CYR61", "KRT17", "NKX3-1", "LMO7", "AR", "KRT8", "KRT18", "KRT19", "ALDH1A3", "PSCA") 
basorlum1$gene <- factor(basorlum1$gene, levels = custom_order)
basorlum1 <- basorlum1[!basorlum1$gene == "KRT19",]
basorlum1 <- basorlum1[!(basorlum1$gene == "ALDH1A3" & basorlum1$cluster == 3),]
basorlum1 <- basorlum1[!(basorlum1$gene == "PSCA" & basorlum1$cluster == 3),]
basorlum1 <- basorlum1[!(basorlum1$gene == "LMO7" & basorlum1$cluster == 3),]
basal_rows1 <- grep("Basal", basorlum1$cell_type, ignore.case = TRUE)

# Plots 

ggplot(basorlum1, aes(x = as.factor(cluster), y = gene)) +
  geom_point(aes(size = pct.1, color = avg_log2FC)) +
  scale_color_gradient(low = "blue", high = "red", name = "Fold Change") +
  scale_size_continuous(range = c(3, 10), name = "% Cells Expressing") +
  labs(
    x = "Cluster",
    y = "Gene",
    title = "Cluster Marker Expression (Basal & Luminal Markers)"
  ) +
  theme_minimal() +
 theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    axis.text.y = element_text(size = 10, face = ifelse(levels(as.factor(basorlum1$gene)) == "KRT19", "bold", "plain"), 
                               color = ifelse(levels(as.factor(basorlum1$gene)) == "KRT19", "red", "black"))
  ) +
  annotate("segment", x = 0.5, xend = length(unique(basorlum1$cluster)) + 0.5, 
           y = 4 + 0.5, yend = 4 + 0.5,
           linetype = "dashed", color = "gray") + 
  annotate("segment", x = 0.5, xend = length(unique(basorlum1$cluster)) + 0.5, 
           y = 9 + 0.5, yend = 9 + 0.5,
           linetype = "dashed", color = "gray")
 

## Cancer analysis 

cancermark1 <- marcado[grep("Proliferation|Cell Cycle|Apoptosis|Signaling|Metastasis|Adhesion|Tumor Suppression|Hypoxia Response|Stress Response|Inflammation", marcado$cell_type, ignore.case = TRUE), ]
clusters_of_interest1 <- c(4, 5, 0, 9, 2)
cancer1 <- marks_fin1[marks_fin1$cluster %in% clusters_of_interest1 &
                                marks_fin1$gene %in% cancermark1$gene, ]

cancermark <- marcado[grep("Proliferation|Cell Cycle|Apoptosis|Signaling|Metastasis|Adhesion|Tumor Suppression|Hypoxia Response|Stress Response|Inflammation", marcado$cell_type, ignore.case = TRUE), ]
clusters_of_interest <- c(1, 0, 7, 6, 3, 5, 11, 13, 14)
cancer <- marks_fin[marks_fin$cluster %in% clusters_of_interest &
                                marks_fin$gene %in% cancermark$gene, ]

## Cluster by cluster analisis
## CHENG

clust_comp <- FindMarkers(matrix1, 
               ident.1 = 2,
               ident.2 = 10,
               logfc.threshold = 0.25, 
               min.pct = 0.25, 
               only.pos = FALSE,
               test.use = 'wilcox',   
               )

clust_comp <- subset(clust_comp, (avg_log2FC < -1.5 | avg_log2FC > 1.5) & p_val_adj < 0.05 & pct.1 > 0.2)
clust_comp$cell_type <- ""
clust_comp$gene <- rownames(clust_comp)
for (n in marcado$gene){
  if (n %in% clust_comp$gene){
cell <- marcado %>% filter(gene == n) %>% pull(cell_type) 
clust_comp$cell_type[clust_comp$gene == n] <- cell
  }
}
 clustvsclust <- ""

clust_comp <- clust_comp[clust_comp$cell_type != "", ]
clustvsclust <- rbind(clustvsclust, clust_comp)

## Cancer cluster by cluster 

cancermarks1 <- marcado[grep("Proliferation|Cell Cycle|Apoptosis|Signaling|Metastasis|Adhesion|Tumor Suppression|Hypoxia Response|Stress Response|Inflammation", marcado$cell_type, ignore.case = TRUE), ]
cancercheng <- clustvsclust[clustvsclust$gene %in% cancermark1$gene, ]

basal_luminal <- marcado[grep("Neuroendocrine", marcado$cell_type, ignore.case = TRUE), ]
clusters_of_interest1 <- c(2)
basorlum1 <- marks_fin1[marks_fin1$cluster %in% clusters_of_interest1 &
                                marks_fin1$gene %in% basal_luminal$gene, ]

basorlum1 <- basorlum1[order(basorlum1$cell_type), ]
basorlum1$gene <- factor(basorlum1$gene, levels = unique(basorlum1$gene))

ggplot(basorlum1, aes(x = as.factor(cluster), y = gene)) +
  geom_point(aes(size = pct.1, color = avg_log2FC)) +
  scale_color_gradient(low = "blue", high = "red", name = "Fold Change") +
  scale_size_continuous(range = c(3, 10), name = "% Cells Expressing") +
  labs(
    x = "Cluster",
    y = "Gene",
    title = "Neuroendocrine markers"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


plot(density(as.vector(expr_tot)), 
     main = "Densidad de Expresión - Datos Originales", 
     xlab = "Valores de Expresión",
     col = "blue", lwd = 2)


