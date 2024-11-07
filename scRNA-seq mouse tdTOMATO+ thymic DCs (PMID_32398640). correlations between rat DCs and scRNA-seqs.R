#https://satijalab.org/seurat/articles/pbmc3k_tutorial
library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)
library(org.Mm.eg.db)
library(limma)
library("pheatmap")

#read matrix----
DCdata1 <- as.data.frame(read_delim("PMID_32398640 Thymic monocyte-derived DCs/DC1_symbol.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE))
DCdata2 <- as.data.frame(read_delim("PMID_32398640 Thymic monocyte-derived DCs/DC2_symbol.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE))

#add unique numbers to duplicated genes and set gene symbols to row names of DCdata1.----
gene_symbols <- DCdata1[,1]

#add unique numbers to duplicated symbols
dup_indices <- duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE)
generate_suffix <- function(num) {
  if(num == 1) {
    return("")
  } else {
    return(paste0("_", num - 1))
  }
}
generate_unique_id <- function(x) {
  if (!any(duplicated(x))) {
    return(x)
  }
  counts <- table(x)
  seq_suffix <- as.numeric(ave(x, x, FUN = seq_along))
  suffix <- sapply(seq_suffix, generate_suffix)
  paste(x, suffix, sep = "")
}

unique_gene_symbols <- gene_symbols
unique_gene_symbols[dup_indices] <- generate_unique_id(gene_symbols[dup_indices])

#finally add symbols to the matrix
rownames(DCdata1) <- unique_gene_symbols
rm(unique_gene_symbols, gene_symbols, na_indices, ensembl_ids, dup_indices)



#add unique numbers to duplicated genes and set gene symbols to row names of DCdata2.----
gene_symbols <- DCdata2[,1]

#add unique numbers to duplicated symbols
dup_indices <- duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE)
generate_suffix <- function(num) {
  if(num == 1) {
    return("")
  } else {
    return(paste0("_", num - 1))
  }
}
generate_unique_id <- function(x) {
  if (!any(duplicated(x))) {
    return(x)
  }
  counts <- table(x)
  seq_suffix <- as.numeric(ave(x, x, FUN = seq_along))
  suffix <- sapply(seq_suffix, generate_suffix)
  paste(x, suffix, sep = "")
}

unique_gene_symbols <- gene_symbols
unique_gene_symbols[dup_indices] <- generate_unique_id(gene_symbols[dup_indices])

#finally add symbols to the matrix
rownames(DCdata2) <- unique_gene_symbols
rm(unique_gene_symbols, gene_symbols, na_indices, ensembl_ids, dup_indices)

#exchange data sets to Seurat objects----
DCdata1 <-DCdata1[,-1]
DCdata2 <-DCdata2[,-1]

DCdata1obj <- CreateSeuratObject(counts = DCdata1, project = "DCdata1")
DCdata2obj <- CreateSeuratObject(counts = DCdata2, project = "DCdata2")
DCdata1obj <- NormalizeData(DCdata1obj)
DCdata2obj <- NormalizeData(DCdata2obj)

integrated_TOMATO_DCs <- merge(DCdata1obj, DCdata2obj, 
                               add.cell.ids =c("DCdata1", "DCdata2"), project= "TOMATO_DCs")
#cut off low quality cells
integrated_TOMATO_DCs[["percent.mt"]] <- PercentageFeatureSet(integrated_TOMATO_DCs, pattern = "^mt-")
head(integrated_TOMATO_DCs@meta.data, 5)
VlnPlot(integrated_TOMATO_DCs, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mean(integrated_TOMATO_DCs@meta.data$nFeature_RNA)

plot1 <- FeatureScatter(integrated_TOMATO_DCs, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(integrated_TOMATO_DCs, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

integrated_TOMATO_DCs <- subset(integrated_TOMATO_DCs, subset = nFeature_RNA < 2000 & percent.mt < 2)

integrated_TOMATO_DCs <- NormalizeData(integrated_TOMATO_DCs)

integrated_TOMATO_DCs <- FindVariableFeatures(integrated_TOMATO_DCs, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(integrated_TOMATO_DCs), 10)

plot1 <- VariableFeaturePlot(integrated_TOMATO_DCs)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(integrated_TOMATO_DCs)
integrated_TOMATO_DCs <- ScaleData(integrated_TOMATO_DCs, features = all.genes)


integrated_TOMATO_DCs <- RunPCA(integrated_TOMATO_DCs, features = VariableFeatures(object = integrated_TOMATO_DCs))
print(integrated_TOMATO_DCs[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(integrated_TOMATO_DCs, dims = 1:2, reduction = "pca")
DimPlot(integrated_TOMATO_DCs, reduction = "pca")

DimHeatmap(integrated_TOMATO_DCs, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(integrated_TOMATO_DCs, dims = 1:15, cells = 500, balanced = TRUE)

integrated_TOMATO_DCs <- JackStraw(integrated_TOMATO_DCs, num.replicate = 100)
integrated_TOMATO_DCs <- ScoreJackStraw(integrated_TOMATO_DCs, dims = 1:15)
JackStrawPlot(integrated_TOMATO_DCs, dims = 1:15)
ElbowPlot(integrated_TOMATO_DCs)

integrated_TOMATO_DCs <- FindNeighbors(integrated_TOMATO_DCs)
integrated_TOMATO_DCs <- FindClusters(integrated_TOMATO_DCs, resolution = 0.7)
head(Idents(integrated_TOMATO_DCs), 5)

integrated_TOMATO_DCs <- RunUMAP (integrated_TOMATO_DCs, dims = 1:15)
DimPlot(integrated_TOMATO_DCs, label = TRUE, reduction = "umap")

DimPlot(integrated_TOMATO_DCs , split.by = "orig.ident")

#identification of clusters
FeaturePlot(integrated_TOMATO_DCs, features = c("Ccr7", "Relb", "Batf3", "Ly75", "Cxcr5", #cDC1a
                                                "Cd8a", "Clec9a", "Cadm1", "Xcr1", "Itgae" #cDC1b
                                                ))

FeaturePlot(integrated_TOMATO_DCs, features = c("Sirpa", "Irf4", #cDC2+moDC
                                                "Mgl2", "Cd209a", #cDC2
                                                "Itgam", "Cx3cr1", "Cxcr4","Lyz2", 
                                                "Mertk", "Fcgr2b", "Csf2ra", "Csf2rb", "Mafb",
                                                "Apoe", "Cd14", "Cebpb", "Fcgr1", "Lyz1", "Fcgr3",
                                                "Fn1", "Ccr2", "Ifitm3", "Ccr1", "Ccr5"))  #moDC

FeaturePlot(integrated_TOMATO_DCs, features = c("Bst2", "Ly6d", "Iglc3", "Siglech", "Klk1", "Ccr6", "Ccr9", "Runx3", "Ly6c1", "Ly6c2")) #pDC


#name the clusters
integrated_TOMATO_DCsNames <- c("cDC2", "cDC1b", "pDC", "cDC1a", "cDC1a")
names(integrated_TOMATO_DCsNames) <- levels(Idents(integrated_TOMATO_DCs))
integrated_TOMATO_DCs <- RenameIdents(object = integrated_TOMATO_DCs, integrated_TOMATO_DCsNames)

Idents(object = integrated_TOMATO_DCs, 
       cells=CellSelector(plot = DimPlot(integrated_TOMATO_DCs, label = TRUE, reduction = "umap")))<- "moDC"

#set the order of clusters
integrated_TOMATO_DCs_order <- c("cDC1a", "cDC1b", "cDC2", "moDC", "pDC")
integrated_TOMATO_DCs <- SetIdent(integrated_TOMATO_DCs, value = factor(Idents(integrated_TOMATO_DCs), levels = integrated_TOMATO_DCs_order))

DimPlot(integrated_TOMATO_DCs, label = TRUE, reduction = "umap", label.box = T, repel = T) + theme(legend.position = "none") #Fig. S5C

saveRDS(integrated_TOMATO_DCs, file = "PMID_32398640 Thymic monocyte-derived DCs/integrated_TOMATO_DCs.rds")
integrated_TOMATO_DCs <- readRDS("PMID_32398640 Thymic monocyte-derived DCs/integrated_TOMATO_DCs.rds")

#find markers----
TOMATO_DC_markers <- FindAllMarkers(JoinLayers(integrated_TOMATO_DCs), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
TOMATO_DC_markers <- TOMATO_DC_markers[order(TOMATO_DC_markers$cluster, TOMATO_DC_markers$avg_log2FC, decreasing = T),]

TOMATO_DC_markers$rat_gene <- convert_orthologs(TOMATO_DC_markers, gene_input = "gene", gene_output = "columns",
                                           input_species = "mouse", output_species = "rat",
                                           drop_nonorths = F, method = "gprofiler",
                                           non121_strategy = "kp")[,2]

TOMATO_DC_markers <- TOMATO_DC_markers[!TOMATO_DC_markers$rat_gene=="N/A", ]

write.table(TOMATO_DC_markers[TOMATO_DC_markers$cluster == "cDC1a",]$rat_gene,
            file = "PMID_32398640 Thymic monocyte-derived DCs/TOMATO_cDC1a_signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

write.table(TOMATO_DC_markers[TOMATO_DC_markers$cluster == "cDC1b",]$rat_gene,
            file = "PMID_32398640 Thymic monocyte-derived DCs/TOMATO_cDC1b_signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

write.table(TOMATO_DC_markers[TOMATO_DC_markers$cluster == "cDC2",]$rat_gene,
            file = "PMID_32398640 Thymic monocyte-derived DCs/TOMATO_cDC2_signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

write.table(TOMATO_DC_markers[TOMATO_DC_markers$cluster == "moDC",]$rat_gene,
            file = "PMID_32398640 Thymic monocyte-derived DCs/TOMATO_moDC_signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

write.table(TOMATO_DC_markers[TOMATO_DC_markers$cluster == "pDC",]$rat_gene,
            file = "PMID_32398640 Thymic monocyte-derived DCs/TOMATO_pDC_signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)


#Inflammatory cDC2/MC marker application----
MouseLungDCs_markers[MouseLungDCs_markers$cluster == "inf-cDC2",]$gene

integrated_TOMATO_DCs<-JoinLayers(integrated_TOMATO_DCs)
integrated_TOMATO_DCs@meta.data$ClusterNames <- integrated_TOMATO_DCs@active.ident

integrated_TOMATO_DCs <- AddModuleScore(object = integrated_TOMATO_DCs,
                                 features = list(MouseLungDCs_markers[MouseLungDCs_markers$cluster == "inf-cDC2",]$gene),
                                 name = "LungInfDC2score")

FeaturePlot(integrated_TOMATO_DCs,"LungInfDC2score1")+ # 6x4inch
  labs(title = expression(paste("Lung Inflammatory cDC2 Score"))) + 
  theme(plot.title = element_text(hjust = 0.5))
VlnPlot(integrated_TOMATO_DCs,"LungInfDC2score1")+theme(legend.position = "none")	+
  labs(title = expression(paste("Lung Inflammatory cDC2 Score")), x = "", y = "") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

integrated_TOMATO_DCs@meta.data %>%
  group_by(ClusterNames) %>%
  summarise(avg_LungInfDC2score1 = mean(LungInfDC2score1, na.rm = TRUE))


MouseLungDCs_markers[MouseLungDCs_markers$cluster == "MC",]$gene

integrated_TOMATO_DCs <- AddModuleScore(object = integrated_TOMATO_DCs,
                                 features = list(MouseLungDCs_markers[MouseLungDCs_markers$cluster == "MC",]$gene),
                                 name = "LungMCscore")

FeaturePlot(integrated_TOMATO_DCs,"LungMCscore1")+ # 6x4inch
  labs(title = expression(paste("Lung MC Score"))) + 
  theme(plot.title = element_text(hjust = 0.5))
VlnPlot(integrated_TOMATO_DCs,"LungMCscore1")+theme(legend.position = "none")	+
  labs(title = expression(paste("Lung MC Score")), x = "", y = "") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

integrated_TOMATO_DCs@meta.data %>%
  group_by(ClusterNames) %>%
  summarise(avg_LungMCscore1 = mean(LungMCscore1, na.rm = TRUE))