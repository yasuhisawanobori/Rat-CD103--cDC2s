#https://satijalab.org/seurat/articles/pbmc3k_tutorial
library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleR)
library(Rtools)
library(loupeR)
library(data.table)
setup()


rm(list=ls())

SplenicDCs <- Read10X(
  "PMID_33997687 mouse spleen CD11c+ scRNA-seq/",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

Splenic_object <- CreateSeuratObject(counts = SplenicDCs$`Gene Expression`, project = "test")
Splenic_object

SplenicDCs[c("CD3D", "TCL1A", "MS4A1")]
dense.size <- object.size(as.matrix(SplenicDCs))
dense.size
sparse.size <- object.size(SplenicDCs)
sparse.size
dense.size/sparse.size

Splenic_object[["percent.mt"]] <- PercentageFeatureSet(Splenic_object, pattern = "^mt-")
head(Splenic_object@meta.data, 5)
VlnPlot(Splenic_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(Splenic_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Splenic_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

Splenic_object <- subset(Splenic_object, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20)
Splenic_object <- NormalizeData(Splenic_object)

Splenic_object <- FindVariableFeatures(Splenic_object, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Splenic_object), 10)

plot1 <- VariableFeaturePlot(Splenic_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(Splenic_object)
Splenic_object <- ScaleData(Splenic_object, features = all.genes)

Splenic_object <- RunPCA(Splenic_object, features = VariableFeatures(object = Splenic_object))
print(Splenic_object[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Splenic_object, dims = 1:2, reduction = "pca")
DimPlot(Splenic_object, reduction = "pca")

DimHeatmap(Splenic_object, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Splenic_object, dims = 1:15, cells = 500, balanced = TRUE)

Splenic_object <- JackStraw(Splenic_object, num.replicate = 100)
Splenic_object <- ScoreJackStraw(Splenic_object, dims = 1:10)
JackStrawPlot(Splenic_object, dims = 1:10)
ElbowPlot(Splenic_object)

Splenic_object <- FindNeighbors(Splenic_object, k.param=16)
Splenic_object <- FindClusters(Splenic_object, resolution = 0.5)
head(Idents(Splenic_object), 5)

Splenic_object <- RunUMAP (Splenic_object, dims = 1:16)
DimPlot(Splenic_object, label = TRUE, reduction = "umap", label.box = T, repel = T)

saveRDS(Splenic_object, file = "PMID_33997687 mouse spleen CD11c+ scRNA-seq/splenic_DCs.rds")
Splenic_object <- readRDS(file = "PMID_33997687 mouse spleen CD11c+ scRNA-seq/splenic_DCs.rds")

#Add clusters from the metadata
Splenic_objectInfo<-read.csv("PMID_33997687 mouse spleen CD11c+ scRNA-seq/all_samples_metadata.csv", row.names = 1)

Splenic_object <- AddMetaData(Splenic_object,Splenic_objectInfo)

Idents(Splenic_object)<-Splenic_object@meta.data$RNA_snn_res.0.5

#0="immature phagocytic cDC2a", 1="mature Il4i1+ cDC2", 2="immature phagocytic cDC2b", 3="immature phagocytic Sox4+ cDC",
#4="immature phagocytic cDC1", 5="mitotic cDC2a-1", 6="immature mitotic cDC1", 7="mitotic cDC2a-2", 8="mature cDC2a",
#9="Mitochondrial cDC2 (mixed)", 10="mitotic cDC2a-3", 11="mature cDC2b", 12="mature cDC1", 13="NK/ILC/T", 14="pDC", 
#15="monocyte", 16="mixed identity", 17="Mitochondrial cDC1", 18="ribosomal cDC1", 19="stem cell"

Splenic_object_order <- c(0:19)
Splenic_object <- SetIdent(Splenic_object , value = factor(Idents(Splenic_object), levels = Splenic_object_order))

Splenic_objectClusterNames <- c("immature cDC2a", "mature Il4i1+ cDC2", "immature cDC2b", 
                                "immature Sox4+ cDC", "immature cDC1", "mitotic cDC2a-1",
                                "mitotic cDC1", "mitotic cDC2a-2", "mature cDC2a", "mitochondrial cDC2",
                                "mitotic cDC2a-3", "mature cDC2b", "mature cDC1", "NK/ILC/T", "pDC", "monocytes", "mixed identity", 
                                "mitochondrial cDC1", "ribosomal cDC1", "stem cells")

names(Splenic_objectClusterNames) <- levels(Idents(Splenic_object))
Splenic_object <- RenameIdents(object = Splenic_object, Splenic_objectClusterNames)

Splenic_object_order <- c("immature cDC1", "ribosomal cDC1", "mitotic cDC1", "mature cDC1", "mitochondrial cDC1",
                          "immature cDC2a", "immature Sox4+ cDC",  "mitotic cDC2a-1", "mitotic cDC2a-2",  "mitotic cDC2a-3",
                          "mitochondrial cDC2", "mature Il4i1+ cDC2", "mature cDC2a",
                          "immature cDC2b", "mature cDC2b", "monocytes", "pDC", "NK/ILC/T", "mixed identity", "stem cells")
Splenic_object <- SetIdent(Splenic_object , value = factor(Idents(Splenic_object), levels = Splenic_object_order))

saveRDS(Splenic_object, file = "PMID_33997687 mouse spleen CD11c+ scRNA-seq/splenic_DCs.rds")
Splenic_object <- readRDS(file = "PMID_33997687 mouse spleen CD11c+ scRNA-seq/splenic_DCs.rds")

DimPlot(Splenic_object, label = TRUE, reduction = "umap", label.box = T, repel = T, label.size = 3)+ theme(legend.position = "none")  #Fig. S5A

#find markers
Splenic_object_simple <- RenameIdents(Splenic_object,  "immature cDC2a"="cDC2a", "mitotic cDC2a-1"="cDC2a",
                                      "mitotic cDC2a-2"="cDC2a",  "mitotic cDC2a-3"="cDC2a", "mature cDC2a"="cDC2a",
                                      "mature cDC2b"="cDC2b", "immature cDC2b"="cDC2b")
write_rds(Splenic_object_simple, "PMID_33997687 mouse spleen CD11c+ scRNA-seq/Splenic_object_simple.rds")

SplenicDC.simple.markers <- FindAllMarkers(Splenic_object_simple, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
SplenicDC.simple.markers <- SplenicDC.simple.markers[order(SplenicDC.simple.markers$cluster, SplenicDC.simple.markers$avg_log2FC, decreasing = T),]

SplenicDC.simple.markers$rat_gene <- convert_orthologs(SplenicDC.simple.markers, gene_input = "gene", gene_output = "columns",
                                              input_species = "mouse", output_species = "rat",
                                              drop_nonorths = F, method = "gprofiler",
                                              non121_strategy = "kp")[,2]

SplenicDC.simple.markers <- SplenicDC.simple.markers[!SplenicDC.simple.markers$rat_gene=="N/A", ]

write.table(SplenicDC.simple.markers[SplenicDC.simple.markers$cluster == "cDC2a",]$rat_gene,
            file = "PMID_33997687 mouse spleen CD11c+ scRNA-seq/Spleen_cDC2a_signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

write.table(SplenicDC.simple.markers[SplenicDC.simple.markers$cluster == "cDC2b",]$rat_gene,
            file = "PMID_33997687 mouse spleen CD11c+ scRNA-seq/Spleen_cDC2b_signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

write.table(SplenicDC.simple.markers[SplenicDC.simple.markers$cluster == "monocytes",]$rat_gene,
            file = "PMID_33997687 mouse spleen CD11c+ scRNA-seq/Spleen_monocytes_signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

write.table(SplenicDC.simple.markers[SplenicDC.simple.markers$cluster == "monocytes",]$rat_gene,
            file = "PMID_33997687 mouse spleen CD11c+ scRNA-seq/Spleen_monocytes_signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

write.table(SplenicDC.simple.markers[SplenicDC.simple.markers$cluster == "mature Il4i1+ cDC2",]$rat_gene,
            file = "PMID_33997687 mouse spleen CD11c+ scRNA-seq/Spleen_mature_Il4i1+_cDC2_signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

#Inflammatory cDC2/MC marker application----
MouseLungDCs_markers[MouseLungDCs_markers$cluster == "inf-cDC2",]$gene

Splenic_object <- AddModuleScore(object = Splenic_object,
                                 features = list(MouseLungDCs_markers[MouseLungDCs_markers$cluster == "inf-cDC2",]$gene),
                                 name = "LungInfDC2score")

FeaturePlot(Splenic_object,"LungInfDC2score1")+ # 6x4inch
  labs(title = expression(paste("Lung Inflammatory cDC2 Score"))) + 
  theme(plot.title = element_text(hjust = 0.5))
VlnPlot(Splenic_object,"LungInfDC2score1")+theme(legend.position = "none")	+
  labs(title = expression(paste("Lung Inflammatory cDC2 Score")), x = "", y = "") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

Splenic_object@meta.data %>%
  group_by(ClusterNames) %>%
  summarise(avg_LungInfDC2score1 = mean(LungInfDC2score1, na.rm = TRUE))


MouseLungDCs_markers[MouseLungDCs_markers$cluster == "MC",]$gene

Splenic_object <- AddModuleScore(object = Splenic_object,
                                 features = list(MouseLungDCs_markers[MouseLungDCs_markers$cluster == "MC",]$gene),
                                 name = "LungMCscore")

FeaturePlot(Splenic_object,"LungMCscore1")+ # 6x4inch
  labs(title = expression(paste("Lung MC Score"))) + 
  theme(plot.title = element_text(hjust = 0.5))
VlnPlot(Splenic_object,"LungMCscore1")+theme(legend.position = "none")	+
  labs(title = expression(paste("Lung MC Score")), x = "", y = "") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

Splenic_object@meta.data %>%
  group_by(ClusterNames) %>%
  summarise(avg_LungMCscore1 = mean(LungMCscore1, na.rm = TRUE))