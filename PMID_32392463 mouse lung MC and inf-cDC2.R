#https://satijalab.org/seurat/articles/pbmc3k_tutorial
library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleR)
library(Rtools)
library(loupeR)
library(data.table)
setup()

#data source: http://bioit2.irc.ugent.be/cdc2/
MouseLungDCsMtx <- Read10X(
  "PMID_32392463 mouse lung MC and inf-cDC2/mm10",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

MouseLungDCs <- CreateSeuratObject(counts = MouseLungDCsMtx)

MouseLungDCs

MouseLungDCs[["percent.mt"]] <- PercentageFeatureSet(MouseLungDCs, pattern = "^mt-")
head(MouseLungDCs@meta.data, 5)
VlnPlot(MouseLungDCs, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(MouseLungDCs, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(MouseLungDCs, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

MouseLungDCs <- subset(MouseLungDCs, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 5)
MouseLungDCs <- NormalizeData(MouseLungDCs)

MouseLungDCs <- FindVariableFeatures(MouseLungDCs, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(MouseLungDCs), 10)

plot1 <- VariableFeaturePlot(MouseLungDCs)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(MouseLungDCs)
MouseLungDCs <- ScaleData(MouseLungDCs, features = all.genes)

MouseLungDCs <- RunPCA(MouseLungDCs, features = VariableFeatures(object = MouseLungDCs))
print(MouseLungDCs[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(MouseLungDCs, dims = 1:2, reduction = "pca")
DimPlot(MouseLungDCs, reduction = "pca")

DimHeatmap(MouseLungDCs, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(MouseLungDCs, dims = 1:15, cells = 500, balanced = TRUE)

MouseLungDCs <- JackStraw(MouseLungDCs, num.replicate = 100)
MouseLungDCs <- ScoreJackStraw(MouseLungDCs, dims = 1:10)
JackStrawPlot(MouseLungDCs, dims = 1:10)
ElbowPlot(MouseLungDCs)

MouseLungDCs <- FindNeighbors(MouseLungDCs, k.param=16)
MouseLungDCs <- FindClusters(MouseLungDCs, resolution = 0.5)
head(Idents(MouseLungDCs), 5)

MouseLungDCs <- RunUMAP (MouseLungDCs, dims = 1:16)
DimPlot(MouseLungDCs, label = TRUE, reduction = "umap", label.box = T, repel = T)

saveRDS(MouseLungDCs, "PMID_32392463 mouse lung MC and inf-cDC2/MouseLungDCs.rds")
MouseLungDCs <- readRDS("PMID_32392463 mouse lung MC and inf-cDC2/MouseLungDCs.rds")

#add metadata
MouseLungDCsInfo<-read.csv("PMID_32392463 mouse lung MC and inf-cDC2/annot_Fig_S5.csv", row.names = 1)

MouseLungDCs@reductions$umap@cell.embeddings[rownames(MouseLungDCs@reductions$umap@cell.embeddings), "umap_1"] <- 
  MouseLungDCsInfo[rownames(MouseLungDCs@reductions$umap@cell.embeddings), "UMAP_1"]
MouseLungDCs@reductions$umap@cell.embeddings[rownames(MouseLungDCs@reductions$umap@cell.embeddings), "umap_2"] <- 
  MouseLungDCsInfo[rownames(MouseLungDCs@reductions$umap@cell.embeddings), "UMAP_2"]

MouseLungDCs<-AddMetaData(MouseLungDCs,MouseLungDCsInfo)

Idents(MouseLungDCs)<-MouseLungDCs@meta.data$cluster

MouseLungDCs_order <- c("Non-migr cDC1", "Migr cDC1", "Non-migr cDC2", "Migr cDC2", "inf-cDC2", "pDC", "Prolif DC", "MC")
MouseLungDCs <- SetIdent(MouseLungDCs , value = factor(Idents(MouseLungDCs), levels = MouseLungDCs_order))

MouseLungDCs <- subset(MouseLungDCs, cluster %in% na.omit(MouseLungDCs$cluster))

MouseLungDCs <- RenameIdents(MouseLungDCs, "Migr cDC1"="Migratory cDC1", "Non-migr cDC1"="Non-migratory cDC1", "Prolif DC"="proliferating DC",
             "Migr cDC2"="Migratory cDC2", "Non-migr cDC2"= "Non-migratory cDC2")

saveRDS(MouseLungDCs, "PMID_32392463 mouse lung MC and inf-cDC2/MouseLungDCs.rds")
MouseLungDCs <- readRDS("PMID_32392463 mouse lung MC and inf-cDC2/MouseLungDCs.rds")

#Fig. S5B, 5.5x4 inches
DimPlot(MouseLungDCs, label = TRUE, reduction = "umap", label.box = T, repel = T) +
  theme(legend.position = "none")

#find markers----
MouseLungDCs_markers <- FindAllMarkers(MouseLungDCs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
MouseLungDCs_markers <- MouseLungDCs_markers[order(MouseLungDCs_markers$cluster, MouseLungDCs_markers$avg_log2FC, decreasing = T),]

MouseLungDCs_markers$rat_gene <- convert_orthologs(MouseLungDCs_markers, gene_input = "gene", gene_output = "columns",
                                                input_species = "mouse", output_species = "rat",
                                                drop_nonorths = F, method = "gprofiler",
                                                non121_strategy = "kp")[,2]

MouseLungDCs_markers <- MouseLungDCs_markers[!MouseLungDCs_markers$rat_gene=="N/A", ]

write.table(MouseLungDCs_markers[MouseLungDCs_markers$cluster == "Non-migratory cDC2",]$rat_gene,
            file = "PMID_32392463 mouse lung MC and inf-cDC2/LungNonMigrDC2signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

write.table(MouseLungDCs_markers[MouseLungDCs_markers$cluster == "Migratory cDC2",]$rat_gene,
            file = "PMID_32392463 mouse lung MC and inf-cDC2/LungMigrDC2signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)


write.table(MouseLungDCs_markers[MouseLungDCs_markers$cluster == "inf-cDC2",]$rat_gene,
            file = "PMID_32392463 mouse lung MC and inf-cDC2/LungInfDC2signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

write.table(MouseLungDCs_markers[MouseLungDCs_markers$cluster == "MC",]$rat_gene,
            file = "PMID_32392463 mouse lung MC and inf-cDC2/LungMCsignature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

#infDC2 vs MC markers
infDC2vsMC_markers <- FindMarkers(MouseLungDCs, ident.1 = "inf-cDC2", ident.2 = "MC", only.pos = F, min.pct = 0.25, logfc.threshold = 1)
infDC2vsMC_markers <- infDC2vsMC_markers[order(infDC2vsMC_markers$avg_log2FC, decreasing = T),]

infDC2vsMC_markers$rat_gene <- convert_orthologs(infDC2vsMC_markers, gene_input = "rownames", gene_output = "columns",
                                                   input_species = "mouse", output_species = "rat",
                                                   drop_nonorths = F, method = "gprofiler",
                                                   non121_strategy = "kp")[,2]

write.table(infDC2vsMC_markers[infDC2vsMC_markers$avg_log2FC>=1,]$rat_gene,
            file = "PMID_32392463 mouse lung MC and inf-cDC2/Lung_infDC2vsMCsignature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

write.table(infDC2vsMC_markers[infDC2vsMC_markers$avg_log2FC<=-1,]$rat_gene,
            file = "PMID_32392463 mouse lung MC and inf-cDC2/Lung_MCvsinfDC2signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)