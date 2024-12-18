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

#GSEA----
#re-clustering
MouseLungDCs_simple <- MouseLungDCs
                                   
Idents(object = MouseLungDCs_simple, 
       cells=intersect(rownames(subset(MouseLungDCs_simple, subset=umap_2>-5)@meta.data), 
                       WhichCells(MouseLungDCs_simple, idents = "proliferating DC")))<- "proliferating cDC2"

Idents(object = MouseLungDCs_simple, 
       cells=intersect(rownames(subset(MouseLungDCs_simple, subset=umap_2 < -5)@meta.data), 
                       WhichCells(MouseLungDCs_simple, idents = "proliferating DC")))<- "proliferating cDC1"

MouseLungDCs_simple <- RenameIdents(MouseLungDCs_simple, "proliferating cDC2"="cDC2a", "Migratory cDC2" ="cDC2a",
                                    "Non-migratory cDC2" ="cDC2a")

#calculate log2FC
MoLungMcDC2aExpression <- FindMarkers(MouseLungDCs_simple, ident.1 = "MC", ident.2 = "cDC2a")
MoLungInfDC2DC2aExpression <- FindMarkers(MouseLungDCs_simple, ident.1 = "inf-cDC2", ident.2 = "cDC2a")

#save calculated log2FC
write.csv(MoLungMcDC2aExpression, "PMID_32392463 mouse lung MC and inf-cDC2/MoLungMcDC2aExpression.csv", row.names = T)
write.csv(MoLungInfDC2DC2aExpression, "PMID_32392463 mouse lung MC and inf-cDC2/MoLungInfDC2DC2aExpression.csv", row.names = T)

MoLungMcDC2aExpression <- read.csv("PMID_32392463 mouse lung MC and inf-cDC2/MoLungMcDC2aExpression.csv", row.names = 1)
MoLungInfDC2DC2aExpression <- read.csv("PMID_32392463 mouse lung MC and inf-cDC2/MoLungInfDC2DC2aExpression.csv", row.names = 1)

#format log2FC data
MoLungMcDC2aExpression <- MoLungMcDC2aExpression[order(MoLungMcDC2aExpression$avg_log2FC, decreasing = T),]
MoLungMcDC2aExpressionLog2FC <- MoLungMcDC2aExpression$avg_log2FC
names(MoLungMcDC2aExpressionLog2FC) <- rownames(MoLungMcDC2aExpression)

MoLungInfDC2DC2aExpression <- MoLungInfDC2DC2aExpression[order(MoLungInfDC2DC2aExpression$avg_log2FC, decreasing = T),]
MoLungInfDC2DC2aExpressionLog2FC <- MoLungInfDC2DC2aExpression$avg_log2FC
names(MoLungInfDC2DC2aExpressionLog2FC) <- rownames(MoLungInfDC2DC2aExpression)

#execute the GSEA
MoLungMcDC2a_GSEA <- gseGO(geneList = MoLungMcDC2aExpressionLog2FC, OrgDb = "org.Mm.eg.db", ont = "BP", keyType = "SYMBOL", pAdjustMethod="none")

MoLungInfDC2DC2a_GSEA <- gseGO(geneList = MoLungInfDC2DC2aExpressionLog2FC, OrgDb = "org.Mm.eg.db", ont = "BP", keyType = "SYMBOL", pAdjustMethod="none")

#divide CD103- equiv and CD103+equiv dominant gene sets
MoLungMcDom <- MoLungMcDC2a_GSEA
MoLungMcDom@result <- subset(MoLungMcDC2a_GSEA, MoLungMcDC2a_GSEA@result$enrichmentScore >0)

MoLungInfDC2Dom <- MoLungInfDC2DC2a_GSEA
MoLungInfDC2Dom@result <- subset(MoLungInfDC2DC2a_GSEA, MoLungInfDC2DC2a_GSEA@result$enrichmentScore >0)

#depict figures
set.seed(123)
emapplot(pairwise_termsim(MoLungMcDom), layout.params = list(layout = "kk"), color="p.adjust", showCategory=30,
         cluster.params = list(cluster = T, legend=T, n=5), cex_label_group = 1.5, node_label = "none")
set.seed(123)
emapplot(pairwise_termsim(MoLungInfDC2Dom), layout.params = list(layout = "kk"), color="p.adjust", showCategory=30,
         cluster.params = list(cluster = T, legend=T, n=5), cex_label_group = 1.5, node_label = "none")

treeplot(pairwise_termsim(MoLungMcDom), showCategory = 30, color = "p.adjust", # 9.65x5 inches
         offset.params = list(bar_tree = rel(5), extend=0.5, hexpand = 0.15), hilight.params = list(hilight = T),
         fontsize = 4, label_format = 20, cluster.params = list(label_format = 10), label_format_tiplab = 100, 
         cex_category = 0.75 )

treeplot(pairwise_termsim(MoLungInfDC2Dom), showCategory = 30, color = "p.adjust", # 9.65x5 inches
         offset.params = list(bar_tree = rel(5), extend=0.5, hexpand = 0.15), hilight.params = list(hilight = T),
         fontsize = 4, label_format = 20, cluster.params = list(label_format = 10), label_format_tiplab = 100, 
         cex_category = 0.75 )

gseaplot2(MoLungMcDC2a_GSEA, 
          geneSetID = c("GO:0006909", "GO:0002253", "GO:0001819",
                        "GO:0006260", "GO:0007059", "GO:0044839"),
          subplots = 1:2, pvalue_table = F, base_size = 14)

gseaplot2(MoLungInfDC2DC2a_GSEA, 
          geneSetID = c("GO:0006909", "GO:0002253", "GO:0001819",
                        "GO:0006260", "GO:0007059", "GO:0044839"),
          subplots = 1:2, pvalue_table = F, base_size = 14)

#save GSEA results
write.xlsx(list("MoLungMcDom"=MoLungMcDom@result, 
                "MoLungInfDC2DomDom"=MoLungInfDC2Dom@result), 
           "PMID_32392463 mouse lung MC and inf-cDC2/MoLung_GSEA.xlsx")
