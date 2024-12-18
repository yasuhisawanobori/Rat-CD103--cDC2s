#https://satijalab.org/seurat/articles/pbmc3k_tutorial
library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleR)
library(loupeR)
library(data.table)
setup()


human_splenic_DC_matrix <-read.delim ("PMID_31668803_human/human_spleen_raw_counts.tsv", row.names=1)
human_splenic_DC_matrix <- t(human_splenic_DC_matrix)

write.csv(human_splenic_DC_matrix, file = "PMID_31668803_human/human_spleen_raw_counts.csv", row.names = TRUE)

human_splenic_DC <- CreateSeuratObject(human_splenic_DC_matrix)
human_splenic_DC

rm(human_splenic_DC_matrix)

human_splenic_DC[["percent.mt"]] <- PercentageFeatureSet(human_splenic_DC, pattern = "^MT.")

VlnPlot(human_splenic_DC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(human_splenic_DC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(human_splenic_DC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

human_splenic_DC <- subset(human_splenic_DC, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 7)
human_splenic_DC <- NormalizeData(human_splenic_DC)

human_splenic_DC <- FindVariableFeatures(human_splenic_DC, selection.method = "vst", nfeatures = 2000)
hsd_top10 <- head(VariableFeatures(human_splenic_DC), 10)

plot1 <- VariableFeaturePlot(human_splenic_DC)
plot2 <- LabelPoints(plot = plot1, points = hsd_top10, repel = TRUE)
plot1 + plot2

hsd.all.genes <- rownames(human_splenic_DC)
human_splenic_DC <- ScaleData(human_splenic_DC, features = hsd.all.genes)

human_splenic_DC <- RunPCA(human_splenic_DC, features = VariableFeatures(object = human_splenic_DC))
print(human_splenic_DC[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(human_splenic_DC, dims = 1:2, reduction = "pca")
DimPlot(human_splenic_DC, reduction = "pca")

DimHeatmap(human_splenic_DC, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(human_splenic_DC, dims = 1:15, cells = 500, balanced = TRUE)

human_splenic_DC <- JackStraw(human_splenic_DC, num.replicate = 100)
human_splenic_DC <- ScoreJackStraw(human_splenic_DC, dims = 1:20)
JackStrawPlot(human_splenic_DC, dims = 1:20)
ElbowPlot(human_splenic_DC)

human_splenic_DC <- FindNeighbors(human_splenic_DC)
human_splenic_DC <- FindClusters(human_splenic_DC, resolution = 0.75)
head(Idents(human_splenic_DC), 5)

human_splenic_DC <- RunUMAP (human_splenic_DC, dims = 1:16)
human_splenic_DC <- RunTSNE(human_splenic_DC, dims = 1:16)

DimPlot(human_splenic_DC, label = TRUE, reduction = "umap")
DimPlot(human_splenic_DC, label = TRUE, reduction = "tsne")

cluster_names <- levels(Idents(human_splenic_DC))

saveRDS(human_splenic_DC, file = "PMID_31668803_human/human_splenic_DC.rds")
human_splenic_DC <- readRDS(file = "PMID_31668803_human/human_splenic_DC.rds")

#import metadata
human_splenic_DCInfo <- read.delim("PMID_31668803_human/human_spleen_cell_metadata_4465x9.tsv", row.names = 1)

human_splenic_DC@reductions$tsne@cell.embeddings[rownames(human_splenic_DC@reductions$tsne@cell.embeddings), "tSNE_1"] <- 
  human_splenic_DCInfo[rownames(human_splenic_DC@reductions$tsne@cell.embeddings), "tsne_x"]
human_splenic_DC@reductions$tsne@cell.embeddings[rownames(human_splenic_DC@reductions$tsne@cell.embeddings), "tSNE_2"] <- 
  human_splenic_DCInfo[rownames(human_splenic_DC@reductions$tsne@cell.embeddings), "tsne_y"]

human_splenic_DC <-AddMetaData (human_splenic_DC, human_splenic_DCInfo)

Idents(human_splenic_DC)<-human_splenic_DC@meta.data$cell_type

human_splenic_DC_order <- c("cDC1", "Mitotic cDC1", "CLEC10A+ cDC2", "CLEC10A- cDC2", "CCR7+ cDC2", "Mitotic cDC2", "AS DC")
human_splenic_DC <- SetIdent(human_splenic_DC, value = factor(Idents(human_splenic_DC), levels = human_splenic_DC_order))

human_splenic_DC <- subset(human_splenic_DC, cell_type %in% na.omit(human_splenic_DC$cell_type))

DimPlot(human_splenic_DC, label = TRUE, reduction = "tsne", label.box = T, repel = T) + theme(legend.position = "none") #Fig. S6A

saveRDS(human_splenic_DC, file = "PMID_31668803_human/human_splenic_DC.rds")
human_splenic_DC <- readRDS(file = "PMID_31668803_human/human_splenic_DC.rds")

#find markers----
human_splenic_DC_markers <- FindAllMarkers(human_splenic_DC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
human_splenic_DC_markers <- human_splenic_DC_markers[order(human_splenic_DC_markers$cluster, human_splenic_DC_markers$avg_log2FC, decreasing = T),]

human_splenic_DC_markers$rat_gene <- convert_orthologs(human_splenic_DC_markers, gene_input = "gene", gene_output = "columns",
                                                   input_species = "human", output_species = "rat",
                                                   drop_nonorths = F, method = "gprofiler",
                                                   non121_strategy = "kp")[,2]

human_splenic_DC_markers <- human_splenic_DC_markers[!human_splenic_DC_markers$rat_gene=="N/A", ]

write.table(human_splenic_DC_markers[human_splenic_DC_markers$cluster == "CLEC10A+ cDC2",]$rat_gene,
            file = "PMID_31668803_human/HumanSpl_CLEC10A+_cDC2_signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

write.table(human_splenic_DC_markers[human_splenic_DC_markers$cluster == "CLEC10A- cDC2",]$rat_gene,
            file = "PMID_31668803_human/HumanSpl_CLEC10A-_cDC2_signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

write.table(human_splenic_DC_markers[human_splenic_DC_markers$cluster == "CCR7+ cDC2",]$rat_gene,
            file = "PMID_31668803_human/HumanSpl_CCR7+_cDC2_signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

write.table(human_splenic_DC_markers[human_splenic_DC_markers$cluster == "Mitotic cDC2",]$rat_gene,
            file = "PMID_31668803_human/HumanSpl_Mitotic_cDC2_signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

write.table(human_splenic_DC_markers[human_splenic_DC_markers$cluster == "AS DC",]$rat_gene,
            file = "PMID_31668803_human/HumanSpl_AS_DC_signature.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

#GSEA----
#calculate log2FC
human_splenic_DC_simple <- RenameIdents(human_splenic_DC, "CLEC10A+ cDC2"="CD103nEq", "CCR7+ cDC2"="CD103nEq",
                                        "CLEC10A- cDC2"="CD103pEq", "Mitotic cDC2"="CD103pEq")
HuSplCD103nCD103pExpression <- FindMarkers(human_splenic_DC_simple, ident.1 = "CD103nEq", ident.2 = "CD103pEq")

#save calculated log2FC
write.csv(HuSplCD103nCD103pExpression, "PMID_31668803_human/HuSplCD103nCD103pExpression.csv", row.names = T)

HuSplCD103nCD103pExpression  <- read.csv("PMID_31668803_human/HuSplCD103nCD103pExpression.csv", row.names = 1)

#format log2FC data
HuSplCD103nCD103pExpression <- HuSplCD103nCD103pExpression[order(HuSplCD103nCD103pExpression$avg_log2FC, decreasing = T),]
HuSplCD103nCD103pExpressionLog2FC <- HuSplCD103nCD103pExpression$avg_log2FC
names(HuSplCD103nCD103pExpressionLog2FC) <- rownames(HuSplCD103nCD103pExpression)

#execute the GSEA
HuSplCD103nCD103p_GSEA <- gseGO(geneList = HuSplCD103nCD103pExpressionLog2FC, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", pAdjustMethod="none")

#divide CD103- equiv and CD103+equiv dominant gene sets
HuSplCD103nDom <- HuSplCD103nCD103p_GSEA
HuSplCD103nDom@result <- subset(HuSplCD103nCD103p_GSEA, HuSplCD103nCD103p_GSEA@result$enrichmentScore >0)

HuSplCD103pDom <- HuSplCD103nCD103p_GSEA
HuSplCD103pDom@result <- subset(HuSplCD103nCD103p_GSEA, HuSplCD103nCD103p_GSEA@result$enrichmentScore <0)

#depict figures
set.seed(123)
emapplot(pairwise_termsim(HuSplCD103nDom), layout.params = list(layout = "kk"), color="p.adjust", showCategory=30,
         cluster.params = list(cluster = T, legend=T, n=5), cex_label_group = 1.5, node_label = "none")
set.seed(123)
emapplot(pairwise_termsim(HuSplCD103pDom), layout.params = list(layout = "kk"), color="p.adjust", showCategory=30,
         cluster.params = list(cluster = T, legend=T, n=5), cex_label_group = 1.5, node_label = "none")

treeplot(pairwise_termsim(HuSplCD103nDom), showCategory = 30, color = "p.adjust", # 9.65x5 inches
         offset.params = list(bar_tree = rel(5), extend=0.5, hexpand = 0.15), hilight.params = list(hilight = T),
         fontsize = 4, label_format = 20, cluster.params = list(label_format = 10), label_format_tiplab = 100, 
         cex_category = 0.75 )

treeplot(pairwise_termsim(HuSplCD103pDom), showCategory = 30, color = "p.adjust", # 9.65x5 inches
         offset.params = list(bar_tree = rel(5), extend=0.5, hexpand = 0.15), hilight.params = list(hilight = T),
         fontsize = 4, label_format = 20, cluster.params = list(label_format = 10), label_format_tiplab = 100, 
         cex_category = 0.75 )

gseaplot2(HuSplCD103nCD103p_GSEA, 
          geneSetID = c("GO:0006909", "GO:0002253", "GO:0001819",
                        "GO:0006260", "GO:0007059", "GO:0044839"),
          subplots = 1:2, pvalue_table = F, base_size = 14)

#save GSEA results
write.xlsx(list("HuSplCD103nDom"=HuSplCD103nDom@result, 
                "HuSplCD103pDom"=HuSplCD103pDom@result), 
           "PMID_31668803_human/HuSplCD103nCD103p_GSEA.xlsx")
