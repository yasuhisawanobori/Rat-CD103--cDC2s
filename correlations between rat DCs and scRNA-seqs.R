library(dplyr)
library(ggplot2)
library(readr)
library(limma)
library("pheatmap")
library(orthogene)
library(stringr)

#extract expression data and mix with rat microarray data for Pearson correlation----
#extract mouse thymic tdTomato-incorpolating DC expression data
integrated_TOMATO_DCs <- readRDS("PMID_32398640 Thymic monocyte-derived DCs/integrated_TOMATO_DCs.rds")
TOMATO_DCs_expression <- log2(as.data.frame(
  AggregateExpression(integrated_TOMATO_DCs, normalization.method = "RC", scale.factor = 1000000)))

TOMATO_DCs_expression$rat_gene <- convert_orthologs(TOMATO_DCs_expression, 
                                                input_species = "mouse", output_species = "rat",
                                                drop_nonorths = F, gene_output = "columns",
                                                non121_strategy = 5)[,2]

TOMATO_DCs_expression <- TOMATO_DCs_expression[!duplicated(TOMATO_DCs_expression$rat_gene),]

rownames(TOMATO_DCs_expression)<-TOMATO_DCs_expression$rat_gene

TOMATO_DCs_expression <- subset(TOMATO_DCs_expression, select= -c(rat_gene))

#extract mouse thymic DC expression data
Thymic_Myeloid <- readRDS(file = "PMID_35637352 mouse thymic CD11bCD11c/Thymic_Myeloid_Analysis.rds")
MoThymicMyeloid_expression <- log2(as.data.frame(
  AggregateExpression(Thymic_Myeloid, normalization.method = "RC", scale.factor = 1000000)))

MoThymicMyeloid_expression$rat_gene <- convert_orthologs(MoThymicMyeloid_expression, 
                                                    input_species = "mouse", output_species = "rat",
                                                    drop_nonorths = F, gene_output = "columns",
                                                    non121_strategy = 5)[,2]

MoThymicMyeloid_expression <- MoThymicMyeloid_expression[!duplicated(MoThymicMyeloid_expression$rat_gene),]

rownames(MoThymicMyeloid_expression)<-MoThymicMyeloid_expression$rat_gene

MoThymicMyeloid_expression <- subset(MoThymicMyeloid_expression, select= -c(rat_gene))

#extract mouse splenic DC expression data
Splenic_object <- readRDS(file = "PMID_33997687 mouse spleen CD11c+ scRNA-seq/splenic_DCs.rds")
MoSplDC_expression <- log2(as.data.frame(
  AggregateExpression(Splenic_object, normalization.method = "RC", scale.factor = 1000000)))

MoSplDC_expression$rat_gene <- convert_orthologs(MoSplDC_expression, 
                                                         input_species = "mouse", output_species = "rat",
                                                         drop_nonorths = F, gene_output = "columns",
                                                         non121_strategy = 5)[,2]

MoSplDC_expression <- MoSplDC_expression[!duplicated(MoSplDC_expression$rat_gene),]

rownames(MoSplDC_expression)<-MoSplDC_expression$rat_gene

MoSplDC_expression <- subset(MoSplDC_expression, select= -c(rat_gene))

#extract mouse lung DC expression data
MouseLungDCs <- readRDS("PMID_32392463 mouse lung MC and inf-cDC2/MouseLungDCs.rds")
MouseLungDCs_expression <- log2(as.data.frame(
  AggregateExpression(MouseLungDCs, normalization.method = "RC", scale.factor = 1000000)))

MouseLungDCs_expression$rat_gene <- convert_orthologs(MouseLungDCs_expression, 
                                                 input_species = "mouse", output_species = "rat",
                                                 drop_nonorths = F, gene_output = "columns",
                                                 non121_strategy = 5)[,2]

MouseLungDCs_expression <- MouseLungDCs_expression[!duplicated(MouseLungDCs_expression$rat_gene),]

rownames(MouseLungDCs_expression)<-MouseLungDCs_expression$rat_gene

MouseLungDCs_expression <- subset(MouseLungDCs_expression, select= -c(rat_gene))

#extract human splenic DC expression data
human_splenic_DC <- readRDS(file = "PMID_31668803_human/human_splenic_DC.rds")
HuSplDC_expression <- log2(as.data.frame(
  AggregateExpression(human_splenic_DC, normalization.method = "RC", scale.factor = 1000000)))

HuSplDC_expression$rat_gene <- convert_orthologs(HuSplDC_expression, 
                                                      input_species = "human", output_species = "rat",
                                                      drop_nonorths = F, gene_output = "columns",
                                                      non121_strategy = 5)[,2]

HuSplDC_expression <- HuSplDC_expression[!duplicated(HuSplDC_expression$rat_gene),]

rownames(HuSplDC_expression)<-HuSplDC_expression$rat_gene

HuSplDC_expression <- subset(HuSplDC_expression, select= -c(rat_gene))

#read rat DC/Mf data
RatDC_Mf_array <- as.data.frame(read_csv("PMID_34777358_35637352 mouse cDC1 cDC2 correlation/rat_DC_Mf_array-TACnormalized.csv"))

gene_symbols <- RatDC_Mf_array$`Gene Symbol` 

dup_indices <- duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE) # Find duplicated values

count <- ave(seq_along(gene_symbols), gene_symbols, FUN = seq_along) # Count the number of occurrences of each value

gene_symbols[dup_indices] <- paste0(gene_symbols[dup_indices], "-", count[dup_indices] - 1) # Append unique numbers to duplicates

gene_symbols <- gsub("-0", "", gene_symbols) #remove "-0"

row.names(RatDC_Mf_array) <- gene_symbols #Set gene symbols to row names

RatDC_Mf_array <- RatDC_Mf_array[, -(1:3)] #remove unnecessary columns

#calculate averages
RatDC_Mf_array$RatThymicDC1_average <- log2(rowMeans(2^RatDC_Mf_array[, 1:2], na.rm = T))
RatDC_Mf_array$RatSplenicDC1_average <- log2(rowMeans(2^RatDC_Mf_array[, 3:4], na.rm = T))
RatDC_Mf_array$RatThymicDC2_average <- log2(rowMeans(2^RatDC_Mf_array[, 5:6], na.rm = T))
RatDC_Mf_array$RatSplenicDC2_average <- log2(rowMeans(2^RatDC_Mf_array[, 7:8], na.rm = T))
RatDC_Mf_array$RatThymicCD103neg_average <- log2(rowMeans(2^RatDC_Mf_array[, 9:10], na.rm = T))
RatDC_Mf_array$RatSplenicCD103neg_average <- log2(rowMeans(2^RatDC_Mf_array[, 11:12], na.rm = T))
RatDC_Mf_array$RatThymicMf_average <- log2(rowMeans(2^RatDC_Mf_array[, 13:14], na.rm = T))
RatDC_Mf_array$RatSplenicMf_average <- log2(rowMeans(2^RatDC_Mf_array[, 15:16], na.rm = T))
RatDC_Mf_array<-RatDC_Mf_array[,-(1:16)]

#Merge the objects by genes (row.names)
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp <- merge(RatDC_Mf_array, TOMATO_DCs_expression, by="row.names", all=T)
row.names(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp) <- RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp$Row.names
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp <-RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp[,-1]

RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp <- merge(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp,
                                                           MoThymicMyeloid_expression, by="row.names", all=T)
row.names(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp) <- RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp$Row.names
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp <-RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp[,-1]

RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp <- merge(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp,
                                                           MoSplDC_expression, by="row.names", all=T)
row.names(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp) <- RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp$Row.names
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp <-RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp[,-1]

RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp <- merge(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp,
                                                           MouseLungDCs_expression, by="row.names", all=T)
row.names(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp) <- RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp$Row.names
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp <-RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp[,-1]

RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp <- merge(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp,
                                                           HuSplDC_expression, by="row.names", all=T)
row.names(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp) <- RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp$Row.names
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp <-RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp[,-1]

#remove NaNs and -Infs
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp[RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp == "NaN" | RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp == "-Inf"] <- "" 

#change class from matrix to data.frame
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp <- as.data.frame(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp) 

#change class of the contents from character to numeric
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp[,1:58] <- lapply(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp[,1:58] , as.numeric) 

#boxplot function
BP_distribution <- function(data, title = "Distribution of expression values", 
                            xlab = "", ylab = "removeBatchEffect-ed Log2 expression") {
  par(mar = c(13, 4, 4, 2) + 0.1)
 # Create the boxplot
  boxplot(data,
          main = title,
          xlab = xlab,
          ylab = ylab,
          las = 2,        # Rotate x-axis labels for better readability
          notch = FALSE, 
          cex.main = 1.5, # Increase the font size of the title
          cex.lab = 1.3)  # Increase the font size of the axis labels
}

BP_distribution(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp)

#removeBatchEffect
batch <- c(rep("RatArray", 8), rep("tdTomato", 5), rep("MoThy", 10), rep("MoSpl", 20), rep("MoLung", 8), rep("HuSpl", 7))
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_exp <- removeBatchEffect(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_exp, batch)

BP_distribution(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_exp)

#removeBatchEffect+standardize
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp <- scale(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_exp)

BP_distribution(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp)

#calculate SDs
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp <-as.data.frame(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp)

sd.p=function(x){sd(x, na.rm = TRUE)*sqrt((length(x)-1)/length(x))}
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp$RatThyCD103npSD <- apply(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[, c(3, 5)], 1, sd.p)
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp$RatSplCD103npSD <- apply(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[, c(4, 6)], 1, sd.p)
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp$RatThySplCD103npSD <- apply(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[, c(3, 4, 5, 6)], 1, sd.p)
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp$HumanDC2SD <- apply(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[, c(55:58)], 1, sd.p)

#export the expression matrix
write.csv(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp, "PMID_32398640 Thymic monocyte-derived DCs/RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp.csv", na ="", row.names = TRUE)

RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp <- read.csv("PMID_32398640 Thymic monocyte-derived DCs/RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp.csv", row.names = 1)

#rat cDCs vs mouse scRNA-seq  heatmaps (rBE+std version)----
#order the matrix
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp <- 
  RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[order(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp$RatThyCD103npSD, decreasing = T),]

#MoThyTomato
MoThyTomatoCorrelation <- data.frame(matrix(nrow=5, ncol=2))
rownames(MoThyTomatoCorrelation) <- c("cDC2", "moDC", "cDC1a", "cDC1b" , "pDC")
colnames(MoThyTomatoCorrelation) <- c("RatThyCD103-", "RatThyCD103+")

row_vars <- c("RNA.cDC2.x", "RNA.moDC",	"RNA.cDC1a", "RNA.cDC1b",	"RNA.pDC.x")
col_vars <- c("RatThymicCD103neg_average", "RatThymicDC2_average")

for (i in seq_len(nrow(MoThyTomatoCorrelation))) {
  for (j in seq_len(ncol(MoThyTomatoCorrelation))) {
    MoThyTomatoCorrelation[i, j] <- (cor.test(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, col_vars[j]], 
                                              RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, row_vars[i]]))$estimate
  }
}

scaled_MoThyTomatoCorrelation <- MoThyTomatoCorrelation
scaled_MoThyTomatoCorrelation <- t(apply(MoThyTomatoCorrelation, 1, function(x) scale(as.numeric(x))))

pheatmap(t(scaled_MoThyTomatoCorrelation),
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F, 
         display_numbers = t(apply(MoThyTomatoCorrelation, 2, function(x) sprintf("%.2f", x))))

pheatmap(t(MoThyTomatoCorrelation), display_numbers = T,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F) #Fig. 5C

pheatmap(t(apply(MoThyTomatoCorrelation[, c(1,2)], 1, sd.p)), display_numbers = T,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F)

#MoThy
MoThyCorrelation <- data.frame(matrix(nrow=7,ncol=2))
rownames(MoThyCorrelation) <- c("cDC2", "mDC2", "cDC1",  "mDC1", "pDC", "Mono/Mac", "Granulocyte")
colnames(MoThyCorrelation) <- c("RatThyCD103-", "RatThyCD103+")

row_vars <- c("RNA.cDC2.y", "RNA.mDC2",	"RNA.cDC1.x",	"RNA.mDC1",	"RNA.pDC.y",	"RNA.Mono.Mac",	"RNA.Granulocytes")
col_vars <- c("RatThymicCD103neg_average", "RatThymicDC2_average")

for (i in seq_len(nrow(MoThyCorrelation))) {
  for (j in seq_len(ncol(MoThyCorrelation))) {
    MoThyCorrelation[i, j] <- (cor.test(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, col_vars[j]], 
                                        RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, row_vars[i]]))$estimate
  }
}

scaled_MoThyCorrelation <- MoThyCorrelation
scaled_MoThyCorrelation <- t(apply(MoThyCorrelation, 1, function(x) scale(as.numeric(x))))

pheatmap(t(scaled_MoThyCorrelation),
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F, 
         display_numbers = apply(t(MoThyCorrelation), 2, function(x) sprintf("%.2f", x)))

pheatmap(t(MoThyCorrelation), display_numbers = T,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F)


pheatmap(t(apply(MoThyCorrelation[, c(1,2)], 1, sd.p)), display_numbers = T,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F)


#MoSpl
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp <- 
  RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[order(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp$RatSplCD103npSD, decreasing = T),]

MoSplCorrelation <- data.frame(matrix(nrow=17,ncol=2))
rownames(MoSplCorrelation) <- c("RNA.mitochondrial.cDC2",	"RNA.mature.Il4i1..cDC2", 
                                "RNA.immature.phagocytic.cDC2a", "RNA.mitotic.cDC2a.1", "RNA.mitotic.cDC2a.2",	"RNA.mitotic.cDC2a.3", "RNA.mature.cDC2a",
                                "RNA.immature.phagocytic.cDC2b",	"RNA.mature.cDC2b", 
                                "RNA.immature.phagocytic.cDC1",	"RNA.ribosomal.cDC1",	"RNA.immature.mitotic.cDC1",	"RNA.mature.cDC1",	"RNA.mitochondrial.cDC1",	"RNA.immature.phagocytic.Sox4..cDC", 
                                "RNA.pDC.x.1", "RNA.monocytes")
colnames(MoSplCorrelation) <- c("RatSplCD103-", "RatSplCD103+")

row_vars <- c("RNA.mitochondrial.cDC2",	"RNA.mature.Il4i1..cDC2", 
              "RNA.immature.cDC2a", "RNA.mitotic.cDC2a.1", "RNA.mitotic.cDC2a.2",	"RNA.mitotic.cDC2a.3", "RNA.mature.cDC2a",
              "RNA.immature.cDC2b",	"RNA.mature.cDC2b", 
              "RNA.immature.cDC1",	"RNA.ribosomal.cDC1",	"RNA.mitotic.cDC1",	"RNA.mature.cDC1",	"RNA.mitochondrial.cDC1",	"RNA.immature.Sox4..cDC", 
              "RNA.pDC.x.1", "RNA.monocytes")
col_vars <- c("RatSplenicCD103neg_average", "RatSplenicDC2_average")

for (i in seq_len(nrow(MoSplCorrelation))) {
  for (j in seq_len(ncol(MoSplCorrelation))) {
    MoSplCorrelation[i, j] <- (cor.test(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, col_vars[j]], 
                                        RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, row_vars[i]]))$estimate
  }
}

scaled_MoSplCorrelation <- MoSplCorrelation
scaled_MoSplCorrelation <- t(apply(MoSplCorrelation, 1, function(x) scale(as.numeric(x))))

pheatmap(t(scaled_MoSplCorrelation),
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F, 
         display_numbers = apply(t(MoSplCorrelation), 2, function(x) sprintf("%.2f", x)))

pheatmap(t(MoSplCorrelation), display_numbers = T,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F) #Fig. 5A

pheatmap(t(apply(MoSplCorrelation[, c(1,2)], 1, sd.p)), display_numbers = T,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F)

#MoLung
MoLungCorrelation <- data.frame(matrix(nrow=8,ncol=2))
rownames(MoLungCorrelation) <- c("RNA.Migratory.cDC2", "RNA.Non.migratory.cDC2",	"RNA.inf.cDC2",	"RNA.MC",
                                 "RNA.Migratory.cDC1",	"RNA.Non.migratory.cDC1", "RNA.pDC.y.1", "RNA.proliferating.DC")
colnames(MoLungCorrelation) <- c("RatSplCD103-", "RatSplCD103+")

row_vars <- c("RNA.Migratory.cDC2", "RNA.Non.migratory.cDC2",	"RNA.inf.cDC2",	"RNA.MC",
              "RNA.Migratory.cDC1",	"RNA.Non.migratory.cDC1", "RNA.pDC.y.1", "RNA.proliferating.DC")
col_vars <- c("RatSplenicCD103neg_average", "RatSplenicDC2_average")

for (i in seq_len(nrow(MoLungCorrelation))) {
  for (j in seq_len(ncol(MoLungCorrelation))) {
    MoLungCorrelation[i, j] <- (cor.test(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, col_vars[j]], 
                                        RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, row_vars[i]]))$estimate
  }
}

scaled_MoLungCorrelation <- MoLungCorrelation
scaled_MoLungCorrelation <- t(apply(MoLungCorrelation, 1, function(x) scale(as.numeric(x))))

pheatmap(t(scaled_MoLungCorrelation),
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F, 
         display_numbers = t(apply(MoLungCorrelation, 2, function(x) sprintf("%.2f", x))))

pheatmap(t(MoLungCorrelation), display_numbers = T,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F) # Fig. 5B


pheatmap(t(apply(MoLungCorrelation[, c(1,2)], 1, sd.p)), display_numbers = T,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F)

#HuSpl
HuSplCorrelation <- data.frame(matrix(nrow=7,ncol=2))
rownames(HuSplCorrelation) <- c("RNA.CLEC10A..cDC2",	"RNA.CLEC10A..cDC2.1",	"RNA.CCR7..cDC2",	"RNA.Mitotic.cDC2",	"RNA.AS.DC", "RNA.cDC1.y",	"RNA.Mitotic.cDC1")
colnames(HuSplCorrelation) <- c("RatSplCD103-", "RatSplCD103+")

row_vars <- c("RNA.CLEC10A..cDC2",	"RNA.CLEC10A..cDC2.1",	"RNA.CCR7..cDC2",	"RNA.Mitotic.cDC2",	"RNA.AS.DC", "RNA.cDC1.y",	"RNA.Mitotic.cDC1")
col_vars <- c("RatSplenicCD103neg_average", "RatSplenicDC2_average")

for (i in seq_len(nrow(HuSplCorrelation))) {
  for (j in seq_len(ncol(HuSplCorrelation))) {
    HuSplCorrelation[i, j] <- (cor.test(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, col_vars[j]], 
                                        RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, row_vars[i]]))$estimate
  }
}

scaled_HuSplCorrelation <- HuSplCorrelation
scaled_HuSplCorrelation <- t(apply(HuSplCorrelation, 1, function(x) scale(as.numeric(x))))

pheatmap(scaled_HuSplCorrelation,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F, 
         display_numbers = apply(HuSplCorrelation, 2, function(x) sprintf("%.2f", x)))

pheatmap(t(HuSplCorrelation), display_numbers = T,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F) #Fig. 6A

pheatmap(t(apply(HuSplCorrelation[, c(1,2)], 1, sd.p)), display_numbers = T,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F)

#mouse spleen vs HuSpl 
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp <- RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[order(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp$HumanDC2SD, decreasing = T),]

HuSplmoSplCorrelation <- data.frame(matrix(nrow=5,ncol=11))
rownames(HuSplmoSplCorrelation) <- c("RNA.CLEC10A..cDC2",	"RNA.CLEC10A..cDC2.1",	"RNA.CCR7..cDC2",	"RNA.Mitotic.cDC2",	"RNA.AS.DC")
colnames(HuSplmoSplCorrelation) <- c("RNA.mitochondrial.cDC2",	"RNA.mature.Il4i1..cDC2", 
                                     "RNA.immature.cDC2a", "RNA.mitotic.cDC2a.1", "RNA.mitotic.cDC2a.2",	"RNA.mitotic.cDC2a.3", "RNA.mature.cDC2a",
                                     "RNA.immature.cDC2b",	"RNA.mature.cDC2b", 
                                     "RNA.immature.Sox4..cDC", "RNA.monocytes")

row_vars <- c("RNA.CLEC10A..cDC2",	"RNA.CLEC10A..cDC2.1",	"RNA.CCR7..cDC2",	"RNA.Mitotic.cDC2",	"RNA.AS.DC")
col_vars <- c("RNA.mitochondrial.cDC2",	"RNA.mature.Il4i1..cDC2", 
              "RNA.immature.cDC2a", "RNA.mitotic.cDC2a.1", "RNA.mitotic.cDC2a.2",	"RNA.mitotic.cDC2a.3", "RNA.mature.cDC2a",
              "RNA.immature.cDC2b",	"RNA.mature.cDC2b", 
              "RNA.immature.Sox4..cDC", "RNA.monocytes")

for (i in seq_len(nrow(HuSplmoSplCorrelation))) {
  for (j in seq_len(ncol(HuSplmoSplCorrelation))) {
    HuSplmoSplCorrelation[i, j] <- (cor.test(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, col_vars[j]], 
                                        RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, row_vars[i]]))$estimate
  }
}

scaled_HuSplmoSplCorrelation <- HuSplmoSplCorrelation
scaled_HuSplmoSplCorrelation <- t(apply(HuSplmoSplCorrelation, 1, function(x) scale(as.numeric(x))))

pheatmap(scaled_HuSplmoSplCorrelation,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F, 
         display_numbers = apply(HuSplmoSplCorrelation, 2, function(x) sprintf("%.2f", x)))

pheatmap(HuSplmoSplCorrelation, display_numbers = T,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F)


#Spleen cDC2b vs Spleen Monocytes vs Lung inf cDC2 vs MC vs thymus moDC----
#extract mouse simplified splenic DC expression data
Splenic_object_simple <- readRDS(file = "PMID_33997687 mouse spleen CD11c+ scRNA-seq/Splenic_object_simple.rds")
MoSplDC_simple_expression <- log2(as.data.frame(
  AggregateExpression(Splenic_object_simple, normalization.method = "RC", scale.factor = 1000000)))

MoSplDC_simple_expression$rat_gene <- convert_orthologs(MoSplDC_simple_expression, 
                                                      input_species = "mouse", output_species = "rat",
                                                      drop_nonorths = F, gene_output = "columns",
                                                      non121_strategy = 5)[,2]

MoSplDC_simple_expression <- MoSplDC_simple_expression[!duplicated(MoSplDC_simple_expression$rat_gene),]

rownames(MoSplDC_simple_expression)<-MoSplDC_simple_expression$rat_gene

MoSplDC_simple_expression <- subset(MoSplDC_simple_expression, select= -c(rat_gene))

#Merge the objects by genes (row.names)
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp <- merge(RatDC_Mf_array, TOMATO_DCs_expression, by="row.names", all=T)
row.names(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp) <- RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp$Row.names
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp <-RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp[,-1]

RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp <- merge(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp,
                                                           MoThymicMyeloid_expression, by="row.names", all=T)
row.names(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp) <- RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp$Row.names
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp <-RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp[,-1]

RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp <- merge(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp,
                                                           MoSplDC_simple_expression, by="row.names", all=T)
row.names(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp) <- RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp$Row.names
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp <-RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp[,-1]

RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp <- merge(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp,
                                                           MouseLungDCs_expression, by="row.names", all=T)
row.names(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp) <- RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp$Row.names
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp <-RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp[,-1]

RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp <- merge(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp,
                                                           HuSplDC_expression, by="row.names", all=T)
row.names(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp) <- RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp$Row.names
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp <-RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp[,-1]

#remove NaNs and -Infs
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp[RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp == "NaN" | RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp == "-Inf"] <- "" 

#change class from matrix to data.frame
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp <- as.data.frame(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp) 

#change class of the contents from character to numeric
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp[,1:53] <- lapply(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp[,1:53] , as.numeric) 

#boxplot function
BP_distribution <- function(data, title = "Distribution of expression values", 
                            xlab = "", ylab = "removeBatchEffect-ed Log2 expression") {
  par(mar = c(13, 4, 4, 2) + 0.1)
  # Create the boxplot
  boxplot(data,
          main = title,
          xlab = xlab,
          ylab = ylab,
          las = 2,        # Rotate x-axis labels for better readability
          notch = FALSE, 
          cex.main = 1.5, # Increase the font size of the title
          cex.lab = 1.3)  # Increase the font size of the axis labels
}

BP_distribution(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp)

#removeBatchEffect
batch <- c(rep("RatArray", 8), rep("tdTomato", 5), rep("MoThy", 10), rep("MoSpl", 15), rep("MoLung", 8), rep("HuSpl", 7))
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_exp <- removeBatchEffect(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_exp, batch)

BP_distribution(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_exp)

#removeBatchEffect+standardize
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp <- scale(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_exp)

BP_distribution(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp)

#calculate SDs
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp <-as.data.frame(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp)

sd.p=function(x){sd(x, na.rm = TRUE)*sqrt((length(x)-1)/length(x))}
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp$RatThyCD103npSD <- apply(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp[, c(3, 5)], 1, sd.p)
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp$RatSplCD103npSD <- apply(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp[, c(4, 6)], 1, sd.p)
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp$RatThySplCD103npSD <- apply(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp[, c(3, 4, 5, 6)], 1, sd.p)

write.csv(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp,
          "PMID_32398640 Thymic monocyte-derived DCs/RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp.csv", na ="", row.names = TRUE)
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp <- 
  read.csv("PMID_32398640 Thymic monocyte-derived DCs/RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp.csv", row.names = 1)

#order the matrix
RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp <- 
  RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp[order(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp$RatThySplCD103npSD, decreasing = T),]

CD103nSimilarSubsetCorrelation <- data.frame(matrix(nrow=6,ncol=6))
rownames(CD103nSimilarSubsetCorrelation) <- c("RNA.cDC2b", "RNA.monocytes", "RNA.Mono.Mac", "RNA.inf.cDC2", "RNA.MC", "RNA.moDC")
colnames(CD103nSimilarSubsetCorrelation) <- c("RNA.cDC2b", "RNA.monocytes", "RNA.Mono.Mac", "RNA.inf.cDC2", "RNA.MC", "RNA.moDC")

row_vars <- c("RNA.cDC2b", "RNA.monocytes", "RNA.Mono.Mac", "RNA.inf.cDC2", "RNA.MC", "RNA.moDC")
col_vars <- c("RNA.cDC2b", "RNA.monocytes", "RNA.Mono.Mac", "RNA.inf.cDC2", "RNA.MC", "RNA.moDC")

for (i in seq_len(nrow(CD103nSimilarSubsetCorrelation))) {
  for (j in seq_len(ncol(CD103nSimilarSubsetCorrelation))) {
    CD103nSimilarSubsetCorrelation[i, j] <- (cor.test(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp[1:1000, col_vars[j]], 
                                        RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp[1:1000, row_vars[i]]))$estimate
  }
}

pheatmap(CD103nSimilarSubsetCorrelation, display_numbers = T,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F)
#the heatmap of cDC2/MC markers----
#curated from PMID: 36642930
pheatmap(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp
         [c("Flt3", "Dpp4", "Ly6i", "RGD1565410; Ly6i", "Adgre1", "Ikzf1", "Gfi1", "Spi1", "Cd209a", "Clec4a1", "Cd14", "Cx3cr1",
            "Fcer1a", "Fcgr1a", "Fcgr1", "Fcgr3a", "Esam", "Clec10a", "Clec12a", "Thbd", "Cd163", "Fcgr2b", "Fcgr3a", "C5ar1",
            "Csf1r", "Mrc1", "Mertk", "Mafb", "Irf4", "Irf8", "Stat4", "Cd86", "Il12b"), 
           c("RatSplenicCD103neg_average", "RatSplenicDC2_average", "RatThymicCD103neg_average", "RatThymicDC2_average", 
             "RNA.cDC2b","RNA.monocytes", "RNA.cDC2a", 
             "RNA.Mono.Mac", "RNA.cDC2.y", "RNA.mDC2", 
             "RNA.inf.cDC2", "RNA.MC", "RNA.Non.migratory.cDC2", "RNA.Migratory.cDC2",
             "RNA.moDC", "RNA.cDC2.x"
             #"RNA.CLEC10A..cDC2",	"RNA.CCR7..cDC2",	"RNA.CLEC10A..cDC2.1", "RNA.Mitotic.cDC2"
             )],
         cellwidth = 25, cellheight = 15, cluster_rows = F, cluster_cols = F, display_numbers = T, fontsize_row = 12, scale = "row")

#curated from PMID: 32392463
pheatmap(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp
         [c("Il12b", "Cd86", "Irf8", "Stat4", "Traf1", "Cxcr5", "Il12rb1", "Il12rb2", "Batf2", "Oas2", "Oas1a", "Tlr7", "Fcgr1a",
            "Ifit1", "Ifit3", "Mnda", 
            "Sct", "Hamp", "LOC100910979", "Apod", "Gpr33-ps1", "B3gnt5", "Amigo3", "Nt5c3a",
            "Ly6i", "Rsad2", "Mnda", "Cxcl9", "Cxcl10", "Il1rn", "Fgl2", "Ifit3", "Ifit1",
            
            "Mertk", "Mafb", "Csf1r", "Sirpa", "Clec4a3",
            "Adgre1", "Mafb", "Ccr2", "Csf1r", "Cx3cr1", "Cebpb"), 
           c("RatSplenicCD103neg_average", "RatSplenicDC2_average", "RatThymicCD103neg_average", "RatThymicDC2_average", 
             "RNA.cDC2b","RNA.monocytes", "RNA.cDC2a", 
             "RNA.Mono.Mac", "RNA.cDC2.y", "RNA.mDC2", 
             "RNA.inf.cDC2", "RNA.MC", "RNA.Non.migratory.cDC2", "RNA.Migratory.cDC2",
             "RNA.moDC", "RNA.cDC2.x"
             #"RNA.CLEC10A..cDC2",	"RNA.CCR7..cDC2",	"RNA.CLEC10A..cDC2.1", "RNA.Mitotic.cDC2"
             )],
         cellwidth = 25, cellheight = 12, cluster_rows = T, cluster_cols = F, display_numbers = T, fontsize_row = 12, scale = "row")

#curated from PMID: 32392463, simplified InfcDC2
pheatmap(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp
         [c("Cd86", "Irf8", "Stat4", "Batf2", "Ifit1", "Phf11b", "Ifit3", "Mnda", 
            "Sct", "Hamp", "Apod", "B3gnt5", "Amigo3", "Nt5c3a",
            "Rsad2", "Cxcl9", "Cxcl10", "Fgl2"), 
           c("RNA.inf.cDC2", "RNA.MC", 
             "RatSplenicCD103neg_average", "RatSplenicDC2_average", "RatThymicCD103neg_average", "RatThymicDC2_average",
             "RNA.cDC2b","RNA.monocytes", "RNA.cDC2a", 
             "RNA.Mono.Mac", "RNA.cDC2.y", "RNA.mDC2", 
             "RNA.Non.migratory.cDC2", "RNA.Migratory.cDC2",
             "RNA.moDC", "RNA.cDC2.x"
             #"RNA.CLEC10A..cDC2",	"RNA.CCR7..cDC2",	"RNA.CLEC10A..cDC2.1", "RNA.Mitotic.cDC2"
             )],
         cellwidth = 25, cellheight = 15, cluster_rows = T, cluster_cols = F, display_numbers = T, fontsize_row = 12, scale = "row")

#curated from PMID: 32392463, simplified MC
pheatmap(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp
         [c( "Fcgr1","Mertk", "Mafb", "Csf1r", 
            "Adgre1", "Ccr2", "Cx3cr1", "Cebpb", "Cd14"), 
           c("RNA.inf.cDC2", "RNA.MC", 
             "RatSplenicCD103neg_average", "RatSplenicDC2_average", "RatThymicCD103neg_average", "RatThymicDC2_average"
             #"RNA.cDC2b","RNA.monocytes", "RNA.cDC2a", 
             #"RNA.Mono.Mac", "RNA.cDC2.y", "RNA.mDC2", 
             #"RNA.Non.migratory.cDC2", "RNA.Migratory.cDC2",
             #"RNA.moDC", "RNA.cDC2.x"
             #"RNA.CLEC10A..cDC2",	"RNA.CCR7..cDC2",	"RNA.CLEC10A..cDC2.1", "RNA.Mitotic.cDC2"
           )],
         cellwidth = 25, cellheight = 15, cluster_rows = T, cluster_cols = F, display_numbers = T, 
         fontsize_row = 12, scale = "row")

#curated from PMID: 32392463, markers expressed by both subsets
pheatmap(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp
         [c( "Fcgr1a", "Fcer1a", "Fcgr3a", "Fcgr2b"), 
           c("RNA.inf.cDC2", "RNA.MC", 
             "RatSplenicCD103neg_average", "RatSplenicDC2_average", "RatThymicCD103neg_average", "RatThymicDC2_average"
             #"RNA.cDC2b","RNA.monocytes", "RNA.cDC2a", 
             #"RNA.Mono.Mac", "RNA.cDC2.y", "RNA.mDC2", 
             #"RNA.Non.migratory.cDC2", "RNA.Migratory.cDC2",
             #"RNA.moDC", "RNA.cDC2.x"
             #"RNA.CLEC10A..cDC2",	"RNA.CCR7..cDC2",	"RNA.CLEC10A..cDC2.1", "RNA.Mitotic.cDC2"
           )],
         cellwidth = 25, cellheight = 15, cluster_rows = F, cluster_cols = F, display_numbers = T, 
         fontsize_row = 12, scale = "row")

#curated from PMID: 32392463, decreased-InfcDC2
pheatmap(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp
         [c("Cd86", "Irf8", "Stat4", "Batf2", "Ifit1", "Phf11b", "Ifit3", "Mnda", 
            "Sct", "Hamp", "Apod", "B3gnt5", "Amigo3", "Nt5c3a",
            "Rsad2", "Cxcl9", "Cxcl10", "Fgl2"), 
           c("RNA.inf.cDC2", "RNA.MC", 
             "RatSplenicCD103neg_average", "RatSplenicDC2_average", "RatThymicCD103neg_average", "RatThymicDC2_average",
             #"RNA.cDC2b","RNA.monocytes", "RNA.cDC2a", 
             #"RNA.Mono.Mac", "RNA.cDC2.y", "RNA.mDC2", 
             #"RNA.Non.migratory.cDC2", "RNA.Migratory.cDC2",
             #"RNA.moDC", "RNA.cDC2.x"
             "RNA.CLEC10A..cDC2",	"RNA.CCR7..cDC2",	"RNA.CLEC10A..cDC2.1", "RNA.Mitotic.cDC2"
           )],
         cellwidth = 25, cellheight = 15, cluster_rows = T, cluster_cols = F, display_numbers = T, fontsize_row = 12, scale = "row")

#curated from PMID: 32392463, decreased-MC
pheatmap(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp
         [c( "Fcgr1","Mertk", "Mafb", "Csf1r", 
             "Adgre1", "Ccr2", "Cx3cr1", "Cebpb", "Cd14"), 
           c("RNA.inf.cDC2", "RNA.MC", 
             "RatSplenicCD103neg_average", "RatSplenicDC2_average", "RatThymicCD103neg_average", "RatThymicDC2_average",
             #"RNA.cDC2b","RNA.monocytes", "RNA.cDC2a", 
             #"RNA.Mono.Mac", "RNA.cDC2.y", "RNA.mDC2", 
             #"RNA.Non.migratory.cDC2", "RNA.Migratory.cDC2",
             #"RNA.moDC", "RNA.cDC2.x"
             "RNA.CLEC10A..cDC2",	"RNA.CCR7..cDC2",	"RNA.CLEC10A..cDC2.1", "RNA.Mitotic.cDC2"
           )],
         cellwidth = 25, cellheight = 15, cluster_rows = T, cluster_cols = F, display_numbers = T, 
         fontsize_row = 12, scale = "row")


SplcDC2bMarkers <- c("Adgre1", "Apobec1", "Arl11", "Blvrb", "Ccdc92; Glul", "Ccl9", "Ccr2", "Cd209a", "Cd300c", "Cd302", "Cebpd", "Celf4", "Clec4a1", "Clec4a3", "Clec12a", "Csf1r", "Ctsc", "Ctss", "Cx3cr1", "Cybb", "Dusp6", "Ehd4", "Emb", "Fcer1g", "Fcgr2b", "Glul", "Hexb", "Hpgd", "Ifitm2", "Ifitm3", "Igsf6", "Il13ra1", "Krt80", "Lamp1", "Lgals3bp", "Lgmn", "LOC691141", "LOC100911469; Ms4a6b", "Lrp1", "Lst1", "Lyl1", "Lyz2", "Mcemp1", "Ms4a6b", "Ms4a6b", "Ms4a6c", "Ncf2", "Nfe2l2", "Nfix", "Nxpe4", "Oasl2", "P2ry6", "Pla2g7", "Ptafr", "Ptpro", "Rnase2", "Sdc3", "Sirpa", "Slfn5", "Tgfbi", "Tmem176a", "Tmem176b", "Tnfrsf21", "Zeb2")
SplMonoMarkers <- c("Abcc5", "Abcd2", "Abhd12", "Abi3", "Acaa2", "Ace3; Ace", "Acer3", "Acot9", "Adam10", "Add3", "Adgre1", "Adgre5", "Adipor1", "Adrb2", "Agpat4", "Ahnak", "Akt1", "Aldh2", "Aldh3a2", "Aldh3b1", "Aldoa", "Anxa1", "Anxa2", "Ap1s2", "Api5", "Apobec1", "Apoe", "Aprt", "Arap1", "Arhgap27", "Arhgap39", "Arhgap39", "Arhgef1", "Arhgef3", "Arhgef3", "Arhgef10l", "Arid3a", "Arrb1", "Asph; LOC100910356; ENSRNOG00000007445", "Atg101", "Atp1a1", "Atp1a3", "Atp1b3", "Atp2b1", "Atp6ap1", "Atp6v1b2", "Atp7a", "Atp8b4", "Atp13a2", "Aup1", "Azi2", "B4galt5", "Bach1", "Bach1", "Bbc3", "Bcl2", "Bcl6", "Blvrb", "Bmpr2", "Bnip3l", "C3", "Camk2d", "Camkk2", "Capn2", "Cast", "Cblb", "Cbx4", "Ccdc88c", "Ccl9", "Ccm2", "Ccr2", "Cd44", "Cd300c", "Cd300e", "Cd300lb", "Cd300lf", "Cd302", "Cdk2ap2", "Cdkn1a", "Cdkn2d", "Ceacam1", "Cebpb", "Cebpd", "Cela1", "Cept1", "Cers5", "Chmp2b", "Chst12", "Clec4a1", "Clec4a3", "Clec12a", "Cmip", "Cmtm7", "Cnppd1", "Cpd", "Cpt1a", "Crk", "Csf1r", "Csf3r", "Ctsb", "Ctsd", "Cx3cr1", "Cxcr4", "Cybb", "Cyfip2", "Cyp4f18", "Cyth3", "Ddx3x", "Degs1", "Dgkd", "Dgkg", "Dgkh", "Dgkz", "Diaph1", "Dok3", "Dram2", "Dstyk", "Dusp3", "Dusp6", "Dusp16", "Ebi3", "Ehd4", "Eif4ebp1", "Elf4", "Elmo2", "Emb", "Emilin2", "Emp3", "Eno3", "Ethe1", "Fabp4", "Fam168a", "Fbxo28", "Fcer1g", "Fcgr2b", "Fcgr3a", "Fcho2", "Fgd4", "Fgr", "Flna", "Fmnl1", "Fndc3a", "Fosl2", "Frmd8", "Fry", "Furin", "Fyn", "G0s2", "G6pd", "Gcnt2", "Gcnt2", "Gcnt2", "Gda", "Ggh", "Ggta1", "Gk", "Gpcpd1", "Gpcpd1", "Gpr141", "Gpx1", "Gramd4", "Grina", "Grk6", "Gsr", "Gyg1", "Hacd4", "Hdac5", "Heg1", "Herc4", "Hes1", "Hip1", "Hipk1", "Hmox1", "Hp", "Hpcal1", "Hpgd", "Hsd17b11", "Idh1", "Ier2", "Ifitm2", "Ifitm3", "Ifitm6", "Igsf6", "Il6r", "Il13ra1", "Il17ra", "Inpp5f", "Inpp5f", "Itgal", "Itgam", "Itgav", "Itgb1", "Itgb2", "Klf2", "Klf3", "Klf4", "Klf10", "Klf13", "Klra2", "Klra2", "Krt80", "L1cam", "Lair1", "Lamp1", "Lamp2", 
                    "Laptm4a", "Lats2", "Ldlr", "Ldlrad3", "Ldlrad3; LOC100909436", "Lgals3", "Lgmn", "Lilra5", "Litaf", "Lmo4", "Lnpep", "Lnpep", "LOC500331", "LOC500331", "LOC691141", "LOC100910497", "LOC100911469; Ms4a6b", "LOC102548265; Cmip", "LOC102548345; Nedd9", "LOC102548350; Cmip", "LOC103690035; Tgfbr1", "Lpcat2", "Lrp1", "Lrrc8d", "Lrrfip1", "Lst1", "Ly6e", "Lyl1", "Lyz2", "Mafb", "Mafg", "Man2a1", "Map3k11", "Mapk3", "Mapkapk2", 
                    "Mapkapk3", "Mbp", "Mef2a", "Metrnl", "Mettl7a", "Mettl9; DREV.0", "Mgst1", "Mier1", "Mnda", "Ms4a6b", "Ms4a6b", "Ms4a6c", "Msn", "Msrb1", "Mxi1", "Myo1c", "Myo1f", "Myo1g", "Myo18a", "N4bp1", "Nab1", "Nadk", "Naip6", "Ncf2", "Ncoa1", "Ncor2", "Ncor2", "Ncstn", "Ndel1", "Nedd9", "Neu1", "Nfam1", "Nfe2l2", "Nfic", "Nfix", "Nfkbiz", "Nhsl2", "Nin", "Ninj1", "Nkiras2", "Nod1", "Notch2", "Nr4a1", "Nxpe4", "Orai1", "Oxr1", "P4ha1", "Pag1", "Pbxip1", "Pdlim5", "Pfkfb4", "Pglyrp1", "Phyh", "Pisd", "Pitpnm1", "Pla2g7", "Plac8", "Plagl2", "Plaur", "Plcg2", "Plec", "Plin2", "Plod3; Mir702", "Pltp", "Plxnb2; Mir349", "Pon2", "Pot1b", "Pou2f2", "Pparg", "Ppm1h", "Ppp1cb", "Ppp2r5a", "Ppp2r5c", "Prdx5", "Prkcd", "Prkch", "Prr13", "Ptger4", "Ptk2b", "Ptpn1", "Ptpn12", "Ptpre", "Ptprj", "Ptpro", "Pxk", "Pygl", "Raf1", "Ralb", "Rap1b", "Rap1gap2", "Rap2a", "Rap2c", "Rara", "Rasa3", "Rasgrp2", "Rassf5", "Rbms1", "Rbpms", "Rcbtb2", "Rftn1", "RGD1304884", "RGD1560687", "RGD1563917", "RGD1565410; Ly6i", "Rin3", "Rnase2", "Rnf149", "Rnf166", "Rras", "S1pr5", "Samsn1", "Scarb1", "Scarb2", "Sec14l1", "Sephs2", "Serinc3", "Sgk1", "Sgms1", "Sh2d1b2", "Sh2d3c", "Sh3bp2", "Sh3kbp1", "Sirpa", "Slc11a1", "Slc12a2", "Slc12a6", "Slc16a3", "Slc29a1", "Slfn2", "Smagp", "Smc6", "Smpdl3a", "Snrk", "Snx18", "Soat1", "Socs3", "Sorl1", "Sowahc", "Spata13", "Spen", "Sppl2a", "Ssh2", "Stat5b", "Stk10", "Stk38", "Stk40", "Stx11", "Susd3", "Svil", "Sypl1", "Tab2", "Tbc1d10b", "Tcf7l2", "Tgfbi", "Tgfbr1", "Tgfbr2", "Thbd", "Tiam2", "Tle4", "Tlr7", "Tm6sf1", "Tm9sf4", "Tmcc1", "Tmem38b", "Tmem50a", "Tmem51", "Tmem164", "Tmem259", "Tnfaip2", "Tnfrsf1a", "Tnfrsf1b", "Tnfrsf21", "Tpd52", "Tppp3", "Tpst2", "Trem3", "Treml4", "Trib1", "Trim8", "Tspan14", "Tubb6", "Tyrobp", "Ubash3b", "Ugp2", "Vps13c; ENSRNOG00000030213", "Wsb1", "Xbp1", "Xdh", "Xkr8; Smpdl3b", "Ypel3", "Zbp1", "Zbtb7b", "Zeb2", "Zfand3", "Zfyve9", "Znhit1", "Zyx")
LungInfDC2Markers <- c("Batf", "Batf2", "Ccl3", "Ccl4", "Ccnd2", "Cd40", "Cd69", "Cxcl9", "Cxcl10", "Daxx", "Ddt", "Ddt", "Dhx58", "Dnajc7", "Dnase1l3", "Dnase1l3", "Eif2ak2", "Fgl2", "Flnb", "Gbp2", "Gbp4", "Gbp5", "Gk", "Gmppb", "Grap2", "Herc6", "Ifi44", "Ifi47", "Ifih1", "Ifit1", "Ifit2", "Ifit3", "Ifitm1", "Ifitm3", "Il1r2", "Il12b", "Irf7", "Irgm", "Isg15", "Isg20", "Kctd14", "Kdr", "Klrb1a; RGD2301395; Klrb1b", "Klrk1", "Kmo", "Kynu", "LOC690000", "LOC100910979; LOC102555392", "LOC100911469; Ms4a6b", "LOC102554102; Trim30c", "Malt1", "Mnda", "Ms4a6a", "Ms4a6b", "Ms4a6b", "Mx1", "Nmi", "Nt5c3a", "Oas1a", "Oas3", "Oasl", "Oasl2", "Parp9", "Parp12", "Parp14", "Phf11b", "Phf11b; Phf11; LOC102551046", "Pml", "Pnp", "Ppa1", "Procr", "Pttg1", "RGD1565410; Ly6i", "Rnf213", "Rsad2", "RT1-DOa", "Rtp4", "Samd9l", "Scimp", "Slamf8", "Slfn1", "Slfn2", "Slfn5", "Slfn13", "Sp100", "Sp100", "Tor3a", "Usp18", "Zbp1")
LungMCMarkers <- c("Abcg1", "Abhd12", "Abi3", "Abr", "Acer3", "Acp2", "Acp5", "Acp5", "Adam17", "Adap2", "Adgre1", "Aga", "Aif1", "Alas1", "Alas1", "Aldh3b1", "Amdhd2", "Anpep", "Anxa1", "Aph1b", "Aph1b", "Aph1b", "Apobec1", "Apoe", "Aprt", "Arhgap25", "Arid5b", "Asph; LOC100910356; ENSRNOG00000007445", "Atp6ap1", "Atp6v1a", "Atp6v1c1", "Atp13a2", "Axl", "Azin1", "B3gnt8", "B4galt1", "B4galt5", "Blnk", "Blvra", "Blvrb", "C1qa", "C1qb", "C1qc", "C3", "C3ar1", "Camk1", "Camk2d", "Capg", "Card19", "Ccdc92; Glul", "Ccl2", "Ccl3", "Ccl4", "Ccl9", "Ccr2", "Ccr5", "Ccrl2", "Cd9", "Cd14", "Cd72", "Cd84", "Cd93", "Cd180", "Cd300lf", "Cd302", "Cebpb", "Cept1", "Cib1", "Clec4a1", "Clec4a3", "Clec4d", "Clec4e", "Clec5a", "Clec6a-ps1", "Clec12a", "Cln3", "Cmtm3", "Cndp2", "Comt", "Comtd1", "Cpd", "Creb5", "Creg1", "Csf1r", "Csf3r", "Cstb", "Ctsa", "Ctsb", "Ctsc", "Ctsd", "Ctss", "Cx3cr1", "Cxcl2", "Cxcr4", "Cybb", "Cyp4f18", "Cyp4v3", "Daglb", "Dgkz", "Dhrs3", "Dnase1l1", "Dnase2", "Dnmt3a", "Dok3", "Dpep2", "Dpp7", "Dram2", "Dtnbp1", "Ehd4", "Emilin2", "Emp1", "Ethe1", "Ets2", "F10", "Fcer1g", "Fcgr1a", "Fcgr2b", "Fcgr3a", "Fos", "Fosl2", "Frrs1", "Fth1", "Fuca2", "Furin", "Gadd45g", "Gatm", "Gbp4", "Gcsh", "Ggh", "Gk", "Gla", "Glrx", "Glul", "Gnaq", "Gns", "Gns", "Gpnmb", "Grina", "Grn", "Hacd4", "Hexa", "Hexb", "Hivep3", "Hivep3", "Hk2", "Hmox1", "Hpgds", "Id3", "Ier3", "Ifnar2", "Igsf6", "Il1b", "Il1rn", "Il6r", "Il10rb", "Il17ra", "Irf2bp2", "Itgb2", "Itgb5", "Itm2b", "Itpkb", "Kcnn4", "Klf4", "Klra2", "Klra2", "Knop1; Gde1", "Lair1", "Lamp2", "Lamtor4", "Lat2", "Lcp2", "Lgals3", "Lgals3bp", "Lgals5", "Lgals8", "Lgmn", "Lipa", "Lmna", "Lmo4", "LOC691141", "LOC100910497", "LOC100911469; Ms4a6b", "LOC103690035; Tgfbr1", "Lpcat2", "Lpxn", "Lrp1", "Lrpap1", "Lrrc25", "Lst1", "Luzp1", "Ly9", "Ly86", "Ly96", "Lyz2", "Mafb", "Mbnl2", "Mcfd2", "Mdfic", "Mdm2", "Metrnl", "Mfsd1", "Milr1", "Mmp14", "Mnda", "Mpeg1", "Ms4a6a", "Ms4a6b", "Ms4a6b", "Ms4a6c", "Ms4a7", "Msrb1", "Myo1f", "Myof", "Naglu", "Naip6", "Nceh1", "Ncf2", "Neu1", "Neurl3", "Ninj1", "Nlrp3", "Nr1h3", "Osm", "P2rx4", "P2ry6", "Pid1", "Pla2g7", "Pla2g15", "Plaur", "Pld3", "Pld4", "Plin2", "Plk3", "Plod3; Mir702", "Pltp", "Plxnb2; Mir349", "Por", "Pou2f2", "Ppt2", "Prdx5", "Psap", "Psen2", "Ptgs2", "Ptpre", "Rab7b", "Rab20", "Rcbtb2", "RGD1307554", "RGD1560687", "RGD1565410; Ly6i", "Rnf130", "Rnf149", "Rnpep", "Rrbp1", "Rtp4", "Sat1", "Scarb2", "Sdc4", "Sdcbp", "Sema4d", "Sgk1", "Sh2d1b2", "Sirpa", "Ski; RGD1565591", "Slamf9", "Slc11a1", "Slc15a3", "Slc16a3", "Slc29a1", "Slc31a2", "Slc35f6", "Slc43a2", "Snx5", "Snx18", "Soat1", "Socs3", "Sod2", "Spp1", "Susd3", "Sypl1", "Tbxas1", "Tcirg1", "Tcn2", "Tex264", "Tgfbi", "Tgfbr1", "Tgif1", "Tgm2", "Tlr2", "Tmed10", "Tmem37", "Tmem50a", "Tmem86a", "Tmem106a", "Tmem205", "Tnfaip2", "Tnfsf12", "Tpd52", "Tpi1", "Tpi1", "Tpp1", "Tppp3", "Tpst2", "Trem2", "Tyrobp", "Uap1l1", "Vav1", "Vegfa", "Ypel3", "Zc3h12a", "Zeb2", "Zfp703")
ThyMoDCMarkers <- c("Aif1", "Anxa5", "Apobec1", "Apoe", "App", "Aprt", "Arl4c", "Atp2b1", "Atp6v1b2", "Axl", "B4galt5", "Bach1", "Bach1", "Bcl3", "C1qa", "C1qb", "C1qc", "C3", "Capg", "Capza2", "Ccl9", "Ccr2", "Ccr5", "Ccrl2", "Cd14", "Cd274", "Cd300c", "Cebpb", "Cfp", "Clec7a", "Cltc", "Cotl1", "Csf1r", "Csf2ra", "Cstb", "Ctrl", "Ctsc", "Ctsd", "Ctss", "Ctsz", "Cxcr4", "Cyba", "Cyth1", "Dnaja2", "Dpysl2", "Fcer1g", "Fcgr2b", "Fgr", "Fth1", "Grina", "Hebp1", "Hexa", "Ifi47", "Ifit3", "Ifitm3", "Ifrd1", "Il1b", "Il1rn", "Il10ra", "Irf7", "Irgm", "Isg15", "Isg20", "Lamp2", "Lgals3", "Lgals3bp", "Lgals5", "Lgmn", "Lmna", "LOC100911469; Ms4a6b", "Lst1", "Lyz2", "Mafb", "Malt1", "Mnda", "Ms4a6a", "Ms4a6b", "Ms4a6b", "Myo1f", "Nadk", "Ncf2", "Nfam1", "Nfe2l2", "Nfkbie", "Oasl2", "Olfm1", "P2rx4", "Parp12", "Pitpna", "Plxnd1", "Ptafr", "RGD1560687", "RGD1565410; Ly6i", "Rnf130", "Rnf213", "Rtp4", "S100a4", "S100a6", "S100a10", "Sat1", "Scimp", "Sema4d", "Slc29a3", "Slfn2", "Slfn13", "Snap23", "Snx2", "Snx20", "Socs3", "Tcirg1", "Tdrd7", "Tgfbi", "Tlr2", "Tmem14c", "Tmem50a", "Tmem176b", "Tnip3", "Vdac2", "Zbp1", "Zeb2")

#Jaccard similarity
JaccardVectors <- list(SplcDC2bMarkers, SplMonoMarkers, LungInfDC2Markers, LungMCMarkers, ThyMoDCMarkers)
JaccardMatrix <- matrix(NA, nrow = length(JaccardVectors), ncol = length(JaccardVectors))
rownames(JaccardMatrix) <- colnames(JaccardMatrix) <- c("Splenic cDC2b", "Splenic monocytes",
                                                        "Lung inf-cDC2", "Lung moDC", "Thymic TdTOM+ moDC")

jaccard_similarity <- function(x, y) {
  length(intersect(x, y)) / length(union(x, y))
}

for (i in 1:length(JaccardVectors)) {
  for (j in 1:length(JaccardVectors)) {
    JaccardMatrix[i, j] <- jaccard_similarity(JaccardVectors[[i]], JaccardVectors[[j]])
  }
}

diag(JaccardMatrix) <- NA

pheatmap(JaccardMatrix, cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F, display_numbers = T) #Fig. S5G

#cDC2b markers
pheatmap(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp[c(SplcDC2bMarkers), 
                                                                         c("RatSplenicCD103neg_average", "RatSplenicDC2_average", "RatThymicCD103neg_average", "RatThymicDC2_average",
                                                                           "RNA.cDC2a", "RNA.cDC2b"
                                                                           #"RNA.monocytes", "RNA.cDC2a", "RNA.cDC2b",
                                                                           #"RNA.Mono.Mac", "RNA.cDC2.y", "RNA.mDC2", 
                                                                           #"RNA.Non.migratory.cDC2", "RNA.Migratory.cDC2", "RNA.inf.cDC2", "RNA.MC", 
                                                                           #"RNA.moDC", "RNA.cDC2.x"
                                                                           #"RNA.CLEC10A..cDC2",	"RNA.CCR7..cDC2",	"RNA.CLEC10A..cDC2.1", "RNA.Mitotic.cDC2"
                                                                         )],
         cellwidth = 25, cellheight = 3, cluster_rows = T, cluster_cols = F, display_numbers = F, fontsize_row = 3, scale = "row")  #Fig. 5F

#Inf-cDC2 markers
pheatmap(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp[c(LungInfDC2Markers), 
           c( 
             "RatSplenicCD103neg_average", "RatSplenicDC2_average", "RatThymicCD103neg_average", "RatThymicDC2_average",
             "RNA.inf.cDC2", "RNA.MC"
             #"RNA.cDC2b" , "RNA.cDC2a", "RNA.monocytes"
             #"RNA.Mono.Mac", "RNA.cDC2.y", "RNA.mDC2", 
             #"RNA.Non.migratory.cDC2", "RNA.Migratory.cDC2",
             #"RNA.moDC", "RNA.cDC2.x"
             #"RNA.CLEC10A..cDC2",	"RNA.CCR7..cDC2",	"RNA.CLEC10A..cDC2.1", "RNA.Mitotic.cDC2"
             )],
         cellwidth = 25, cellheight = 3, cluster_rows = T, cluster_cols = F, display_numbers = F, fontsize_row = 3, scale = "row") #Fig. 5G

#Lung moDC markers
pheatmap(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp[c(LungMCMarkers), 
                                                                         c("RatSplenicCD103neg_average", "RatSplenicDC2_average", "RatThymicCD103neg_average", "RatThymicDC2_average",
                                                                           "RNA.inf.cDC2", "RNA.MC"
                                                                           #"RNA.cDC2b" , "RNA.cDC2a", "RNA.monocytes"
                                                                           #"RNA.Mono.Mac", "RNA.cDC2.y", "RNA.mDC2", 
                                                                           #"RNA.Non.migratory.cDC2", "RNA.Migratory.cDC2",
                                                                           #"RNA.moDC", "RNA.cDC2.x"
                                                                           #"RNA.CLEC10A..cDC2",	"RNA.CCR7..cDC2",	"RNA.CLEC10A..cDC2.1", "RNA.Mitotic.cDC2"
                                                                         )],
         cellwidth = 25, cellheight = 1, cluster_rows = T, cluster_cols = F, display_numbers = F, fontsize_row = 1, scale = "row") #Fig. 5H

#thymus moDC markers
pheatmap(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp[c(ThyMoDCMarkers), 
                                                                         c("RNA.inf.cDC2", "RNA.MC", 
                                                                           "RatSplenicCD103neg_average", "RatSplenicDC2_average", "RatThymicCD103neg_average", "RatThymicDC2_average"
                                                                           #"RNA.cDC2b" , "RNA.cDC2a", "RNA.monocytes"
                                                                           #"RNA.Mono.Mac", "RNA.cDC2.y", "RNA.mDC2", 
                                                                           #"RNA.Non.migratory.cDC2", "RNA.Migratory.cDC2",
                                                                           #"RNA.moDC", "RNA.cDC2.x"
                                                                           #"RNA.CLEC10A..cDC2",	"RNA.CCR7..cDC2",	"RNA.CLEC10A..cDC2.1", "RNA.Mitotic.cDC2"
                                                                         )],
         cellwidth = 25, cellheight = 3, cluster_rows = T, cluster_cols = F, display_numbers = F, fontsize_row = 3, scale = "row")

#Human CLEC10A+ cDC2 markers
pheatmap(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp
         [c("Anxa1", "Il1r2", "Il18", "Sap30", "Lpar6", "Ifitm1", "Card9", "Smap2", "Il18", "Hcst", "Ptafr", "Traf3ip3", "Calhm2", "Cacna2d3", "Capg", "Pid1", "Lgals3", "Lilra5", "Stfa3l1", "Tnfsf13b", "Casp1", "Pdk4", "Tyrobp", "Serpina1", "Atf3", "Crebl2", "Dhrs9", "Rcsd1", "Ccdc92; Glul", "Tnfrsf1a", "Timp1", "Tnfaip3", "Klf2", "Msrb2", "Sectm1a", "Plaur", "Fcer1a", "Capn2", "Rb1", "Camk1", "Mad1l1", "Rnf130", "Gimap7", "Dok2", "Glul", "Il13ra1", "Tspo", "Kctd12", "Cd1d1", "Gimap4", "Phactr1", "Ethe1", "Fcer1g", "Gpat3", "Il1b", "Zfand5", "Cmtm3", "Cfd", "Siglec5", "Sgk1", "Csf1r", "Tnfrsf1b", "Rassf4", "Klf9", "Arrdc3", "S100a6", "Aplp2", "Clec10a", "Ets2", "Iqgap2; LOC100360623", "Mnda", "Oas1a", "Slc7a7", "Cited2", "Ms4a6a", "Fcgr2b", "Samsn1", "Cebpd", "Igsf6", "C1qa", "Zeb2", "Nlrp3", "Ms4a7", "Clec12a", "Cd302", "C1qb", "Csf3r", "Cebpb", "Gbp2", "S100a4"),
           c(#"RNA.inf.cDC2", "RNA.MC", 
             "RatSplenicCD103neg_average", "RatSplenicDC2_average", "RatThymicCD103neg_average", "RatThymicDC2_average",
             "RNA.CLEC10A..cDC2","RNA.CLEC10A..cDC2.1",	"RNA.CCR7..cDC2"
             #"RNA.cDC2b","RNA.monocytes", "RNA.cDC2a", 
             #"RNA.Mono.Mac", "RNA.cDC2.y", "RNA.mDC2", 
             #"RNA.Non.migratory.cDC2", "RNA.Migratory.cDC2",
             #"RNA.moDC", "RNA.cDC2.x"
             # "RNA.Mitotic.cDC2"
           )],
         cellwidth = 25, cellheight = 3, cluster_rows = T, cluster_cols = F, 
         display_numbers = F, fontsize_row = 2, scale = "row", show_rownames = T) #Fig. 6D

#Human CCR7+ cDC2 markers
pheatmap(RatDC_TomatoDC_MoThyDC_MoSimpleSplDC_HuSplDC_merged_rBE_std_exp
         [c("Vamp5", "Ap1s3", "Ifitm1", "Gpr137b", "Snn", "Pttg1", "Snx8", "Cfp", "Arhgap22", "Poglut1", "Gpr132", "Orai1", "Cul1", "Socs1", "Cd63", "Relb", "Ube2z", "Trim69", "Fscn1", "Dnajb11", "Samd9l", "Lgals3", "Tifab", "Lilra5", "Ccdc28b", "Ccl19", "Stfa3l1", "Ccr7", "Slc50a1", "Tnfsf13b", "Pim3", "Dusp5", "Acat2l1", "Serpina1", "Atf3", "Mt2A", "Traf1", "Spag9", "Crebl2", "Arl8b", "Bhlhe40", "Tymp", "Il27ra", "Timp1", "Rel", "Atf5", "Tnfaip3", "Birc3", "Cd83", "Qpct", "Lap3", "Plaur", "Hspb1", "Icam1", "Tbc1d8", "Pdlim5", "Gpx4", "Nfkb2", "Gch1", "Bcl2a1", "Hapln3", "Marcks; LOC681252", "Atp6v1a", "Isg20", "Casp7", "Ube2e1; LOC680961; ENSRNOG00000008794", "Nfkb1", "Mgll", "Tnfsf10", "Grina", "Creld2", "Nfkbie", "Ido1", "Serpinb9d", "Plek", "Gadd45b", "Bcl3", "Irf1", "Tnfaip2", "Ebi3", "Epsti1", "Ptpn1", "Ogfrl1", "Rab9a", "Il1b", "Lrrfip2", "MARCKSL1; Marcksl1; LOC100911158; LOC100911122", "Filip1l", "Ptgir", "Cd40", "Pea15", "Glipr2", "Stat4", "Ccdc167", "G0s2", "Il4i1", "Marcksl1", "Gbp4", "Rassf4", "Parp14", "Tnip2", "Ier3", "Marcksl1; LOC102550530", "Nampt", "Cdkn1a", "Ninj1", "Isg15", "Parvb", "Cflar", "C1qa", "Cst7", "Socs3", "C1qb", "Sod2", "Irf7", "Cxcl10", "Mx1", "Gbp2", "Prdm1", "Cxcl9", "Cd44")
           ,c(#"RNA.inf.cDC2", "RNA.MC",
             "RatSplenicCD103neg_average", "RatSplenicDC2_average", "RatThymicCD103neg_average", "RatThymicDC2_average",
             "RNA.CLEC10A..cDC2","RNA.CLEC10A..cDC2.1",	"RNA.CCR7..cDC2"
             #"RNA.cDC2b","RNA.monocytes", "RNA.cDC2a", 
             #"RNA.Mono.Mac", "RNA.cDC2.y", "RNA.mDC2", 
             #"RNA.Non.migratory.cDC2", "RNA.Migratory.cDC2",
             #"RNA.moDC", "RNA.cDC2.x"
             #"RNA.CLEC10A..cDC2.1", "RNA.Mitotic.cDC2"
           )],
         cellwidth = 25, cellheight = 3, cluster_rows = T, cluster_cols = F, 
         display_numbers = F, fontsize_row = 3, scale = "row", show_rownames = T) #Fig. 6E

#mouse lung MC vs other tissues heatmaps (rBE+std version)----
#calculate SD
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp$MouseLungSD <- apply(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[, c(46, 47, 48, 51)], 1, sd.p)

#order the matrix
RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp <- 
  RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[order(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp$MouseLungSD , decreasing = T),]

#MoThyTomato
MoThyTomatoCorrelation <- data.frame(matrix(ncol=4, nrow=2))
rownames(MoThyTomatoCorrelation) <- c("cDC2", "moDC")
colnames(MoThyTomatoCorrelation) <- c("Non-migr cDC2", "Migr cDC2", "inf-cDC2", "MC")

row_vars <- c("RNA.cDC2.x", "RNA.moDC")
col_vars <- c("RNA.Non.migr.cDC2",	"RNA.Migr.cDC2",	"RNA.inf.cDC2",	"RNA.MC")

for (i in seq_len(nrow(MoThyTomatoCorrelation))) {
  for (j in seq_len(ncol(MoThyTomatoCorrelation))) {
    MoThyTomatoCorrelation[i, j] <- (cor.test(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, col_vars[j]], 
                                              RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, row_vars[i]]))$estimate
  }
}

scaled_MoThyTomatoCorrelation <- MoThyTomatoCorrelation
scaled_MoThyTomatoCorrelation <- t(apply(MoThyTomatoCorrelation, 1, function(x) scale(as.numeric(x))))

pheatmap(scaled_MoThyTomatoCorrelation,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F, 
         display_numbers = apply(MoThyTomatoCorrelation, 2, function(x) sprintf("%.2f", x)))

pheatmap(MoThyTomatoCorrelation, display_numbers = T,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F)

#MoThy
MoThyCorrelation <- data.frame(matrix(nrow=10,ncol=4))
rownames(MoThyCorrelation) <- c("cDC2", "mDC2", "cDC1",  "mDC1", "pDC", "Mono/Mac", "Granulocyte", "T cell", "B cell", "NK cell")
colnames(MoThyCorrelation) <- c("Non-migr cDC2", "Migr cDC2", "inf-cDC2", "MC")

row_vars <- c("RNA.cDC2.y", "RNA.mDC2",	"RNA.cDC1.x",	"RNA.mDC1",	"RNA.pDC.y",	"RNA.Mono.Mac",	"RNA.Granulocytes",	"RNA.T.cells",	"RNA.B.cells",	"RNA.NK.cells")
col_vars <- c("RNA.Non.migr.cDC2",	"RNA.Migr.cDC2",	"RNA.inf.cDC2",	"RNA.MC")

for (i in seq_len(nrow(MoThyCorrelation))) {
  for (j in seq_len(ncol(MoThyCorrelation))) {
    MoThyCorrelation[i, j] <- (cor.test(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, col_vars[j]], 
                                        RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, row_vars[i]]))$estimate
  }
}

scaled_MoThyCorrelation <- MoThyCorrelation
scaled_MoThyCorrelation <- t(apply(MoThyCorrelation, 1, function(x) scale(as.numeric(x))))

pheatmap(scaled_MoThyCorrelation,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F, 
         display_numbers = apply(MoThyCorrelation, 2, function(x) sprintf("%.2f", x)))

pheatmap(MoThyCorrelation, display_numbers = T,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F)

#MoSpl
MoSplCorrelation <- data.frame(matrix(nrow=20,ncol=4))
rownames(MoSplCorrelation) <- c("RNA.immature.phagocytic.cDC1",	"RNA.ribosomal.cDC1",	"RNA.immature.mitotic.cDC1",	"RNA.mature.cDC1",	"RNA.mitochondrial.cDC1",	"RNA.immature.phagocytic.cDC2a",
                                "RNA.immature.phagocytic.Sox4..cDC",	"RNA.mitotic.cDC2a.1", "RNA.mitotic.cDC2a.2",	"RNA.mitotic.cDC2a.3",	"RNA.mitochondrial.cDC2..mixed.",	"RNA.mature.Il4i1..cDC2",
                                "RNA.mature.cDC2a",	"RNA.immature.phagocytic.cDC2b",	"RNA.mature.cDC2b",	"RNA.monocyte",	"RNA.pDC.x.1",	"RNA.NK.ILC.T",	"RNA.mixed.identity",	"RNA.stem.cell")
colnames(MoSplCorrelation) <- c("Non-migr cDC2", "Migr cDC2", "inf-cDC2", "MC")

row_vars <- c("RNA.immature.phagocytic.cDC1",	"RNA.ribosomal.cDC1",	"RNA.immature.mitotic.cDC1",	"RNA.mature.cDC1",	"RNA.mitochondrial.cDC1",	"RNA.immature.phagocytic.cDC2a",
              "RNA.immature.phagocytic.Sox4..cDC",	"RNA.mitotic.cDC2a.1", "RNA.mitotic.cDC2a.2",	"RNA.mitotic.cDC2a.3",	"RNA.mitochondrial.cDC2..mixed.",	"RNA.mature.Il4i1..cDC2",
              "RNA.mature.cDC2a",	"RNA.immature.phagocytic.cDC2b",	"RNA.mature.cDC2b",	"RNA.monocyte",	"RNA.pDC.x.1",	"RNA.NK.ILC.T",	"RNA.mixed.identity",	"RNA.stem.cell")
col_vars <- c("RNA.Non.migr.cDC2",	"RNA.Migr.cDC2",	"RNA.inf.cDC2",	"RNA.MC")

for (i in seq_len(nrow(MoSplCorrelation))) {
  for (j in seq_len(ncol(MoSplCorrelation))) {
    MoSplCorrelation[i, j] <- (cor.test(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, col_vars[j]], 
                                        RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, row_vars[i]]))$estimate
  }
}

scaled_MoSplCorrelation <- MoSplCorrelation
scaled_MoSplCorrelation <- t(apply(MoSplCorrelation, 1, function(x) scale(as.numeric(x))))

pheatmap(scaled_MoSplCorrelation,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F, 
         display_numbers = apply(MoSplCorrelation, 2, function(x) sprintf("%.2f", x)))

pheatmap(MoSplCorrelation, display_numbers = T,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F)

#HuSpl
HuSplCorrelation <- data.frame(matrix(nrow=7,ncol=4))
rownames(HuSplCorrelation) <- c("RNA.cDC1.y",	"RNA.Mitotic.cDC1",	"RNA.CLEC10A..cDC2",	"RNA.CLEC10A..cDC2.1",	"RNA.CCR7..cDC2",	"RNA.Mitotic.cDC2",	"RNA.AS.DC")
colnames(HuSplCorrelation) <- c("Non-migr cDC2", "Migr cDC2", "inf-cDC2", "MC")

row_vars <- c("RNA.cDC1.y",	"RNA.Mitotic.cDC1",	"RNA.CLEC10A..cDC2",	"RNA.CLEC10A..cDC2.1",	"RNA.CCR7..cDC2",	"RNA.Mitotic.cDC2",	"RNA.AS.DC")
col_vars <- c("RNA.Non.migr.cDC2",	"RNA.Migr.cDC2",	"RNA.inf.cDC2",	"RNA.MC")

for (i in seq_len(nrow(HuSplCorrelation))) {
  for (j in seq_len(ncol(HuSplCorrelation))) {
    HuSplCorrelation[i, j] <- (cor.test(RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, col_vars[j]], 
                                        RatDC_TomatoDC_MoThyDC_MoSplDC_HuSplDC_merged_rBE_std_exp[1:1000, row_vars[i]]))$estimate
  }
}

scaled_HuSplCorrelation <- HuSplCorrelation
scaled_HuSplCorrelation <- t(apply(HuSplCorrelation, 1, function(x) scale(as.numeric(x))))

pheatmap(scaled_HuSplCorrelation,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F, 
         display_numbers = apply(HuSplCorrelation, 2, function(x) sprintf("%.2f", x)))

pheatmap(HuSplCorrelation, display_numbers = T,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F)