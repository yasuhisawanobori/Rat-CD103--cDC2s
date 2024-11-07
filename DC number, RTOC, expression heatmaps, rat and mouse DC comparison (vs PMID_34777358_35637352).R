library(readr)
library(org.Rn.eg.db)
library(limma)
library(clusterProfiler)
library(openxlsx)
library(enrichplot)
library(ggplot2)
library("pheatmap")
library(orthogene)

#figures of cell numbers (Fig. 1C)----
# read data of DC numbers
DC_count.matrix <- read.delim("PMID_34777358_35637352 mouse cDC1 cDC2 correlation/20231031 CD103- cell number.txt")

# Define the order of groups
group_order <- c("CD103-",  "cDC2", "cDC1", "noDC")

# Convert the 'group' column to a factor with specified order
DC_count.matrix$group <- factor(DC_count.matrix$group, levels = group_order)

# Custom colors for points
point_colors <- c("noDC" = "#000000", "cDC1" = "#E60012", "cDC2" = "#0000FF", "CD103-" = "#009944")

# Create a dot plot for DC number data (4.95x3.55 inches)
#thymic DC% in anti-MHCII magnetic sorted cells
ggplot(subset(DC_count.matrix, organ == "thymus"), aes(x = group, y = `X..in.concentrated.cells`, fill = group)) +
  geom_point(size = 5, shape = 21, color = "transparent") +  # Set color to "transparent"
  stat_summary(fun = "mean", geom = "crossbar", width = 0.5) +
  labs(title = "", y = "", x = "") +
  scale_fill_manual(values = point_colors) +
  guides(fill = FALSE) +  # Remove the fill legend
  theme(panel.background = element_rect(fill = "white"),  # Set plot area color to white
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Set axis line color
        axis.text = element_text(size = 18, color = "black"),  # Set axis text size and color
        axis.ticks = element_line(color = "black", size = 1),  # Set axis ticks color and size
        axis.title.y = element_text(vjust = 2), # Adjust Y-axis title position
        axis.text.x = element_blank()) +  # Remove X-axis text
  ylim(0, 10)

ggplot(subset(DC_count.matrix, organ == "thymus"), aes(x = group, y = `cells..x10.4cells.100mg.tissue.`, fill = group)) +
  geom_point(size = 5, shape = 21, color = "transparent") +  # Set color to "transparent"
  stat_summary(fun = "mean", geom = "crossbar", width = 0.5) +
  labs(title = "", y = "", x = "") +
  scale_fill_manual(values = point_colors) +
  guides(fill = FALSE) +  # Remove the fill legend
  theme(panel.background = element_rect(fill = "white"),  # Set plot area color to white
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Set axis line color
        axis.text = element_text(size = 18, color = "black"),  # Set axis text size and color
        axis.ticks = element_line(color = "black", size = 1),  # Set axis ticks color and size
        axis.title.y = element_text(vjust = 2), # Adjust Y-axis title position
        axis.text.x = element_blank()) +  # Remove X-axis text
  ylim(0, 13)

ggplot(subset(DC_count.matrix, organ == "spleen"), aes(x = group, y = `X..in.concentrated.cells`, fill = group)) +
  geom_point(size = 5, shape = 21, color = "transparent") +  # Set color to "transparent"
  stat_summary(fun = "mean", geom = "crossbar", width = 0.5) +
  labs(title = "", y = "", x = "") +
  scale_fill_manual(values = point_colors) +
  guides(fill = FALSE) +  # Remove the fill legend
  theme(panel.background = element_rect(fill = "white"),  # Set plot area color to white
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Set axis line color
        axis.text = element_text(size = 18, color = "black"),  # Set axis text size and color
        axis.ticks = element_line(color = "black", size = 1),  # Set axis ticks color and size
        axis.title.y = element_text(vjust = 2), # Adjust Y-axis title position
        axis.text.x = element_blank()) +  # Remove X-axis text
  ylim(0, 13)

ggplot(subset(DC_count.matrix, organ == "spleen"), aes(x = group, y = `cells..x10.4cells.100mg.tissue.`, fill = group)) +
  geom_point(size = 5, shape = 21, color = "transparent") +  # Set color to "transparent"
  stat_summary(fun = "mean", geom = "crossbar", width = 0.5) +
  labs(title = "", y = "", x = "") +
  scale_fill_manual(values = point_colors) +
  guides(fill = FALSE) +  # Remove the fill legend
  theme(panel.background = element_rect(fill = "white"),  # Set plot area color to white
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(color = "black"),  # Set axis line color
        axis.text = element_text(size = 18, color = "black"),  # Set axis text size and color
        axis.ticks = element_line(color = "black", size = 1),  # Set axis ticks color and size
        axis.title.y = element_text(vjust = 2), # Adjust Y-axis title position
        axis.text.x = element_blank()) +  # Remove X-axis text
  ylim(0, 55)  # Set the y-axis limits  # Force Y-axis to start from 0

#RTOC1-1 + 1-2 + 1-3 (Fig. 4B-D)----
RTOC1_1to3matrix <- read_delim("20231225 RTOC1-1~3/20240819 RTOC1-1+1-2+1-3 matrix.txt", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)

# Define the order of groups
group_order <- c("CD103-cDC2", "CD103+cDC2", "cDC1", "noDC")

# Convert the 'group' column to a factor with specified order
RTOC1_1to3matrix$`Sample type` <- factor(RTOC1_1to3matrix$`Sample type`, levels = group_order)

# Custom colors for points
point_colors <- c("noDC" = "#000000", "cDC1" = "#E60012", "CD103+cDC2" = "#0000FF", "CD103-cDC2" = "#009944")
experiment_colors <- c("RTOC1-1" = "red", "RTOC1-2" = "blue", "RTOC1-3" = "green")

# Create a dot plot with crossbars for mean, and custom point colors
RTOC_plot <- function(data, y_column) {
  ggplot(data, aes(x = `Sample type`, y = .data[[y_column]], fill = `Sample type`)) +
    geom_jitter(size = 5, shape = 21, color = "transparent") +  # Set color to "transparent"
    stat_summary(fun = "mean", geom = "crossbar", width = 0.5) +
    labs(title = y_column, y = "", x = "") +
    scale_fill_manual(values = point_colors) +
    guides(fill = 'none') +  # Remove the fill legend
    theme(panel.background = element_rect(fill = "white"),  # Set plot area color to white
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          axis.line = element_line(color = "black"),  # Set axis line color
          axis.text = element_text(size = 18, color = "black"),  # Set axis text size and color
          axis.ticks = element_line(color = "black", size = 1),  # Set axis ticks color and size
          axis.title.y = element_text(vjust = 2),   # Adjust Y-axis title position
          axis.text.x = element_blank())+  # Remove x-axis labels
    expand_limits(y = c(0, 2))  # Set Y-axis limits
}

RTOC_plot2 <- function(data, y_column) {
  ggplot(data, aes(x = `Sample type`, y = .data[[y_column]], fill = `Sample type`)) +
    geom_jitter(size = 5, shape = 21, color = "transparent") +  # Set color to "transparent"
    stat_summary(fun = "mean", geom = "crossbar", width = 0.5) +
    labs(title = y_column, y = "", x = "") +
    scale_fill_manual(values = point_colors) +
    guides(fill = 'none') +  # Remove the fill legend
    theme(panel.background = element_rect(fill = "white"),  # Set plot area color to white
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          axis.line = element_line(color = "black"),  # Set axis line color
          axis.text = element_text(size = 18, color = "black"),  # Set axis text size and color
          axis.ticks = element_line(color = "black", size = 1),  # Set axis ticks color and size
          axis.title.y = element_text(vjust = 2),   # Adjust Y-axis title position
          axis.text.x = element_blank())+  # Remove x-axis labels
    expand_limits(y = c(0.8, 1.2))  # Set Y-axis limits
}

#6.5x4.5 inches
RTOC_plot(RTOC1_1to3matrix, "CD45+ N")
RTOC_plot(RTOC1_1to3matrix, "DP %")

RTOC_plot(RTOC1_1to3matrix, "CD4SP %")
RTOC_plot(RTOC1_1to3matrix, "CD8SP %") 
RTOC_plot(RTOC1_1to3matrix, "DP N")
RTOC_plot(RTOC1_1to3matrix, "CD4SP N")
RTOC_plot(RTOC1_1to3matrix, "CD8SP N")
RTOC_plot2(RTOC1_1to3matrix, "CD4SPCD62L %")
RTOC_plot2(RTOC1_1to3matrix, "CD8SPCD62L %")
RTOC_plot(RTOC1_1to3matrix, "CD4TCRhi p")
RTOC_plot(RTOC1_1to3matrix, "CD25FOXP3%")
RTOC_plot(RTOC1_1to3matrix, "% of CD25+Foxp3+")


#read mouse thymic DC data----
MoThyDC_Seq<- as.data.frame(read_delim("PMID_34777358_35637352 mouse cDC1 cDC2 correlation/GSE198789_subread_counts_gene_symbol.txt.gz", 
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE))

rownames(MoThyDC_Seq) <- MoThyDC_Seq$GeneName #set values of GeneName as row names

MoThyDC_Seq <- MoThyDC_Seq[, -1] #remove first (GeneName) column

MoThyDC_Seq <- log2(MoThyDC_Seq-1) # convert values to log2

MoThyDC_Seq_conv <- convert_orthologs(MoThyDC_Seq, input_species = "mouse", output_species = "rat", 
                                 drop_nonorths = F, non121_strategy = "1")

#read mouse splenic DC data----
MoSplDC_seq <- as.data.frame(read_csv("PMID_34777358_35637352 mouse cDC1 cDC2 correlation/GSE181475_gene_count_matrix.csv.gz"))

#convert rownames to symbols. If there is no symbols for ENSEMBL IDs, leave the IDs intact.
ensembl_ids <- MoSplDC_seq[,1]
gene_symbols <- mapIds(org.Mm.eg.db, keys = ensembl_ids, column = 'SYMBOL', keytype = 'ENSEMBL')
na_indices <- is.na(gene_symbols)
gene_symbols[na_indices] <- ensembl_ids[na_indices]

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
rownames(MoSplDC_seq) <- unique_gene_symbols
rm(unique_gene_symbols, gene_symbols, na_indices, ensembl_ids, dup_indices)

MoSplDC_seq <- MoSplDC_seq[,-1] #remove first column
colnames(MoSplDC_seq) <- c("TLI_DC1a-1", "TLI_DC1b-1", "TLI_DC3-1", "TLI_DC2-1", #1-4
                           "TLI_DC1a-2", "TLI_DC1b-2", "TLI_DC3-2", #5-7
                           "DC1a-1", "DC1b-1", "DC3-1", "DC2-1", #9-12
                           "DC1a-2", "DC1b-2", "DC3-2", "DC2-2", #13-16
                           "TLI_DC1a-3", "TLI_DC1b-3", "TLI_DC3-3", "TLI_DC2-3", #17-20
                           "TLI_DC1a-4", "TLI_DC1b-4", "TLI_DC3-4", "TLI_DC2-4") #21-24

MoSplDC_seq <- log2(MoSplDC_seq-1) # convert values to log2

MoSplDC_seq_conv <- convert_orthologs(MoSplDC_seq, input_species = "mouse", output_species = "rat", 
                                 drop_nonorths = F, non121_strategy = "1")

#read rat DC/Mf data----
RatDC_Mf_array <- as.data.frame(read_csv("PMID_34777358_35637352 mouse cDC1 cDC2 correlation/rat_DC_Mf_array-TACnormalized.csv"))

gene_symbols <- RatDC_Mf_array$`Gene Symbol` 

dup_indices <- duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE) # Find duplicated values

count <- ave(seq_along(gene_symbols), gene_symbols, FUN = seq_along) # Count the number of occurrences of each value

gene_symbols[dup_indices] <- paste0(gene_symbols[dup_indices], "-", count[dup_indices] - 1) # Append unique numbers to duplicates

gene_symbols <- gsub("-0", "", gene_symbols) #remove "-0"

row.names(RatDC_Mf_array) <- gene_symbols #Set gene symbols to row names

RatDC_Mf_array <- RatDC_Mf_array[, -(1:3)] #remove unnecessary columns

#merge the converted objects----
merged_seq <- merge(RatDC_Mf_array, MoThyDC_Seq_conv, by="row.names", all=T)
row.names(merged_seq) <- merged_seq$Row.names
merged_seq<-merged_seq[,-1]

merged_seq <- merge(merged_seq, MoSplDC_seq_conv, by="row.names", all=T)
row.names(merged_seq) <- merged_seq$Row.names
merged_seq<-merged_seq[,-1]

#remove batch effects
batch <- c(rep("RatArray", 16), rep("MoThy", 9), rep("MoSpl", 23))
merged_seq <- removeBatchEffect(merged_seq, batch)

#remove NaNs and -Infs
merged_seq[merged_seq  == "NaN" | merged_seq  == "-Inf"] <- "" 

#calculate averages
merged_seq <- as.data.frame(merged_seq) #change class from matrix to data.frame
merged_seq[,1:48] <- lapply(merged_seq[,1:48], as.numeric) #change class of the contents from character to numeric

merged_seq$RatThymicDC1_average <- log2(rowMeans(2^merged_seq[, 1:2], na.rm = T))
merged_seq$RatSplenicDC1_average <- log2(rowMeans(2^merged_seq[, 3:4], na.rm = T))
merged_seq$RatThymicDC2_average <- log2(rowMeans(2^merged_seq[, 5:6], na.rm = T))
merged_seq$RatSplenicDC2_average <- log2(rowMeans(2^merged_seq[, 7:8], na.rm = T))
merged_seq$RatThymicCD103neg_average <- log2(rowMeans(2^merged_seq[, 9:10], na.rm = T))
merged_seq$RatSplenicCD103neg_average <- log2(rowMeans(2^merged_seq[, 11:12], na.rm = T))
merged_seq$RatThymicMf_average <- log2(rowMeans(2^merged_seq[, 13:14], na.rm = T))
merged_seq$RatSplenicMf_average <- log2(rowMeans(2^merged_seq[, 15:16], na.rm = T))

merged_seq$MoThymicDC1_average <- log2(rowMeans(2^merged_seq[, 23:25], na.rm = T))
merged_seq$MoThymicCD301bposDC2_average <- log2(rowMeans(2^merged_seq[, 17:19], na.rm = T))
merged_seq$MoThymicCD301bnegDC2_average <- log2(rowMeans(2^merged_seq[, 20:22], na.rm = T))
merged_seq$MoThymicDC2_average <- as.numeric(ifelse(
  is.na(merged_seq[, "MoThymicCD301bposDC2_average"]) & is.na(merged_seq[, "MoThymicCD301bnegDC2_average"]),"NaN",
  ifelse(
    is.na(merged_seq[, "MoThymicCD301bposDC2_average"]),
    log2(0.54*2^merged_seq[, "MoThymicCD301bnegDC2_average"]),
    ifelse(
      is.na(merged_seq[, "MoThymicCD301bnegDC2_average"]),
      log2(0.46*2^merged_seq[, "MoThymicCD301bposDC2_average"]),
      log2(0.46 * 2^merged_seq[, "MoThymicCD301bposDC2_average"] + 0.54 * 2^merged_seq[, "MoThymicCD301bnegDC2_average"])
    ))))

merged_seq$MoSplenicDC1a_average <- log2(rowMeans(2^merged_seq[, c(33, 37)], na.rm = T))
merged_seq$MoSplenicDC1b_average <- log2(rowMeans(2^merged_seq[, c(34, 38)], na.rm = T))
merged_seq$MoSplenicDC2_average <- log2(rowMeans(2^merged_seq[, c(36, 40)], na.rm = T))

#scaling
merged_seq <- as.data.frame(scale(merged_seq))

#calculate and add SD columns
sd.p=function(x){sd(x, na.rm = TRUE)*sqrt((length(x)-1)/length(x))}
merged_seq$RatThyDcSD <- as.numeric(apply(merged_seq[, c("RatThymicDC1_average", "RatThymicDC2_average", "RatThymicCD103neg_average")], 1, sd.p))
merged_seq$RatSplDcSD <- as.numeric(apply(merged_seq[, c("RatSplenicDC1_average", "RatSplenicDC2_average", "RatSplenicCD103neg_average")], 1, sd.p))

merged_seq <- merged_seq[order(merged_seq$RatThyDcSD, decreasing = T),]

#export the expression matrix
write.csv(merged_seq, "PMID_34777358_35637352 mouse cDC1 cDC2 correlation/merged_seq.csv", na ="", row.names = TRUE)

merged_seq <- read.csv("PMID_34777358_35637352 mouse cDC1 cDC2 correlation/merged_seq.csv", row.names = 1)

#calculate correlations and depict heatmaps----
RatDcMoDcCorrelation <- data.frame(matrix(nrow=4, ncol=6))
rownames(RatDcMoDcCorrelation) <- c("MoSplDC2", "MoSplDC1a", "MoThyDC2", "MoThyDC1")
colnames(RatDcMoDcCorrelation) <- c("RatSplCD103-", "RatSplCD103+", "RatSplDC1", "RatThyCD103-", "RatThyCD103+", "RatThyDC1")

merged_seq <- merged_seq[order(merged_seq$RatThyDcSD, decreasing = T),]

RatDcMoDcCorrelation["MoThyDC1", "RatThyDC1"] <-(cor.test(merged_seq[1:1000, "MoThymicDC1_average"], 
                                                          merged_seq[1:1000, "RatThymicDC1_average"]))$estimate
RatDcMoDcCorrelation["MoThyDC2", "RatThyDC1"] <-(
  cor.test(merged_seq[1:1000, "MoThymicDC2_average"], 
                                                          merged_seq[1:1000, "RatThymicDC1_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC1a", "RatThyDC1"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC1a_average"], 
                                                          merged_seq[1:1000, "RatThymicDC1_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC2", "RatThyDC1"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC2_average"], 
                                                           merged_seq[1:1000, "RatThymicDC1_average"]))$estimate

RatDcMoDcCorrelation["MoThyDC1", "RatThyCD103+"] <-(cor.test(merged_seq[1:1000, "MoThymicDC1_average"], 
                                                          merged_seq[1:1000, "RatThymicDC2_average"]))$estimate
RatDcMoDcCorrelation["MoThyDC2", "RatThyCD103+"] <-(cor.test(merged_seq[1:1000, "MoThymicDC2_average"], 
                                                          merged_seq[1:1000, "RatThymicDC2_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC1a", "RatThyCD103+"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC1a_average"], 
                                                           merged_seq[1:1000, "RatThymicDC2_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC2", "RatThyCD103+"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC2_average"], 
                                                          merged_seq[1:1000, "RatThymicDC2_average"]))$estimate

RatDcMoDcCorrelation["MoThyDC1", "RatThyCD103-"] <-(cor.test(merged_seq[1:1000, "MoThymicDC1_average"], 
                                                             merged_seq[1:1000, "RatThymicCD103neg_average"]))$estimate
RatDcMoDcCorrelation["MoThyDC2", "RatThyCD103-"] <-(cor.test(merged_seq[1:1000, "MoThymicDC2_average"], 
                                                             merged_seq[1:1000, "RatThymicCD103neg_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC1a", "RatThyCD103-"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC1a_average"], 
                                                              merged_seq[1:1000, "RatThymicCD103neg_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC2", "RatThyCD103-"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC2_average"], 
                                                             merged_seq[1:1000, "RatThymicCD103neg_average"]))$estimate

merged_seq <- merged_seq[order(merged_seq$RatSplDcSD, decreasing = T),]

RatDcMoDcCorrelation["MoThyDC1", "RatSplDC1"] <-(cor.test(merged_seq[1:1000, "MoThymicDC1_average"], 
                                                          merged_seq[1:1000, "RatSplenicDC1_average"]))$estimate
RatDcMoDcCorrelation["MoThyDC2", "RatSplDC1"] <-(cor.test(merged_seq[1:1000, "MoThymicDC2_average"], 
                                                          merged_seq[1:1000, "RatSplenicDC1_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC1a", "RatSplDC1"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC1a_average"], 
                                                           merged_seq[1:1000, "RatSplenicDC1_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC2", "RatSplDC1"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC2_average"], 
                                                          merged_seq[1:1000, "RatSplenicDC1_average"]))$estimate

RatDcMoDcCorrelation["MoThyDC1", "RatSplCD103+"] <-(cor.test(merged_seq[1:1000, "MoThymicDC1_average"], 
                                                          merged_seq[1:1000, "RatSplenicDC2_average"]))$estimate
RatDcMoDcCorrelation["MoThyDC2", "RatSplCD103+"] <-(cor.test(merged_seq[1:1000, "MoThymicDC2_average"], 
                                                          merged_seq[1:1000, "RatSplenicDC2_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC1a", "RatSplCD103+"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC1a_average"], 
                                                           merged_seq[1:1000, "RatSplenicDC2_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC2", "RatSplCD103+"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC2_average"], 
                                                          merged_seq[1:1000, "RatSplenicDC2_average"]))$estimate

RatDcMoDcCorrelation["MoThyDC1", "RatSplCD103-"] <-(cor.test(merged_seq[1:1000, "MoThymicDC1_average"], 
                                                             merged_seq[1:1000, "RatSplenicCD103neg_average"]))$estimate
RatDcMoDcCorrelation["MoThyDC2", "RatSplCD103-"] <-(cor.test(merged_seq[1:1000, "MoThymicDC2_average"], 
                                                             merged_seq[1:1000, "RatSplenicCD103neg_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC1a", "RatSplCD103-"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC1a_average"], 
                                                              merged_seq[1:1000, "RatSplenicCD103neg_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC2", "RatSplCD103-"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC2_average"], 
                                                             merged_seq[1:1000, "RatSplenicCD103neg_average"]))$estimate

scaled_RatDcMoDcCorrelation <- RatDcMoDcCorrelation
scaled_RatDcMoDcCorrelation["MoThyDC1",] <- scale(as.numeric(RatDcMoDcCorrelation["MoThyDC1",]))
scaled_RatDcMoDcCorrelation["MoThyDC2",] <- scale(as.numeric(RatDcMoDcCorrelation["MoThyDC2",]))
scaled_RatDcMoDcCorrelation["MoSplDC1a",] <- scale(as.numeric(RatDcMoDcCorrelation["MoSplDC1a",]))
scaled_RatDcMoDcCorrelation["MoSplDC2",] <- scale(as.numeric(RatDcMoDcCorrelation["MoSplDC2",]))

#Fig. 2D
pheatmap(scaled_RatDcMoDcCorrelation,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F, 
         display_numbers = apply(RatDcMoDcCorrelation, 2, function(x) sprintf("%.2f", x)))

pheatmap(RatDcMoDcCorrelation,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F, 
         display_numbers = T)

#calculate correlations and depict heatmaps (CD301b+/-)----
RatDcMoDcCorrelation <- data.frame(matrix(nrow=6, ncol=6))
rownames(RatDcMoDcCorrelation) <- c("MoSplDC2", "MoSplDC1a", "MoThymicCD301bnegDC2", "MoThymicCD301bposDC2", "MoThyDC2", "MoThyDC1")
colnames(RatDcMoDcCorrelation) <- c("RatSplCD103-", "RatSplCD103+", "RatSplDC1", "RatThyCD103-", "RatThyCD103+", "RatThyDC1")

merged_seq <- merged_seq[order(merged_seq$RatThyDcSD, decreasing = T),]

RatDcMoDcCorrelation["MoThyDC1", "RatThyDC1"] <-(cor.test(merged_seq[1:1000, "MoThymicDC1_average"], 
                                                          merged_seq[1:1000, "RatThymicDC1_average"]))$estimate
RatDcMoDcCorrelation["MoThyDC2", "RatThyDC1"] <-(cor.test(merged_seq[1:1000, "MoThymicDC2_average"], 
                                                          merged_seq[1:1000, "RatThymicDC1_average"]))$estimate
RatDcMoDcCorrelation["MoThymicCD301bnegDC2", "RatThyDC1"] <-(cor.test(merged_seq[1:1000, "MoThymicCD301bnegDC2_average"], 
                                                                      merged_seq[1:1000, "RatThymicDC1_average"]))$estimate
RatDcMoDcCorrelation["MoThymicCD301bposDC2", "RatThyDC1"] <-(cor.test(merged_seq[1:1000, "MoThymicCD301bposDC2_average"], 
                                                                      merged_seq[1:1000, "RatThymicDC1_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC1a", "RatThyDC1"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC1a_average"], 
                                                          merged_seq[1:1000, "RatThymicDC1_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC2", "RatThyDC1"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC2_average"], 
                                                           merged_seq[1:1000, "RatThymicDC1_average"]))$estimate

RatDcMoDcCorrelation["MoThyDC1", "RatThyCD103+"] <-(cor.test(merged_seq[1:1000, "MoThymicDC1_average"], 
                                                          merged_seq[1:1000, "RatThymicDC2_average"]))$estimate
RatDcMoDcCorrelation["MoThyDC2", "RatThyCD103+"] <-(cor.test(merged_seq[1:1000, "MoThymicDC2_average"], 
                                                          merged_seq[1:1000, "RatThymicDC2_average"]))$estimate
RatDcMoDcCorrelation["MoThymicCD301bnegDC2", "RatThyCD103+"] <-(cor.test(merged_seq[1:1000, "MoThymicCD301bnegDC2_average"], 
                                                                      merged_seq[1:1000, "RatThymicDC2_average"]))$estimate
RatDcMoDcCorrelation["MoThymicCD301bposDC2", "RatThyCD103+"] <-(cor.test(merged_seq[1:1000, "MoThymicCD301bposDC2_average"], 
                                                                      merged_seq[1:1000, "RatThymicDC2_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC1a", "RatThyCD103+"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC1a_average"], 
                                                           merged_seq[1:1000, "RatThymicDC2_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC2", "RatThyCD103+"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC2_average"], 
                                                          merged_seq[1:1000, "RatThymicDC2_average"]))$estimate

RatDcMoDcCorrelation["MoThyDC1", "RatThyCD103-"] <-(cor.test(merged_seq[1:1000, "MoThymicDC1_average"], 
                                                             merged_seq[1:1000, "RatThymicCD103neg_average"]))$estimate
RatDcMoDcCorrelation["MoThyDC2", "RatThyCD103-"] <-(cor.test(merged_seq[1:1000, "MoThymicDC2_average"], 
                                                             merged_seq[1:1000, "RatThymicCD103neg_average"]))$estimate
RatDcMoDcCorrelation["MoThymicCD301bnegDC2", "RatThyCD103-"] <-(cor.test(merged_seq[1:1000, "MoThymicCD301bnegDC2_average"], 
                                                                         merged_seq[1:1000, "RatThymicCD103neg_average"]))$estimate
RatDcMoDcCorrelation["MoThymicCD301bposDC2", "RatThyCD103-"] <-(cor.test(merged_seq[1:1000, "MoThymicCD301bposDC2_average"], 
                                                                         merged_seq[1:1000, "RatThymicCD103neg_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC1a", "RatThyCD103-"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC1a_average"], 
                                                              merged_seq[1:1000, "RatThymicCD103neg_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC2", "RatThyCD103-"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC2_average"], 
                                                             merged_seq[1:1000, "RatThymicCD103neg_average"]))$estimate

merged_seq <- merged_seq[order(merged_seq$RatSplDcSD, decreasing = T),]

RatDcMoDcCorrelation["MoThymicCD301bnegDC2", "RatSplDC1"] <-(cor.test(merged_seq[1:1000, "MoThymicCD301bnegDC2_average"], 
                                                                         merged_seq[1:1000, "RatSplenicDC1_average"]))$estimate
RatDcMoDcCorrelation["MoThymicCD301bposDC2", "RatSplDC1"] <-(cor.test(merged_seq[1:1000, "MoThymicCD301bposDC2_average"], 
                                                                         merged_seq[1:1000, "RatSplenicDC1_average"]))$estimate
RatDcMoDcCorrelation["MoThyDC1", "RatSplDC1"] <-(cor.test(merged_seq[1:1000, "MoThymicDC1_average"], 
                                                          merged_seq[1:1000, "RatSplenicDC1_average"]))$estimate
RatDcMoDcCorrelation["MoThyDC2", "RatSplDC1"] <-(cor.test(merged_seq[1:1000, "MoThymicDC2_average"], 
                                                          merged_seq[1:1000, "RatSplenicDC1_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC1a", "RatSplDC1"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC1a_average"], 
                                                           merged_seq[1:1000, "RatSplenicDC1_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC2", "RatSplDC1"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC2_average"], 
                                                          merged_seq[1:1000, "RatSplenicDC1_average"]))$estimate

RatDcMoDcCorrelation["MoThymicCD301bnegDC2", "RatSplCD103+"] <-(cor.test(merged_seq[1:1000, "MoThymicCD301bnegDC2_average"], 
                                                                      merged_seq[1:1000, "RatSplenicDC2_average"]))$estimate
RatDcMoDcCorrelation["MoThymicCD301bposDC2", "RatSplCD103+"] <-(cor.test(merged_seq[1:1000, "MoThymicCD301bposDC2_average"], 
                                                                      merged_seq[1:1000, "RatSplenicDC2_average"]))$estimate
RatDcMoDcCorrelation["MoThyDC1", "RatSplCD103+"] <-(cor.test(merged_seq[1:1000, "MoThymicDC1_average"], 
                                                          merged_seq[1:1000, "RatSplenicDC2_average"]))$estimate
RatDcMoDcCorrelation["MoThyDC2", "RatSplCD103+"] <-(cor.test(merged_seq[1:1000, "MoThymicDC2_average"], 
                                                          merged_seq[1:1000, "RatSplenicDC2_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC1a", "RatSplCD103+"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC1a_average"], 
                                                           merged_seq[1:1000, "RatSplenicDC2_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC2", "RatSplCD103+"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC2_average"], 
                                                          merged_seq[1:1000, "RatSplenicDC2_average"]))$estimate

RatDcMoDcCorrelation["MoThymicCD301bnegDC2", "RatSplCD103-"] <-(cor.test(merged_seq[1:1000, "MoThymicCD301bnegDC2_average"], 
                                                                         merged_seq[1:1000, "RatSplenicCD103neg_average"]))$estimate
RatDcMoDcCorrelation["MoThymicCD301bposDC2", "RatSplCD103-"] <-(cor.test(merged_seq[1:1000, "MoThymicCD301bposDC2_average"], 
                                                                         merged_seq[1:1000, "RatSplenicCD103neg_average"]))$estimate
RatDcMoDcCorrelation["MoThyDC1", "RatSplCD103-"] <-(cor.test(merged_seq[1:1000, "MoThymicDC1_average"], 
                                                             merged_seq[1:1000, "RatSplenicCD103neg_average"]))$estimate
RatDcMoDcCorrelation["MoThyDC2", "RatSplCD103-"] <-(cor.test(merged_seq[1:1000, "MoThymicDC2_average"], 
                                                             merged_seq[1:1000, "RatSplenicCD103neg_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC1a", "RatSplCD103-"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC1a_average"], 
                                                              merged_seq[1:1000, "RatSplenicCD103neg_average"]))$estimate
RatDcMoDcCorrelation["MoSplDC2", "RatSplCD103-"] <-(cor.test(merged_seq[1:1000, "MoSplenicDC2_average"], 
                                                             merged_seq[1:1000, "RatSplenicCD103neg_average"]))$estimate

scaled_RatDcMoDcCorrelation <- RatDcMoDcCorrelation
scaled_RatDcMoDcCorrelation["MoThyDC1",] <- scale(as.numeric(RatDcMoDcCorrelation["MoThyDC1",]))
scaled_RatDcMoDcCorrelation["MoThyDC2",] <- scale(as.numeric(RatDcMoDcCorrelation["MoThyDC2",]))
scaled_RatDcMoDcCorrelation["MoThymicCD301bnegDC2",] <- scale(as.numeric(RatDcMoDcCorrelation["MoThymicCD301bnegDC2",]))
scaled_RatDcMoDcCorrelation["MoThymicCD301bposDC2",] <- scale(as.numeric(RatDcMoDcCorrelation["MoThymicCD301bposDC2",]))
scaled_RatDcMoDcCorrelation["MoSplDC1a",] <- scale(as.numeric(RatDcMoDcCorrelation["MoSplDC1a",]))
scaled_RatDcMoDcCorrelation["MoSplDC2",] <- scale(as.numeric(RatDcMoDcCorrelation["MoSplDC2",]))


pheatmap(scaled_RatDcMoDcCorrelation,
         cellwidth = 25, cellheight = 25, cluster_rows = F, cluster_cols = F, 
         display_numbers = apply(RatDcMoDcCorrelation, 2, function(x) sprintf("%.2f", x)))

#merge the objects with all rows----
merged_seq_all <- merge(RatDC_Mf_array, MoThyDC_Seq, by="row.names", all=T)
row.names(merged_seq_all) <- merged_seq_all$Row.names
merged_seq_all<-merged_seq_all[,-1]

merged_seq_all <- merge(merged_seq_all, MoSplDC_seq, by="row.names", all=T)
row.names(merged_seq_all) <- merged_seq_all$Row.names
merged_seq_all<-merged_seq_all[,-1]

#remove NaNs and -Infs
merged_seq_all[merged_seq_all  == "NaN" | merged_seq_all  == "-Inf"] <- "" 
merged_seq_all[] <- lapply(merged_seq_all, function(x) as.numeric(as.character(x)))


#remove batch effects
batch <- c(rep("RatArray", 16), rep("MoThy", 9), rep("MoSpl", 23))
merged_seq_all <- removeBatchEffect(merged_seq_all, batch)

#scaling
merged_seq_all <- scale(merged_seq_all)

#calculate averages
merged_seq_all <- as.data.frame(merged_seq_all) #change class from matrix to data.frame
merged_seq_all[,1:48] <- lapply(merged_seq_all[,1:48], as.numeric) #change class of the contents from character to numeric

merged_seq_all$RatThymicDC1_average <- log2(rowMeans(2^merged_seq_all[, 1:2], na.rm = T))
merged_seq_all$RatSplenicDC1_average <- log2(rowMeans(2^merged_seq_all[, 3:4], na.rm = T))
merged_seq_all$RatThymicDC2_average <- log2(rowMeans(2^merged_seq_all[, 5:6], na.rm = T))
merged_seq_all$RatSplenicDC2_average <- log2(rowMeans(2^merged_seq_all[, 7:8], na.rm = T))
merged_seq_all$RatThymicCD103neg_average <- log2(rowMeans(2^merged_seq_all[, 9:10], na.rm = T))
merged_seq_all$RatSplenicCD103neg_average <- log2(rowMeans(2^merged_seq_all[, 11:12], na.rm = T))
merged_seq_all$RatThymicMf_average <- log2(rowMeans(2^merged_seq_all[, 13:14], na.rm = T))
merged_seq_all$RatSplenicMf_average <- log2(rowMeans(2^merged_seq_all[, 15:16], na.rm = T))

merged_seq_all$MoThymicDC1_average <- log2(rowMeans(2^merged_seq_all[, 23:25], na.rm = T))
merged_seq_all$MoThymicCD301bposDC2_average <- log2(rowMeans(2^merged_seq_all[, 17:19], na.rm = T))
merged_seq_all$MoThymicCD301bnegDC2_average <- log2(rowMeans(2^merged_seq_all[, 20:22], na.rm = T))
merged_seq_all$MoThymicDC2_average <- as.numeric(ifelse(
  is.na(merged_seq_all[, "MoThymicCD301bposDC2_average"]) & is.na(merged_seq_all[, "MoThymicCD301bnegDC2_average"]),"NaN",
  ifelse(
    is.na(merged_seq_all[, "MoThymicCD301bposDC2_average"]),
    log2(0.54*2^merged_seq_all[, "MoThymicCD301bnegDC2_average"]),
    ifelse(
      is.na(merged_seq_all[, "MoThymicCD301bnegDC2_average"]),
      log2(0.46*2^merged_seq_all[, "MoThymicCD301bposDC2_average"]),
      log2(0.46 * 2^merged_seq_all[, "MoThymicCD301bposDC2_average"] + 0.54 * 2^merged_seq_all[, "MoThymicCD301bnegDC2_average"])
    ))))


merged_seq_all$MoSplenicDC1a_average <- log2(rowMeans(2^merged_seq_all[, c(33, 37)], na.rm = T))
merged_seq_all$MoSplenicDC1b_average <- log2(rowMeans(2^merged_seq_all[, c(34, 38)], na.rm = T))
merged_seq_all$MoSplenicDC2_average <- log2(rowMeans(2^merged_seq_all[, c(36, 40)], na.rm = T))

#export the expression matrix
write.csv(merged_seq_all, "PMID_34777358_35637352 mouse cDC1 cDC2 correlation/merged_seq_all.csv", na ="", row.names = TRUE)

merged_seq_all <- read.csv("PMID_34777358_35637352 mouse cDC1 cDC2 correlation/merged_seq_all.csv", row.names = 1, na ="")

#the boxplot of distribution, 4.86 x 8.4 inches.
par(mar = c(14, 4, 4, 2), mgp = c(2, 0.8, 0))  # Increase bottom margin (first value in c()) to fit x-axis labels
boxplot(merged_seq_all[,c("RatThymicDC1_average", "RatThymicDC2_average", "RatThymicCD103neg_average", "RatThymicMf_average",
                          "RatSplenicDC1_average", "RatSplenicDC2_average", "RatSplenicCD103neg_average", "RatSplenicMf_average",
                          "MoThymicDC1_average", "MoThymicDC2_average", "MoSplenicDC1a_average", "MoSplenicDC2_average")],
        main = "Distribution of expression values",
        xlab = "",
        ylab = "removeBatchEffect-ed Log2 expression",
        las = 2,  # Rotate x-axis labels for better readability
        notch = FALSE, 
        cex.main = 1.5,  # Increase the font size of the title
        cex.lab = 1.3)

#depict heatmaps
#DC/macrophage/moDC genes (Fig. 2F)
italic_row_labels <- as.expression(lapply(rownames(merged_seq_all[c("Flt3", "Itgax", "Zbtb46", "Klrb1a; RGD2301395; Klrb1b",
                                                                    "Cd14", "Mafb", "C5ar1", "Csf1r", "Adgre1","Fcgr1a", "Fcgr3a", "Fcgr2a","Fcgr2b","Mertk", "Nos2", "Cx3cr1", "Cd163"),]), function(x) bquote(italic(.(x)))))
pheatmap(merged_seq_all[c("Flt3", "Itgax", "Zbtb46", "Klrb1a; RGD2301395; Klrb1b",
                          "Cd14", "Mafb", "C5ar1", "Csf1r", "Adgre1","Fcgr1a", "Fcgr3a", "Fcgr2a","Fcgr2b","Mertk", "Nos2", "Cx3cr1", "Cd163"), 
                          c("RatSplenicCD103neg_average", "RatSplenicDC2_average", "RatSplenicDC1_average", "RatSplenicMf_average",
                            "RatThymicCD103neg_average", "RatThymicDC2_average", "RatThymicDC1_average", "RatThymicMf_average")],
         cellwidth = 30, cellheight = 20, cluster_rows = F, cluster_cols = F, display_numbers = T,  scale = "row", labels_row = italic_row_labels)

#cDC2+cDC1 genes (Fig. 2E)
italic_row_labels <- as.expression(lapply(rownames(merged_seq_all[c("Itgam", "Clec4a2", "Cx3cr1", "Irf4", "Klf4", "Notch2", "Tbx21", 
                                                                    "Esam", "Clec9a", "Xcr1", "Irf8", "Batf3", "Itgae"
                                                                    ),]), function(x) bquote(italic(.(x)))))

pheatmap(merged_seq_all[c("Itgam", "Clec4a2", "Cx3cr1", "Irf4", "Klf4", "Notch2", "Tbx21", "Esam", "Clec9a", "Xcr1", "Irf8", "Batf3", "Itgae"), 
                        c("RatSplenicCD103neg_average", "RatSplenicDC2_average", "RatSplenicDC1_average",
                          "RatThymicCD103neg_average", "RatThymicDC2_average", "RatThymicDC1_average",
                          "MoSplenicDC2_average", "MoSplenicDC1a_average", "MoThymicDC2_average", "MoThymicDC1_average")],
         cellwidth = 30, cellheight = 15, cluster_rows = F, cluster_cols = F, display_numbers = T, scale="row", labels_row = italic_row_labels)

#monocyte-derived cell genes
pheatmap(merged_seq_all[c("Dpp4", "C5ar1", "Fcgr1a", "Fcgr1", "Fcgr3a", "Fcgr4"), 
                        c("RatThymicDC1_average", "RatThymicDC2_average", "RatThymicCD103neg_average",
                          "RatSplenicDC1_average", "RatSplenicDC2_average", "RatSplenicCD103neg_average",
                          "MoThymicDC1_average", "MoThymicDC2_average", "MoSplenicDC1a_average", "MoSplenicDC2_average")],
         cellwidth = 30, cellheight = 15, cluster_rows = F, cluster_cols = F, display_numbers = T, fontsize_row = 12)

#GSEA analysis between CD103-DC vs cDC2----
cDC1_vs_cDC2_vs_thymus_and_spleen <- as.data.frame(read_delim("PMID_34777358_35637352 mouse cDC1 cDC2 correlation/TAC export/Comparison interaction cDC1 vs cDC2 vs thymus and spleen.txt", 
                                                                       delim = "\t", escape_double = FALSE, trim_ws = TRUE))
CD103_and_cDC2_vs_thymus_and_spleen <- as.data.frame(read_delim("PMID_34777358_35637352 mouse cDC1 cDC2 correlation/TAC export/Comparison interaction CD103- and cDC2 vs thymus and spleen.txt", 
                                                                         delim = "\t", escape_double = FALSE, trim_ws = TRUE))

rownames(cDC1_vs_cDC2_vs_thymus_and_spleen)<-cDC1_vs_cDC2_vs_thymus_and_spleen[,1]
rownames(CD103_and_cDC2_vs_thymus_and_spleen)<-CD103_and_cDC2_vs_thymus_and_spleen[,1]

RatDC_Mf_array_GESA<- merge(cDC1_vs_cDC2_vs_thymus_and_spleen[, 2:6], CD103_and_cDC2_vs_thymus_and_spleen[,c(3,5)], by = 'row.names', all = TRUE)

gene_symbols <- RatDC_Mf_array_GESA$`Gene Symbol` 

dup_indices <- duplicated(gene_symbols) | duplicated(gene_symbols, fromLast = TRUE) # Find duplicated values

count <- ave(seq_along(gene_symbols), gene_symbols, FUN = seq_along) # Count the number of occurrences of each value

gene_symbols[dup_indices] <- paste0(gene_symbols[dup_indices], "-", count[dup_indices] - 1) # Append unique numbers to duplicates

gene_symbols <- gsub("-0", "", gene_symbols) #remove "-0"

row.names(RatDC_Mf_array_GESA) <- gene_symbols #Set gene symbols to row names

RatDC_Mf_array_GESA <- RatDC_Mf_array_GESA[, -(1)] #remove unnecessary columns

#calculate Log2FC values
Rat_Thy_CD103neg_vs_cDC2 <-(RatDC_Mf_array_GESA$`thymus: CD103- DC-like Avg (log2)` - RatDC_Mf_array_GESA$`thymus: cDC2 Avg (log2)`)
names(Rat_Thy_CD103neg_vs_cDC2) <- row.names(RatDC_Mf_array_GESA)
Rat_Thy_CD103neg_vs_cDC2<- sort(Rat_Thy_CD103neg_vs_cDC2, decreasing = T)

Rat_Spl_CD103neg_vs_cDC2 <-(RatDC_Mf_array_GESA$`spleen: CD103- DC-like Avg (log2)` - RatDC_Mf_array_GESA$`spleen: cDC2 Avg (log2)`)
names(Rat_Spl_CD103neg_vs_cDC2) <- row.names(RatDC_Mf_array_GESA)
Rat_Spl_CD103neg_vs_cDC2<- sort(Rat_Spl_CD103neg_vs_cDC2, decreasing = T)

Rat_ThyCD103neg_vs_SplCD103neg <-(RatDC_Mf_array_GESA$`thymus: CD103- DC-like Avg (log2)` - RatDC_Mf_array_GESA$`spleen: CD103- DC-like Avg (log2)`)
names(Rat_ThyCD103neg_vs_SplCD103neg) <- row.names(RatDC_Mf_array_GESA)
Rat_ThyCD103neg_vs_SplCD103neg<- sort(Rat_ThyCD103neg_vs_SplCD103neg, decreasing = T)

#analyse the Log2FC values with GSEA
Rat_Thy_CD103neg_vs_cDC2_GSEA <- gseGO(geneList = Rat_Thy_CD103neg_vs_cDC2, OrgDb = "org.Rn.eg.db", ont = "BP", keyType = "SYMBOL", pAdjustMethod="none")
Rat_Spl_CD103neg_vs_cDC2_GSEA <- gseGO(geneList = Rat_Spl_CD103neg_vs_cDC2, OrgDb = "org.Rn.eg.db", ont = "BP", keyType = "SYMBOL", pAdjustMethod="none")
Rat_ThyCD103neg_vs_SplCD103neg_GSEA <- gseGO(geneList = Rat_ThyCD103neg_vs_SplCD103neg, OrgDb = "org.Rn.eg.db", ont = "BP", keyType = "SYMBOL", pAdjustMethod="none")

#split GSEA results to +/-enriched gene sets
Rat_Thy_CD103negDom <- Rat_Thy_CD103neg_vs_cDC2_GSEA
Rat_Thy_CD103negDom@result <- Rat_Thy_CD103negDom@result[Rat_Thy_CD103negDom@result$enrichmentScore >0, ]

Rat_Thy_cDC2Dom <- Rat_Thy_CD103neg_vs_cDC2_GSEA
Rat_Thy_cDC2Dom@result <- subset(Rat_Thy_CD103neg_vs_cDC2_GSEA, Rat_Thy_CD103neg_vs_cDC2_GSEA@result$enrichmentScore <0)

Rat_Spl_CD103negDom <- Rat_Spl_CD103neg_vs_cDC2_GSEA
Rat_Spl_CD103negDom@result <- subset(Rat_Spl_CD103neg_vs_cDC2_GSEA, Rat_Spl_CD103neg_vs_cDC2_GSEA@result$enrichmentScore >0)

Rat_Spl_cDC2Dom  <-Rat_Spl_CD103neg_vs_cDC2_GSEA
Rat_Spl_cDC2Dom@result <- subset(Rat_Spl_CD103neg_vs_cDC2_GSEA, Rat_Spl_CD103neg_vs_cDC2_GSEA@result$enrichmentScore <0)

ThyCD103negDom <- Rat_ThyCD103neg_vs_SplCD103neg_GSEA
ThyCD103negDom@result <- subset(ThyCD103negDom, ThyCD103negDom@result$enrichmentScore>0)

SplCD103negDom <- Rat_ThyCD103neg_vs_SplCD103neg_GSEA
SplCD103negDom@result <- subset(SplCD103negDom, SplCD103negDom@result$enrichmentScore<0)

#split CD103neg GSEA results to common and unique gene sets between thymus and spleen
Rat_Thy_CD103negDomCommon <- Rat_Thy_CD103negDom
Rat_Thy_CD103negDomCommon@result <- Rat_Thy_CD103negDomCommon@result[Rat_Thy_CD103negDomCommon@result$ID %in%
                                                                       intersect(Rat_Thy_CD103negDom@result$ID, Rat_Spl_CD103negDom@result$ID),]

Rat_Thy_CD103negDomUnique <- Rat_Thy_CD103negDom
Rat_Thy_CD103negDomUnique@result <- Rat_Thy_CD103negDomUnique@result[Rat_Thy_CD103negDomUnique@result$ID %in%
                                                                       setdiff(Rat_Thy_CD103negDom@result$ID, Rat_Spl_CD103negDom@result$ID),]

Rat_Spl_CD103negDomCommon <- Rat_Spl_CD103negDom
Rat_Spl_CD103negDomCommon@result <- Rat_Spl_CD103negDomCommon@result[Rat_Spl_CD103negDomCommon@result$ID %in%
                                                                       intersect(Rat_Thy_CD103negDom@result$ID, Rat_Spl_CD103negDom@result$ID),]

Rat_Spl_CD103negDomUnique <- Rat_Spl_CD103negDom
Rat_Spl_CD103negDomUnique@result <- Rat_Spl_CD103negDomUnique@result[Rat_Spl_CD103negDomUnique@result$ID %in%
                                                                       setdiff(Rat_Spl_CD103negDom@result$ID, Rat_Thy_CD103negDom@result$ID),]

#Treeplots
treeplot(pairwise_termsim(Rat_Thy_CD103negDom), showCategory = 30, color = "p.adjust", #Fig. 3C, 9.65x5 inches
         offset.params = list(bar_tree = rel(5), extend=0.5, hexpand = 0.15), hilight.params = list(hilight = T),
         fontsize = 4, label_format = 20, cluster.params = list(label_format = 10), label_format_tiplab = 100, 
        cex_category = 0.75 )

treeplot(pairwise_termsim(Rat_Thy_cDC2Dom),  showCategory = 30, color = "p.adjust", #Fig. 3C, 9.65x5 inches
         offset.params = list(bar_tree = rel(5), extend=0.5, hexpand = 0.15), hilight.params = list(hilight = T),
         fontsize = 4, label_format = 25, cluster.params = list(label_format = 10), label_format_tiplab = 100, 
         cex_category = 0.75 )

treeplot(pairwise_termsim(Rat_Spl_CD103negDom),  showCategory = 30, color = "p.adjust", #Fig. 3A, 9.65x5 inches
         offset.params = list(bar_tree = rel(5), extend=0.5, hexpand = 0.15), hilight.params = list(hilight = T),
         fontsize = 4, label_format = 20, cluster.params = list(label_format = 10), label_format_tiplab = 100, 
         cex_category = 0.75 )

treeplot(pairwise_termsim(Rat_Spl_cDC2Dom),  showCategory = 30, color = "p.adjust", #Fig. 3A, 9.65x5 inches
         offset.params = list(bar_tree = rel(5), extend=0.5, hexpand = 0.15), hilight.params = list(hilight = T),
         fontsize = 4, label_format = 25, cluster.params = list(label_format = 10), label_format_tiplab = 100, 
         cex_category = 0.75 )

treeplot(pairwise_termsim(ThyCD103negDom),  showCategory = 30, color = "p.adjust", # 9.65x5 inches
         offset.params = list(bar_tree = rel(5), extend=0.5, hexpand = 0.3), hilight.params = list(hilight = T),
         fontsize = 4, label_format = 30, cluster.params = list(label_format = 10), label_format_tiplab = 100, 
         cex_category = 0.75 )

treeplot(pairwise_termsim(SplCD103negDom), showCategory = 30, color = "p.adjust", # 9.65x5 inches
         offset.params = list(bar_tree = rel(7), extend=0.5, hexpand = 0.2), hilight.params = list(hilight = T),
         fontsize = 4, label_format = 25, cluster.params = list(label_format = 10), label_format_tiplab = 100, 
         cex_category = 0.75 )

#GSEA plots
gseaplot2(Rat_Thy_CD103neg_vs_cDC2_GSEA, #Fig. 3D, 8.2x3.5 inches
          geneSetID = c("GO:0006909", "GO:0002253", "GO:0001819"),
          subplots = 1:2, pvalue_table = F, base_size = 14)

gseaplot2(Rat_Thy_CD103neg_vs_cDC2_GSEA, #Fig. 3D, 8x4 inches
          geneSetID = c("GO:0006260", "GO:0007059", "GO:0044839"),
          subplots = 1:2, pvalue_table = F, base_size = 14)

gseaplot2(Rat_Spl_CD103neg_vs_cDC2_GSEA, #Fig. 3B, 8x4 inches
          geneSetID = c("GO:0006909", "GO:0002253", "GO:0001819"),
          subplots = 1:2, pvalue_table = F, base_size = 14)

gseaplot2(Rat_Spl_CD103neg_vs_cDC2_GSEA, #Fig. 3B, 8x4 inches
          geneSetID = c("GO:0006260", "GO:0007059", "GO:0044839"),
          subplots = 1:2, pvalue_table = F, base_size = 14)

gseaplot2(Rat_ThyCD103neg_vs_SplCD103neg_GSEA,  # 8x4 inches
          geneSetID = c("GO:0009060", "GO:0022900", "GO:0032963"),
          subplots = 1:2, pvalue_table = F, base_size = 14)

gseaplot2(Rat_ThyCD103neg_vs_SplCD103neg_GSEA,  # 8x4 inches
          geneSetID = c("GO:0045088", "GO:0002221", "GO:0002833"),
          subplots = 1:2, pvalue_table = F, base_size = 14)

#save gene_list
write.table(c(Rat_Thy_CD103neg_vs_cDC2_GSEA@geneSets$`GO:0006909`,Rat_Thy_CD103neg_vs_cDC2_GSEA@geneSets$`GO:0002253`,
              Rat_Thy_CD103neg_vs_cDC2_GSEA@geneSets$`GO:0001819`),
            file = "PMID_34777358_35637352 mouse cDC1 cDC2 correlation/immune_associated.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

write.table(c(Rat_Thy_CD103neg_vs_cDC2_GSEA@geneSets$`GO:0006260`,Rat_Thy_CD103neg_vs_cDC2_GSEA@geneSets$`GO:0007059`,
              Rat_Thy_CD103neg_vs_cDC2_GSEA@geneSets$`GO:0044839`),
            file = "PMID_34777358_35637352 mouse cDC1 cDC2 correlation/cycle_associated.gene_list.txt",
            col.names="gene_symbol", row.names=F, quote = F)

#save GSEA results
write.xlsx(list("Rat_Thy_CD103neg_vs_cDC2"=Rat_Thy_CD103neg_vs_cDC2_GSEA, "Rat_Spl_CD103neg_vs_cDC2"=Rat_Spl_CD103neg_vs_cDC2_GSEA,
                "Rat_ThyCD103neg_vs_SplCD103neg"=Rat_ThyCD103neg_vs_SplCD103neg_GSEA), 
           "PMID_34777358_35637352 mouse cDC1 cDC2 correlation/CD103neg_vs_cDC2_GSEA.xlsx")