#####
#COVID-19 snRNAseq data analysis
library(dplyr)
library(patchwork)
library("Seurat") 
library(openxlsx)
#load the data
data <- readRDS(file = "~/CP_COVID.rds")
#Chose the "RNA" assay
DefaultAssay(data) <- "RNA" 
#Exclude the cells from a flu patient
Idents(data) <- data$orig.ident
table(data$orig.ident)
data <- subset(data, cells = WhichCells(data, idents = c("ct_15_12_CP","ct_15_17_CP","ct_16_18_CP"
                                                         ,"ct_16_23_CP","ct_16_24_CP","ct_16_25_CP",
                                                         "cv_72_CP","cv_73_CP","cv_75_CP","cv_76_CP",
                                                         "cv_78_CP","cv_90_CP","cv_91_CP")))
#Select the cells from healthy control
Control <- subset(data, cells = WhichCells(data, idents = c("ct_15_12_CP","ct_15_17_CP","ct_16_18_CP"
                                                            ,"ct_16_23_CP","ct_16_24_CP","ct_16_25_CP")))
#Show the UMAP of CYP46A1 expressing cells in healthy control
FeaturePlot(Control, features = c("CYP46A1"),pt.size = 0.1, order = TRUE,
            cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) ))
#Show the UMAP of each cell type
Idents(Control) <- Control$cellID
DimPlot(Control, reduction = "umap", label = FALSE,pt.size = 0.1)
#Dot plot of each cell type without scaling
DotPlot(Control, features = c("CYP46A1"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
DotPlot(Control, features = c("NR1H2","NR1H3"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
DotPlot(Control, features = c("TNFRSF1A","TNFRSF1B"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
#Chose the Epthelial cells and make the Dot plot in Epthelial cells type without scaling
Epithelial <- subset(Control, cells = WhichCells(Control, idents = c("Epithelial")))
DotPlot(Epithelial, features = c("CYP7A1","CYP46A1","CYP27A1","CYP11A1","CH25H"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
#Find the DEGs
Idents(data) <- data$Biogroup
DEG <- FindMarkers(data, ident.1 = "Case",logfc.threshold = 0,test.use = "MAST", 
                   latent.vars = c("Sex","Batch"))
#Export to Excel file
write.xlsx(DEG, col_names = TRUE,format_headers = TRUE, rowNames = TRUE, "~/DEG.xlsx")

#####
#Mouse aging snRNAseq data analysis
library(dplyr)
library(patchwork)
library(ggplot2)
library(Seurat)
library(openxlsx)
library("RColorBrewer")
#load the data
CP <- Read10X(data.dir = "~/Desktop/Mouse CP")
colnames(CP)
colname <- colnames(CP)
write.csv(colname, file="~/Desktop/colname.csv")
rownames(CP)
#select the cells from Adult and Aged cells
adult <- CP[,14560:50912]
aged <- CP[,50913:83040]
colnames(adult)
rownames(adult)
colnames(aged)
rownames(aged)
Adult <- CreateSeuratObject(counts = adult, project = "Adult", min.cells = 3, min.features = 200)
Aged <- CreateSeuratObject(counts = aged, project = "Aged", min.cells = 3, min.features = 200)
#merge files
data <- merge(Adult, y = Aged, add.cell.ids = c("Adult", "Aged"), project = "data")
#Making new Metadata
sapply(X = strsplit(colnames(data), split = "_"), FUN = "[", 1)
new.cluster.ids <- sapply(X = strsplit(colnames(data), split = "_"), FUN = "[", 1)
data <- AddMetaData(object = data,metadata = new.cluster.ids,col.name = 'Sample')
head(x = data[[]])
#Removing low-quality cells
data.genes <- grep(pattern = "^mt-", x = rownames(GetAssayData(data)), value = TRUE)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
# Normalizing the data
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
# Identification of highly variable features (feature selection)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data), 10)
# Scaling the data
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
# Perform linear dimensional reduction(PCA)
data <- RunPCA(data, features = VariableFeatures(object = data))
# Examine and visualize PCA results a few different ways
print(data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(data, dims = 1:2, reduction = "pca")
DimPlot(data, reduction = "pca")
# Determine the ‘dimensionality’ of the dataset
ElbowPlot(data)

# Cluster the cells
data <- FindNeighbors(data, dims = 1:12)
data <- FindClusters(data, resolution = 0.5)
data <- RunUMAP(data, dims = 1:12)
DimPlot(data, reduction = "umap", label = FALSE,pt.size = 0.1)
table(data@active.ident)
#Show the UMAP grouped by CP regions
DimPlot(data, reduction = "umap", label = FALSE, group.by = "orig.ident",
        cols = c(rgb(202,67,54,max=255), rgb(93,117,162,max=255), rgb(108,162,77,max=255)),pt.size = 0.1)
#Show the UMAP grouped by "Adult" or "Aged"
DimPlot(data, reduction = "umap", label = FALSE, group.by = "Sample",
        cols = c(rgb(202,67,54,max=255), rgb(93,117,162,max=255)),pt.size = 0.1) 
# Determined the cluster based on the original paper, Dani et. al.
Idents(data) <- data$seurat_clusters
new.cluster.ids <- c("Epithelial", "Epithelial", "Epithelial", "Epithelial","Epithelial", "Mesenchymal","Epithelial",
                     "Epithelial","Epithelial","Epithelial","Mesenchymal","Endothelial","Immune",
                     "Neuronal/Glial","Epithelial")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
#Show the UMAP grouped by cell types
DimPlot(data, reduction = "umap", label = FALSE,pt.size = 0.1)

#Show the UMAP of CYP46A1 expressing cells in healthy control
FeaturePlot(data, features = c("Cyp46a1"),pt.size = 1, order = TRUE,
            cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) ))

#Dot plot of each cell type without scaling
DotPlot(data, features = c("Cyp46a1"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
DotPlot(data, features = c("Nr1h2","Nr1h3"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
DotPlot(data, features = c("Tnfrsf1a","Tnfrsf1b"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
Epithelial <- subset(data, cells = WhichCells(data, idents = c("Epithelial")))
DotPlot(Epithelial, features = c("Cyp7a1","Cyp46a1","Cyp27a1","Cyp11a1","Ch25h"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
#Find the DEGs
Idents(data) <- data$Sample
CP_DEG <- FindMarkers(data, ident.1 = "Aged",ident.2 = "Adult",logfc.threshold = 0,test.use = "MAST")
write.xlsx(CP_DEG, col_names = TRUE,format_headers = TRUE, rowNames = TRUE, "~/CP_DEG.xlsx")

#####
#Dot plots for RNAseq datasets
library("ggpubr")
library(ggplot2)
library("ggrepel")
library(tidyverse)
data <- as.matrix(read.csv(file = "~/data.csv",header=T, row.names=1))
data2 <- data.frame(data,face = c("italic"))
p <- as.matrix(read.csv(file = "~/data with p-value.csv",header=T, row.names=1))
p2 <- data.frame(p,face = c("italic"))
#select the genes which we want to show
names <- rownames(data2)
selected_names <- c("SERPINA3")
names[!names %in% selected_names] <- ""
#Select genes which have significance
data2 <- p2 %>%  
  mutate(
    Expression = case_when(p < 0.05 ~ "Significant",
                           p > 0.05 ~ "Non-Significant"))
#Make the plot
ggplot(data2, aes(x, y, label = names, fontface = face)) +
  geom_point(aes(color = "black"), size = 1.5) +
  xlab(expression("24OH (log"[2]*"FC)")) + 
  ylab(expression("Disease (log"[2]*"FC)")) +
  scale_color_manual(values = c("black", "black", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1))) +
  theme_minimal() +
  geom_text_repel() +
  theme(panel.background = element_rect(fill = "transparent", colour = "black", size = 1.2),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title =element_text(size=15,  family="sans", color="black"),
        axis.text.x =element_text(size=12,  family="sans", color="black"),
        axis.text.y =element_text(size=12,  family="sans", color="black"),
        axis.ticks = element_line(colour = "black", size = 0.7)) 
#####
#Volcano plot
library(tidyverse)
library("ggplot2")
library("ggrepel")
data <- as.matrix.data.frame(read.csv(file = "~/volcano.csv",header=T,row.names=1))
data2 <- data.frame(data,face = c("italic"))
#Select genes which have significance
data2 <- data2 %>% 
  mutate(
    Expression = case_when(x >= 1 & y >= 1.3 ~ "Up-regulated",
                           x <= -1 & y >= 1.3 ~ "Down-regulated",
                           TRUE ~ "No-significant")
  )
names <- rownames(data2)
#Make the plot
ggplot(data2, aes(x, y, label = names, fontface = face)) +
  geom_point(aes(color = Expression), size = 1.5) +
  xlab(expression("24OH (log"[2]*"FC)")) + 
  ylab(expression("-log"[10]*italic("P")*"adj")) +
  scale_color_manual(values = c("#FFA500", "gray80", "#8B0000")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  theme_minimal()+
  geom_text_repel()+ 
  theme(panel.background = element_rect(fill = "transparent", colour = "black", size = 1.2),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title =element_text(size=15,  family="sans", color="black"),
        axis.text.x =element_text(size=12,  family="sans", color="black"),
        axis.text.y =element_text(size=12,  family="sans", color="black"),
        axis.ticks = element_line(colour = "black", size = 0.7))
