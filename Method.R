#####
#COVID-19 snRNAseq data analysis
library(dplyr)
library(patchwork)
library("Seurat") 
library(openxlsx)
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
#Select the samples from healty control
Control <- subset(data, cells = WhichCells(data, idents = c("ct_15_12_CP","ct_15_17_CP","ct_16_18_CP"
                                                            ,"ct_16_23_CP","ct_16_24_CP","ct_16_25_CP")))

#Separete to each cell type
Idents(data) <- data$cellID
Epithelial <- subset(data, cells = WhichCells(data, idents = c("Epithelial")))
#Find the DEGs
Idents(data) <- data$Biogroup
DEG <- FindMarkers(data, ident.1 = "Case",logfc.threshold = 0,test.use = "MAST", 
                   latent.vars = c("Sex","Batch"))
#Export to Excel file
write.xlsx(DEG, col_names = TRUE,format_headers = TRUE, rowNames = TRUE, "~/DEG.xlsx")
#Plot the gene distributions
FeaturePlot(data, features = c("CYP46A1"),pt.size = 1, order = TRUE,
            cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) ))
FeaturePlot(data, features = c("IFNAR1","IFNAR2"),pt.size = 1, order = TRUE,
            cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) ))
FeaturePlot(data, features = c("IFNGR1","IFNGR2"),pt.size = 1, order = TRUE)
#Violin plot of each cell type
VlnPlot(data, features = c("CYP46A1"), split.by = "Biogroup")
#Bubble plot of each cell type without scaling
Idents(Control) <- Control$cellID
DotPlot(Control, features = c("CYP46A1"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
DotPlot(Control, features = c("IFNAR1","IFNAR2"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
DotPlot(Control, features = c("IFNGR1","IFNGR2")
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
DotPlot(Control, features = c("NR1H2","NR1H3"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
DotPlot(Control, features = c("TNFRSF1A","TNFRSF1B"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
Epithelial <- subset(Control, cells = WhichCells(Control, idents = c("Epithelial")))
DotPlot(Epithelial, features = c("CYP7A1","CYP46A1","CYP27A1","CYP11A1","CH25H"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")

#####
#Mouse aging snRNAseq data analysis
library(dplyr)
library(patchwork)
library(ggplot2)
library(Seurat)
library(openxlsx)
library("RColorBrewer")
CP <- Read10X(data.dir = "~/Desktop/Mouse CP")
colnames(CP)
rownames(CP)
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
sapply(X = strsplit(colnames(data), split = "_"), FUN = "[", 1)#barcodeの一番最初の部分のみ抽出
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
# a standard pre-processing step prior to dimensional reduction techniques like PCA.
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
# Perform linear dimensional reduction(PCA)
data <- RunPCA(data, features = VariableFeatures(object = data))
# Examine and visualize PCA results a few different ways
print(data[["pca"]], dims = 1:5, nfeatures = 5)
# 3 ways to choose PC number for clustering
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
DimPlot(data, reduction = "umap", label = FALSE, group.by = "orig.ident",
        cols = c(rgb(202,67,54,max=255), rgb(93,117,162,max=255), rgb(108,162,77,max=255)),pt.size = 0.1)
DimPlot(data, reduction = "umap", label = FALSE, group.by = "Sample",
        cols = c(rgb(202,67,54,max=255), rgb(93,117,162,max=255)),pt.size = 0.1) 
# Determined the cluster based on the original paper, Dani et. al.
Idents(data) <- data$seurat_clusters
new.cluster.ids <- c("Epithelial", "Epithelial", "Epithelial", "Epithelial","Epithelial", "Mesenchymal","Epithelial",
                     "Epithelial","Epithelial","Epithelial","Mesenchymal","Endothelial","Immune",
                     "Neuronal/Glial","Epithelial")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
#Plot the gene distributions
FeaturePlot(data, features = c("Cyp46a1"),pt.size = 1, order = TRUE,
            cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) ))
FeaturePlot(data, features = c("Ifnar1","Ifnar2"),pt.size = 1, order = TRUE,
            cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) ))
FeaturePlot(data, features = c("Ifngr1","Ifngr2"),pt.size = 1, order = TRUE)
#Bubble plot of Adult sample
Idents(Adult) <- Adult$seurat_clusters
new.cluster.ids <- c("Epithelial", "Epithelial", "Epithelial", "Epithelial","Epithelial", "Mesenchymal","Epithelial",
                     "Epithelial","Epithelial","Epithelial","Mesenchymal","Endothelial","Immune",
                     "Neuronal/Glial","Epithelial")
names(new.cluster.ids) <- levels(Adult)
Adult <- RenameIdents(Adult, new.cluster.ids)
DotPlot(Adult, features = c("Cyp46a1"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
DotPlot(Adult, features = c("Ifnar1","Ifnar2"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
DotPlot(Adult, features = c("Ifngr1","Ifngr2")
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
DotPlot(Adult, features = c("Nr1h2","Nr1h3"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
DotPlot(Adult, features = c("Tnfrsf1a","Tnfrsf1b"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
Epithelial <- subset(Adult, cells = WhichCells(Adult, idents = c("Epithelial")))
DotPlot(Epithelial, features = c("Cyp7a1","Cyp46a1","Cyp27a1","Cyp11a1","Ch25h"),cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE, scale.by = "radius")
#Violin plot of each sample
VlnPlot(data, features = c("Cyp46a1"), split.by = "Sample")
#Find the DEGs
Idents(data) <- data$Sample
CP_DEG <- FindMarkers(data, ident.1 = "Aged",ident.2 = "Adult",logfc.threshold = 0,test.use = "MAST")
write.xlsx(CP_DEG, col_names = TRUE,format_headers = TRUE, rowNames = TRUE, "~/CP_DEG.xlsx")
#Extract epithelial cells
Idents(data) <- data@active.ident
Epi <- subset(data, cells = WhichCells(data,idents = c("Epithelial")))
DotPlot(Epi, features = c("Cyp7a1","Cyp46a1","Cyp27a1","Cyp11a1","Ch25h"),
        cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )
        ,dot.min = 0,scale.min = 0,col.min=0,scale = FALSE)
VlnPlot(Epi, features = c("Cyp46a1"), split.by = "Sample")
Idents(Epi) <- Epi$Sample
Epi_DEG <- FindMarkers(Epi, ident.1 = "Aged",ident.2 = "Adult",logfc.threshold = 0,test.use = "MAST")
write.xlsx(Epi_DEG, col_names = TRUE,format_headers = TRUE, rowNames = TRUE, "~/Desktop/Epi_DEG.xlsx")
#plot UMAP
DimPlot(data, reduction = "umap", label = FALSE, group.by = "orig.ident",
        cols = c(rgb(110,142,247,max=255), rgb(202,58,126,max=255), rgb(244,179,62,max=255)),pt.size = 0.1)
DimPlot(data, reduction = "umap", label = FALSE, group.by = "Sample",
        cols = c(rgb(110,142,247,max=255), rgb(244,179,62,max=255)),pt.size = 0.1)
DimPlot(data, reduction = "umap", label = FALSE, group.by = "orig.ident",
        cols = c(rgb(110,142,247,max=255), rgb(202,58,126,max=255), rgb(244,179,62,max=255)),pt.size = 0.1)
DimPlot(data, reduction = "umap", label = FALSE,
        cols = c(rgb(112,178,228,max=255), rgb(211,162,56,max=255), rgb(70,156,118,max=255), 
                 rgb(48,112,173,max=255), rgb(193,125,165,max=255)),
        pt.size = 0.1)

#####
#Venn plot
library(VennDiagram)
library(tidyverse)
AD <- as.matrix(read.csv(file="~/AD.csv",header=F))
FTD <- as.matrix(read.csv(file="~/FTD.csv",header=F))
HD <- as.matrix(read.csv(file="~/HD.csv",header=F))
MS <- as.matrix(read.csv(file="~/PMS.csv",header=F))
COVID <- as.matrix(read.csv(file="~/COVID.csv",header=F))
Aging <- as.matrix(read.csv(file="~/mAging.csv",header=F))
SCZ <- as.matrix(read.csv(file="~/SCZ.csv",header=F))
A <- intersect(intersect(AD, FTD), HD)
#Venn plot for human and mouse merged data
venn.diagram(x=list(A,Aging), 
             filename='~/Human+Mouse.png',
             category.names=c("AD-FTD-HD","mAging"),
             output=TRUE,col=c("#440154ff", '#21908dff'),
             fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
             cex = 2,
             fontfamily = "sans",
             cat.cex = 0)
#Venn plot for human data
venn.diagram(x=list(AD,MS,COVID,SCZ), 
             filename='~/Human.png',
             category.names=c("AD-FTD-HD","MS","COVID","SCZ"),
             output=TRUE,col=c("#440154ff", '#21908dff', '#fde725ff','#5dc963ff'),
             fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3), alpha('#5dc963ff',0.3)),
             cex = 2,
             fontfamily = "sans",
             cat.cex = 0)

#####
#Heatmap
library(pheatmap)
library("RColorBrewer")
data <- as.matrix(read.csv(file="~/heatmap.csv",header=T, row.names=1))
out <- pheatmap(data,color = colorRampPalette(brewer.pal(10, "YlOrBr"))(256),
                show_rownames=F, cluster_cols=T, cluster_rows=T,scale="row",
                clustering_distance_rows="euclidean",
                clustering_distance_cols="euclidean", 
                clustering_method="complete", border_color=FALSE)

#####
#Bubble plot of GO terms from Metascape
library("BiocParallel")
library(ggplot2)
library(Hmisc)
library(stringr)
library(ggpubr)
library("RColorBrewer")
#Prepare data for plotting GO terms from metascape results sheet
eval_vec=Vectorize(eval, vectorize.args="expr") 
data=read.delim("~/data.txt",header = T,stringsAsFactors = F)
data=data[grep("*Summary*",data$GroupID),]
data$LogP_pos=abs(data$LogP)
data$ratio=eval_vec(parse(text = data$InTerm_InList))
data$counts=as.numeric(lapply(data$InTerm_InList,function(x) {unlist(strsplit(x,"/"))[1]}))
data$Description=capitalize(data$Description)
#Make the plot
ggplot(data=data[1:15,], aes(x=LogP_pos, y=reorder(Description,LogP_pos),size=ratio,color=counts))+geom_point()+
  scale_size(range = c(0.5, 3.5), name="Gene Ratio",breaks = c(0.1,0.15,0.2))+
  theme(axis.text.x=element_text(size=6.5,  family="sans",color="black"),axis.text.y=element_text(size=6.5,  family="sans",color="black"), 
        axis.title =element_text(size=6.5,  family="sans", color="black"),axis.line = element_line(colour = 'black', size = 0.46),axis.ticks.length=unit(.15, "cm"),
        axis.ticks =element_line(color="black",size=0.46),axis.title.y=element_blank(),
        panel.grid.major=element_line(color = "grey90",size=0.46),panel.background = element_rect(fill = "transparent"),
        legend.title = element_text(size=6.5),legend.key.width = unit(0.5,"cm"),legend.text =element_text(size=6),legend.key.height = unit(0.1,"cm"),legend.position = "bottom",legend.direction = "horizontal",legend.box="horizontal",legend.spacing.x = unit(0.08,"cm"),legend.box.margin=margin(-10, 0, 0, -100))+
  labs(col=expression("Counts"))+xlab(expression("-log"[10]*italic("P")))+
  scale_color_gradient(low="#F4B570", high = "#F06B37")+guides(color=guide_colorbar(order=1),size=guide_legend(order=2))

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

#####
#DESeq2
library( "DESeq2" )
library(ggplot2)
#Before DESeq, I deleted duplicated genes from data. Chose the most abundant one.
countData <- read.csv('~/countdata.csv', header = TRUE)
metaData <- read.csv('~/metadata.csv', header = TRUE)
dds <- DESeqDataSetFromMatrix(countData=countData,colData=metaData,design=~Condition, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
write.csv(res, file="~/result.csv")
