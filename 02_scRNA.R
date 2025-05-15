rm(list=ls())
library(GEOquery)
 
if (!dir.exists("01_scRNA")) {dir.create("01_scRNA")}
 
options(stringsAsFactors = F)
gc()
library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(cowplot)
library(dplyr)
library(data.table)
library(SeuratDisk)
file <- list.files('../00_Rawdata/single_cell/')  
file   
metadata <- read.table('../00_Rawdata/02.metadata.txt',header = T,sep = '\t',
                       check.names = F)
scRNAlist <- list()   
for(i in 1:length(file)){
  rt=fread(paste("../00_Rawdata/single_cell/",file[i],sep='/'),header = TRUE) 
  length(rt$Gene)   
  length(unique(rt$Gene))
  hg=rt$Gene
  dat=rt[,2:ncol(rt)]  
  rownames(dat)=hg  
  hg[grepl('^MT-',hg)]  
  colnames(dat)  
  rownames(dat)  
  class(dat)
  project <- str_split((str_split(file[i],'[.]')[[1]]),'/')[[1]]   
  
  
  
  
  scRNAlist[[i]] <- CreateSeuratObject(counts = dat,
                                       min.cells = 3, min.features = 200,project=project)   
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = project)
  
  meta <- metadata[i,]
  n_cell <- length(colnames(scRNAlist[[i]]))  
  Type <- rep(x =meta[1,2],n_cell)
  Patients_ID <- rep(x =meta[1,1],n_cell)
  scRNAlist[[i]] $Type <- Type
  scRNAlist[[i]] $Patients_ID <- Patients_ID
  
  if(T){    
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") 
  }
  
  if(T){
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]")
  }
}
samples_name <- sapply(str_split(file,'_'),'[',1)   
names(scRNAlist) <- samples_name
system.time(save(scRNAlist, file = "01.scRNAlist.Rdata"))   
scRNA=merge(x=scRNAlist[[1]],
            scRNAlist[2:length(scRNAlist)])
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")  
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[SL]")
head(scRNA@meta.data)
system.time(save(scRNA, file = "01.scRNA_orig.Rdata")) 
rm(list = ls())
if (! dir.exists("./01_QC")){
  dir.create("./01_QC")
}
load('../01.scRNA_orig.Rdata')
Screening_res <- data.frame()
temp1 <- data.frame(cell_number = length(colnames(scRNA)),
                    gene_number = length(rownames(scRNA)))  
Screening_res <- rbind(Screening_res, temp1)
rownames(Screening_res) <- 'before_quality'
theme.set = theme(
  
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 12,  family = "Times"),
  axis.text.y = element_text(size = 14,  family = "Times"),
  legend.text = element_text(size = 14, family = "Times"),
  legend.title = element_text(size = 16,face='bold',family = "Times"),
  text = element_text(family = "Times"))
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")   
group = "Patients_ID"
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=1,ncol = 3)    
ggsave("01.vlnplot_before_qc.pdf", plot = violin, width = 18, height = 6)
ggsave("01.vlnplot_before_qc.png", plot = violin, width = 18, height = 6)
length(colnames(scRNA))   
length(rownames(scRNA))  
library(showtext)
pctMT=10
minGene=200
maxGene=3000
maxUMI=10000
minUMI=500
scRNA <- subset(scRNA, subset =
                  nCount_RNA > minUMI &
                  nCount_RNA < maxUMI & 
                  nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene& 
                  percent.mt < pctMT)
length(colnames(scRNA))   
length(rownames(scRNA))  
table(scRNA$orig.ident)
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")
theme.set = theme(
  
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 12,  family = "Times"),
  axis.text.y = element_text(size = 14,  family = "Times"),
  legend.text = element_text(size = 14, family = "Times"),
  legend.title = element_text(size = 16,face='bold',family = "Times"),
  text = element_text(family = "Times"))
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")   
group = "Patients_ID"
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=1,ncol = 3)  
ggsave(filename = "02.vlnplot_after_qc.pdf", plot = violin, width = 18, height = 6)
ggsave(filename = "02.vlnplot_after_qc.png", plot = violin, width = 18, height = 6)
system.time(save(scRNA, file = "scRNA_qc.Rdata"))
rm(list = ls())
if (! dir.exists("./02_Integrate")){
  dir.create("./02_Integrate")
}
load('../01_QC/scRNA_qc.Rdata')
library(Seurat)
library(tidyverse)
combined <- NormalizeData(scRNA)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)
library(ggplot2)
library(Seurat)
theme_set(theme_minimal(base_family = "Times"))
pdf("01.before_intergrated_PCA.pdf", width = 7, height = 5)
DimPlot(combined, reduction = "pca", group.by = 'orig.ident') +
  theme(
    text = element_text(family = "Times")
  )
dev.off()
png("01.before_intergrated_PCA.png", width = 7, height = 5, units = 'in', res = 600)
DimPlot(combined, reduction = "pca", group.by = 'orig.ident') +
  theme(
    text = element_text(family = "Times")
  )
dev.off()
system.time(save(combined , file = "unintergrated.Rdata"))
load("unintergrated.Rdata")
combined <- IntegrateLayers(
  object = combined,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  k.weight = 60,
  verbose = FALSE)
scRNA <- combined
scRNA <- JoinLayers(scRNA)
system.time(save(scRNA , file = "scRNA_intergrated.Rdata"))
rm(list = ls())
if (! dir.exists("./03_PCA")){
  dir.create("./03_PCA")
}
load('../02_Integrate/scRNA_intergrated.Rdata')
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(scRNA), 10)
p <- VariableFeaturePlot(scRNA)
p <- LabelPoints(plot = p, 
                 points = top10, 
                 repel = T)
pdf(file="01.feature_selection.pdf",width=10,height=6,family='Times')              
theme.set= theme(
  
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 12,  family = "Times"),
  axis.text.y = element_text(size = 14,  family = "Times"),
  legend.text = element_text(size = 14, family = "Times"),
  legend.title = element_text(size = 16,face='bold',family = "Times"),
  text = element_text(family = "Times"))
plot1 <- VariableFeaturePlot(object = scRNA)+theme.set
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot2))
dev.off()
png(file="01.feature_selection.png",width=10,height=6,family='Times',units='in',res=600)              
theme.set= theme(axis.title.x=element_blank(),
                 axis.title = element_text(size = 20, face = "bold", family = "Times"),
                 axis.text.x = element_text(size = 12,  family = "Times"),
                 axis.text.y = element_text(size = 14,  family = "Times"),
                 legend.text = element_text(size = 14, family = "Times"),
                 legend.title = element_text(size = 16,face='bold',family = "Times"),
                 text = element_text(family = "Times"))
plot1 <- VariableFeaturePlot(object = scRNA)+theme.set
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot2))
dev.off()
scRNA.nor.sca <- ScaleData(scRNA)
scRNA.norm.pca <- RunPCA(scRNA.nor.sca, features = VariableFeatures(object = scRNA.nor.sca), npcs = 50)
print(scRNA.norm.pca[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(scRNA.norm.pca, dims = 1:2, reduction = "pca")
DimPlot(scRNA.norm.pca, reduction = "pca",group.by = 'orig.ident')
pdf("02.PCA.pdf",w=8,h=5)
DimPlot(scRNA.norm.pca, reduction = "pca",group.by = 'orig.ident')+theme.set
dev.off()
png("02.PCA.png",w=8,h=5,units = 'in',res = 600)
DimPlot(scRNA.norm.pca, reduction = "pca",group.by = 'orig.ident')+theme.set
dev.off()
pdf("03.PCA_heatmap.pdf",w=9,h=9)
DimHeatmap(object = scRNA, dims = 1:9, cells = 500,balanced = TRUE)+theme.set
dev.off()
png("03.PCA_heatmap.png",w=9,h=9,units = 'in',res = 600)
DimHeatmap(scRNA.norm.pca, dims = 1:9, cells = 500, balanced = TRUE)+theme.set
dev.off()
scRNA.norm.pca <- JackStraw(scRNA.norm.pca, num.replicate = 100, dims = 50) 
scRNA.norm.pca <- ScoreJackStraw(scRNA.norm.pca, dims = 1:50)
system.time(save(scRNA.norm.pca, file = "scRNA.norm.pca.Jack.Rdata"))
plot_pca <- JackStrawPlot(scRNA.norm.pca, dims = 1:50)
plot_pca
ggsave(filename = '04.pca_cluster.pdf',plot_pca,w=13,h=8)
ggsave(filename = '04.pca_cluster.png',plot_pca,w=13,h=8,dpi = 600)
plot_elbow <- ElbowPlot(scRNA.norm.pca, ndims = 50)
plot_elbow
ggsave(filename = '05.pca_sd.pdf',plot_elbow,w=6,h=5)
ggsave(filename = '05.pca_sd.png',plot_elbow,w=6,h=5,dpi = 600)
load('./scRNA.norm.pca.Jack.Rdata')
scRNA.norm.pca.clu <- FindNeighbors(scRNA.norm.pca, dims = 1:30)
scRNA.norm.pca.clu <- FindClusters(scRNA.norm.pca.clu, resolution = 0.4) 
UMAP <- RunUMAP(object = scRNA.norm.pca.clu, dims = 1:20)                      
system.time(save(UMAP, file = "UMAP.Rdata"))
theme.set <-  theme(
  axis.title = element_text(size = 16, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 10,  family = "Times"),
  axis.text.y = element_text(size = 14,  family = "Times"),
  legend.position = 'right',legend.direction = 'vertical',
  legend.text = element_text(size = 14, family = "Times"),
  legend.title = element_text(size = 16,face='bold',family = "Times"),
  text = element_text(family = "Times"))
pdf(file = paste0("04.UMAP.pdf"),width = 8,height = 6)
a <- dev.cur()   
png(file = paste0("04.UMAP.png"),width= 8, height= 6, units="in", res=300)
dev.control("enable")
par(mar = c(2,2,2,2),cex=1.5,family="Times")
DimPlot(UMAP,reduction = 'umap',group.by = 'orig.ident')+
  labs(x = "UMAP1", y = "UMAP2",title = "UMAP") +theme.set
dev.copy(which = a) 
dev.off()
dev.off()
pdf(file = paste0("04.UMAP1.pdf"),width = 8,height = 6)
a <- dev.cur()   
png(file = paste0("04.UMAP1.png"),width= 8, height= 6, units="in", res=300)
dev.control("enable")
par(mar = c(2,2,2,2),cex=1.5,family="Times")
DimPlot(UMAP,reduction = 'umap',label = T)+
  labs(x = "UMAP1", y = "UMAP2",title = "UMAP") +theme.set
dev.copy(which = a) 
dev.off()
dev.off()
pdf(file = paste0("04.UMAP_sample.pdf"),width = 6,height = 8)
a <- dev.cur()   
png(file = paste0("04.UMAP_sample.png"),width= 6, height= 8, units="in", res=300)
dev.control("enable")
par(mar = c(2,2,2,2),cex=1.5,family="Times")
DimPlot(UMAP,reduction = 'umap',split.by = 'orig.ident',ncol = 3)+
  labs(x = "UMAP1", y = "UMAP2",title = "UMAP") +theme.set
dev.copy(which = a) 
dev.off()
dev.off()
rm(list = ls())
options(stringsAsFactors = F)
if(!dir.exists("04_Marker")){dir.create("04_Marker")}
load('../03_PCA/UMAP.Rdata')
library(tidyverse)
all.markers<-FindAllMarkers(UMAP,only.pos = TRUE,min.pct = 0.25,
                            logfc.threshold = 0.5,test.use = 'wilcox',return.thresh = 0.01)
write.table(all.markers,file = '04.AllMarkers.xls',sep = '\t',row.names = F,quote = F)
all.markers <- read.table('04.AllMarkers.xls',header = T)
library(tidyverse)
all.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
top<-all.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
write.table(top,file = '04.topMarkers_10.xls',sep = '\t',row.names = F,quote = F)
top<-all.markers%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
write.table(top,file = '04.topMarkers_5.xls',sep = '\t',row.names = F,quote = F)
pbmc.markers <- read.table('04.AllMarkers.xls',header = T)
pbmc.markers <- pbmc.markers%>%group_by(cluster)%>%top_n(n=30,wt=avg_log2FC)
markerKU<-readxl::read_xlsx('./Cell_marker_Human.xlsx')
index<-which(markerKU$tissue_type == 'Bone marrow') 
markerKU<-markerKU[index,c('cell_name','marker','PMID','year')]
markerKU<-markerKU %>% na.omit()
the_marker<-markerKU[which(markerKU$marker %in% pbmc.markers$gene),] 
markerlist<-merge(pbmc.markers,the_marker,by.x='gene',by.y= 'marker')
markerlist<-markerlist[order(markerlist$cluster),] 
clus_list<-markerlist$cluster %>% as.factor() %>% levels()
theme.set = theme(
  
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 18,  face = "bold", family = "Times"),
  axis.text.y = element_text(size = 18,  face = "bold", family = "Times"),
  legend.text = element_text(size = 18, face = "bold", family = "Times"),
  legend.title = element_text(size = 20,face='bold',family = "Times"),
  text = element_text(family = "Times"))
for (i in  1:length(clus_list)) {
  the_clu<-clus_list[i]
  my_gene<-markerlist$gene[which(markerlist$cluster==the_clu)] %>% unique()
  pdf(file=paste(the_clu,"_DotPlot.pdf"),width=(length(my_gene)/2 +5) ,height=9,family='Times')
  
  p<-DotPlot(UMAP, features = my_gene, assay = 'RNA') + 
    RotatedAxis() + 
    theme.set+
    scale_color_gradientn(colours = c('#b2e7cb','#87CEFA','#66CC66','#FFCC33'))
  p %>% print()
  dev.off()
} 
theme.set <-  theme(
  axis.title = element_text(size = 16, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 10,  family = "Times"),
  axis.text.y = element_text(size = 14,  family = "Times"),
  legend.position = 'right',legend.direction = 'vertical',
  legend.text = element_text(size = 14, family = "Times"),
  legend.title = element_text(size = 16,face='bold',family = "Times"),
  text = element_text(family = "Times"))
pdf('05.topheat.pdf',w=20,h=25)
DoHeatmap(UMAP,features = top$gene,label = F)+
  labs(title="", x="Cells separated by clusters", y = "",size=40)+theme.set
dev.off()
png('05.topheat.png',w=20,h=25,units = 'in',res = 600)
DoHeatmap(UMAP,features = top$gene,label = F)+
  labs(title="", x="Cells separated by clusters", y = "",size=40)+theme.set
dev.off()
rm(list = ls())
options(stringsAsFactors = F)
if(!dir.exists("05_SingleR")){dir.create("05_SingleR")}
library(Seurat)
library(dplyr)
library(Matrix)
library(scales)
library(harmony)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(SingleR)
load("../03_PCA/UMAP.Rdata")
scRNA_singleR <- GetAssayData(UMAP, slot = "data")
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
clusters=UMAP@meta.data$seurat_clusters
pred.hesc <- SingleR(test = scRNA_singleR,
                     ref = hpca.se,
                     labels = hpca.se$label.main,
                     method = "cluster", 
                     clusters = clusters, 
                     assay.type.test = "logcounts", 
                     assay.type.ref = "logcounts")
table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), 
                      celltype=pred.hesc$labels, 
                      stringsAsFactors = F) 
UMAP@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
theme.set= theme(
  
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 12,  family = "Times"),
  axis.text.y = element_text(size = 14,  family = "Times"),
  legend.text = element_text(size = 14, family = "Times"),
  legend.title = element_text(size = 16,face='bold',family = "Times"),
  text = element_text(family = "Times",size=14,face='bold'))
pdf('1000.cell_type_singleR.pdf',w=10,h=8,family='Times')
p1 <- DimPlot(UMAP, reduction = "umap", group.by = "singleR", label = TRUE, 
              label.size = 5, cols = brewer.pal(9, "Set3")) + theme.set
p1
dev.off()
png('1000.cell_type_singleR.png',w=10,h=8,family='Times',units='in',res=600)
p1 <- DimPlot(UMAP, reduction = "umap", group.by = "singleR", label = TRUE, 
              label.size = 5, cols = brewer.pal(9, "Set3")) + theme.set
p1
dev.off()
system.time(save(UMAP, file = "1000.UMAP_singleR.Rdata")) 
load('../03_PCA/UMAP.Rdata')
theme.set <-  theme(
  axis.title = element_text(size = 16, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 10, family = "Times"),
  axis.text.y = element_text(size = 14, family = "Times"),
  legend.position = 'right',legend.direction = 'vertical',
  legend.text = element_text(size = 14, family = "Times"),
  legend.title = element_text(size = 16,face='bold', family = "Times"),
  text = element_text(family = "Times"))
gene<-c("TIMP1","SELENOP","SELL","IGLL1","FUT4","RETN","PRTN3","CLEC5A", 
        "CD3E","CD2","IL32","IL7R","LILRB2",'MAFB','CPA3','CD3D','FUS','TUBB','CD8A','BCL11B',
        'AHSP','HBA1','HBB','HBD','G0S2','CD163','CD93','FGR','IKZF2','IL2RA','MEIS1','MMRN1',
        'SDC1','TNFRSF17','MS4A1','PAX5','LILRA4','GZMB')
pdf(file="01.DotPlot.pdf",width=9,height=8,family='Times')
DotPlot(UMAP,features = gene)+theme.set+ coord_flip()+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  theme(legend.text=element_text(size=15))+
  theme(legend.title=element_text(size=15))+
  theme(axis.title.x =element_text(size=15), axis.title.y=element_text(size=15))+
  
  theme(axis.text = element_text(size = 15))+
  theme(axis.text.x = element_text(angle = 0, color = "black", size = 15, face = 2))
dev.off()
png(file="01.DotPlot.png",width=9,height=8,family='Times',units='in',res=600)
DotPlot(UMAP,features = gene)+theme.set+ coord_flip()+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  theme(legend.text=element_text(size=15))+
  theme(legend.title=element_text(size=15))+
  theme(axis.title.x =element_text(size=15), axis.title.y=element_text(size=15))+
  
  theme(axis.text = element_text(size = 15))+
  theme(axis.text.x = element_text(angle = 0, color = "black", size = 15, face = 2))
dev.off()
UMAP_sub <- UMAP
UMAP_sub<-subset(UMAP_sub,seurat_clusters %in% c(0:19))
new.cluster.ids = c(
  "0" = "Macrophage", 
  "1" = "Progenitor cell",
  "2" = "Myeloid cell",  
  "3" = "GMP",  
  "4" = "Granulocyte", 
  "5" = "Myeloid cell",   
  "6" = "T cell",  
  "7" = "Monocyte", 
  "8" = "Granulocyte-monocyte progenitor", 
  '9' = 'T cell', 
  '10' = 'Natural killer T (NKT) cell',
  '11' = "T cell", 
  '12' = "Erythroid lineage cell", 
  '13' = 'Monocyte', 
  '14' = 'Monocyte', 
  '15' = 'Regulatory T(Treg) cell',
  '16' = 'Hematopoietic stem cell\nHematopoietic progenitor cell' ,
  '17' = 'Plasma cell' ,
  '18' = 'B cell',
  '19' = 'Plasmacytoid dendritic cell' 
)  
UMAP2name <- RenameIdents(UMAP_sub, new.cluster.ids)
UMAP2name @meta.data$seurat_clusters<-UMAP2name @active.ident
UMAP2name @meta.data$seurat_clusters %>% as.factor() %>% levels()
UMAP2name $celltype = Idents(UMAP2name )
save(UMAP2name,file = "./UMAP2name.Rdata")
pdf(file="02.Cell_marker_DotPlot.pdf",width=9,height=8,family='Times')
DotPlot(UMAP2name,features = gene)+theme.set+ coord_flip()+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  theme(legend.text=element_text(size=15))+
  theme(legend.title=element_text(size=15))+
  theme(axis.title.x =element_text(size=15), axis.title.y=element_text(size=15))+
  
  theme(axis.text = element_text(size = 15))+
  theme(axis.text.x = element_text(angle = 0, color = "black", size = 15, face = 2))
dev.off()
pdf(file="02.cellmarker_DotPlot.pdf",width=(length(gene)/2+3) ,height=5,family='Times')
p<-DotPlot(UMAP2name, features = gene, assay = 'RNA') + 
  RotatedAxis() + 
  theme.set+
  scale_color_gradientn(colours = c('#b2e7cb','#87CEFA','#66CC66','#FFCC33'))
p %>% print()
dev.off()
pdf(file="02.cellmarker_DotPlot.pdf",width=(length(gene)/2+3) ,height=5,family='Times')
p<-DotPlot(UMAP2name, features = gene, assay = 'RNA') + 
  RotatedAxis() + 
  theme.set+
  scale_color_gradientn(colours = c('#b2e7cb','#87CEFA','#66CC66','#FFCC33'))
p %>% print()
dev.off()
table(UMAP2name$celltype)
system.time(save(UMAP_sub, file = "UMAP_sub_before_UMAP2name.Rdata"))
load('./UMAP_sub_before_UMAP2name.Rdata')
theme.set = theme(
  
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 18,  face = "bold", family = "Times"),
  axis.text.y = element_text(size = 18,  face = "bold", family = "Times"),
  legend.text = element_text(size = 18, face = "bold", family = "Times"),
  legend.title = element_text(size = 20,face='bold',family = "Times"),
  text = element_text(family = "Times"))
pdf(file="03.UMAP_2th_select.pdf",width=20,height=15,family='Times')
DimPlot(UMAP2name  ,reduction = "umap",label.size = 6,label = TRUE,pt.size=0.8)+theme.set
dev.off()
png(file="03.UMAP_2th_select.png",width=20,height=15,units="in", res=300,family='Times')
DimPlot(UMAP2name  ,reduction = "umap",label.size = 6,label = TRUE,pt.size=0.8)+theme.set
dev.off()
pdf(file="03.UMAP_2th_select_Group.pdf",width=30,height=15,family='Times')
DimPlot(UMAP2name ,reduction = "umap",label.size = 6,label = T,pt.size=0.8, split.by= "Type",dims=c(2,1))+theme.set
dev.off()
png(file="03.UMAP_2th_select_Group.png",width=30,height=15,units="in", res=300,family='Times')
DimPlot(UMAP2name ,reduction = "umap",label.size = 6,label = T,pt.size=0.8, split.by= "Type",dims=c(2,1))+theme.set
dev.off()
pdf(file = paste0("04.UMAP.pdf"),width = 8,height = 6)
a <- dev.cur()   
png(file = paste0("04.UMAP.png"),width= 8, height= 6, units="in", res=300)
dev.control("enable")
par(mar = c(2,2,2,2),cex=1.5,family="Times")
DimPlot(UMAP,reduction = 'umap',label =T)+
  labs(x = "UMAP1", y = "UMAP2",title = "UMAP") +theme.set
dev.copy(which = a) 
dev.off()
dev.off()
anno_plot2 <- DimPlot(UMAP2name, reduction = "umap", split.by = "orig.ident",label = T)+  
  labs(x = "UMAP1", y = "UMAP2",title = "UMAP") +
  theme_bw()+
  theme(legend.position = "top")+
  theme(axis.title.x =element_text(size=20,family = "Times", face = "bold"),
        axis.text.x =element_text(angle=0,size=20,hjust = 0.5,family = "Times", face = "bold"),
        axis.title.y =element_text(size=20,family = "Times", face = "bold"),
        axis.text.y=element_text(size=16,family = "Times", face = "bold"))
ggsave(filename = '02.celltype.group1.pdf',anno_plot2,w=40,h=6)
ggsave(filename = '02.celltype.group1.png',anno_plot2,w=40,h=6)
load("../../05_SingleR/UMAP2name.Rdata")
UMAP<-UMAP2name
UMAP@meta.data$seurat_clusters
scRNA_case <- subset(UMAP, Type == 'AML')
phe <- scRNA_case@meta.data
table(phe$celltype)
cell_type_case <- as.data.frame(sort(table(phe$celltype)))
colnames(cell_type_case) <- c('Celltype', 'AML')
scRNA_Normal <- subset(UMAP,Type == 'Healthy')
phe2 <- scRNA_Normal@meta.data
table(phe$celltype)
cell_type_Normal <- as.data.frame(sort(table(phe2$celltype)))
colnames(cell_type_Normal) <- c('Celltype', 'Healthy')
cell_type <- merge(cell_type_Normal, cell_type_case, by='Celltype')
write.csv(cell_type, 'cell_type.csv')
library(reshape2)
library(plyr)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(dplyr)
mydata <- melt(cell_type, id.vars = "Celltype", variable.name = "Group", value.name = "variable")
custom_colors <- c("#FF9999", "#66B3FF", "#99FF99", "#FFCC99", "#FF6699",  
                   "#CCCCFF", "#FFB3FF", "#FFFF99", "#99FFFF", "#FF99CC",  
                   "#CC99FF", "#C0C0C0", "#FF6666", "#66FFFF", "#99FFCC",  
                   "#FFCC66", "#66CCFF", "#CCCCFF", "#FFFFCC")  
p <- ggplot(mydata, aes(x=Group, y=variable, fill=Celltype)) +  
  geom_bar(position = "fill", stat="identity", alpha=0.7) +  
  theme_bw() +  
  scale_fill_manual(values=custom_colors) +  
  theme(axis.title.x = element_text(size=22, color='black', face = "bold", family='Times'),  
        axis.text.x = element_text(size=18, face = "bold", family='Times'),  
        axis.title.y = element_text(size=22, color='black', face = "bold", family='Times'),  
        axis.text.y = element_text(size=18, face = "bold", family='Times'),  
        legend.title = element_text(size=20, color='black', face = "bold", family='Times'),  
        legend.text = element_text(size=18, color='black', face = "bold", family='Times'),  
        title = element_text(size=20, color='black', face = "bold", family='Times'),  
        strip.text = element_text(size = 14, family = "Times", face = "bold")) +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  labs(x="Celltype", y="Proportion", fill="")  
ggsave('04.cell_proportion.pdf', width = 10, height = 12, plot = p)  
ggsave('04.cell_proportion.png', w=10, h=12, plot = p)
UMAP@meta.data$seurat_clusters %>% as.factor() %>% table()
mastUMAP<-UMAP
Idents(mastUMAP)<-mastUMAP@meta.data$Type
levels(mastUMAP)
cell_ls<-mastUMAP@meta.data$celltype %>% levels()
if (! dir.exists("./mast")){
  dir.create("./mast")
}
library(Seurat)
for (i in cell_ls) {
  the_name <- NULL
  mastUMAP <- subset(UMAP, celltype == i)
  Idents(mastUMAP) <- mastUMAP@meta.data$Type
  levels(mastUMAP)
  
  
  MAST_Result <- tryCatch(
    {
      Seurat::FindMarkers(mastUMAP,
                          ident.1 = "AML",
                          ident.2 = "Healthy",
                          
                          min.pct = 0.1,
                          logfc.threshold = log(1),
                          min.diff.pct = 0.1)
    },
    error = function(e) {
      message(paste("Error in cell type:", i, "-", e$message, "Skipping to next..."))
      return(NULL)  
    }
  )
  
  
  if (is.null(MAST_Result)) {
    next
  }
  
  the_name <- gsub('\n', ' ', i)
  write.csv(MAST_Result, file = paste(the_name, 'wilcox.csv'), row.names = TRUE)
}
list.files()
load("../../05_SingleR/UMAP2name.Rdata")
UMAP_Tcell <- subset(UMAP2name, seurat_clusters  %in%c('T cell')) 
save(UMAP_Tcell,file = "UMAP_Tcell.Rdata")
UMAP_Tcell@assays
sc_DEGs <- FindAllMarkers(object = UMAP2name, only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)
label_data_df <-sce.markers%>%group_by(cluster)%>%top_n(n=5,wt=abs(avg_log2FC))
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
sc_DEGs_order <- sc_DEGs[order(sc_DEGs$cluster,sc_DEGs$avg_log2FC),]
cut_point <- 0.25
sc_DEGs$regulated1 <- ifelse(sc_DEGs$avg_log2FC < -cut_point & sc_DEGs$p_val<0.05 ,'down',(
  ifelse(sc_DEGs$avg_log2FC > cut_point & sc_DEGs$p_val<0.05,'up','no sig')
))
table(sc_DEGs$cluster,sc_DEGs$regulated1)
ggplot(sc_DEGs) +
  
  geom_jitter(aes(cluster, avg_log2FC, color = regulated1),
              size = 1, width = 0.4, alpha = 0.8,family = "Times") +
  
  geom_tile(aes(cluster, 0, fill = cluster),
            height = 0.4,
            color = "black",
            alpha = 0.5,
            show.legend = F,
            width = 0.8,family = "Times") +
  
  
  
  
  
  
  geom_text_repel(
    data = label_data_df,
    aes(cluster, avg_log2FC, label = gene),
    size = 3, max.overlaps = 100,family = "Times"
  ) +
  xlab("") +
  ylab("avg_log2FC") +
  
  scale_fill_hue() +
  
  scale_color_hue(name = "Regulate") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5,family = "Times"),
        legend.position = "top")
ggsave("manhadun.pdf", height = 8, width = 12)
ggsave("manhadun.png", height = 8, width = 12)
rm(list = ls())
options(stringsAsFactors = F)
if(!dir.exists("06_Single_DEG")){dir.create("06_Single_DEG")}
load("../05_SingleR/UMAP2name.Rdata")
UMAP2name@meta.data$celltype <- Idents(UMAP2name)
table(UMAP2name@meta.data$celltype,UMAP2name@meta.data$Patients_ID)
df_bar<-as.data.frame(prop.table(table(UMAP2name@meta.data$celltype,UMAP2name@meta.data$Patients_ID)))
colnames(df_bar)<-c("cell","patient","proportion")
df_bar$group <- ifelse(grepl("^AML",df_bar$patient),"AML","Healthy")
df_bar$proportion<-df_bar$proportion*100
pdf("01.cell_type_by_sample_bar_plot.pdf",height = 6,width = 12)
ggplot(df_bar,aes(x=patient,y=proportion,fill=cell))+
  geom_bar(stat = "identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length = unit(0.5,"cm"),axis.text = element_text(size = 10,color = "black"),
        axis.text.x =element_text(angle = 50, vjust = 1, hjust = 1,size=18, face = "bold",family='Times',color='black'),
        axis.text.y=element_text(size=18,  face = "bold",family='Times',color='black'),
        axis.title.y =element_text(size=22, face = "bold",family='Times',color='black'))+
  guides(fill=guide_legend(title=NULL))
dev.off()  
png("01.cell_type_by_sample_bar_plot.png",height = 6,width = 12,units='in',res = 600)
ggplot(df_bar,aes(x=patient,y=proportion,fill=cell))+
  geom_bar(stat = "identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length = unit(0.5,"cm"),axis.text = element_text(size = 10,color = "black"),
        axis.text.x =element_text(angle = 50, vjust = 1, hjust = 1,size=18, face = "bold",family='Times',color='black'),
        axis.text.y=element_text(size=18,  face = "bold",family='Times',color='black'),
        axis.title.y =element_text(size=22, face = "bold",family='Times',color='black'))+
  guides(fill=guide_legend(title=NULL))
dev.off()  
write.csv(df_bar,file = "./df_bar.csv")
pdf("./02.cell_type_by_sample_bar_plot_group.pdf",height = 10,width = 7)
ggplot(df_bar,aes(x=cell,y=proportion,fill=group))+
  geom_bar(stat = "identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length = unit(0.5,"cm"),axis.text = element_text(size = 10,color = "black"),
        axis.text.x =element_text(angle = 50, vjust = 1, hjust = 1,size=18, face = "bold",family='Times',color='black'),
        axis.text.y=element_text(size=18,  face = "bold",family='Times',color='black'),
        axis.title.y =element_text(size=22, face = "bold",family='Times',color='black'))+
  guides(fill=guide_legend(title=NULL))
dev.off()  
library(ggpubr)
pdf("./03.pdf",height = 10,width = 12)
ggplot(df_bar,aes(x=cell,y=proportion,fill=group))+
  geom_boxplot(width=0.5,position = position_dodge(0.5))+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length = unit(0.5,"cm"),axis.text = element_text(size = 10,color = "black"))+
  guides(fill=guide_legend(title=NULL))+
  stat_compare_means(aes(group=group),
                     method = "wilcox.test",
                     label = "p.format",
                     symnum.args = list(cutpoints = c(0,0.001,0.01,0.05,1),
                                        symbol = c("***","**","*",""))
  )
dev.off()
pdf("./04.cell_type_by_sample_bar_plot_group_deg.pdf",height = 6,width = 7)
ggplot(df_bar,aes(x=cell,y=proportion,fill=group)) +
  scale_fill_manual(values = c("red","#20B2AA")) + 
  ylim(0,0.2)+
  geom_boxplot()+
  stat_compare_means(aes(group = group),method = 'wilcox.test',label = "p.signif",cex=4)+
  xlab("cell") + 
  ylab("proportion") + 
  theme(panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA))+
  theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))
dev.off()
library(ggplot2)
library(ggpubr)
pdf("./03.cell_type_by_sample_box.pdf", height = 6, width = 9, family = "serif") 
ggplot(df_bar, aes(x = cell, y = proportion, fill = group)) +  
  geom_boxplot(width = 0.5, position = position_dodge(0.5)) +  
  ggtitle("") +  
  theme_bw() +  
  theme(axis.ticks.length = unit(0.5, "cm"),  
        axis.text = element_text(size = 10, color = "black", family = "serif"), 
        legend.text = element_text(family = "serif"), 
        legend.title = element_text(family = "serif"), 
        axis.title.x = element_text(size = 15, face = "bold", family = "serif"), 
        axis.title.y = element_text(size = 15, face = "bold", family = "serif")) + 
  guides(fill = guide_legend(title = NULL)) +  
  stat_compare_means(aes(group = group),  
                     method = "wilcox.test",  
                     label = "p.format",  
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),  
                                        symbol = c("***", "**", "*", ""))) +  
  scale_y_continuous(limits = c(0, 0.2)) 
dev.off()
pdf("./04.cell_type_by_sample_bar_plot_group_deg.pdf", height = 6, width = 7, family = "serif")  
ggplot(df_bar, aes(x = cell, y = proportion, fill = group)) +  
  scale_fill_manual(values = c("red", "#20B2AA")) +  
  ylim(0, 0.2) +  
  geom_boxplot() +  
  stat_compare_means(aes(group = group), method = 'wilcox.test', label = "p.signif", cex = 4) +  
  xlab("cell") +   
  ylab("proportion") +   
  theme(panel.background = element_blank(),  
        panel.border = element_rect(linetype = "solid", fill = NA),  
        axis.text.x = element_text(face = "bold", color = "gray50", angle = 70, vjust = 1, hjust = 1, family = "serif"), 
        axis.text.y = element_text(family = "serif"), 
        axis.title.x = element_text(family = "serif"), 
        axis.title.y = element_text(family = "serif")) 
dev.off()
load("../02_Integrate/scRNA_intergrated.Rdata")
s_cmb <- UMAP2name
s_cmb@meta.data$seurat_clusters<-factor(s_cmb@meta.data$seurat_clusters)
s_cmb@assays$RNA$data %>% head(10)
LayerData(s_cmb, assay="RNA", layer='data')
GetAssayData(s_cmb, assay="RNA", slot='data')
s_cmb[["RNA"]]$scale.data[c(1:4),c(1:4)]
gene<-read.csv("/data/nas1/anbingzhuang/project/17_BJTC-596-3/04.Cox_LASSO/03.cox.csv", header = T)
gene<-gene$x
dat <- s_cmb[["RNA"]]$scale.data[c(1:2000),c(1:15714)]
dat <- dat %>%  t(.)%>% as.data.frame %>% tibble::rownames_to_column(var = "sample")
gene %in% colnames(dat)
group<-data.frame(sample=colnames(s_cmb),Group=s_cmb$Type)
dat <- merge(dat, group, by = "sample")
seurat_clusters<-data.frame(sample=colnames(s_cmb),group=s_cmb$seurat_clusters)
identical(seurat_clusters$sample,dat$sample )
dat2 <- merge(dat,seurat_clusters,by='sample')
for (i in c(1,2,3)){
  dat3 <- dat2[,c("sample",gene[i],"group","Group")] %>% as.data.frame()
  colnames(dat3)<-c('sample','Expression','seurat_clusters','group')
  dat3<-dat3[order(dat3$seurat_clusters,decreasing = F),]
  library(ggpubr)
  p <- ggboxplot(dat3 , x = "seurat_clusters", y = "Expression",
                 fill= "group",
                 palette = c("#4DAF4A","#984EA3"),
                 
                 bxp.errorbar = T,
                 
                 
                 
  ) +
    stat_compare_means(method = "wilcox.test",label = "p.signif",cex=5,color='black',family="Times",face = "bold",aes(group = group))+
    
    theme_bw()+
    labs(x = "", y =paste("Expression of",gene[i],sep=' '), color = "") +
    theme(axis.title.x = element_text(size = 20, face = "bold", family = "Times"),
          axis.title.y = element_text(size = 20, face = "bold", family = "Times"),
          axis.text.x = element_text(size = 13, color='black',face = "bold", family = "Times", angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 13, face = "bold", family = "Times"),
          legend.text = element_text(size = 15, family = "Times",face = "bold"),
          legend.title = element_text(size = 18, family = "Times",face = "bold"),
          
          plot.title = element_text(size = 20, color='black',family = "Times",face = "bold",hjust=0.5),
          text = element_text(family = "Times"),
          panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    theme(legend.position = "top")
  p_name <- paste('p',c(1,2),sep='')   
  assign(p_name[i], p) 
}
library(ggpubr)
library(cowplot)
p1
ggsave('05.CALR_boxplot.pdf',w=18,h=12)
ggsave('05.CALR_boxplot.png',w=18,h=12)
p2
ggsave('05.CDK6_boxplot.pdf',w=18,h=12)
ggsave('05.CDK6_boxplot.png',w=18,h=12)
pdf("05.gene.pdf",width = 12,height = 6)
DotPlot(s_cmb,features = gene)+RotatedAxis()
dev.off()
if(!dir.exists("07_cellchat")){dir.create("07_cellchat")}
load("../05_SingleR/UMAP2name.Rdata")
library(CellChat) 
library(Seurat)
library(ggplot2)
UMAP<-UMAP2name
data.input <- UMAP@assays$RNA$counts 
meta = UMAP@meta.data 
table(meta$Type)
cell.use = rownames(meta)[meta$Type == "Healthy"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
table(meta$seurat_clusters)
unique(meta$seurat_clusters)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype",)
cellchat <- setIdent(cellchat, ident.use = "celltype") 
levels(cellchat@idents)
unique(cellchat@idents)
cellchat@idents <- droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents), unique(cellchat@idents)))
groupSize <- as.numeric(table(cellchat@idents)) 
groupSize
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB) 
library(tidyverse)
dplyr::glimpse(CellChatDB$interaction)  
CellChatDB_interaction <- CellChatDB$interaction
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 
cellchat <- subsetData(cellchat, features = NULL) 
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat) 
cellchat <- projectData(cellchat, PPI.human) 
cellchat <- computeCommunProb(cellchat) 
df.net <- subsetCommunication(cellchat) 
if (! dir.exists("./control")){
  dir.create("./control")
}
write.csv(df.net, 'df.net.csv')
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
pdf('04.number_of_interactions.pdf',w=10,h=10,family='serif')
netVisual_heatmap(cellchat)
dev.off()
png('04.number_of_interactions.png',w=1000,h=1000,family='serif')
netVisual_heatmap(cellchat)
dev.off()
pdf('01.net_circle_number.pdf',w=7,h=9,family='serif')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
png('01.net_circle_number.png',w=7,h=9,units='in',res=600,family='serif')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf('02.net_circle_weight.pdf',w=7,h=9,family='serif')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()
png('02.net_circle_weight.png',w=7,h=9,units='in',res=600,family='serif')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()
levels(cellchat@idents)
cellchat@data.signaling
pdf('03.single_circle.pdf',w=30,h=10)
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()
png('03.single_circle.png',w=30,h=10,units = 'in',res = 600)
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()
if(!dir.exists("07_cellchat")){dir.create("07_cellchat")}
library(CellChat) 
library(Seurat)
library(ggplot2)
data.input <- UMAP@assays$RNA$counts 
meta = UMAP@meta.data 
table(meta$Type)
cell.use = rownames(meta)[meta$Type == "AML"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
table(meta$celltype)
unique(meta$celltype)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype",)
cellchat <- setIdent(cellchat, ident.use = "celltype") 
levels(cellchat@idents)
unique(cellchat@idents)
cellchat@idents <- droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents), unique(cellchat@idents)))
groupSize <- as.numeric(table(cellchat@idents)) 
groupSize
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB) 
library(tidyverse)
dplyr::glimpse(CellChatDB$interaction)  
CellChatDB_interaction <- CellChatDB$interaction
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 
cellchat <- subsetData(cellchat, features = NULL) 
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat) 
cellchat <- projectData(cellchat, PPI.human) 
cellchat <- computeCommunProb(cellchat) 
df.net <- subsetCommunication(cellchat) 
if (! dir.exists("./disease")){
  dir.create("./disease")
}
write.csv(df.net, 'df.net.csv')
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
pdf('04.number_of_interactions.pdf',w=10,h=10,family='serif')
netVisual_heatmap(cellchat)
dev.off()
png('04.number_of_interactions.png',w=600,h=600,family='serif')
netVisual_heatmap(cellchat)
dev.off()
pdf('01.net_circle_number.pdf',w=10,h=12,family='serif')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
png('01.net_circle_number.png',w=10,h=12,units='in',res=600,family='serif')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf('02.net_circle_weight.pdf',w=7,h=9,family='serif')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()
png('02.net_circle_weight.png',w=7,h=9,units='in',res=600,family='serif')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()
levels(cellchat@idents)
cellchat@data.signaling
pdf('03.single_circle.pdf',w=30,h=10)
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()
png('03.single_circle.png',w=30,h=10,units = 'in',res = 600)
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()
rm(list = ls())
library(Seurat)
devtools::load_all('/data/nas1/luchunlin/project/YQQL0106-11/monocle')
library(monocle)
if (! dir.exists("./08_Trajectory")){
  dir.create("./08_Trajectory")
}
load('../05_SingleR/UMAP2name.Rdata') 
if (! dir.exists("./T_cell")){
  dir.create("./T_cell")
}
UMAP <-  subset(UMAP2name, seurat_clusters  %in%c('T cell') )
markers.gene1 <- read.csv('../../05_SingleR/mast/T cell wilcox.csv')[,1]
markers.gene<-c(markers.gene1) %>% unique()
data <- as(as.matrix(UMAP@assays$RNA$counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = UMAP@meta.data)
fData <- data.frame(gene_short_name = row.names(data), 
                    row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
mycds <- detectGenes(mycds, min_expr = 0.1)
mycds <- setOrderingFilter(mycds,markers.gene)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)
save(mycds, file = "mycds_order.t2d.RData")
load('./mycds_order.t2d.RData')
p1 <- plot_cell_trajectory(mycds, 
                           color_by = "Pseudotime",
                           cell_size=0.5,
                           show_backbone=TRUE,
                           show_tree=TRUE, 
                           backbone_color="black", 
                           markers=NULL, 
                           use_color_gradient = FALSE,
                           markers_linear = FALSE,
                           show_cell_names=FALSE,
                           show_state_number = FALSE,
                           cell_link_size=0.75,
                           cell_name_size=2,
                           state_number_size = 2.9,
                           show_branch_points=TRUE,
                           theta = 0) +
  theme(legend.position = "right") +
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p1
ggsave(paste0('01.', '3 clus Endothelial', '_Trajectory.png'), p1, width = 6, height = 6, units = "in", dpi = 300)
ggsave(paste0('01.', '3 clus Endothelial', '_Trajectory.pdf'), p1, width = 6, height = 6, units = "in", dpi = 300)
p2 <- plot_cell_trajectory(mycds, color_by = "State",cell_size=0.5,show_backbone=TRUE) +
  theme(legend.position = "right") +
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2
ggsave(paste0('02.', '3 clus Endothelial', '_State_Trajectory.png'), p2, width = 6, height = 6, units = "in", dpi = 300)
ggsave(paste0('02.', '3 clus Endothelial', '_State_Trajectory.pdf'), p2, width = 6, height = 6, units = "in", dpi = 300)
p3 <- plot_cell_trajectory(mycds, color_by = "State",cell_size=0.5,show_backbone=TRUE) +
  theme(legend.position = "right") +
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  facet_wrap(~Type,nrow = 1)
p3
ggsave(filename = paste0('03.', '3 clus Endothelial', '_State_Trajectory(group).png'),p3,w=10,h=6, units = "in", dpi = 300)
ggsave(filename = paste0('03.', '3 clus Endothelial', '_State_Trajectory(group).pdf'),p3,w=10,h=6, units = "in", dpi = 300)
p4 <- plot_cell_trajectory(mycds, color_by = "Type",cell_size=0.5,show_backbone=TRUE) +
  theme(legend.position = "right") +
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p4
ggsave(filename = paste0('04.', '3 clus Endothelial', '_State_Trajectory(group_merge).png'),p4,w=10,h=6, units = "in", dpi = 300)
ggsave(filename = paste0('04.', '3 clus Endothelial', '_State_Trajectory(group_merge).pdf'),p4,w=10,h=6, units = "in", dpi = 300)
mycds@phenoData@data$seurat_clusters<-gsub('\\\n',' & ', mycds@phenoData@data$seurat_clusters)
p8 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters",cell_size=0.5,show_backbone=TRUE) +
  theme(legend.position = "right") +
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p8
ggsave(filename = paste0('04.', '3 clus Endothelial', '_State_Trajectory(group_merge).png'),p8,w=10,h=6, units = "in", dpi = 300)
ggsave(filename = paste0('04.', '3 clus Endothelial', '_State_Trajectory(group_merge).pdf'),p8,w=10,h=6, units = "in", dpi = 300)
hubgene <- read.csv('/data/nas1/anbingzhuang/project/17_BJTC-596-3/04.Cox_LASSO/03.cox.csv')
hubgene <- hubgene$x
cds_subset <- mycds[hubgene,]
p5 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime") +
  theme(legend.position = "right") +
  theme(axis.title.x = element_text(size = 22, color = 'black', face = "bold", family = 'Times'),
        axis.text.x = element_text(size = 18, family = 'Times'),
        axis.title.y = element_text(size = 22, color = 'black', face = "bold", family = 'Times'),
        axis.text.y = element_text(size = 18, family = 'Times'),
        legend.title = element_text(size = 20, color = 'black', face = "bold", family = 'Times'),
        legend.text = element_text(size = 18, face = "bold", family = 'Times'),
        plot.title = element_text(size = 30, color = 'black', face = "bold", family = 'Times'),
        strip.text = element_text(size = 24, face = "bold", family = 'Times')) +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p5
p6<-plot_genes_in_pseudotime(cds_subset ,color_by = "State")+
  theme(legend.position = "right")+
  theme(axis.title.x =element_text(size=40,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 24, face = "bold", family = 'Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p7<-plot_genes_in_pseudotime(cds_subset ,color_by = "Type")+
  theme(legend.position = "right")+
  theme(axis.title.x =element_text(size=40,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 24, face = "bold", family = 'Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p9<-plot_genes_in_pseudotime(cds_subset ,color_by = "seurat_clusters")+
  theme(legend.position = "right")+
  theme(axis.title.x =element_text(size=40,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 24, face = "bold", family = 'Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
plotc <- p5|p6|p7|p9
plotc
ggsave(filename = paste0('05.', '3 clus Endothelial', 'Trajectory_gene.png'),plotc,w=32,h=20, units = "in", dpi = 300)
ggsave(filename = paste0('05.', '3 clus Endothelial', 'Trajectory_gene.pdf'),plotc,w=32,h=20, units = "in", dpi = 300)
library(patchwork)
patch <- p1 + p2 + p3 + p4 + p8
patch
ggsave(filename = paste0('06.', '3 clus Endothelial', 'merge_Trajectory.png'),patch,w=20,h=12, units = "in", dpi = 300)
ggsave(filename = paste0('06.', '3 clus Endothelial', 'merge_Trajectory.pdf'),patch,w=20,h=12, units = "in", dpi = 300)
if (!dir.exists("subcell")){
  dir.create('subcell')
}
sce.car <-  subset(UMAP2name, seurat_clusters  %in%c('T cell') )
sce.car <- ScaleData(sce.car)
sce.car <- RunPCA(sce.car, features = VariableFeatures(object = sce.car))
ElbowPlot(sce.car, ndims = 50)
sce.car <- FindNeighbors(sce.car, dims = 1:30)
sce.car <- FindClusters(sce.car, resolution = 0.5)
table(sce.car@meta.data$seurat_clusters)
sce.car <- RunUMAP(sce.car, dims = 1:30)
p <- DimPlot(sce.car, reduction = "umap", group.by = "seurat_clusters", label = T) + ggtitle(paste0('Subtype of ', 'T cell' ))
pdf(file = paste0('01.T_cell', '.subtype.pdf'),w=6,h=5)
print(p)
dev.off()
png(file = paste0('01.T_cell', '.subtype.png'),w=500,h=400)
print(p)
dev.off()
all.markers <- FindAllMarkers(sce.car, only.pos = TRUE,min.pct = 0.1,logfc.threshold = 0.5,test.use = 'wilcox',return.thresh = 0.05)
marker_gene <- all.markers %>% group_by(cluster) %>% top_n(n=3,wt=avg_log2FC)
write.csv(marker_gene, file = paste0("T_cell", '_marker(subtype_car).csv'))
marker_gene <- marker_gene[!duplicated(marker_gene$gene), ]
p2 <- DotPlot(sce.car, features = marker_gene$gene,group.by = 'seurat_clusters') + RotatedAxis()+
  labs(title="", y="Subtype", x = "",size=40)
pdf(file = paste0('02.T_cell', '.marker.exp.pdf'),w=8,h=6,family = "Times")
print(p2)
dev.off()
png(file = paste0('02.T_cell', '.marker.exp.png'),w=600,h=450,family = "Times")
print(p2)
dev.off()
Mono_tj <- sce.car
Mono_matrix = GetAssayData(Mono_tj, slot = "count", assay = "RNA")
feature_ann = data.frame(gene_id=rownames(Mono_matrix),gene_short_name=rownames(Mono_matrix))
rownames(feature_ann) = rownames(Mono_matrix)
Mono_fd = new("AnnotatedDataFrame", data = feature_ann)
sample_ann = Mono_tj@meta.data
cell.type = Idents(Mono_tj) %>% data.frame
sample_ann = cbind(sample_ann, cell.type)
colnames(sample_ann)[ncol(sample_ann)] = "Subtype"
Mono_pd = new("AnnotatedDataFrame", data =sample_ann)
Mono.cds = newCellDataSet(Mono_matrix, phenoData =Mono_pd, featureData=Mono_fd, expressionFamily=negbinomial.size())
rm(Mono_tj)
Mono.cds = estimateSizeFactors(Mono.cds)
Mono.cds = estimateDispersions(Mono.cds)
WGCNA::collectGarbage()
Mono.cds = reduceDimension(Mono.cds, max_components = 2, verbose = T,
                           norm_method = "log"
)
Mono.cds = orderCells(Mono.cds)
save(Mono.cds, file = "mycds_order_mono.t2d.RData")
load("mycds_order_mono.t2d.RData")
table(Mono.cds$Subtype)
Mono.cds$Subtype = factor(Mono.cds$Subtype, levels = c('0','1','2','3','4','5'))
Mono.cds$Subtype
p = plot_cell_trajectory(Mono.cds,color_by="Pseudotime",cell_size = 1, theta = 180,
                         size=1, show_backbone=TRUE, show_branch_points = F) +
  
  theme(text = element_text(size = 14,family = "Times"))
p
ggsave("03.PSD_Trajectory_mono.png", width = 6, height = 6, units = "in", dpi = 300)
ggsave("03.PSD_Trajectory_mono.pdf", width = 6, height = 6, units = "in", dpi = 300)
p2 = plot_cell_trajectory(Mono.cds,color_by="Subtype",cell_size = 0.5, theta = 180,
                          size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','#89da6e')) +
  theme(text = element_text(size = 14,family = "Times"))
p2
ggsave("04.cell_type_Trajectory_mono.png", width = 6, height = 6, units = "in", dpi = 300)
ggsave("04.cell_type_Trajectory_mono.pdf", width = 6, height = 6, units = "in", dpi = 300)
p3 = plot_cell_trajectory(Mono.cds,color_by="State",cell_size = 1, theta = 180,
                          size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(text = element_text(size = 14,family = "Times"))
p3
ggsave("05.State_Trajectory_mono.png", width = 6, height = 6, units = "in", dpi = 300)
ggsave("05.State_Trajectory_mono.pdf", width = 6, height = 6, units = "in", dpi = 300)
hubgene <- read.csv('/data/nas1/anbingzhuang/project/17_BJTC-596-3/04.Cox_LASSO/03.cox.csv')
hubgene <- hubgene$x
p4 <- plot_genes_in_pseudotime(Mono.cds[hubgene,],
                               color_by = "State",
                               ncol = 2)+
  theme(axis.title = element_text(face = "bold",family = "Times"),
        strip.text = element_text(face = "bold",family = "Times"))
p4
ggsave(filename = '04.sdheat_hubgene_mono.pdf',p4,w=8,h=4)
ggsave(filename = '04.sdheat_hubgene_mono.png',p4,w=8,h=4,units = 'in',dpi = 300)
rm(list = ls())
options(stringsAsFactors = F)
if(!dir.exists("09_DEGs")){dir.create("09_DEGs")}
load('../03_PCA/UMAP.Rdata')
table(UMAP$seurat_clusters)
UMAP <- subset(UMAP, seurat_clusters %in% c(7,13,14))
UMAP$celltype.group <- paste(UMAP$seurat_clusters, UMAP$Type, sep = "_")
UMAP$celltype <- Idents(UMAP)
Idents(UMAP) <- "celltype.group"
table(UMAP$celltype)
table(UMAP$celltype.group)
i <- 1
DE_gene <- data.frame()
DEGs_res_all <- data.frame()
for (i in c(13)) {
  temp1 <- FindMarkers(UMAP,
                       ident.1=paste0(i,'_AML'),
                       ident.2 = paste0(i,'_Healthy'),
                       test.use = 'wilcox',
                       min.pct = 0.1,
                       verbose = FALSE)
  write.csv(temp1, paste0(i,"_DE.csv"))
  temp1$cell <- i
  DEGs_res_all <- rbind(DEGs_res_all, temp1)
  
  temp1$symbol <- rownames(temp1)
  temp1 <- dplyr::select(temp1, symbol)
  DE_gene <- rbind(DE_gene, temp1)
}
write.csv(DEGs_res_all, file = './DEGs_res_all.csv')
DE_gene <- DEGs_res_all[DEGs_res_all$p_val < 0.05 & abs(DEGs_res_all$avg_log2FC) > 1, ]
DE_gene$symbol <- rownames(DE_gene)
DE_gene <- DE_gene[!duplicated(DE_gene$symbol), ]
write.csv(DE_gene, file = './DE_gene_res_13.csv')
DE_gene <- read.csv('./DE_gene_res_13.csv', row.names = 1)
colnames(DE_gene)[7] <- 'symbol'
all.markers <- FindAllMarkers(UMAP, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5, test.use = 'wilcox', return.thresh = 0.01)
write.csv(all.markers,file = 'TSNE_AllMarkers_13.csv')
TSNE_AllMarkers <- read.csv('./TSNE_AllMarkers_13.csv', row.names = 1)
malignant_gene <- TSNE_AllMarkers[TSNE_AllMarkers$cluster == 'malignant cells', ]
malignant_gene <- dplyr::select(malignant_gene, gene)
colnames(malignant_gene) <- 'symbol'
finally_DE_gene <- rbind(malignant_gene, DE_gene)
finally_DE_gene <- finally_DE_gene[!duplicated(finally_DE_gene$symbol), ] %>% as.data.frame()
colnames(finally_DE_gene) <- 'symbol'
write.csv(finally_DE_gene, file = './finally_DE_gene_13.csv')
DEGs_res_all$change <- as.factor(
  ifelse(DEGs_res_all$p_val < 0.05 & abs(DEGs_res_all$avg_log2FC) > 1,
         ifelse(DEGs_res_all$avg_log2FC > 1,'Up','Down'),'Not'))
table(DEGs_res_all$change)
DEGs_res_all
p1 <- ggplot(DEGs_res_all, aes(x = -log(p_val), y = avg_log2FC, color = change), repel = T)+
  geom_point()+
  theme_classic2()+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "#909090",color = "black"),
        strip.text = element_text(face="bold",color = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=17,face = "bold"))+
  scale_color_manual(values = c("darkgreen","grey","firebrick"))+
  facet_wrap(~cell,ncol = 5,nrow = 2)+
  geom_hline(yintercept = c(1,-1))
p1
ggsave('./01.13_volcano.png', p1,width = 12, height = 8)
ggsave('./01.13_volcano.pdf', p1,width = 12, height = 8,dpi = 600)
rm(list = ls())
options(stringsAsFactors = F)
if(!dir.exists("09_DEGs")){dir.create("09_DEGs")}
load('../03_PCA/UMAP.Rdata')
table(UMAP$seurat_clusters)
UMAP <- subset(UMAP, seurat_clusters %in% c(7,13,14))
UMAP$celltype.group <- paste(UMAP$seurat_clusters, UMAP$Type, sep = "_")
UMAP$celltype <- Idents(UMAP)
Idents(UMAP) <- "celltype.group"
table(UMAP$celltype)
table(UMAP$celltype.group)
i <- 1
DE_gene <- data.frame()
DEGs_res_all <- data.frame()
for (i in c(7)) {
  temp1 <- FindMarkers(UMAP,
                       ident.1=paste0(i,'_AML'),
                       ident.2 = paste0(i,'_Healthy'),
                       test.use = 'wilcox',
                       min.pct = 0.1,
                       verbose = FALSE)
  write.csv(temp1, paste0(i,"_DE.csv"))
  temp1$cell <- i
  DEGs_res_all <- rbind(DEGs_res_all, temp1)
  
  temp1$symbol <- rownames(temp1)
  temp1 <- dplyr::select(temp1, symbol)
  DE_gene <- rbind(DE_gene, temp1)
}
write.csv(DEGs_res_all, file = './DEGs_res_7.csv')
DE_gene <- DEGs_res_all[DEGs_res_all$p_val < 0.05 & abs(DEGs_res_all$avg_log2FC) > 1, ]
DE_gene$symbol <- rownames(DE_gene)
DE_gene <- DE_gene[!duplicated(DE_gene$symbol), ]
write.csv(DE_gene, file = './DE_gene_res_7.csv')
DE_gene <- read.csv('./DE_gene_res_7.csv', row.names = 1)
colnames(DE_gene)[7] <- 'symbol'
all.markers <- FindAllMarkers(UMAP, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5, test.use = 'wilcox', return.thresh = 0.01)
write.csv(all.markers,file = 'TSNE_AllMarkers_7.csv')
TSNE_AllMarkers <- read.csv('./TSNE_AllMarkers_7.csv', row.names = 1)
malignant_gene <- TSNE_AllMarkers[TSNE_AllMarkers$cluster == 'malignant cells', ]
malignant_gene <- dplyr::select(malignant_gene, gene)
colnames(malignant_gene) <- 'symbol'
finally_DE_gene <- rbind(malignant_gene, DE_gene)
finally_DE_gene <- finally_DE_gene[!duplicated(finally_DE_gene$symbol), ] %>% as.data.frame()
colnames(finally_DE_gene) <- 'symbol'
write.csv(finally_DE_gene, file = './finally_DE_gene_7.csv')
DEGs_res_all$change <- as.factor(
  ifelse(DEGs_res_all$p_val < 0.05 & abs(DEGs_res_all$avg_log2FC) > 1,
         ifelse(DEGs_res_all$avg_log2FC > 1,'Up','Down'),'Not'))
table(DEGs_res_all$change)
DEGs_res_all
p1 <- ggplot(DEGs_res_all, aes(x = -log(p_val), y = avg_log2FC, color = change), repel = T)+
  geom_point()+
  theme_classic2()+
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "#909090",color = "black"),
        strip.text = element_text(face="bold",color = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=17,face = "bold"))+
  scale_color_manual(values = c("darkgreen","grey","firebrick"))+
  facet_wrap(~cell,ncol = 5,nrow = 2)+
  geom_hline(yintercept = c(1,-1))
p1
ggsave('./01.7_volcano.png', p1,width = 12, height = 8)
ggsave('./01.7_volcano.pdf', p1,width = 12, height = 8,dpi = 600)
rm(list = ls())
library(Seurat)
library(monocle)
library(dplyr)
if (! dir.exists("./10_subcell")){
  dir.create("./10_subcell")
}
load('../05_SingleR/UMAP2name.Rdata') 
library(monocle)
if (! dir.exists("./T cell")){
  dir.create("./T cell")
}
sce.car <-  subset(UMAP2name, seurat_clusters  %in%c('T cell') )
sce.car <- ScaleData(sce.car)
sce.car <- RunPCA(sce.car, features = VariableFeatures(object = sce.car))
ElbowPlot(sce.car, ndims = 50)
s_cmb<-FindNeighbors(sce.car,dims = 1:20)
library(clustree)
resolutions <- seq(0.1,1,0.1)
for (res in resolutions) {
  s_cmb <- FindClusters(s_cmb, resolution = res, algorithm = 1, verbose = FALSE)
}
clustree(s_cmb@meta.data,prefix='RNA_snn_res.')
ggsave("00.clustreePlot20.png",width = 12,height = 9,dpi = 300)
ggsave("00.clustreePlot20.pdf",width = 12,height = 9)
sce.car <- FindNeighbors(sce.car, dims = 1:20)
sce.car <- FindClusters(sce.car, resolution = 0.6)
table(sce.car@meta.data$seurat_clusters)
sce.car <- RunUMAP(sce.car, dims = 1:30)
p <- DimPlot(sce.car, reduction = "umap", group.by = "seurat_clusters", label = T) + ggtitle(paste0('Subtype of ', 'CD4+Tcell' ))
pdf(file = paste0('01.CD4+Tcell', '.subtype.pdf'),w=6,h=5)
print(p)
dev.off()
png(file = paste0('01.CD4+Tcell', '.subtype.png'),w=500,h=400)
print(p)
dev.off()
all.markers <- FindAllMarkers(sce.car, only.pos = TRUE,min.pct = 0.1,logfc.threshold = 0.5,test.use = 'wilcox',return.thresh = 0.05)
marker_gene <- all.markers %>% group_by(cluster) %>% top_n(n=3,wt=avg_log2FC)
write.csv(marker_gene, file = paste0("T_cell", '_marker(subtype_car).csv'))
marker_gene <- marker_gene[!duplicated(marker_gene$gene), ]
gene <- c("CCR7","CCR6","GZMA",
          "CD4","CD3E","IL7R",
          "CD8A","CD8B","GZMB", "PRF1",
          "TCRD","NKG7" 
) 
gene <- c(
          "CD3D","CD3E",
          
          "CD4","CCR7","SELL","TCF7","IL2","TNF","IL4R","FOXP3","IL2RA","CTLA4",
          "CCL5","GZMB","GZMK","IFNG","CCL4","CCL3","PRF1","NKG7")
p2 <- DotPlot(sce.car, features = gene,group.by = 'seurat_clusters') + RotatedAxis()+
  labs(title="", y="Subtype", x = "",size=40)
p2
pdf(file = paste0('02.Tcell', '.marker.exp.pdf'),w=8,h=6,family = "Times")
print(p2)
dev.off()
png(file = paste0('02.Tcell', '.marker.exp.png'),w=600,h=450,family = "Times")
print(p2)
dev.off()
new.cluster.id = c(
  "0" = "CD4+ T cell",
  "1" = "Other T cell",
  "2" = "CD8+ T cell",
  "3" = "Other T cell",
  "4" = "CD8+ T cell",
  "5" = "CD4+ T cell",
  "6" = "CD4+/CD8+ T cell"
)
Idents(sce.car) <- new.cluster.id[as.character(Idents(sce.car))]
sce.car$celltype <- Idents(sce.car) 
p <- DimPlot(sce.car, reduction = "umap", group.by = "celltype", label = T) + ggtitle(paste0('Subtype of ', 'T cell' ))
pdf(file = paste0('03.Tcell', '.celltype.pdf'),w=8,h=6,family = "Times")
print(p)
dev.off()
png(file = paste0('03.Tcell', '.celltype.png'),w=8,h=6,units = "in",res = 600,family = "Times")
print(p)
dev.off()
gene<- c("CDK6","CALR","HOXA9","PARP1")
p <- FeaturePlot(sce.car, features = gene, ncol = 2,repel = 'celltype',col=c('grey','red')
                 
)
p
png("04.expression_signature.png",w=8,h=8,units = "in",res = 600,family = 'Times')
p
dev.off()
pdf("04.expression_signature.pdf",w=8,h=8,family = 'Times')
p
dev.off()
Idents(sce.car)<-'celltype'
p <- VlnPlot(sce.car, features = gene, ncol = 4, pt.size = 1)
p
png("05.expression_signature.png",w=16,h=8,units = "in",res = 600,family = 'Times')
p
dev.off()
pdf("05.expression_signature.pdf",w=16,h=8,family = 'Times')
p
dev.off()
theme.set = theme(
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 14,  face = "bold", family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 16, face = "bold", family = "Times"),
  legend.title = element_blank(),
  text = element_text(family = "Times"))
p <- DotPlot(UMAP2name,features = gene)+theme.set+ 
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
p
png("06.expression_signature.png",w=14,h=8,units = "in",res = 600,family = 'Times')
p
dev.off()
pdf("06.expression_signature.pdf",w=14,h=8,family = 'Times')
p
dev.off()
table(UMAP2name$celltype)
dat.tpm <- as.data.frame(UMAP2name@assays$RNA$data) %>% apply(2,function(x){x/sum(x) * 10000})
dat.tpm <- log2(dat.tpm+1)
save.image(file = "dat.tpm.RData")
load("dat.tpm.RData")
i <- 1
library(ggpubr)
DE_cell <- data.frame()
for (i in c(1 : length(gene))) {
  
  var <- gene[i]
  temp1 <- dat.tpm[rownames(dat.tpm) %in% gene[i], ] %>% as.data.frame()
  
  
  names(temp1) <- 'expr'
  temp1$sample <- rownames(temp1)
  group <- data.frame(sample = colnames(UMAP2name), Group = UMAP2name@meta.data$group)
  celltype <- data.frame(sample = colnames(UMAP2name), celltype = UMAP2name$celltype)
  dat <- merge(temp1, group, by = "sample")
  dat <- merge(dat, celltype, by = "sample")
  table(dat$expr)  
  p1 <- ggplot(dat, aes(x = Group,
                        y = expr,
                        fill = Group)) +
    geom_point()+
    geom_violin(trim=T,color="black") + 
    
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,fill='white')+ 
    scale_fill_manual(values= c("#20B2AA","#F08080"), name = "Group")+
    labs(title="", x="", y = "expr",size=20) +
    stat_compare_means(data = dat,
                       mapping = aes(group = Group),
                       label ="p.signif",
                       method = 'wilcox.test',
                       paired = F,label.x = 1.4) +
    theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
          axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10),
          axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12),
          axis.title.x=element_text(size=16,face="bold"),
          axis.title.y=element_text(size=16,face="bold"),
          legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=8),
          legend.title = element_text(face = "bold", size = 10),
          legend.position = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+facet_wrap(~celltype,scales = "free",nrow = 2) +
    guides(fill='none')
  p1
  
  
  
  ggsave(paste0('0', i+3, gene[i], '_cell_DE.pdf'),p1,w=14,h=8)
  ggsave(paste0('0', i+3, gene[i], '_cell_DE.png'),p1,w=14,h=8)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}
