rm(list=ls())
library(GEOquery)
 
if (!dir.exists("01_DEG_GSE114868")) {dir.create("01_DEG_GSE114868")}
 
library(magrittr)
train_expr <- read.csv("../00_Rawdata/GSE114868_expr_use.csv", header = T, sep = ",",row.names = 1)
group <- read.csv("../00_Rawdata/GSE114868_group.csv", header = T, sep = ",",row.names = 1)
library(limma)
group_list <- group
expr <- as.data.frame(lapply(train_expr, as.numeric))
rownames(expr) <- rownames(train_expr)
group_list$group <- factor(group_list$group)
design <- model.matrix(~0 + group_list$group)
rownames(design) <- group_list$sample
colnames(design) <- levels(group_list$group)
table(group$group)
compare <- makeContrasts("AML-Control", levels = design)
compare
fit <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
res <- topTable(fit3, coef = 1, number = Inf, sort.by = "logFC", 
                adjust.method = "BH")
dim(res)
logfc <- 1
res$change <- as.factor(
  ifelse(res$P.Value < 0.05 & abs(res$logFC) > logfc,
         ifelse(res$logFC > logfc,'Up','Down'),'Not')
)
table(res$change)
sig_res <- topTable(fit3, coef = 1, number = Inf, sort.by = "logFC",
                    adjust.method = "none", p.value = 0.05, lfc = logfc)
dim(sig_res)
sig_res$change <- as.factor(
  ifelse(sig_res$P.Value < 0.05 & abs(sig_res$logFC) > logfc,
         ifelse(sig_res$logFC > logfc,'Up','Down'),'Not')
)
table(sig_res$change)
write.table(res, file = "01.GSE114868_DEG_all_res_CA.csv", quote = F, sep = ",")
write.table(sig_res, file = "02.GSE114868_DEG_sig_res_CA.csv", quote = F, sep = ",")
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
dat_rep <- res[rownames(res)%in%
                 rownames(rbind(head(sig_res[order(sig_res$logFC,decreasing = T),],10),
                                head(sig_res[order(sig_res$logFC,decreasing = F),],10))),]
volcano_plot <- ggplot(data = res, 
                       aes(x = logFC,
                           y = -log10(P.Value), 
                           color =change)) +
  scale_color_manual(values = c("#20B2AA", "darkgray","red")) +
  scale_x_continuous(breaks = c(-6,-3,-1,0,1,3,6)) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-1,1),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 15),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 15)) +
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 4,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log2 (Fold Change)",
       y = "-log10 (padj)")
volcano_plot
ggsave('03.GSE114868_volcano.png', volcano_plot,width = 8, height = 6)
ggsave('03.GSE114868_volcano.pdf', volcano_plot,width = 8, height = 6)
library(ComplexHeatmap)
rt <- train_expr
rt <- rt[rownames(sig_res),]
heat<-rt[rownames(rt)%in%
           c(head(rownames(subset(sig_res,sig_res$logFC>0.5)),10),head(rownames(subset(sig_res,sig_res$logFC< -0.5)),10)),]
heat[,1:77] <- as.data.frame(lapply(heat[,1:77],as.numeric))
x<-log2(heat+1)
mat <- t(scale(t(x)))
df1 <- as.data.frame(mat)
mat[mat < (-2)] <- (-2)
pdf('04.GSE114868_heatmap.pdf',  w=6,h=6,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(6, "cm")) %v%
  HeatmapAnnotation(Group = group$group, col = list(Group = c("AML" = "#B72230", "Control" = "#104680"))) %v%
  Heatmap(mat, 
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          
          name = "expression", 
          
          cluster_rows = T,
          height = unit(6, "cm"),
          
          
          col = colorRampPalette(c("#20B2AA", "white","red"))(100))
dev.off()
png('04.GSE114868_heatmap.png',w=6,h=6,units='in',res=600,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(6, "cm")) %v%
  HeatmapAnnotation(Group = group$group, col = list(Group = c("AML" = "#B72230", "Control" = "#104680"))) %v%
  Heatmap(mat, 
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          
          name = "expression", 
          
          cluster_rows = T,
          height = unit(6, "cm"),
          
          
          col = colorRampPalette(c("#20B2AA", "white","red"))(100))
dev.off()
