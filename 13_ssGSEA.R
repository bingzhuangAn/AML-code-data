rm(list = ls())
if (!dir.exists("13.CIBERSORT")) {dir.create("13.CIBERSORT")}
library(tidyverse)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(immunedeconv)
library(GSVA)
expr <- read.csv('../04.Cox_LASSO/dat_fpkm.csv',row.names = 1,check.names = F)
colnames(expr) <- gsub("\\.","-",colnames(expr))
riskScore <- read.table("../05_riskModel/risk.csv", header = T,sep = ",",row.names = 1)
train_group <-riskScore[,c("sample","risk")]
gene_set <- read.table("./mmc3.txt",header = T,sep ="\t")
dat.final <- as.matrix(expr)
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])
ssgsea_score = gsva(dat.final, gene_list, 
                    method = "ssgsea", 
                    ssgsea.norm = TRUE, 
                    verbose = TRUE)
write.csv(ssgsea_score,
          file = "01.ssgsea_result_cell.csv",
          quote = F)
colnames(train_group) <- c('sample', 'group')
train_group<-train_group[order(train_group$group),]
ssgsea_score<-read.csv('./01.ssgsea_result_cell.csv',check.names = F, row.names = 1)
ssgsea_score<-ssgsea_score[,train_group$sample]
annotation_col<-as.data.frame(train_group$group)
colnames(annotation_col)='Group'
rownames(annotation_col)=colnames(ssgsea_score)
library(pheatmap)
color.key<-c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
ann_colors<-list(Group = c('Low'="#45a9b8",'High'="#f76a56"))
p <- pheatmap(
  ssgsea_score,
  color = colorRampPalette(color.key)(50),
  border_color = 'darkgrey',
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  labels_row = NULL,
  clustering_method = 'ward.D2',
  show_rownames = T,
  show_colnames = F,
  fontsize_col = 5,
  cluster_cols = F,
  cluster_rows = T)
p
png(filename = "01.heatmap.png", height = 5, width = 8, units = 'in',res = 600,family='Times')
p
dev.off()
pdf(file = "01.heatmap.pdf", height = 5,width = 8,family='Times')
p
dev.off()
colnames(train_group) <- c('sample', 'group')
group_case<-train_group[train_group$group=='High',]
group_case<-as.character(group_case$sample)
group_control<-train_group[train_group$group=='Low',]
group_control<-as.character(group_control$sample)
tiics_result <- read.csv('./01.ssgsea_result_cell.csv',check.names = F, row.names = 1) 
tiics_result <- tiics_result[, train_group$sample] %>% as.matrix()
pvalue = padj <- matrix(0, nrow(tiics_result), 1)
for (i in 1:nrow(tiics_result)){
  pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, group_case],
                                       tiics_result[i, group_control])$p.value
  
  
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(
  pvalue, 
  padj,
  row.names = rownames(tiics_result))
rTable$immune_cell <- rownames(rTable)
rTable$sig <- ifelse(rTable$pvalue < 0.05,
                     ifelse(rTable$pvalue < 0.01, 
                            ifelse(rTable$pvalue < 0.001,
                                   ifelse(rTable$pvalue < 0.0001,
                                          paste(rTable$immune_cell, "****",  sep = ""),
                                          paste(rTable$immune_cell, "***", sep = "")),
                                   paste(rTable$immune_cell, "**", sep = "")),
                            paste(rTable$immune_cell, "*",  sep = "")), 
                     rTable$immune_cell)
diff_Table<-rTable[which(rTable$pvalue<0.05),]
dim(diff_Table)
write.csv(rTable,
          file = "02.cell_all_wilcox_test.csv",
          quote = F,
          row.names = F)
write.csv(diff_Table,
          file = "02.cell_diff_wilcox_test.csv",
          quote = F,
          row.names = F)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
cell.data <- data.frame(Immune_Cell=rownames(tiics_result), 
                        tiics_result, 
                        pvalue=rTable$pvalue)
plot.cell <- cell.data[which(cell.data$pvalue<0.05),]
diff_tiics <- rownames(plot.cell)
violin_dat <- gather(plot.cell, key=Group, value=score, -c("Immune_Cell","pvalue"))
violin_dat$Group <- ifelse(gsub("\\.","-",violin_dat$Group) %in% group_control,
                           "Low risk", "High risk") 
violin_dat$Group <- factor(violin_dat$Group, levels = c("Low risk", "High risk"))
violin_dat <- violin_dat[,-2]
head(violin_dat)
boxplot_diff_TIICs <- ggplot(violin_dat, aes(x=Immune_Cell, 
                                             y=score,
                                             fill=Group)) +
  
  
  
  
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = 21,
               outlier.fill = "black",
               outlier.size = 0.5,
  )+ 
  
  
  
  scale_fill_manual(values= c("#45a9b8","#f76a56"))+ 
  labs(title="", x="", y = "Score",size=20) +
  stat_compare_means(data = violin_dat,
                     method = "wilcox.test", 
                     mapping = aes(group = Group),
                     label ="p.signif",
                     hide.ns = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18,family = "Times"), 
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10,family = "Times"), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=12,family = "Times"), 
        axis.title.x=element_text(size=16,face="bold",family = "Times"),
        axis.title.y=element_text(size=14,face="bold",family = "Times"), 
        legend.text=element_text(face="bold", hjust = 0.5,colour="black", size=11,family = "Times"), 
        legend.title=element_text(face="bold", colour="black", size=11,family = "Times"),
        text=element_text(family = 'Times'),
        
        
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
boxplot_diff_TIICs
ggsave('03.ssgsea_Box.pdf',boxplot_diff_TIICs,w=15,h=6)
ggsave('03.ssgsea_Box.png',boxplot_diff_TIICs,w=15,h=6)
tiics_result <- read.csv('./01.ssgsea_result_cell.csv',check.names = F, row.names = 1)
tiics_result <- t(tiics_result) %>% as.data.frame()
diff <- read.csv('./02.cell_diff_wilcox_test.csv')
tiics_result <- tiics_result[,diff$immune_cell]
res3 <- tiics_result
library(ggcorrplot)
cor_data <- cor(res3,method="spearman")
corp <- cor_pmat(res3)
write.csv(cor_data,'04.cell_cor_r.csv',quote=F)
write.csv(corp,'04.cell_cor_p.csv',quote=F)
corp <-data.frame(corp)
rownames(corp) <- corp$rowname
corp <- corp[,-1]
colnames(corp) <-gsub("."," ",colnames(corp))
corp <-as.matrix(corp)
library(ggcorrplot)
library(corrplot)
pdf("04.cell_corHeatmap.pdf",height=14,width=14,family='Times')
col1 <- colorRampPalette(c("#3300CC", "white", "#CC0000"))
corr_plot<-corrplot(cor_data,
                    tl.cex=1.5,
                    method = "square", 
                    is.corr = T,
                    type = "full", 
                    p.mat = corp,
                    insig = "blank",
                    outline = "white",
                    addCoef.col ="black",
                    col = col1(200),
                    tl.col = 'black',
                    tl.offset = 0.4,
                    number.font = 2,
                    cl.cex = 1.5,
                    number.cex = 0.7
)
dev.off()
png("04.cell_corHeatmap.png",height=14,width=14,family='Times',units='in',res=600)
col1 <- colorRampPalette(c("#3300CC", "white", "#CC0000"))
corr_plot<-corrplot(cor_data,
                    tl.cex=1.5,
                    method = "square", 
                    is.corr = T,
                    type = "full", 
                    p.mat = corp,
                    insig = "blank",
                    outline = "white",
                    addCoef.col ="black",
                    col = col1(200),
                    tl.col = 'black',
                    tl.offset = 0.4,
                    number.font = 2,
                    cl.cex = 1.5,
                    number.cex = 0.7
)
dev.off()
expr1<-expr
model_coef <- read.table("../04.Cox_LASSO/03.cox.csv", header = T,sep = ',')
model_gene <- model_coef$x
expr1 <- expr1[model_gene,] %>% t %>% as.data.frame()
expr1$sample<-rownames(expr1)
expr1 <-merge(expr1,riskScore,id="sample")
rownames(expr1)<-expr1$sample
expr1 <-expr1[,!colnames(expr1) %in% c("sample","risk","OS","OS.time","riskScore","risk")]
tiics_result <- read.csv('./01.ssgsea_result_cell.csv',check.names = F, row.names = 1) 
tiics_result <- tiics_result[, train_group$sample] %>% as.matrix()
tiics_result <- t(tiics_result) %>% as.data.frame()
tiics_result <- tiics_result[rownames(expr1),]
diff <- read.csv('./02.cell_diff_wilcox_test.csv')
tiics_result <- tiics_result[,diff$immune_cell]
identical(rownames(expr1), rownames(tiics_result))
expr1<- as.data.frame(lapply(expr1,as.numeric))
cor_r <- cor(expr1,tiics_result,method = "spearman") 
library(psych)
d <- corr.test(expr1,tiics_result,use="complete",method = 'spearman')
cor_p <- d$p
cor_r2 <- cor_r %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>%
  tidyr::gather(., cell,Correlation,-gene)
cor_p2 <- cor_p %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>% 
  tidyr::gather(., cell, Pvalue, -gene)
cor_dat <- cbind(cor_r2, cor_p2)[,c("gene","cell","Correlation","Pvalue")]
write.csv(cor_dat,"05.correlation_cor.csv")
data <- read.csv('05.correlation_cor.csv',row.names = 1,check.names = F)
data <- data %>%
  mutate(text = case_when( 
    Pvalue <= 0.001 ~ "\n***", 
    between(Pvalue, 0.001, 0.01) ~ "\n**", 
    between(Pvalue, 0.01, 0.05) ~ "\n*",  
    T ~ ""))
p <- 
  ggplot(data, aes(gene, cell)) + 
  geom_tile(aes(fill = Correlation), colour = "white", size = 1, family = "Times")+
  scale_fill_gradient2(low = "#5C5DAF",mid = "white",high = "#EA2E2D") +
  geom_text(aes(label = text),col ="black",size = 3, family = "Times") +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.x = element_text(hjust = 1, size = 9, face = "bold", family = "Times"), 
        axis.text.y = element_text(size = 10, face = "bold", family = "Times")) + 
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"), family = "Times") +   
  scale_x_discrete(position = "top") 
p
ggsave(file=paste0('05.correlation_cor.pdf'), height = 6, width = 10, p, family = "Times")
ggsave(file=paste0('05.correlation_cor.png'), height = 6, width = 10, p)
if (! dir.exists("./00_ESTIMATE")){
  dir.create("./00_ESTIMATE")
}
library(magrittr)
library(dplyr)
library(plyr)
dat <- expr
group<- riskScore[,c('sample','riskScore','risk')]
colnames(group)<-c('sample','riskScore','label')
Low.sample<-group$sample[which(group$label=='Low')]
High.sample<-group$sample[which(group$label=='High')]
dat<-dat[,group$sample]
group_estimate<-group$label%>%as.factor()
design<-model.matrix(~0 + group_estimate)
rownames(design)<-group$sample
colnames(design)<-levels(group_estimate)
design<-as.data.frame(design)
Low<-rownames(design)[which(design$`Low`==1)]
High<-rownames(design)[which(design$`High`==1)]
length(Low)
length(High)
library(estimate)
expr_train <- dat
write.csv(expr_train, 
          '01.expr.csv', 
          col.names = T, 
          row.names = T, 
          quote = F, sep="\t")
write.table(expr_train, 
            'expr.txt', 
            col.names = T, 
            row.names = T, 
            quote = F, sep="\t")
filterCommonGenes(input.f = './expr.txt', 
                  output.f = 'expr_train.gct', 
                  id = 'GeneSymbol')
estimateScore('expr_train.gct', 'train_purity.gct', platform="affymetrix")
es_score <- read.table('train_purity.gct', skip = 2, header = T, check.names = F)
immu_score <- es_score[,3:length(es_score)]
rownames(immu_score) <- es_score$NAME
write.csv(es_score,
          file = "02.es_score.csv",
          quote = F,
          row.names = F)
violin_dat <- data.frame(t(immu_score))
rownames(violin_dat)<-group$sample
violin_dat$sample <- rownames(violin_dat)
violin_dat$group <- ifelse(violin_dat$sample %in% Low.sample,"Low Risk", "High Risk")
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
{
  p1 <- ggplot(violin_dat, aes(x=group,y=StromalScore, fill=group))+
    geom_violin() + 
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ 
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#FF6347","#6495ED"), name = "Group") + 
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ 
    
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), 
          axis.text.y=element_text(family="Times",size=12,face="bold"), 
          axis.title.y=element_text(family="Times",size = 15,face="bold"), 
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Stromal Score", x="", y="Stromal Score")+
    guides(fill='none')
  p1
  p2 <- ggplot(violin_dat, aes(x=group,y=ImmuneScore, fill=group))+
    geom_violin() + 
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ 
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#FF6347","#6495ED"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ 
    
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), 
          axis.text.y=element_text(family="Times",size=12,face="bold"), 
          axis.title.y=element_text(family="Times",size = 15,face="bold"), 
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),     
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Immune Score", x="", y="Immune Score")+
    guides(fill='none')
  p2
  p3 <- ggplot(violin_dat, aes(x=group, y=ESTIMATEScore, fill=group))+
    geom_violin() + 
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ 
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#FF6347","#6495ED"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ 
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), 
          axis.text.y=element_text(family="Times",size=12,face="bold"), 
          axis.title.y=element_text(family="Times",size = 15,face="bold"), 
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="ESTIMATE Score", x="", y="ESTIMATE Score")+
    guides(fill='none')
  p3
  p5 <- cowplot::plot_grid(p1,p2,p3,
                           nrow = 1, 
                           align = 'h', 
                           vjust = -0.3)
  p5
}
ggsave(filename = '03.estimate.all.pdf',p5,w=10,h=6)
ggsave(filename = '03.estimate.all.png',p5,w=10,h=6)
