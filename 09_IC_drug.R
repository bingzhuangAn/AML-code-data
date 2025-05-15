rm(list = ls())
if (!dir.exists("09_IC_drug")) {dir.create("09_IC_drug")}
library(tidyverse)
library(lance)
library(oncoPredict)
library(ggplot2)
data <- read.csv("../04.Cox_LASSO/dat_fpkm.csv",row.names = 1)
colnames(data) <- gsub("\\.","-",colnames(data))
FPKM_log <- data
sub_group <- read.csv("../05_riskModel/risk.csv",header = T)
sub_group <- sub_group[,c("sample","risk")]
colnames(sub_group) <- c("sample","group")
GDSC2_expr <- readRDS('/data/nas1/anbingzhuang/project/all_drrug/database/Drug/GDSC2_Expr (RMA Normalized and Log Transformed).rds')
GDSC2_res <- readRDS('/data/nas1/anbingzhuang/project/all_drrug/database/Drug/GDSC2_Res.rds')
train_data <- FPKM_log
train_data <- train_data[, colnames(train_data) %in% sub_group$sample]
IC50<-read.csv('calcPhenotype_Output/DrugPredictions.csv',row.names=1)
IC50$sample <- rownames(IC50)
sub_group$group <- as.factor(sub_group$group)
dat.IC50 <- merge(IC50, sub_group, by = "sample")
dat.IC50_2 <- dat.IC50 %>% 
  pivot_longer(
    cols = -c("sample", "group"),
    names_to = "drug",
    values_to = "Score"
  )
colnames(dat.IC50_2)
library(rstatix)
stat_res <- dat.IC50_2 %>% 
  group_by(drug) %>%
  wilcox_test(Score ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  
  add_significance("p")
DE.res <- stat_res[which(stat_res$p < 0.05), ] 
write.csv(stat_res, file = 'stat.IC50.csv')
write.csv(DE.res, file = 'DE.IC50.csv')
library(reshape)
write.csv(stat_res,'02.stat_res.csv')
dif_stat_res<-stat_res[stat_res$p<0.05,]
dif_stat_res<-dif_stat_res[order(dif_stat_res$p,decreasing = F),]
write.csv(dif_stat_res,'03.dif_stat_res.csv')
hubgene <- read.csv("../04.Cox_LASSO/03.cox.csv")
hubgene <- hubgene$x
dif_drug<-dif_stat_res$drug
library(stringr )
hub_eset<-train_data[hubgene,]
dif_IC50<-subset(IC50,select=dif_drug)
colnames(dif_IC50)<-str_split_fixed(colnames(dif_IC50),pattern = "_",2)[,1]
dif_drug <-str_split_fixed(dif_drug,pattern = "_",2)[,1]
dif_IC50<-dif_IC50[,order(colnames(dif_IC50),decreasing = F)]
identical(rownames(dif_IC50),colnames(hub_eset))
library(psych)
library(reshape2)
y <- hub_eset%>%t(.)
x <-dif_IC50
d <- corr.test(x,y,use="complete",method = 'spearman')
r <- data.frame(t(d$r))
r$Type<-rownames(r)
r<-melt(r,id='Type')
p <- data.frame(t(d$p))
p$Type<-rownames(p)
p<-melt(p,id='Type')
data<-cbind(r,p)[,c(-4,-5)]
colnames(data)<-c('Type','cell','correlation','p')
data$cell <- as.character(data$cell)  
data$cell <- factor(data$cell)
write.csv(data,file="04.hubgene_drug_cor.csv",quote=F,row.names=F)
data1 <- data[abs(data$correlation) > 0.3 & data$p < 0.05,]
data1$cell <- as.character(data1$cell)  
data1$cell <- factor(data1$cell)
write.csv(data1,file="04.hubgene_drug_cor0.3_P0.05.csv",quote=F,row.names=F)
data$corr <- ''
data$corr[which(data$correlation >0)]='postive'
data$corr[which(data$correlation <0)]='negative'
data$adjp <- ''
data$adjp[which(data$p>0.05)]='>0.05'
data$adjp[which(data$p <0.05)]='<0.05'
data$adjp[which(data$p <0.01)]='<0.01'
data$adjp[which(data$p <0.001)]='<0.001'
data$adjp[which(data$p <0.0001)]='<0.0001'
data$sper <- abs(data$correlation)
library(ggplot2)
library(extrafont)
ggplot() +
  geom_point(data = data, aes(x = cell, y = Type, size = sper, fill = adjp), color = "black", shape = 21) + 
  geom_point(data = data[which(data$corr == 'postive'),], aes(x = cell, y = Type, size = abs(correlation), color = adjp), shape = 16) + 
  scale_fill_manual(values = c("#333366", "#336699", "#3399CC", "#66CCFF", "#CCCCCC")) + 
  scale_color_manual(values = c("#FF6666", "#FF9999", "#FFCCCC", "#FFCCCC", "#CCCCCC")) + 
  labs(x = NULL, y = NULL) + 
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12, family = "Times"),
    axis.text = element_text(size = 10, color = "black", family = "Times"),
    axis.text.x = element_text(size = 8, angle = 90, vjust = .5, hjust = 1, family = "Times"),
    axis.ticks = element_blank(),
    legend.margin = margin(0.1, unit = "pt"),
    legend.text = element_text(family = "Times"),
    legend.title = element_text(family = "Times")
  ) +
  guides(
    size = guide_legend(title = "Spearman p", order = 1),
    fill = guide_legend(title = "Negative\ncorrelation\np-value", order = 2, override.aes = list(size = 5)),
    color = guide_legend(title = "Postive\ncorrelation\np-value", order = 3, override.aes = list(size = 5))
  )
ggsave('05.drug_hubgene.png', width = 15, height = 7)
ggsave('05.drug_hubgene.pdf', width = 15, height = 7)
sub_group$group<-ifelse(sub_group$group== "Low",'Low risk','High risk')
dif_IC50$sample <- rownames(dif_IC50)
data<-merge(sub_group,dif_IC50,by='sample')
library(gghalves)
mycolor <- c('#7599bc','#2c648a')
j <- dif_drug[1:5]
colnames(data)
library(rlang)
library(ggplot2)
library(Cairo)
library(dplyr)
format_p_value <- function(p) {
  ifelse(p < 0.001, sprintf("p < %.0e", p), sprintf("p = %.3f", p))
}
for(i in j){
  data1 <- data[,c("group",i)]
  
  median_by_group <- tapply(data1[[i]], data1$group, median)
  
  
  df_median <- data.frame(
    group = names(median_by_group),
    median_v = median_by_group
  )
  
  p_value <- wilcox.test(data[[i]] ~ data$group)$p.value
  formatted_p <- format_p_value(p_value)
  
  p <- ggplot(data, aes(x = group, y = !!sym(i), fill = group)) +
    scale_color_manual(values = rev(mycolor)) +
    scale_fill_manual(values = rev(mycolor)) +
    geom_half_violin(position = position_nudge(x = 0.4, y = 0),
                     side = 'R', adjust = 1.2, trim = FALSE, color = NA, alpha = 1) +
    geom_point(aes(x = group, y = !!sym(i), color = group),
               position = position_jitter(width = 0.04), size = 0.2, shape = 20) +
    geom_boxplot(outlier.shape = NA, 
                 width = 0.1, alpha = 0.7, position = position_nudge(x = 0.2, y = 0)) +
    geom_point(data = df_median,aes(x=group,y=median_v,group=1),color="red",size=3,position = position_nudge(x = 0.2, y = 0))+
    geom_line(data = df_median,
              mapping = aes(x = group, y = median_v, group=1),
              color="red",size = 1,position = position_nudge(x = 0.2, y = 0))+
    theme_classic() +
    coord_flip() +
    labs(y = paste0(i, ' '), x = '') +
    theme(text = element_text(family = "Times", color = "gray30", size = 12),
          axis.line = element_line(size = 0.6, color = "gray30"),
          axis.ticks = element_line(size = 0.6, color = "gray30"),
          axis.ticks.length = unit(1.5, units = "mm"),
          axis.title.y = element_text(face = "bold", size = 24,family = "Times"), 
          axis.title.x = element_text(face = "bold", size = 24,family = "Times"), 
          axis.text.y = element_text(face = "bold", size = 24,family = "Times"), 
          axis.text.x = element_text(face = "bold", size = 24,family = "Times"),
          legend.position = 'none',
          legend.text = element_text(size = 24),
          plot.margin = unit(x = c(0.2, 0.2, 0.2, 0.2), units = "inches"),
          panel.grid.major.x = element_line(color = "gray80", size = 0.5)) +
    annotate("text", x = 1.9, y = (max(data[[i]]))*1.01, label = formatted_p, size = 10,family = "Times")  
  
  pdf(paste0(i, '_expression_boxplot_sequence.pdf'), width = 10, height = 7)
  print(p)
  dev.off()
  
  png(paste0(i, '_expression_boxplot_sequence.png'), width = 10, height = 7, units = 'in', res = 600)
  print(p)
  dev.off()
}
cor1 <- read.csv("./04.hubgene_drug_cor0.3_P0.05.csv")
FDA_drug <- c("Dasatinib","Erlotinib","Lapatinib","Savolitinib")
cor1$cell %>% unique() 
cor2 <- cor1[abs(cor1$correlation) > 0.3,]
cell_counts <- as.data.frame(table(cor2$cell))
colnames(cell_counts) <- c("cell", "count")
drugs_with_count_3 <- subset(cell_counts, count == 3)
dif_drug <- c("Lapatinib","Gefitinib","Erlotinib","Docetaxel","Afatinib")
hubgene <- c("AHNAK","CSPG4","NCAM1")
dif_stat_res$drug[dif_stat_res$drug == "Nutlin.3a...._1047"] <- "Nutlin-3a (-)_1047"
colnames(IC50)[colnames(IC50) == "Nutlin.3a...._1047"] <- "Nutlin-3a (-)_1047"
library(stringr )
hub_eset<-train_data[hubgene,]
colnames(IC50)<-str_split_fixed(colnames(IC50),pattern = "_",2)[,1]
dif_IC50 <- IC50[,dif_drug]
dif_IC50<-dif_IC50[,order(colnames(dif_IC50),decreasing = F)]
identical(rownames(dif_IC50),colnames(hub_eset))
library(psych)
library(reshape2)
y <- hub_eset%>%t(.)
x <-dif_IC50
d <- corr.test(x,y,use="complete",method = 'spearman')
r <- data.frame(t(d$r))
r$Type<-rownames(r)
r<-melt(r,id='Type')
p <- data.frame(t(d$p))
p$Type<-rownames(p)
p<-melt(p,id='Type')
data<-cbind(r,p)[,c(-4,-5)]
colnames(data)<-c('Type','cell','correlation','p')
data$cell <- as.character(data$cell)  
data$cell[data$cell == "Nutlin.3a...."] <- "Nutlin-3a (-)"  
data$cell <- factor(data$cell)
write.csv(data,file="04.hubgene_drug_cor.csv",quote=F,row.names=F)
data1 <- data[abs(data$correlation) > 0 & data$p < 100,]
data1$cell <- as.character(data1$cell)  
data1$cell[data1$cell == "Nutlin.3a...."] <- "Nutlin-3a (-)"  
data1$cell <- factor(data1$cell)
write.csv(data1,file="05.hubgene_FDAdrug.csv",quote=F,row.names=F)
h <- c("Afatinib","Lapatinib","Erlotinib","Docetaxel","AZD4547","PD173074","Gefitinib")
colnames(IC50)<-str_split_fixed(colnames(IC50),pattern = "_",2)[,1]
colnames(IC50)[order(colnames(IC50))]
dif_IC50 <- IC50[,h]
dif_IC50$sample <- rownames(dif_IC50)
IC50_FE <- merge(dif_IC50,sub_group,by = "sample")
IC50_FE <- gather(data = IC50_FE,key = gene,value = Expression,-c("sample","group"))
IC50_FE$Expression <- as.numeric(IC50_FE$Expression)
IC50_FE$group <- factor(IC50_FE$group)
IC50_FE <- IC50_FE %>% drop_na(gene, Expression)
y_limits <- c(min(IC50_FE$Expression, na.rm = TRUE), max(IC50_FE$Expression, na.rm = TRUE))
data <- IC50_FE[IC50_FE$gene == "Afatinib",]
pdf("Afatinib_box.pdf", height = 6, width = 5)
ggplot(data, aes(x = factor(group), y = Expression)) +
  stat_boxplot(geom = "errorbar", width = 0.3, size = 1, color = "black") + 
  geom_boxplot(aes(fill = factor(group)), 
               color = "black", 
               width = 0.3,
               outlier.shape = 21,
               outlier.alpha = 0,
               size = 0.5) + 
  
  labs(title = NULL, x = "", y = "IC50") +
  theme_classic() +
  stat_compare_means(label = "p.format", aes(group = group), size = 8, label.x = 1.4) +
  theme(
    panel.grid = element_blank(), 
    axis.title.y = element_text(face = "bold", size = 20), 
    axis.title.x = element_text(face = "bold", size = 20), 
    axis.text.y = element_text(face = "bold", size = 20), 
    axis.text.x = element_text(face = "bold", size = 20), 
    legend.position = "none"
  )
dev.off()
png("Afatinib_box.png", height = 600, width = 500)
ggplot(data, aes(x = factor(group), y = Expression)) +
  stat_boxplot(geom = "errorbar", width = 0.3, size = 1, color = "black") + 
  geom_boxplot(aes(fill = factor(group)), 
               color = "black", 
               width = 0.3,
               outlier.shape = 21,
               outlier.alpha = 0,
               size = 0.5) + 
  
  labs(title = NULL, x = "", y = "IC50") +
  theme_classic() +
  stat_compare_means(label = "p.format", aes(group = group), size = 8, label.x = 1.4) +
  theme(
    panel.grid = element_blank(), 
    axis.title.y = element_text(face = "bold", size = 20), 
    axis.title.x = element_text(face = "bold", size = 20), 
    axis.text.y = element_text(face = "bold", size = 20), 
    axis.text.x = element_text(face = "bold", size = 20), 
    legend.position = "none"
  )
dev.off()
