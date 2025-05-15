rm(list = ls())
if(!dir.exists("05_riskModel")){dir.create("05_riskModel")}
library(lance)
library(readr)
library(tidyverse)
survival<-read.csv(file = '../04.Cox_LASSO/survival_data_fpkm.csv',header = T,check.names = F,row.names = 1)
dat.tcga <- read.csv("../04.Cox_LASSO/dat_fpkm.csv",header = T,check.names = F,row.names = 1)
dat.tcga$sample <- rownames(dat.tcga)
inter <- read.csv("../04.Cox_LASSO/03.cox.csv",sep = ",", header = T,row.names = 1)
group <- read.csv("../00_Rawdata/02.TCGA_Group_all.csv",header = T,check.names = F,row.names = 1) 
tumor.sample<-rownames(group)[which(group$Group=='AML')]
DEDRG<-inter
rownames(survival) <- survival$sample
survival<-survival[colnames(dat.tcga),]
survival <- na.omit(survival)
train_dat<-dat.tcga[,colnames(dat.tcga) %in% tumor.sample]
train_dat<-t(train_dat) 
train_dat<-as.data.frame(train_dat)
train_dat<-train_dat[,DEDRG$x]
train_dat$sample<-rownames(train_dat)
train_dat<-merge(survival,train_dat,by='sample')
train_dat<-column_to_rownames(train_dat,var = 'sample')
train_data<-train_dat
rt <- train_data
write.table(train_data,file = '01.survival.csv',sep = '\t',quote = F)
diff_expr_clinical <-train_data
step_Cox <- readRDS("../04.Cox_LASSO/step_Cox.RDS")
riskScore <- predict(step_Cox,type="lp",newdata=diff_expr_clinical)
write.table(data.frame(riskScore) %>% tibble::rownames_to_column(var = "sample"),
            file = "10.train_riskScore.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
riskScore1 <- data.frame(riskScore)
diff_expr_clinical$riskScore <- riskScore1$riskScore
OS_dat <- survival
colnames(OS_dat)[1] <- "sample"
diff_expr_clinical$sample <- rownames(diff_expr_clinical)
score_dat <- diff_expr_clinical
OS_dat <- merge(OS_dat,score_dat,by = "sample",.keep_all = T) %>% na.omit()
OS.cut <- survminer::surv_cutpoint(OS_dat,
                                   time = "OS.time.y",
                                   event = "OS.y",
                                   variables = "riskScore")
res <- data.frame(summary(OS.cut))
res
diff_expr_clinical$risk <- ifelse(diff_expr_clinical$riskScore > median(diff_expr_clinical$riskScore), "High", "Low")
write.csv(diff_expr_clinical, "risk.csv")
kmfit <- survfit(Surv(OS.time, OS) ~ risk, data = diff_expr_clinical)
library(survminer)
customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}
train_km <- ggsurvplot(kmfit,
                       pval = TRUE, 
                       pval.method = T,
                       conf.int = F,
                       legend.labs=c("High risk","Low risk" ),
                       legend.title="Risk Score",
                       xlab = "Overall Survival (days)",
                       title="Train KM",
                       font.main = c(15,"bold"),
                       risk.table = TRUE, 
                       risk.table.col = "strata", 
                       linetype = "strata", 
                       
                       font.family = "Times",
                       risk.table.y.text.col = T,
                       risk.table.y.text = T,
                       risk.table.height = 0.35,
                       ggtheme = theme_bw(), 
                       palette = c("#A73030FF", "#0073C2FF"))
train_km
train_km$plot <- train_km$plot + labs(
  title    = "Survival Curves",
  subtitle = "Based on Kaplan-Meier Estimates"
)
train_km$table <- train_km$table + labs(
  caption  = "Created with Train Data"
)
train_km <- customize_labels(
  train_km,
  font.title    = c(16, "bold"),
  font.subtitle = c(15, "bold.italic"),
  font.caption  = c(14, "plain", "orange"),
  font.x        = c(14, "bold.italic"),
  font.y        = c(14, "bold.italic"),
  font.xtickslab = c(12, "plain")
)
pvalue <- stringr::str_extract(train_km[["plot"]][["layers"]][[4]][["aes_params"]][["label"]],
                               "\\d.*")
train_km[["plot"]][["layers"]][[4]][["aes_params"]][["label"]] <- as.expression(bquote(italic('p')==.(pvalue)))
train_km
pdf(file = "1.train_km.pdf", family = "Times", height = 6, width = 6, onefile = F)
print(train_km)
dev.off()
png(filename = "1.train_km.png", family = "Times", height = 6, width = 6, units = "in", res = 600)
print(train_km)
dev.off()
rt = subset(diff_expr_clinical, select = c(OS, OS.time, riskScore))
rt$OS.time <- rt$OS.time / 365
colnames(rt) <- c("fustat","futime","riskscore")
ROC <- rt
cutoff_1 <- 1
cutoff_2 <- 3
cutoff_3 <- 5
year_1= survivalROC(Stime=ROC$futime,
                    status=ROC$fustat,
                    marker = ROC$riskscore, 
                    predict.time = cutoff_1,
                    method = 'KM')
year_2= survivalROC(Stime=ROC$futime,
                    status=ROC$fustat,
                    marker = ROC$riskscore, 
                    predict.time = cutoff_2,
                    method = 'KM')
year_3= survivalROC(Stime=ROC$futime,
                    status=ROC$fustat,
                    marker = ROC$riskscore, 
                    predict.time = cutoff_3,
                    method = 'KM')
if(T){
  pdf(file = paste0("2.ModelROC_train.pdf"),width = 8,height = 8, family = "Times")
  a <- dev.cur()   
  png(file = paste0("2.ModelROC_train.png"),width= 8, height= 8, units="in", res=300, family = "Times")
  dev.control("enable")
  par(mar = c(5,5,5,2))
  plot(year_1$FP, year_1$TP,
       type="l",col="red",xlim=c(0,1), ylim=c(0,1),
       xlab="False Positive Fraction",
       ylab="True Positive Fraction",
       main="OS-train, Method = KM\n Year = 1,3,5",
       cex.lab = 1.5,
       cex.axis = 1.5,
       cex.main = 1.5
  )
  abline(0,1,col="gray",lty=2)
  lines(year_2$FP, year_2$TP, type="l",col="#EB4B17",xlim=c(0,1), ylim=c(0,1))
  lines(year_1$FP, year_1$TP, type="l",col="#2775AB",xlim=c(0,1), ylim=c(0,1))
  lines(year_3$FP, year_3$TP, type="l",col="#4C8045",xlim=c(0,1), ylim=c(0,1))
  legend(0.6,0.2,c(paste("AUC of 1 year =",round(year_1$AUC,2)),
                   paste("AUC of 3 year =",round(year_2$AUC,2)),
                   paste("AUC of 5 year =",round(year_3$AUC,2))),
         x.intersp=1, y.intersp=0.8,
         lty= 1 ,lwd= 2,col=c("#2775AB","#EB4B17",'#4C8045'),
         bty = "n",
         seg.len=1,cex=1.2)
  dev.copy(which = a)  
  dev.off()
  dev.off()
}
library(ggplot2)
library(ggthemes)
median(diff_expr_clinical$riskScore)
risk_dis <- ggplot(diff_expr_clinical, aes(x=reorder(rownames(diff_expr_clinical), riskScore), y=riskScore, color = risk)) +
  geom_point() +
  scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
  scale_x_discrete(breaks = rownames(diff_expr_clinical)[order(diff_expr_clinical$riskScore)][c(1,50,100,150,200,250,300,350,400)],
                   labels = c(1,50,100,150,200,250,300,350,400),
                   expand = c(0.02,0)) +
  geom_vline(xintercept = nrow(diff_expr_clinical[which(diff_expr_clinical$risk=="Low"),]) + 0.5, lty = 2) +
  geom_hline(yintercept = median(diff_expr_clinical$riskScore), lty =2) +
  labs(x = "Patients(increasing risk score)",
       y = "Risk Score",
       title = "Train Dataset Risk Score Distribution") + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        panel.grid = element_blank(),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15),
        text = element_text(family = "Times", face = "bold"))
risk_dis
ggsave(filename = "3.riskScore_dis.pdf", height = 3, width = 5, risk_dis)
ggsave(filename = "3.riskScore_dis.png", height = 3, width = 5, risk_dis)
surv_stat <- ggplot(diff_expr_clinical, aes(x=reorder(rownames(diff_expr_clinical), riskScore),
                                            y=OS.time,
                                            color = factor(OS,
                                                           levels = c(0,1),
                                                           labels = c("Alive", 
                                                                      "Dead")))) +
  geom_point() +
  scale_color_manual(values = c("#0073C2FF","#A73030FF")) +
  scale_x_discrete(breaks = rownames(diff_expr_clinical)[order(diff_expr_clinical$riskScore)][c(1,50,100,150,200,250,300,350,400)],
                   labels = c(1,50,100,150,200,250,300,350,400),
                   expand = c(0.02,0)) +
  geom_vline(xintercept = nrow(diff_expr_clinical[which(diff_expr_clinical$risk=="High"),]) + 0.5, lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Overall Survival (days)",
       title = "Train Dataset Overall Survival (days) Distribution") + 
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        panel.grid = element_blank(),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15),
        text = element_text(family = "Times", face = "bold"))
surv_stat
ggsave(filename = "4.OS_dis.pdf", height = 5, width = 7, surv_stat)
ggsave(filename = "4.OS_dis.png", height = 5, width = 7, surv_stat)
library(magrittr)
model_gene <-inter$x
group <- rownames(group)[which(group$Group == c("AML"))]  
train_dat<- dat.tcga[, colnames(dat.tcga) %in% group]
riskScore <- riskScore %>% as.data.frame()
colnames(riskScore) <- "riskScore"
riskScore$risk <- ifelse(riskScore$riskScore> median(diff_expr_clinical$riskScore), "High Risk", "Low Risk")
riskScore$risk <- factor(riskScore$risk)
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
mat <- train_dat[model_gene,]
mat <- t(scale(t(mat)))
riskScore$sample <- rownames(riskScore)
mat <- mat[,riskScore[order(riskScore$risk),]$sample]
newnames <- lapply(
  rownames(mat),
  function(x) bquote(italic(.(x))))
annotation_col <- data.frame(
  Group = riskScore$risk,
  row.names = riskScore$sample
)
ann_colors <- list(
  Group = c(`High Risk`="#A73030FF", `Low Risk`="#0073C2FF")
)
pheatplot <- pheatmap(mat=mat,
                      
                      annotation_col = annotation_col,
                      annotation_colors = ann_colors,
                      color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                      fontsize = 16,
                      fontsize_row = 16,
                      
                      fontfamily = "Times",
                      labels_row = as.expression(newnames),
                      show_colnames = FALSE,
                      cluster_cols = F,
                      cluster_rows = T,
                      annotation_names_row = T
)
pheatplot
png(filename = "5.modelGene_pheatmap.png", height = 3, width = 8, units = 'in',res = 600,family='Times')
pheatplot
dev.off()
pdf(file = "5.modelGene_pheatmap.pdf", height = 3,width = 8,family='Times')
pheatplot
dev.off()
dat <- train_dat[model_gene,] %>% t %>% as.data.frame
dat[,1:4] <- as.data.frame(lapply(dat[,1:3],as.numeric))
cor <- stats::cor
library(corrplot)
library(ggcorrplot)
library(RColorBrewer)
corr <- round(cor(dat), 3)
p.mat <- round(ggcorrplot::cor_pmat(dat), 3)
write.table(corr, file = "6.modelgene_corr.txt", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(p.mat, file = "7.modelgene_pmat.txt", sep = "\t", quote = F, row.names = T, col.names = T)
corr_plot <- ggcorrplot(corr, 
                        method = "circle",
                        hc.order = TRUE,
                        hc.method = "ward.D",
                        outline.color = "white",
                        ggtheme = theme_bw(),
                        type = "upper",
                        colors = c("deepskyblue4", "white", "darkorange4"),
                        lab = TRUE,
                        lab_size = 2,
                        p.mat = p.mat,
                        insig = "blank")+
  theme(plot.margin=margin(t = 0.5, r = 0.1, b = 0.5, l = 0.1, unit = "cm"))
corr_plot
corr_plot2 <- ggcorrplot(corr, 
                         method = "circle",
                         hc.order = TRUE,
                         hc.method = "ward.D",
                         outline.color = "white",
                         ggtheme = theme_bw(),
                         type = "upper",
                         colors = c("deepskyblue4", "white", "darkorange4"),
                         lab = TRUE,
                         lab_size = 2) +
  theme(plot.margin=margin(t = 0.5, r = 0.1, b = 0.5, l = 0.1, unit = "cm"))
corr_plot2
col2 = colorRampPalette(rev(c('#67001F', '#B2182B', '#D6604D', '#F4A582',
                              '#FDDBC7', '#FFFFFF', '#D1E5F0', '#92C5DE',
                              '#4393C3', '#2166AC', '#053061')))
corrplot(corr, type = "upper", addCoef.col="grey",
         col = col2(200),
         tl.cex = 1.5,
         tl.col = "black",
         tl.srt = 45
)
pdf(file = "8.modelgene_corr.pdf", family = "Times", height = 8, width = 12)
corrplot(corr,type="upper",
         tl.pos="tp",tl.col="black",tl.srt=45,tl.cex=1.5,
         p.mat = p.mat, insig = "label_sig",
         sig.level = c(.01, .05),pch.cex=1,pch.col = "black",
         order = "AOE",
         col = col2(200),
         mar = c(0,3,1,0)
)
corrplot(corr = corr,type="lower",add=TRUE,method="number",
         tl.pos="n",
         col="black",diag=FALSE, cl.pos="n",pch.col = "black",
         number.cex = 1.2,order = "AOE",
         mar = c(0,3,1,0)
)
dev.off()
png(filename = "8.modelgene_corr.png", family = "Times", height = 8, width = 12, units = "in", res = 300)
corrplot(corr,type="upper",
         tl.pos="tp",tl.col="black",tl.srt=45,tl.cex=1.5,
         p.mat = p.mat, insig = "label_sig",
         sig.level = c(.01, .05),pch.cex=1,pch.col = "black",
         order = "AOE",
         col = col2(200),
         mar = c(0,3,1,0)
)
corrplot(corr = corr,type="lower",add=TRUE,method="number",
         tl.pos="n",
         col="black",diag=FALSE, cl.pos="n",pch.col = "black",
         number.cex = 1.2,order = "AOE",
         mar = c(0,3,1,0)
)
dev.off()
