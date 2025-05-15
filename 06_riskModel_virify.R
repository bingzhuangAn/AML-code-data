rm(list = ls())
if (!dir.exists("06_riskModel_verify")) {dir.create("06_riskModel_verify")}
library(magrittr)
fpkm <- read.csv('../00_Rawdata/GSE71014_expr_use.csv',sep=",",row.names = 1)
fpkm <- na.omit(fpkm)
survival_dat <- read.csv('../00_Rawdata/OS_dat_GSE71014.csv',row.names = 1)
str(survival_dat)
colnames(survival_dat) <- c("id","OS","OS.time")
survival_dat <-na.omit(survival_dat)
survival_dat$OS.time<-round(survival_dat$OS.time)
model_coef <- read.csv("../04.Cox_LASSO/03.step_Cox_coefficients.csv", header = T,row.names = 1)
model_gene <- model_coef %>% rownames()
fpkm_1 <- fpkm[model_gene,] %>%t %>% as.data.frame() %>% rownames_to_column(var = "id")
fpkm_tem <- merge(fpkm_1,survival_dat,by = "id")
rownames(fpkm_tem) <- fpkm_tem$id
fpkm_tem <- fpkm_tem[,-1]
df_merge <- fpkm_tem
write.table(df_merge, file = "01.survival.csv", sep = "\t", quote = F)
library(survival)
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
coef <- as.numeric(model_coef$coef)
df.exp <- df_merge %>% subset(select = -c(OS, OS.time))
colnames(df.exp) <- gsub("-", "_", colnames(df.exp))
df.exp <- apply(df.exp, c(1,2), as.numeric) %>% as.data.frame
riskScore_out <- as.matrix(apply(df.exp, 1, function(x){x %*% coef}))
write.table(riskScore_out, file = "02.verify_riskScore.txt", sep = "\t", row.names = T, col.names = T, quote = F)
identical(rownames(riskScore_out), rownames(df_merge))
df_merge$riskScore_out <- riskScore_out[,1]
df_merge$risk <- ifelse(df_merge$riskScore_out > median(df_merge$riskScore), "High risk", "Low risk")
table(df_merge$risk)
kmfit_out <- survfit(Surv(OS.time, OS) ~ risk, data = df_merge)
verify_km <-  ggsurvplot(kmfit_out,
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
verify_km
verify_km$plot <- verify_km$plot + labs(
  title    = "Survival Curves",
  subtitle = "Based on Kaplan-Meier Estimates"
)
verify_km$table <- verify_km$table + labs(
  caption  = "Created with Verify Data"
)
verify_km <- customize_labels(
  verify_km,
  font.title    = c(16, "bold"),
  font.subtitle = c(15, "bold.italic"),
  font.caption  = c(14, "plain", "orange"),
  font.x        = c(14, "bold.italic"),
  font.y        = c(14, "bold.italic"),
  font.xtickslab = c(12, "plain")
)
pvalue <- stringr::str_extract(verify_km[["plot"]][["layers"]][[4]][["aes_params"]][["label"]],
                               "\\d.*")
verify_km[["plot"]][["layers"]][[4]][["aes_params"]][["label"]] <- as.expression(bquote(italic('p')==.(pvalue)))
verify_km
pdf(file = "02.verify_km.pdf", family = "Times", height = 6, width = 6, onefile = F)
print(verify_km)
dev.off()
png(filename = "02.verify_km.png", family = "Times", height = 6, width = 6, units = "in", res = 600)
print(verify_km)
dev.off()
rt = subset(df_merge, select = c(OS, OS.time, riskScore_out))
rt$OS.time <- rt$OS.time / 12
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
year_1
year_2
year_3
if(T){
  pdf(file = paste0("3.ModelROC_train.pdf"),width = 8,height = 8, family = "Times")
  a <- dev.cur()   
  png(file = paste0("3.ModelROC_train.png"),width= 8, height= 8, units="in", res=300, family = "Times")
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
risk_dis <- ggplot(df_merge, aes(x=reorder(rownames(df_merge), riskScore_out), y=riskScore_out, color = risk)) +
  geom_point() +
  scale_color_manual(values = c("#A73030FF", "#0073C2FF")) +
  scale_x_discrete(breaks = rownames(df_merge)[order(df_merge$riskScore_out)][c(1,50,100,150,200)],
                   labels = c(1,50,100,150,200),
                   expand = c(0.02,0)) +
  geom_vline(xintercept = nrow(df_merge[which(df_merge$risk=="Low risk"),]) + 0.5, lty = 2) +
  geom_hline(yintercept = median(df_merge$riskScore_out), lty =2) +
  labs(x = "Patients(increasing risk score)",
       y = "Risk Score",
       title = "Verify Dataset Risk Score Distribution") +
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
ggsave(filename = "04.riskScore_dis.pdf", height = 3, width = 5, risk_dis)
ggsave(filename = "04.riskScore_dis.png", height = 3, width = 5, risk_dis)
surv_stat <- ggplot(df_merge, aes(x=reorder(rownames(df_merge), riskScore_out),
                                  y=OS.time,
                                  color = factor(OS,
                                                 levels = c(0,1),
                                                 labels = c("Alive",
                                                            "Dead")))) +
  geom_point() +
  scale_color_manual(values = c("#0073C2FF","#A73030FF")) +
  scale_x_discrete(breaks = rownames(df_merge)[order(df_merge$riskScore_out)][c(1,30,60,90,120,150)],
                   labels = c(1,30,60,90,120,150),
                   expand = c(0.02,0)) +
  geom_vline(xintercept = nrow(df_merge[which(df_merge$risk=="High risk"),]) + 0.5, lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Overall Survival (days)",
       title = "Verify Dataset Overall Survival (days) Distribution") +
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
ggsave(filename = "05.OS_dis.pdf", height = 5, width = 7, surv_stat)
ggsave(filename = "05.OS_dis.png", height = 5, width = 7, surv_stat)
library(magrittr)
train_dat <- fpkm
riskScore <- as.data.frame(riskScore_out)
riskScore$risk <- ifelse(median(df_merge$riskScore_out), "High Risk", "Low Risk")
riskScore$risk <- factor(riskScore$risk)
colnames(riskScore) <- c("riskScore","risk")
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
mat <- train_dat[model_gene,]
mat <- t(scale(t(mat)))
riskScore$sample <- rownames(riskScore)
mat <- mat[,riskScore[order(riskScore$risk),]$sample]
newnames <- lapply(
  rownames(mat),
  function(x) bquote(italic(.(x))))
str(riskScore)
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
png(filename = "06.modelGene_pheatmap.png", height = 3, width = 8, units = 'in',res = 600,family='Times')
pheatplot
dev.off()
pdf(file = "06.modelGene_pheatmap.pdf", height = 3,width = 8,family='Times')
pheatplot
dev.off()
dat <- train_dat[model_gene,] %>% t %>% as.data.frame
cor <- stats::cor
library(corrplot)
library(ggcorrplot)
library(RColorBrewer)
corr <- round(cor(dat), 3)
p.mat <- round(ggcorrplot::cor_pmat(dat), 3)
write.table(corr, file = "7.modelgene_corr.txt", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(p.mat, file = "8.modelgene_pmat.txt", sep = "\t", quote = F, row.names = T, col.names = T)
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
pdf(file = "9.modelgene_corr.pdf", family = "Times", height = 8, width = 12)
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
png(filename = "9.modelgene_corr.png", family = "Times", height = 8, width = 12, units = "in", res = 300)
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
