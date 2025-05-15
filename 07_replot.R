rm(list = ls())
if (!dir.exists("07_replot")) {dir.create("07_replot")}
library(rms)
library(dplyr)
risk <- read.table('../05_riskModel/10.train_riskScore.txt', header = T,row.names = 1, check.names = F)
risk$id <- rownames(risk)
train_dat <- read.csv('../04.Cox_LASSO/dat_fpkm.csv', row.names = 1, check.names = F)
max(train_dat)
train_dat <- train_dat[, colnames(train_dat) %in%risk$id]
phenoType <- read.csv("../05_riskModel/01.survival.csv",sep = "\t")
phenoType$sample <- rownames(phenoType)
phenoType <- phenoType[,c("OS","OS.time","sample")]
survival <- read.table("../00_Rawdata/TCGA-LAML.GDC_phenotype.tsv", row.names = 1,sep = "\t",header = T)
survival <- survival[,c(2,40)]
survival$sample <- rownames(survival)
survival <- survival[survival$sample %in% phenoType$sample,]
colnames(survival) <- c("Age","Type","sample")
train_phenoType <- merge(survival, phenoType, by = 'sample')
write.csv(train_phenoType, file = 'phenoType.csv')
train_phenoType$OS <- as.numeric(train_phenoType$OS)
train_phenoType$OS.time <- as.numeric(train_phenoType$OS.time)
train_phenoType2 <- train_phenoType
table(train_phenoType2$Type)
train_phenoType2$Type <- gsub('M0 Undifferentiated','not reported',train_phenoType2$Type)
train_phenoType2$Type <- gsub(c("M2|M1|M4|M5|M6|M7"),'M_other',train_phenoType2$Type)
train_phenoType2 <- subset(train_phenoType2, Type != "not reported")
table(train_phenoType2$Type)
train_phenoType2$Age <- ifelse(train_phenoType2$Age>60 , 1, 0)
table(train_phenoType2$Age)
colnames(train_phenoType2)
colnames(train_phenoType2) <- c('id', "Age","Type","OS","OS.time")
risk <- read.table('../05_riskModel/10.train_riskScore.txt', header = T,row.names = 1, check.names = F) %>% lc.tableToNum()
risk$id <- rownames(risk)
sub_risk <- subset(risk, select = c(id, riskScore))
train_risk_clinical <- merge(train_phenoType2, sub_risk, by = "id")
rownames(train_risk_clinical) <- train_risk_clinical$id
train_risk_clinical <- subset(train_risk_clinical, select = -c(id))
dim(train_risk_clinical)
colnames_train <- colnames(train_risk_clinical)
covariates_train <- colnames_train[-which(colnames_train %in% c("OS", "OS.time"))]
train_risk_clinical$Type <- factor(train_risk_clinical$Type)
library(survival)
res.risk <- coxph(Surv(time = OS.time, event = OS) ~ riskScore, data = train_risk_clinical) %>% summary
res.risk_ph <- cox.zph(coxph(Surv(time = OS.time, event = OS) ~ riskScore, data = train_risk_clinical))
res.risk_ph <- data.frame(res.risk_ph$table)
if(res.risk_ph$p[1] < 0.05){
  print('ph假设检验未通过')
}else{
  print('ph假设检验通过')
}
res.risk <- c(res.risk$conf.int[-2], res.risk$coefficients[5])
res.age <- coxph(Surv(time = OS.time, event = OS) ~ Age, data = train_risk_clinical) %>% summary
res.age_ph <- cox.zph(coxph(Surv(time = OS.time, event = OS) ~ Age, data = train_risk_clinical))
res.age_ph <- data.frame(res.age_ph$table)
if(res.age_ph$p[1] < 0.05){
  print('ph假设检验未通过')
}else{
  print('ph假设检验通过')
}
res.age <- c(res.age$conf.int[-2], res.age$coefficients[5])
res.Type <- coxph(Surv(time = OS.time, event = OS) ~ Type, data = train_risk_clinical) %>% summary
res.Type_ph <- cox.zph(coxph(Surv(time = OS.time, event = OS) ~ Type, data = train_risk_clinical))
res.Type_ph <- data.frame(res.Type_ph$table)
if(res.Type_ph$p[1] < 0.05){
  print('ph假设检验未通过')
}else{
  print('ph假设检验通过')
}
res.Type <- c(res.Type$conf.int[,-2], res.Type$coefficients[,5])
res.ref <- c(1,1,1,NA)
res.ref1 <- c(1,1,NA)
res.ref2 <- c(1,NA)
res <- rbind(res.risk, res.age, res.Type) %>% as.data.frame()
rownames(res)
res$Indicators <- c("riskScore", "Age", "Type" )
colnames(res) <- c("hr","low","up","pv","Indicator")
res$p <- signif(res$pv, 2) %>% paste0("p = ", .)
res$p[is.na(res$pv)] <- NA
res$Indicator <- factor(res$Indicator, levels = rev(res$Indicator))
rownames(res) <- res$Indicator
res2 <- data.frame(p.value=res$pv,
                   HR=res$hr,
                   HR.95L=res$low,
                   HR.95H=res$up,
                   Indicator=res$Indicator)
rownames(res2) <- res2$Indicator
write.table(res2, file = "univariate_cox_prog_forest.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))
library(tidyr)
hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.0001, "< 0.0001", round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
library(forestplot)
{p <- forestplot(labeltext = tabletext, 
                 graph.pos=4,  
                 is.summary = c(FALSE,FALSE,rep(FALSE, 7)),
                 col=fpColors(box="red", lines="royalblue", zero = "gray50"),
                 mean=c(NA,res2$HR),
                 lower=c(NA,res2$HR.95L), 
                 upper=c(NA,res2$HR.95H), 
                 boxsize=0.1,lwd.ci=3,   
                 ci.vertices.height = 0.08,ci.vertices=TRUE, 
                 zero=1,lwd.zero=0.5,    
                 colgap=unit(5,"mm"),    
                 xticks = c(0,1,2,3,4,5,6,7,8,9), 
                 lwd.xaxis=2,            
                 lineheight = unit(1.2,"cm"), 
                 graphwidth = unit(.6,"npc"), 
                 cex=0.9, fn.ci_norm = fpDrawCircleCI, 
                 hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
                 
                 
                 
                 
                 
                 txt_gp=fpTxtGp(label=gpar(cex=1),
                                ticks=gpar(cex=0.8, fontface = "bold",family='Times'),
                                xlab=gpar(cex = 1, fontface = "bold",family='Times'),
                                title=gpar(cex = 1.25, fontface = "bold",family='Times')),
                 xlab="Hazard Ratio",
                 grid = T,
                 title = "Univariate") 
}
p
pdf(file = '01.uncoxforest.pdf',height = 4, width =10, onefile = F,family='Times')
print(p)
dev.off()
png(file = '01.uncoxforest.png',height = 300, width =700,family='Times')
print(p)
dev.off()
cox_more <- coxph(Surv(time = OS.time, event = OS) ~ riskScore + Type  , data = train_risk_clinical)
cox_zph <- cox.zph(cox_more)
cox_table <- cox_zph$table[-nrow(cox_zph$table),]
cox_table
res.mul <- coxph(Surv(time = OS.time, event = OS) ~ riskScore+ Type + Age, data = train_risk_clinical) %>% summary
res.mul  <-  cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5]) %>% as.data.frame()
rownames(res.mul)
res.mul$Indicators = c('riskScore',
                       "Type","Age")
colnames(res.mul) = c("hr","low","up","pv","Indicator")
multi_res <- data.frame(p.value=res.mul$pv,
                        HR=res.mul$hr,
                        HR.95L=res.mul$low,
                        HR.95H=res.mul$up,
                        Indicator=res.mul$Indicator)
rownames(multi_res) <- multi_res$Indicator
multi_res
multi_res <- subset(multi_res, select = -c(Indicator))
write.csv(multi_res,
          file = "03.multivariate_cox_prog_result.csv",
          quote = F,
          row.names = T)
library(tidyr)
hz <- paste(round(multi_res$HR,3),
            "(",
            round(multi_res$HR.95L,3),
            "-"
            ,round(multi_res$HR.95H,3),
            ")",
            sep = "")
hz
tabletext <- cbind(c(NA,rownames(multi_res)),
                   c("P value",ifelse(multi_res$p.value<0.05,
                                      "< 0.05",
                                      round(multi_res$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
pdf(file = '03.mulcoxforest.pdf',height = 4, width = 8, onefile = F,family='Times')
forestplot(labeltext=tabletext, 
           graph.pos=4,  
           is.summary = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), 
           upper=c(NA,multi_res$HR.95H), 
           boxsize=0.2,lwd.ci=3,   
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=0.5,    
           colgap=unit(5,"mm"),    
           xticks = c(0,1,2,3,4,5,6,7), 
           lwd.xaxis=2,            
           lineheight = unit(1.2,"cm"), 
           graphwidth = unit(.6,"npc"), 
           cex=1.2, fn.ci_norm = fpDrawCircleCI, 
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           
           
           
           
           
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate Cox") 
dev.off()
png(file = '03.mulcoxforest.png',height = 4, width = 8, units = 'in',res = 600,family='Times')
forestplot(labeltext=tabletext, 
           graph.pos=4,  
           is.summary = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), 
           upper=c(NA,multi_res$HR.95H), 
           boxsize=0.2,lwd.ci=3,   
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=0.5,    
           colgap=unit(5,"mm"),    
           xticks = c(0,1,2,3,4,5,6,7), 
           lwd.xaxis=2,            
           lineheight = unit(1.2,"cm"), 
           graphwidth = unit(.6,"npc"), 
           cex=1.2, fn.ci_norm = fpDrawCircleCI, 
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           
           
           
           
           
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate Cox") 
dev.off()
multi_cov <- c("riskScore",  "Age")
cox_data_prog <- as.formula(paste0('Surv(OS.time, OS)~', paste(multi_cov, sep = '', collapse = '+')))
cox_more_prog <- coxph(formula = cox_data_prog,
                       data = as.data.frame(train_risk_clinical))
library(rms)
ddist <- datadist(train_risk_clinical)
options(datadist= 'ddist')
train_risk_clinical$Age <- ifelse(train_risk_clinical$Age == 1, '>60', '<=60')
res.cox <- psm(cox_data_prog,
               data = train_risk_clinical, dist = 'lognormal')
surv <- Survival(res.cox) 
function(x) surv(365, x) 
function(x) surv(730, x) 
function(x) surv(1095, x) 
function(x) surv(1825, x) 
nom.cox <- nomogram(res.cox,
                    fun = list(function(x) surv(365, x),
                               function(x) surv(1095, x),
                               function(x) surv(1825, x)),
                    funlabel=c("1-year Survival Probability", "3-year Survival Probability", "5-year Survival Probability"),
                    maxscale = 10,
                    fun.at = c(0.01,seq(0.1,0.9,by=0.2),0.95,0.99),
                    lp=F)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
png(filename = "03.nomogram_line_points.png", height = 800, width = 1700)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
pdf(file = "03.nomogram_line_points.pdf", height = 8, width = 20)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
train_risk_clinical <- merge(train_phenoType2, sub_risk, by = "id")
res.mul = coxph(Surv(time = OS.time, event = OS) ~ riskScore + Type + Age, data = train_risk_clinical)
train_risk_clinical$riskScore_prog = predict(res.mul, newdata = train_risk_clinical, Type = "lp")
train_risk_clinical$risk_prog = ifelse(train_risk_clinical$riskScore_prog > median(train_risk_clinical$riskScore_prog, na.rm = T), "high", "low")
train_risk_clinical$risk_prog = factor(train_risk_clinical$risk_prog, levels = c("high", "low"), labels = c("High risk", "Low risk"))
train_risk_clinical2 <- train_risk_clinical
library(survival)
library(survminer)
multi_ROC_out <- function(time_vector, risk_score_table){
  library(survivalROC)
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore_prog,
                           predict.time=single_time,method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
for_multi_ROC <- multi_ROC_out(time_vector = c(365*c(1,3,5)),
                               risk_score_table = train_risk_clinical2)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)
library(scales)
library(geomROC)
library(plotROC)
auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y2 <- round(for_multi_ROC[which(for_multi_ROC$Time==730),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)
ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive,
                                 label=Cut_values,
                                 color=Time)) +
  scale_color_manual(breaks = c("365","1095","1825"),
                     labels = c("1 year", "3 year","5 year"),
                     values = c("#4682B4", "#FF4040", "#20B2AA")) +
  geom_roc(labels = F, stat = 'identity') +
  style_roc() +
  geom_abline(slope = 1, intercept = 0, color = 'gray', lineType=2) +
  theme_bw() +
  labs(title = "Nomogram ROC") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
  
  
  
  annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 year =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 year =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 year =', format(auc_y5,nsmall=2))))
ROC
ggsave('05.Prog_ROC.png', ROC,width = 5, height = 4)
ggsave('05.Prog_ROC.pdf', ROC,width = 5, height = 4)
