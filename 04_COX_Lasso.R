rm(list = ls())
if (! dir.exists("./04.Cox_LASSO")){dir.create("./04.Cox_LASSO")}
library(GEOquery)
library(Biobase)
library(tidyverse)
library(ggsci)
group <- read.csv("../00_Rawdata/02.TCGA_Group_all.csv",header=T)
survival_data <- read_tsv("../00_Rawdata/TCGA-LAML.survival.tsv")
survival_data <- survival_data[survival_data$sample %in% group$X, c(1,2,4)]
write.csv(survival_data,"survival_data_fpkm.csv")
hub_gene <- read.csv("../03_common_gene/hub_gene.csv", header=T,row.names = 1)
dat1 <- read.csv("../00_Rawdata/03.TCGA_fpkms_all.csv", row.names = 1)
colnames(dat1) <- gsub('\\.','-',colnames(dat1))
dat <- dat1[, colnames(dat1) %in% survival_data$sample]
write.csv(dat,"dat_fpkm.csv")
max(dat1) 
train_dat <- dat1[rownames(dat1) %in% hub_gene$x, colnames(dat1) %in% survival_data$sample]
write.csv(train_dat,"train_dat_fpkm.csv")
train_dat <- t(train_dat) %>% as.data.frame()
train_dat$sample <- rownames(train_dat)
train_dat <- merge(survival_data, train_dat, by='sample')
train_dat <- column_to_rownames(train_dat, var = 'sample')
train_dat <- train_dat[rownames(train_dat) %in% group$X, ]
train_data <- train_dat 
library(survival)
library(survminer)
pfilter <- 0.05  
uniresult <- data.frame()   
ph_res <- data.frame()
for(i in colnames(train_data[, 3:ncol(train_data)])){   
  unicox <- coxph(Surv(time = OS.time, event = OS) ~ train_data[, i], data = train_data)  
  cox_zph <- cox.zph(unicox)
  cox_table <- data.frame(cox_zph$table)
  if (cox_table$p[1] < 0.05) {
    next
  } else{
    rownames(cox_table) <- c(i, 'GLOBAL')
    ph_res <- rbind(ph_res, cox_table)
    unisum <- summary(unicox)   
    pvalue <- round(unisum$coefficients[, 5], 3) 
    if(pvalue < pfilter){ 
      uniresult <- rbind(uniresult,
                         cbind(gene = i,
                               HR = unisum$coefficients[, 2], 
                               pvalue = unisum$coefficients[, 5], 
                               L95CI = unisum$conf.int[, 3], 
                               H95CI = unisum$conf.int[, 4]
                         ))
    }
  }
  
}   
ph_res$symbol <- rownames(ph_res)
ph_res <- ph_res[ph_res$symbol %in% uniresult$gene, ]
uniresult$gene
uniresult <- column_to_rownames(uniresult, var = "gene")
uniresult[,1:4] <- as.data.frame(lapply(uniresult[,1:4],as.numeric))
uniresult <- signif(uniresult, digits = 4)
res_results_0.05 <- as.data.frame(cbind(rownames(uniresult), uniresult$pvalue, paste0(uniresult$HR, " (", uniresult$L95CI, "-", uniresult$H95CI, ")")))
res_results_0.05[,2] <- as.data.frame(sapply(res_results_0.05[,2],as.numeric))
names(res_results_0.05) <- c("gene", "p.value", "HR (95% CI for HR)")
res_results_0.05 <- column_to_rownames(res_results_0.05, var = "gene")
write.csv(uniresult, file = "01.uniresult.csv")
res_results_0.05_2 <- separate(res_results_0.05, "HR (95% CI for HR)",
                               into = c("HR", "HR.95L", "HR.95H"),
                               sep = " ")
res_results_0.05_2 <- separate(res_results_0.05_2, "HR.95L",
                               into = c("HR.95L", "HR.95H"),
                               sep = "\\-")
res_results_0.05_2$HR.95L <- gsub("\\(", "", res_results_0.05_2$HR.95L)
res_results_0.05_2$HR.95H <- gsub("\\)", "", res_results_0.05_2$HR.95H)
res_results_0.05_2[, 1:ncol(res_results_0.05_2)] <- as.numeric(unlist(res_results_0.05_2[,1:ncol(res_results_0.05_2)]))
res_results_0.05_2 <- res_results_0.05_2[order(res_results_0.05_2$HR),]
hz <- paste(round(res_results_0.05_2$HR, 4),
            "(", round(res_results_0.05_2$HR.95L, 4),
            "-", round(res_results_0.05_2$HR.95H, 4),")", sep = "")
tabletext <- cbind(c(NA, "Gene", rownames(res_results_0.05_2)),
                   c(NA, "P value", ifelse(res_results_0.05_2$p.value<0.001,
                                           "< 0.001",
                                           round(res_results_0.05_2$p.value,3))),
                   c(NA, "Hazard Ratio(95% CI)", hz))
gene_list <- rownames(res_results_0.05)
cox_data <- as.formula(paste0('Surv(OS.time, OS)~', paste(gene_list, sep = '', collapse = '+')))
cox_more <- coxph(cox_data, data = train_data)
cox_zph <- cox.zph(cox_more)
ggcoxzph(cox_zph)
pdf("01.PH_plot.pdf",width = 40,height = 30,family = "Times")
ggcoxzph(cox_zph)
dev.off()
png("01.PH_plot.png",width = 2500,height = 1500,family = "Times")
ggcoxzph(cox_zph)
dev.off()
library(forestplot)
p <- forestplot(labeltext = tabletext,
                mean = c(NA, NA, res_results_0.05_2$HR),
                lower = c(NA, NA, res_results_0.05_2$HR.95L), 
                upper = c(NA, NA, res_results_0.05_2$HR.95H), 
                col = fpColors(box = "#FF0000", lines = "#4169E1", zero = "gray50"),
                graph.pos = 3,  
                is.summary = c(TRUE, TRUE, rep(FALSE, 70)),
                boxsize = 0.06, 
                lwd.ci = 3,   
                ci.vertices.height = 0.08, 
                ci.vertices = TRUE,
                zero = 1, 
                lwd.zero = 0.5,    
                colgap = unit(5, "mm"),    
                xticks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), 
                lwd.xaxis = 2,  
                lineheight = unit(1.2, "cm"), 
                graphwidth = unit(.5, "npc"), 
                cex = 0.9, 
                fn.ci_norm = fpDrawCircleCI, 
                hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
                txt_gp = fpTxtGp(label = gpar(cex = 1, fontfamily = "Times"),
                                 ticks = gpar(cex = 0.8, fontface = "bold", fontfamily = "Times"),
                                 xlab = gpar(cex = 1, fontface = "bold", fontfamily = "Times"),
                                 title = gpar(cex = 1.25, fontface = "bold", fontfamily = "Times")),
                xlab = "Hazard Ratio",
                grid = T 
)
pdf(file = "01.Univariate_cox_forest.pdf", family = "Times", height = 6, width = 8, onefile = F)
print(p)
dev.off()
png(filename = "01.Univariate_cox_forest.png", family = "Times", height = 6, width = 8, units = 'in', res = 600)
print(p)
dev.off()
library(survival)
set.seed(10)
gene_list <- uniresult %>% rownames()
diff_expr_clinical <-train_data
cox_data <- as.formula(paste0('Surv(OS.time, OS)~', paste(gene_list, sep = '', collapse = '+')))
cox_more_2 <- coxph(cox_data, data = diff_expr_clinical)
mul_cox_result <- summary(cox_more_2)$coefficients
mul_cox_gene <- rownames(mul_cox_result)
mul_cox_gene
write.table(mul_cox_result, file = "02.multicox_coefficients.csv", sep = "\t", quote = F)
library(survminer)
train_HRforest <- ggforest(model = cox_more_2,
                           data = diff_expr_clinical,
                           main = "Hazard ratio of Candidate Genes",
                           fontsize = 1) +
  theme(plot.title = element_text(face = "bold", size = 10))
train_HRforest
ggsave(filename = "02.multicox_forest.png", height = 7, width = 10, train_HRforest)
ggsave(filename = "02.multicox_forest.pdf", height = 7, width = 10, train_HRforest)
tdmultiCoxSum=summary(cox_more_2)
outResult=data.frame()
outResult=cbind(
  HR=tdmultiCoxSum$conf.int[,"exp(coef)"],
  pvalue=tdmultiCoxSum$coefficients[,"Pr(>|z|)"],
  L95CI=tdmultiCoxSum$conf.int[,"lower .95"],
  H95CI=tdmultiCoxSum$conf.int[,"upper .95"]
  ) %>% as.data.frame()
outResult[,1:4] <- as.data.frame(lapply(outResult[,1:4],as.numeric))
outResult <- signif(outResult, digits = 4)
res_results_0.05 <- as.data.frame(cbind(rownames(outResult), outResult$pvalue, paste0(outResult$HR, " (", outResult$L95CI, "-", outResult$H95CI, ")")))
res_results_0.05[,2] <- as.data.frame(sapply(res_results_0.05[,2],as.numeric))
names(res_results_0.05) <- c("gene", "p.value", "HR (95% CI for HR)")
res_results_0.05 <- column_to_rownames(res_results_0.05, var = "gene")
write.csv(uniresult, file = "01.uniresult.csv")
outResult2 <- separate(res_results_0.05, "HR (95% CI for HR)",
                               into = c("HR", "HR.95L", "HR.95H"),
                               sep = " ")
outResult2 <- separate(outResult2, "HR.95L",
                               into = c("HR.95L", "HR.95H"),
                               sep = "\\-")
outResult2$HR.95L <- gsub("\\(", "", outResult2$HR.95L)
outResult2$HR.95H <- gsub("\\)", "", outResult2$HR.95H)
outResult2[, 1:ncol(outResult2)] <- as.numeric(unlist(outResult2[,1:ncol(outResult2)]))
outResult2 <- outResult2[order(outResult2$HR),]
hz <- paste(round(outResult2$HR, 4),
            "(", round(outResult2$HR.95L, 4),
            "-", round(outResult2$HR.95H, 4),")", sep = "")
tabletext <- cbind(c(NA, "Gene", rownames(outResult2)),
                   c(NA, "P value", ifelse(outResult2$p.value<0.001,
                                           "< 0.001",
                                           round(outResult2$p.value,3))),
                   c(NA, "Hazard Ratio(95% CI)", hz))
library(forestplot)
p <- forestplot(labeltext = tabletext,
                mean = c(NA, NA, outResult2$HR),
                lower = c(NA, NA, outResult2$HR.95L), 
                upper = c(NA, NA, outResult2$HR.95H), 
                col = fpColors(box = "#FF0000", lines = "#4169E1", zero = "gray50"),
                graph.pos = 3,  
                is.summary = c(TRUE, TRUE, rep(FALSE, 70)),
                boxsize = 0.06, 
                lwd.ci = 3,   
                ci.vertices.height = 0.08, 
                ci.vertices = TRUE,
                zero = 1, 
                lwd.zero = 0.5,    
                colgap = unit(5, "mm"),    
                xticks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), 
                lwd.xaxis = 2,  
                lineheight = unit(1.2, "cm"), 
                graphwidth = unit(.5, "npc"), 
                cex = 0.9, 
                fn.ci_norm = fpDrawCircleCI, 
                hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
                txt_gp = fpTxtGp(label = gpar(cex = 1, fontfamily = "Times"),
                                 ticks = gpar(cex = 0.8, fontface = "bold", fontfamily = "Times"),
                                 xlab = gpar(cex = 1, fontface = "bold", fontfamily = "Times"),
                                 title = gpar(cex = 1.25, fontface = "bold", fontfamily = "Times")),
                xlab = "Hazard Ratio",
                grid = T 
)
pdf(file = "02.multicox_forest.pdf", family = "Times", height = 6, width = 8, onefile = F)
print(p)
dev.off()
png(filename = "02.multicox_forest.png", family = "Times", height = 6, width = 8, units = 'in', res = 600)
print(p)
dev.off()
step_Cox <- step(cox_more_2 ,direction = "backward") 
saveRDS(step_Cox,"step_Cox.RDS")
cox_summary <- summary(step_Cox)
cox_result <- summary(step_Cox)$coefficients
cox_result <- cox_result[,]
cox_gene <- rownames(cox_result)
cox_gene
write.csv(cox_result, file = "03.step_Cox_coefficients.csv", quote = F)
write.csv(cox_gene, file = "03.cox.csv", quote = F)
train_HRforest <- ggforest(model = step_Cox,
                           data = diff_expr_clinical,
                           main = "Hazard ratio of Candidate Genes",
                           fontsize = 1) +
  theme(plot.title = element_text(face = "bold", size = 10))
train_HRforest
ggsave(filename = "03.step_Cox_forest.png", height = 7, width = 10, train_HRforest)
ggsave(filename = "03.step_Cox_forest.pdf", height = 7, width = 10, train_HRforest)
tdmultiCoxSum=summary(step_Cox)
outResult=data.frame()
outResult=cbind(
  HR=tdmultiCoxSum$conf.int[,"exp(coef)"],
  pvalue=tdmultiCoxSum$coefficients[,"Pr(>|z|)"],
  L95CI=tdmultiCoxSum$conf.int[,"lower .95"],
  H95CI=tdmultiCoxSum$conf.int[,"upper .95"]
) %>% as.data.frame()
outResult[,1:4] <- as.data.frame(lapply(outResult[,1:4],as.numeric))
outResult <- signif(outResult, digits = 4)
res_results_0.05 <- as.data.frame(cbind(rownames(outResult), outResult$pvalue, paste0(outResult$HR, " (", outResult$L95CI, "-", outResult$H95CI, ")")))
res_results_0.05[,2] <- as.data.frame(sapply(res_results_0.05[,2],as.numeric))
names(res_results_0.05) <- c("gene", "p.value", "HR (95% CI for HR)")
res_results_0.05 <- column_to_rownames(res_results_0.05, var = "gene")
write.csv(uniresult, file = "01.uniresult.csv")
outResult2 <- separate(res_results_0.05, "HR (95% CI for HR)",
                       into = c("HR", "HR.95L", "HR.95H"),
                       sep = " ")
outResult2 <- separate(outResult2, "HR.95L",
                       into = c("HR.95L", "HR.95H"),
                       sep = "\\-")
outResult2$HR.95L <- gsub("\\(", "", outResult2$HR.95L)
outResult2$HR.95H <- gsub("\\)", "", outResult2$HR.95H)
outResult2[, 1:ncol(outResult2)] <- as.numeric(unlist(outResult2[,1:ncol(outResult2)]))
outResult2 <- outResult2[order(outResult2$HR),]
hz <- paste(round(outResult2$HR, 4),
            "(", round(outResult2$HR.95L, 4),
            "-", round(outResult2$HR.95H, 4),")", sep = "")
tabletext <- cbind(c(NA, "Gene", rownames(outResult2)),
                   c(NA, "P value", ifelse(outResult2$p.value<0.001,
                                           "< 0.001",
                                           round(outResult2$p.value,3))),
                   c(NA, "Hazard Ratio(95% CI)", hz))
library(forestplot)
p <- forestplot(labeltext = tabletext,
                mean = c(NA, NA, outResult2$HR),
                lower = c(NA, NA, outResult2$HR.95L), 
                upper = c(NA, NA, outResult2$HR.95H), 
                col = fpColors(box = "#FF0000", lines = "#4169E1", zero = "gray50"),
                graph.pos = 3,  
                is.summary = c(TRUE, TRUE, rep(FALSE, 70)),
                boxsize = 0.06, 
                lwd.ci = 3,   
                ci.vertices.height = 0.08, 
                ci.vertices = TRUE,
                zero = 1, 
                lwd.zero = 0.5,    
                colgap = unit(5, "mm"),    
                xticks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), 
                lwd.xaxis = 2,  
                lineheight = unit(1.2, "cm"), 
                graphwidth = unit(.5, "npc"), 
                cex = 0.9, 
                fn.ci_norm = fpDrawCircleCI, 
                hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
                txt_gp = fpTxtGp(label = gpar(cex = 1, fontfamily = "Times"),
                                 ticks = gpar(cex = 0.8, fontface = "bold", fontfamily = "Times"),
                                 xlab = gpar(cex = 1, fontface = "bold", fontfamily = "Times"),
                                 title = gpar(cex = 1.25, fontface = "bold", fontfamily = "Times")),
                xlab = "Hazard Ratio",
                grid = T 
)
p
pdf(file = "03.step_forest.pdf", family = "Times", height = 6, width = 8, onefile = F)
print(p)
dev.off()
png(filename = "03.step_forest.png", family = "Times", height = 6, width = 8, units = 'in', res = 600)
print(p)
dev.off()
