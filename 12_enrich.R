rm(list = ls())
if (!dir.exists("12_enrich")) {dir.create("12_enrich")}
library(magrittr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
hub_gene <- read.csv("../03_common_gene/hub_gene.csv", sep = ",",header = T)
hub_gene_symbol <- hub_gene$x
length(hub_gene_symbol)
gene_transform <- bitr(hub_gene_symbol,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID, 
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID", 
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable = TRUE)
write.table(ego@result,file = "01.GO.xls",sep = "\t",quote = F,row.names = F)
dim(ego@result)
table(ego@result$ONTOLOGY)
str(ego@result$ONTOLOGY)
library(enrichplot)
go_dot <- dotplot(ego, showCategory=10, split="ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free", space = "free") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width= 60)) +
  theme(text = element_text(family = "Times"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        legend.title = element_text(face = "bold"),
        strip.text = element_text(size = 15, family = "Times"))
go_dot
ggsave(filename = "02.GO_dot.pdf", width = 12, height = 10, go_dot)
ggsave(filename = "02.GO_dot.png", width = 12, height = 10, go_dot)
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.9, 
                 qvalueCutoff = 1
)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "04.KEGG.xls",sep = "\t",quote = F,row.names = F)
dim(kk@result)
significant_results <- kk@result[kk@result$pvalue < 0.05, ]
slot(kk, "result") <- significant_results 
kk_dot <- dotplot(kk, showCategory=10) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 80)) +
  theme(text = element_text(family = "Times"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        legend.title = element_text(face = "bold"))
kk_dot
ggsave(filename = "03.KEGG_dot.pdf", width = 12, height = 5, kk_dot)
ggsave(filename = "03.KEGG_dot.png", width = 12, height = 5, kk_dot)
