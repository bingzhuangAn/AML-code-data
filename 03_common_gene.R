rm(list = ls())
library(ggvenn)
DEGs <- read.csv("../01_DEG_GSE114868/02.GSE114868_DEG_sig_res_CA.csv")
senes <- read.csv("../00_Rawdata/03.CSRGs.csv",header = T)
Tcell <- read.csv("../01_scRNA/05_SingleR/mast/T cell wilcox.csv")
Tcell_p <- Tcell[Tcell$p_val <- 0.05 & abs(Tcell$avg_log2FC) > 0.1 ,]$X
common <- intersect(rownames(DEGs),senes$Gene)
common2 <- intersect(common,Tcell_p)
write.csv(common2,"hub_gene.csv")
p2<-ggvenn(list( "DEGs" = rownames(DEGs), "CSRGs" = senes$Gene, "T cell" = Tcell_p),
           show_percentage = T,
           stroke_color = "white",
           stroke_size = 0.5,
           fill_color = c("#FF8C00","#1E90FF","purple"),
           set_name_color =c("#FF8C00","#1E90FF","purple"),
           set_name_size = 4,text_size=3)
ggsave("./venn.pdf",width=6,height=6,p2)
ggsave("./venn.png",width=6,height=6,p2)
