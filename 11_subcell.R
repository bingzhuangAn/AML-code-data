rm(list = ls())
if (!dir.exists("./11_suncell")) {
  dir.create("./11_suncell")
}
subcell <- read.csv("suncell.csv",header = T)
library(ggplot2)
custom_colors <- c("#FF9999", "#66B3FF", "#99FF99", "#FFCC99", "#FF6699",  
                   "#CCCCFF", "#FFB3FF", "#FFFF99", "#99FFFF", "#FF99CC",  
                   "#CC99FF", "#C0C0C0", "#FF6666", "#66FFFF", "#99FFCC",  
                   "#FFCC66", "#66CCFF", "#CCCCFF", "#FFFFCC")  
p <- ggplot(subcell, aes(x = gene, y = Confidence, fill = Compartment)) +  
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
        strip.text = element_text(size = 14,color='black', family = "Times", face = "bold")) +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  labs(x="Gene", y="Confidence", fill="Compartment")  
ggsave('subcell.pdf', width = 12, height = 8, plot = p)  
ggsave('subcell.png', w=12, h=8, plot = p)
png('subcell.png',w = 1200, h=800,family ="Times")
print(p)
dev.off()
