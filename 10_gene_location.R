rm(list = ls())
if (!dir.exists("./10_gene_location")) {
  dir.create("./10_gene_location")
}
library(AnnoProbe)
library(circlize)
library(RCircos)
library(maftools)
library(tidyverse)
model_coef <- read.csv("../04.Cox_LASSO/03.cox.csv", row.names = 1)
hub_gene <- model_coef$x
deg_mrna <- hub_gene
ag <- annoGene(deg_mrna, ID_type = 'SYMBOL', species = 'human')
data <- ag[, c(4, 5, 6, 1)]
colnames(data) <- c("Chromosome", "chromStart", "chromEnd", "Gene")
data(RCircos.Scatter.Data);
data(UCSC.HG19.Human.CytoBandIdeogram)
chr.exclude <- NULL 
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram; 
tracks.inside <- 10 
tracks.outside <- 0 
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside);
RCircos.List.Plot.Parameters()
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot();
data(RCircos.Gene.Label.Data);
side <- "in";
track.num <- 1;
RCircos.Gene.Connector.Plot(data, track.num,side);
name.col <- 4;
track.num <- 2;
RCircos.Gene.Name.Plot(data, name.col,track.num, side);
dev.off()
params <- RCircos.Get.Plot.Parameters()
params$text.size <- 1.5  
RCircos.Reset.Plot.Parameters(params)
params$label.distance <- 1.2  
RCircos.Reset.Plot.Parameters(params)
pdf("circos_hub_gene3_adjusted.pdf", width = 8, height = 8)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
connector_track_num <- 2
RCircos.Gene.Connector.Plot(data, connector_track_num, side)
name_track_num <- 3
RCircos.Gene.Name.Plot(data, name.col, name_track_num, side)
dev.off()
png("circos_hub_gene3_adjusted.png", width = 800, height = 800)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
connector_track_num <- 2
RCircos.Gene.Connector.Plot(data, connector_track_num, side)
name_track_num <- 3
RCircos.Gene.Name.Plot(data, name.col, name_track_num, side)
dev.off()
