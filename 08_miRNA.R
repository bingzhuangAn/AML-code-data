rm(list = ls())
if (!dir.exists("08_miRNA")) {dir.create("08_miRNA")}
p <- 5
library(multiMiR)
library(dplyr)
hub_gene <- read.csv("../04.Cox_LASSO/03.cox.csv")
hub_gene <- as.vector(hub_gene$x)
combined_list <- list()
combined_list1 <- list()
for (an in hub_gene) {
  example <- tryCatch({
    print(an)
    get_multimir(org = "hsa",                  
                 target = an,                   
                 table = "predicted",          
                 summary = TRUE,
                 predicted.cutoff = p,        
                 predicted.cutoff.type = "p",  
                 predicted.site = "all") 
  }, error = function(e) {
    cat("Error in get_multimir for target:", an, "\n", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(example)) {
    next
  }
  
  
  example_result <- example@data
  data_list <- list(
    targetscan = example_result[example_result$database == "targetscan",]$mature_mirna_id,
    diana_microt = example_result[example_result$database == "diana_microt",]$mature_mirna_id,
    miranda = example_result[example_result$database == "miranda",]$mature_mirna_id,
    mirdb = example_result[example_result$database == "mirdb",]$mature_mirna_id,
    pictar = example_result[example_result$database == "pictar",]$mature_mirna_id,
    pita = example_result[example_result$database == "pita",]$mature_mirna_id
  )
  
  
  n <- length(data_list)
  intersection_matrix <- matrix(0, n, n, dimnames = list(names(data_list), names(data_list)))
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        intersection_matrix[i, j] <- length(intersect(data_list[[i]], data_list[[j]]))
      }
    }
  }
  print(intersection_matrix)
  combined_list1[[an]] <- intersection_matrix
  
  intersection_df <- as.data.frame(as.table(intersection_matrix))
  intersection_df$Gene <- an
  
  
  combined_list[[an]] <- intersection_df
}
combined_list1
combined_df <- do.call(rbind, combined_list)
for (j in hub_gene) {
  example <- tryCatch({
    print(j)
    get_multimir(org = "hsa",                  
                 target = j,                   
                 table = "predicted",          
                 summary = TRUE,
                 predicted.cutoff = p,        
                 predicted.cutoff.type = "p",  
                 predicted.site = "all") 
  }, error = function(e) {
    
    cat("Error in get_multimir for target:", j, "\n", e$message, "\n")
    return(NULL)  
  })
  
  
  if (is.null(example)) {
    next
  }
  
  
  table(example@data$type)  
  example_result <- example@data
  head(example_result)
  
  
  targetscan = example_result[example_result$database == "targetscan",]$mature_mirna_id
  diana_microt = example_result[example_result$database == "diana_microt",]$mature_mirna_id
  miranda = example_result[example_result$database == "miranda",]$mature_mirna_id
  mirdb = example_result[example_result$database == "mirdb",]$mature_mirna_id
  pictar = example_result[example_result$database == "pictar",]$mature_mirna_id
  pita = example_result[example_result$database == "pita",]$mature_mirna_id
  
  
  miRNA <- intersect(targetscan, pita)
  miRNA1 <- intersect(miRNA,diana_microt)
  
  
  if (length(miRNA) == 0) {
    next
  }
  
  
  write.csv(miRNA, paste0(j, "_miRNA.csv"))
}
file_list <- list.files(pattern = "miRNA.csv$", full.names = F)
merged_data <- data.frame()
for (file in file_list) {
  file_data <- read.csv(file, header = TRUE, row.names = 1, sep = ",")
  if (nrow(file_data) == 0) {
    next
  }
  file_data$ID <- strsplit(file, "_")[[1]][1]
  selected_columns <- file_data[, c("x", "ID"), drop = FALSE]
  merged_data <- rbind(merged_data, selected_columns)
}
head(merged_data)
colnames(merged_data)<-c("miRNA","mRNA")
write.csv(merged_data, "merged_data.csv", row.names = FALSE)
starbase <- read.table("/data/nas1/anbingzhuang/project/database/lncRNA_miRNA_interaction.txt",header = T,sep = "\t")
merge_gene <- read.csv("merged_data.csv")
lnc_data <- starbase[starbase$geneType == "lincRNA",]
all_lncRNA <- lnc_data[lnc_data$miRNAname %in% merge_gene$miRNA,c("miRNAname","geneName")]
colnames(all_lncRNA)<-c("miRNA","lncRNA")
write.csv(all_lncRNA, "all_lncRNA.csv", row.names = FALSE)
colnames(merge_gene) <- c("RNA1","RNA2")
colnames(all_lncRNA) <- c("RNA1","RNA2")
all_m_lnc_mi <- rbind(merge_gene, all_lncRNA)
write.csv(all_m_lnc_mi,"all_net.csv")
