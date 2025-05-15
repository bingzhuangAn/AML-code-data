rm(list=ls())
library(GEOquery)
 
if (!dir.exists("00_Rawdata")) {dir.create("00_Rawdata")}
 
GSEID <- "GSE114868"
library(GEOquery)
library(Biobase)
library(tidyverse)
library(GEOquery)
library(Biobase)
library(tidyverse)
gset<-getGEO(GSEID,
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
gset$GSE114868_series_matrix.txt.gz@annotation
expr<-as.data.frame(exprs(gset[[1]]))
gpls = gset[[1]]@annotation  
options(stringsAsFactors = F)
gpl <- getGEO(gpls,destdir = ".")
gpl <- Table(gpl)
a=gset[[1]]
colnames(gpl)
gpl$GENE_SYMBOL<-data.frame(sapply(gpl$gene_assignment,function(x)unlist(strsplit(x," // "))[2]),
                            stringsAsFactors=F)[,1]
probe2symbol<-gpl %>%
  
  dplyr::select('ID','GENE_SYMBOL')%>%
  filter('GENE_SYMBOL'!='')%>%
  na.omit()
colnames(probe2symbol)[1] <- "ID"
colnames(probe2symbol)[2] <- 'geneSymbol'
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID <- as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     
  dplyr::select(geneSymbol,everything())%>%     
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    
  arrange(desc(rowMean))%>%       
  distinct(geneSymbol,.keep_all = T)%>%      
  dplyr::select(-rowMean)%>%     
  tibble::column_to_rownames(colnames(.)[1])   
pd<-pData(a)
colnames(pd)
group<-data.frame(sample=pd$geo_accession,group=pd$title)
group$group <- ifelse(grepl("^AML", group$group),"AML","Control")
table(group$group)
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.csv(dat,file =  paste0(GSEID,'_expr_use.csv'))
save(dat,group,file = paste0(GSEID,'.RData'))
write.csv(group,file =  paste0(GSEID,'_group.csv'))
count<-read.table('./TCGA-LAML.htseq_counts.tsv',header=T,check.names=FALSE)
rownames(count)<-count[,1]
count<-count[-1]
library(tinyarray) 
library(openxlsx)
data<- trans_exp(count)   
data<-2^data-1
counts<- as.data.frame(apply(data, 2, function(x) as.integer(x)))
rownames(counts)<-rownames(data)
library(readr)
phenotype<-read_delim('./TCGA-LAML.GDC_phenotype.tsv')%>%data.frame(.)
clinical<-subset(phenotype,sample_type.samples=='Primary Blood Derived Cancer - Peripheral Blood'|sample_type.samples=='Solid Tissue Normal')
clinical<-clinical[clinical$submitter_id.samples%in%colnames(counts),]
clinical<-clinical[order(clinical$sample_type.samples,decreasing = T),]
write.xlsx(clinical,'clinical.xlsx') 
Group<-data.frame(sample=clinical$submitter_id.samples,group=clinical$sample_type.samples)
Group$group<-ifelse(Group$group=='Solid Tissue Normal','Normal','AML')
Group<-data.frame(row.names=Group$sample ,Group=Group$group )
counts<-subset(counts,select=rownames(Group)) 
identical(colnames(counts),rownames(Group))
table(Group$Group)   
counts_all<-counts
Group_all<-Group
write.csv(counts_all,'01.TCGA_counts_all.csv',quote=F)
write.csv(Group_all,'02.TCGA_Group_all.csv',quote=F)
mycounts<-counts
mycounts$Gene<-rownames(mycounts)
gene_length<-read.table('./All_hg19gene_len.txt',header=T)
mycounts<-merge(mycounts,gene_length,by='Gene')
rownames(mycounts)<-mycounts$Gene
mycounts<-mycounts[,-1]
kb <- mycounts$Length / 1000
kb
countdata <- mycounts[,1:152]
rpk <- countdata / kb
rpk
fpkm <- t(t(rpk)/colSums(countdata) * 10^6) 
head(fpkm)
fpkms_all<-log2(fpkm+1)%>%data.frame()
colnames(fpkms_all)<-gsub('[.]','-',colnames(fpkms_all))
write.csv(fpkms_all,'03.TCGA_fpkms_all.csv',quote=F) 
fpkm <- read.table("./TCGA-LAML.htseq_fpkm.tsv",header=T,check.names=FALSE)
rownames(fpkm)<-fpkm[,1]
fpkm<-fpkm[-1]
data<- trans_exp(fpkm)   
group<-data.frame(sample=rownames(Group),Group=Group$Group)
AML <- subset(group,Group=='AML')
data <- subset(data,select=AML$sample)   
survival<-read.table('./TCGA-LAML.survival.tsv',header=T)
survival<-data.frame(row.names=survival$sample,OS.time=survival$OS.time,OS=survival$OS)
fpkms<-data[,colnames(data)%in%rownames(survival)]
survival<-survival[colnames(fpkms),]
survival<-na.omit(survival)
fpkms<-subset(fpkms,select=rownames(survival))
identical(colnames(fpkms),rownames(survival))
counts_AML<-subset(counts_all,select=rownames(survival))
survival_AML<-survival
fpkms_AML<-fpkms
identical(colnames(fpkms_AML),rownames(survival_AML))
write.csv(counts_AML,'04.counts_AML.csv',quote=F) 
write.csv(fpkms_AML,'05.fpkms_AML.csv',quote=F) 
write.csv(survival_AML,'06.survival_AML.csv',quote=F)
GSEID <- "GSE30029"
library(GEOquery)
library(Biobase)
library(tidyverse)
gset<-getGEO(GSEID,
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
gset$GSE30029_series_matrix.txt.gz@annotation
expr<-as.data.frame(exprs(gset[[1]]))
gpls = gset[[1]]@annotation  
options(stringsAsFactors = F)
gpl <- parseGEO(fname = "./GPL6947.annot.gz", GSElimits = NULL,AnnotGPL = F, 
                getGPL = F)
gpl <- Table(gpl)
a=gset[[1]]
colnames(gpl)
probe2symbol<-gpl %>%
  
  dplyr::select('ID','Gene symbol')%>%
  filter('Gene symbol'!='')%>%
  na.omit()
colnames(probe2symbol)[1] <- "ID"
colnames(probe2symbol)[2] <- 'geneSymbol'
id <- colnames(expr)
dat <- expr %>%
  t() %>%                    
  as.data.frame() %>%        
  lapply(as.numeric) %>%     
  as.data.frame()            
rownames(dat) <- id
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID <- as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     
  dplyr::select(geneSymbol,everything())%>%     
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    
  arrange(desc(rowMean))%>%       
  distinct(geneSymbol,.keep_all = T)%>%      
  dplyr::select(-rowMean)%>%     
  tibble::column_to_rownames(colnames(.)[1])   
pd<-pData(a)
colnames(pd)
group<-data.frame(sample=pd$geo_accession,group=pd$title)
group$group <- ifelse(grepl("^AML", group$group),"AML","Control")
table(group$group)
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.csv(dat,file =  paste0(GSEID,'_expr_use.csv'))
save(dat,group,file = paste0(GSEID,'.RData'))
write.csv(group,file =  paste0(GSEID,'_group.csv'))
GSEID <- "GSE12417"
library(GEOquery)
library(Biobase)
library(tidyverse)
gset<-getGEO(GSEID,
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
gset$GSE12417_series_matrix.txt.gz@annotation
expr<-as.data.frame(exprs(gset[[1]]))
gpls = gset[[1]]@annotation  
options(stringsAsFactors = F)
gpl <- parseGEO(fname = "./GPL10558.annot.gz", GSElimits = NULL,AnnotGPL = F, 
                getGPL = F)
gpl <- Table(gpl)
a=gset[[1]]
colnames(gpl)
probe2symbol<-gpl %>%
  
  dplyr::select('ID','Gene symbol')%>%
  filter('Gene symbol'!='')%>%
  na.omit()
colnames(probe2symbol)[1] <- "ID"
colnames(probe2symbol)[2] <- 'geneSymbol'
id <- colnames(expr)
dat <- expr %>%
  t() %>%                    
  as.data.frame() %>%        
  lapply(as.numeric) %>%     
  as.data.frame()            
rownames(dat) <- id
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID <- as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     
  dplyr::select(geneSymbol,everything())%>%     
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    
  arrange(desc(rowMean))%>%       
  distinct(geneSymbol,.keep_all = T)%>%      
  dplyr::select(-rowMean)%>%     
  tibble::column_to_rownames(colnames(.)[1])   
pd<-pData(a)
colnames(pd)
group<-data.frame(sample=pd$geo_accession,group=pd$title)
group$group <- ifelse(grepl("AML$", group$group),"AML","Control")
table(group$group)
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.csv(dat,file =  paste0(GSEID,'_expr_use.csv'))
save(dat,group,file = paste0(GSEID,'.RData'))
write.csv(group,file =  paste0(GSEID,'_group.csv'))
colnames(pd)
write.csv(pd,"pd.csv")
OS_dat_GSE71014 <- pd[,c(35,36)]
library(tidyverse)
OS_dat_GSE71014 <- OS_dat_GSE71014 %>%
  separate(`event (1:ch1`, into = c("column1", "column2", "OS"), sep = ":", fill = "right") %>%
  mutate(`event (1:ch1` = OS) %>%  
  select(-c("column1","column2","event (1:ch1"))  
colnames(OS_dat_GSE71014) <- c("OS","OS.time")
write.csv(OS_dat_GSE71014,"OS_dat_GSE71014.csv")
GSEID <- "GSE71014"
library(GEOquery)
library(Biobase)
library(tidyverse)
gset<-getGEO(GSEID,
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
gset$GSE71014_series_matrix.txt.gz@annotation
expr<-as.data.frame(exprs(gset[[1]]))
gpls = gset[[1]]@annotation  
options(stringsAsFactors = F)
gpl <- parseGEO(fname = "./GPL10558.annot.gz", GSElimits = NULL,AnnotGPL = F, 
                getGPL = F)
gpl <- Table(gpl)
a=gset[[1]]
colnames(gpl)
probe2symbol<-gpl %>%
  
  dplyr::select('ID','Gene symbol')%>%
  filter('Gene symbol'!='')%>%
  na.omit()
colnames(probe2symbol)[1] <- "ID"
colnames(probe2symbol)[2] <- 'geneSymbol'
id <- colnames(expr)
dat <- expr %>%
  t() %>%                    
  as.data.frame() %>%        
  lapply(as.numeric) %>%     
  as.data.frame()            
rownames(dat) <- id
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID <- as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     
  dplyr::select(geneSymbol,everything())%>%     
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    
  arrange(desc(rowMean))%>%       
  distinct(geneSymbol,.keep_all = T)%>%      
  dplyr::select(-rowMean)%>%     
  tibble::column_to_rownames(colnames(.)[1])   
pd<-pData(a)
colnames(pd)
group<-data.frame(sample=pd$geo_accession,group=pd$title)
group$group <- ifelse(grepl("AML$", group$group),"AML","Control")
table(group$group)
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.csv(dat,file =  paste0(GSEID,'_expr_use.csv'))
save(dat,group,file = paste0(GSEID,'.RData'))
write.csv(group,file =  paste0(GSEID,'_group.csv'))
colnames(pd)
write.csv(pd,"pd.csv")
OS_dat_GSE71014 <- pd[,c(35,36)]
library(tidyverse)
OS_dat_GSE71014 <- OS_dat_GSE71014 %>%
  separate(`event (1:ch1`, into = c("column1", "column2", "OS"), sep = ":", fill = "right") %>%
  mutate(`event (1:ch1` = OS) %>%  
  select(-c("column1","column2","event (1:ch1"))  
colnames(OS_dat_GSE71014) <- c("OS","OS.time")
write.csv(OS_dat_GSE71014,"OS_dat_GSE71014.csv")
GSEID <- "GSE12417"
library(GEOquery)
library(Biobase)
library(tidyverse)
gset <- getGEO(filename = "GSE12417-GPL97_series_matrix.txt.gz",getGPL = F)
expr <- gset@assayData$exprs %>% as.data.frame()
options(stringsAsFactors = F)
gpl <- parseGEO(fname = "./GPL97.annot.gz", GSElimits = NULL,AnnotGPL = F, 
                getGPL = F)
gpl <- Table(gpl)
colnames(gpl)
gpl$GENE_SYMBOL<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),
                            stringsAsFactors=F)[,1]
probe2symbol<-gpl %>%
  
  dplyr::select('ID','GENE_SYMBOL')%>%
  filter('GENE_SYMBOL'!='')%>%
  na.omit()
colnames(probe2symbol)[1] <- "ID"
colnames(probe2symbol)[2] <- 'geneSymbol'
dat <- expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID <- as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     
  dplyr::select(geneSymbol,everything())%>%     
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    
  arrange(desc(rowMean))%>%       
  distinct(geneSymbol,.keep_all = T)%>%      
  dplyr::select(-rowMean)%>%     
  tibble::column_to_rownames(colnames(.)[1])   
library(tidyr)
pd<-pData(gset)
OS <- pd[,c(2,10)]
OS$characteristics_ch1
OS_split <- separate(OS, col = characteristics_ch1,
                     into = c("stage","Age","OS.time","OS"), sep = ";")
OS <- separate(OS_split, col = OS.time,
               into = c("a","b","d","OS.time","c"),sep = " ")
OS_final <- separate(OS, col = OS,
                     into = c("n","OS"),sep = ":")
colnames(OS_final)
OS <- OS_final[,c(1,10,7)]
write.csv(dat,file =  paste0(GSEID,'97_expr_use.csv'))
write.csv(OS,"OS_dat_GSE12417_97.csv")
GSEID <- "GSE37642"
library(GEOquery)
library(Biobase)
library(tidyverse)
gset <- getGEO(filename = "GSE37642-GPL96_series_matrix.txt.gz",getGPL = F)
expr <- gset@assayData$exprs %>% as.data.frame()
options(stringsAsFactors = F)
gpl <- parseGEO(fname = "./GPL96.annot.gz", GSElimits = NULL,AnnotGPL = F, 
                getGPL = F)
gpl <- Table(gpl)
colnames(gpl)
gpl$GENE_SYMBOL<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),
                            stringsAsFactors=F)[,1]
probe2symbol<-gpl %>%
  
  dplyr::select('ID','GENE_SYMBOL')%>%
  filter('GENE_SYMBOL'!='')%>%
  na.omit()
colnames(probe2symbol)[1] <- "ID"
colnames(probe2symbol)[2] <- 'geneSymbol'
dat <- expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID <- as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     
  dplyr::select(geneSymbol,everything())%>%     
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    
  arrange(desc(rowMean))%>%       
  distinct(geneSymbol,.keep_all = T)%>%      
  dplyr::select(-rowMean)%>%     
  tibble::column_to_rownames(colnames(.)[1])   
library(tidyr)
pd<-pData(gset)
OS <- pd[,c(2,40,41)]
colnames(OS) <- c("id","OS","OS.time")
OS$OS <- ifelse(OS$OS == "dead",1,0)
write.csv(dat,file =  paste0(GSEID,'_96_expr_use.csv'))
write.csv(OS,paste0(GSEID,"OS_dat_96.csv"))
