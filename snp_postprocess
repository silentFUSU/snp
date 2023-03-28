rm(list=ls())
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-redhat-linux-gnu-library/4.2","/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R/lib/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/")
#source("code/library.R")
# library(AnnoProbe)
library(vegan)
library(tibble)
#library(infercnv)
library(stringr)
library(ggplot2)
library(vcfR)
library(VennDiagram)
# library(Rsamtools)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(tidyr)
load("data/merge_data/before_doubletfinder/snp_Rdata/snp.rdata")
pri_da<-as.data.frame(pri_se$seurat_clusters)
rownames(pri_da)<-paste0("pri_",rownames(pri_da))
pri_da$rownames<-rownames(pri_da)
pri_da$`pri_se$seurat_clusters`<-paste0("pri_",pri_da$`pri_se$seurat_clusters`)
pri_da<-as.data.frame(pri_da[which(rownames(pri_da)%in%pri_se_name),])
colnames(pri_da)[1]<-"seurat_cluster"

rec_da<-as.data.frame(rec_se$seurat_clusters)
rownames(rec_da)<-paste0("rec_",rownames(rec_da))
rec_da$rownames<-rownames(rec_da)
rec_da$`rec_se$seurat_clusters`<-paste0("rec_",rec_da$`rec_se$seurat_clusters`)
rec_da<-as.data.frame(rec_da[which(rownames(rec_da)%in%rec_se_name),])
colnames(rec_da)[1]<-"seurat_cluster"
da<-rbind(pri_da,rec_da)
df.empty <- data.frame(matrix(ncol = length(overlap), nrow = (length(pri_se_name)+length(rec_se_name))))
rownames(df.empty)<-rownames(da)
df.empty <- cbind(seurat_cluster=da[,1],df.empty)
colnames(df.empty)[2:ncol(df.empty)]<-overlap
df.empty[1:5,1:5]
for(i in 1:length(overlap)){
  for(j in 1:2){
    if(file.info(paste0("result/merge_data/before_doubletfinder/snp/snp_s04_pri/",overlap[i],"-pri1.txt"))$size>0 &
       file.info(paste0("result/merge_data/before_doubletfinder/snp/snp_s04_pri/",overlap[i],"-pri2.txt"))$size>0){
      barcode<- read.table(paste0("result/merge_data/before_doubletfinder/snp/snp_s04_pri/",overlap[i],"-pri",j,".txt"))
      barcode$V1<-paste0("pri_",barcode$V1)
      df.empty[which(rownames(df.empty)%in%barcode$V1),i+1]<-(3-j)
    }
  }
}

for(i in 1:length(overlap)){
  for(j in 1:2){
    if(file.info(paste0("result/merge_data/before_doubletfinder/snp/snp_s04_rec/",overlap[i],"-rec1.txt"))$size>0 &
       file.info(paste0("result/merge_data/before_doubletfinder/snp/snp_s04_rec/",overlap[i],"-rec2.txt"))$size>0){
      barcode<- read.table(paste0("result/merge_data/before_doubletfinder/snp/snp_s04_rec/",overlap[i],"-rec",j,".txt"))
      barcode$V1<-paste0("rec_",barcode$V1)
      df.empty[which(rownames(df.empty)%in%barcode$V1),i+1]<-(3-j)
    }
  }
}
df.empty[is.na(df.empty)]<-0

saveRDS(df.empty,"result/merge_data/before_doubletfinder/snp/s04_19037.rds")
df.empty<-readRDS("result/merge_data/before_doubletfinder/snp/s04_19037.rds")
df.empty2<-df.empty[which(rowSums(df.empty[,2:ncol(df.empty)]) > 0),]
#df.empty$group<-"pri"
#df.empty$group[(length(pri_se_name)+1):(length(pri_se_name)+length(rec_se_name))]<-"rec"

pca1 <- prcomp(df.empty[,2:ncol(df.empty)],center = F,scale = F)
saveRDS(df.empty,"result/merge_data/before_doubletfinder/snp/s04_19037_pca.rds")
df1 <- pca1$x
df1 <- as.data.frame(df1)
summ1 <- summary(pca1)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")
df1$group<-df.empty$group
ggplot(data = df1,aes(x = PC1,y = PC2,color = group))+
 # 添加置信椭圆
  geom_point(size = 3.5, alpha = 0.05)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_colour_manual(values = c("purple","orange","pink"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
ggsave(filename = "result/merge_data/before_doubletfinder/snp/s04_19037_PCA.png",width = 10,height = 10)


annotation_col<-as.data.frame(df.empty[,19038])
colnames(annotation_col)[1]<-"condition"
rownames(annotation_col)<-rownames(df.empty)
df.empty<-as.data.frame(t(df.empty))
df.empty<-as.matrix(df.empty)
df.empty2<-as.data.frame(df.empty[1:200,])
df.empty2[1:200,]<-as.numeric(df.empty2[1:200,])
df.empty2<-as.matrix(df.empty2)
df.empty<-df.empty[-19038,]
df.empty[1:19037,]<-as.numeric(df.empty[1:19037,])
pdf("result/merge_data/before_doubletfinder/snp/s04_19037.pdf", width = 30,height = 50)
Heatmap(
  df.empty,cluster_columns = T,cluster_rows = T,col = c("blue","red","black","white")
)
dev.off()
se <-readRDS("data/merge_data/before_doubletfinder/annotation_rds/add_all_normal/s04_pri_rec_merge.rds")
rownames(se@meta.data)
rownames<-as.data.frame(rownames(df.empty))
rownames<-separate(rownames,`rownames(df.empty)`,c("1","2"),"_")
rownames<-paste0(rownames$`1`,"_cancer_",rownames$`2`)
DimPlot(se,group.by = "orig.ident",cells.highlight = rownames)

rownames <- as.data.frame(rownames(se@meta.data))
colnames(rownames)[1]<-"barcode"
rownames <- as.data.frame(rownames[which(str_detect(rownames$barcode,"cancer")),])
colnames(rownames)[1]<-"barcode"
rownames<-separate(rownames,barcode,c("c1","c2","c3"),"_")
rownames$barcode<- paste0(rownames$c1,"_",rownames$c3)
rownames<-as.data.frame(rownames$barcode)
write.table(rownames,"result/merge_data/before_doubletfinder/gatk/barcode_cancer/S04.txt",row.names = F,col.names = F,quote = FALSE)
