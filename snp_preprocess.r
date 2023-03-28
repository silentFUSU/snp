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
pri<-read.vcfR('/storage/zhangyanxiaoLab/suzhuojie/software/varCA/out_gatk_varscan/classify/S04_pri/final.vcf.gz', verbose = FALSE )
rec<-read.vcfR('/storage/zhangyanxiaoLab/suzhuojie/software/varCA/out_gatk_varscan/classify/S04_rec/final.vcf.gz', verbose = FALSE )

filter <- read.table(file='/storage/zhangyanxiaoLab/suzhuojie/software/varCA/out_gatk_varscan/classify/S04_pri/filter.tsv.gz',header = T)

pri_bed <- as.data.frame(read.table("/storage/zhangyanxiaoLab/suzhuojie/software/varCA/out_gatk_varscan/classify/S04_pri/final_filter.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

rec_bed <- as.data.frame(read.table("/storage/zhangyanxiaoLab/suzhuojie/software/varCA/out_gatk_varscan/classify/S04_rec/final_filter.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

pri_snp<-as.data.frame(pri@fix)
pri_snp$label<-paste0(pri_snp$CHROM,"-",pri_snp$POS)
pri_bed$label<-paste0(pri_bed$V1,"-",pri_bed$V3)
filter$label<-paste0(filter$CHROM,"-",filter$POS)
pri_snp<-pri_snp[which(pri_snp$label %in% pri_bed$label),]
pri_filter<-filter[which(filter$label %in% pri_snp$label),]
pri_filter<-pri_filter[-which(pri_filter$varscan.snp.RD==0),]
pri_filter<-pri_filter[-which(pri_filter$varscan.snp.AD==0),]
pri_filter<-pri_filter[-which(pri_filter$varscan.snp.RD<10),]
pri_snp<-pri_snp[which(pri_snp$label%in%pri_filter$label),]

filter2<-read.table(file='/storage/zhangyanxiaoLab/suzhuojie/software/varCA/out_gatk_varscan/classify/S04_rec/filter.tsv.gz',header = T)
rec_snp<-as.data.frame(rec@fix)
rec_snp$label<-paste0(rec_snp$CHROM,"-",rec_snp$POS)
rec_bed$label<-paste0(rec_bed$V1,"-",rec_bed$V3)
rec_snp<-rec_snp[which(rec_snp$POS %in% rec_bed$V3),]
filter2$label<-paste0(filter2$CHROM,"-",filter2$POS)
rec_filter<-filter2[which(filter2$label %in% rec_snp$label),]
rec_filter<-rec_filter[-which(rec_filter$varscan.snp.RD==0),]
rec_filter<-rec_filter[-which(rec_filter$varscan.snp.AD==0),]
rec_filter<-rec_filter[-which(rec_filter$varscan.snp.RD<10),]
rec_snp<-rec_snp[which(rec_snp$label%in%rec_filter$label),]

rec_snp$symbol<-paste0(rec_snp$label,"-",rec_snp$REF,"-",rec_snp$ALT)
pri_snp$symbol<-paste0(pri_snp$label,"-",pri_snp$REF,"-",pri_snp$ALT)
overlap <- intersect(pri_snp$symbol,rec_snp$symbol)
pri_snp2<-pri_snp[which(pri_snp$symbol %in% overlap),]
rec_snp2<-rec_snp[which(rec_snp$symbol %in% overlap),]

overlap<-read.table("result/merge_data/before_doubletfinder/snp/s04_pri_snp_position.txt",header = F,row.names = NULL)
pri_se<-readRDS("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/s04_pri.rds")
rec_se<-readRDS("data/merge_data/before_doubletfinder/annotation_rds/clean_annotation/s04_rec.rds")
pri_se_name<-rownames(pri_se@meta.data)[which(!pri_se$annotation %in% c("Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T cell"))]
pri_se_name<-paste0("pri_",pri_se_name)
rec_se_name<-rownames(rec_se@meta.data)[which(!rec_se$annotation %in% c("Astrocyte","Endothelial","Fibroblast","Microglia","Oligodendrocyte","T cell"))]
rec_se_name<-paste0("rec_",rec_se_name)
overlap<-overlap[-which(!overlap$V4%in%c("A","T","C","G") ),]
overlap<-paste0(overlap$V1,"-",overlap$V2,"-",overlap$V3,"-",overlap$V4)
save.image("data/merge_data/before_doubletfinder/snp_Rdata/snp.rdata")

