library("optparse")
option_list <- list(
  make_option(c("-i","--input"),type="character"),
  make_option(c("-n","--name"),type = "character")
)

args <- parse_args(OptionParser(option_list=option_list))

data_dir = args$input
object_name = args$name

# 工作路径设置-----
setwd("./")
# 读取数据路径-----
type = "scRNA-seq"
wkdir  =  getwd()
Batch_list <- list.files(data_dir,pattern = "\\.rds",full.names = T)
rds_Path = file.path(wkdir,Batch_list)
#------------------------------------------------------------

library(dplyr)
library(Seurat)

#华大测序数据读取需要加：，gene.column = 1
# min.cells = 3 代表筛选至少在三个细胞中表达的基因
# project="R002" 写入这个Seurat对象的初始标识
# min.features = 300 代表筛选至少表达300个基因的细胞，低质量的细胞或空飞沫通常很少有基因
rds_list = list()
for (i in 1:length(rds_Path)) {
  Batch <- readRDS(rds_Path[i])  
  #Batch  <- RenameCells(Batch ,add.cell.id = i)#不同样本之间会存在重名
  rds_list[[i]] = Batch
}
rm(Batch)

# normalize and identify variable features for each dataset independently
rds_list  <- lapply(X = rds_list , FUN = function(x) {
  x@active.assay = "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = rds_list )

# Step2 数据整合 --------------------------------------------------------------
anchors <- FindIntegrationAnchors(object.list = rds_list, anchor.features = features)
# this command creates an 'integrated' data assay
stereo.combined <- IntegrateData(anchorset = anchors)
#-----------------------
require(Spat)
load_spat_env()
time = Sys.Date()
object_name = "merge"
# 文件夹设置---------------
parameter_workdir = paste(object_name,"report",time,sep = "_")
dir.create(parameter_workdir)
QC_dir = paste("./",parameter_workdir,"/","1.QC",sep = "")
dir.create(QC_dir)
dd_dir = paste("./",parameter_workdir,"/","2.remove_doublet",sep = "")
dir.create(dd_dir)
cell_dir = paste("./",parameter_workdir,"/","2.cell cluster",sep = "")
dir.create(cell_dir)
cell_dir_cluster = paste(cell_dir,"/","4.split every cluster in umap",sep = "")
dir.create(cell_dir_cluster)
filter_dir = paste(parameter_workdir,"/","3.cluster filter",sep = "")
dir.create(filter_dir)
#---
diff_dir = paste(parameter_workdir,"/","4.Diff expression",sep = "")
dir.create(diff_dir)
#---#---
diff_dir_vlnplot = paste(diff_dir,"/","1.top20_diff_gene_vlnplot",sep = "")
dir.create(diff_dir_vlnplot)
diff_dir_sheatmap = paste(diff_dir,"/","2.top20_diff_gene_exp_in_umap",sep = "")
dir.create(diff_dir_sheatmap)
diff_dir_enrich = paste(diff_dir,"/","3.unique_DEG_enrichment",sep = "")
dir.create(diff_dir_enrich)
#---
annot_dir = paste(parameter_workdir,"/","5.Marker_gene_display",sep = "")
dir.create(annot_dir)
# Read data----------------
stereo = stereo.combined
rm(anchors,rds_list,stereo.combined)
DefaultAssay(stereo) = "RNA"
stereo$orig.ident = factor(stereo$orig.ident,levels = unique(stereo$orig.ident))
stereo@active.ident = stereo$orig.ident
# qc过滤-------------------
stereo[["percent.mt"]] <- PercentageFeatureSet(stereo, pattern = "^MT-")
#写入核糖体基因
stereo[["percent.rb"]] <- PercentageFeatureSet(stereo, pattern = "^RP[SL]")
#写入红细胞基因
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(stereo))
stereo[["percent.HB"]] <- PercentageFeatureSet(stereo, features=HB.genes)
#rm(HB.genes)

#绘制质控前QC图-----------
#可选："percent.rb","percent.HB"
VlnPlot(stereo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"),sort = "increasing",pt.size = 0,ncol = 5)*
  geom_boxplot(colour = "#002333",position=position_dodge(0.6),size=0.5,width=0.2)
ggplot2::ggsave(paste(QC_dir, "/","1.1_cell_QC-VlnPlot_before_filter.pdf",sep = ""),width = 18,height = 8,dpi = 300)
#  qc输出
meta_data = as.data.frame(stereo@meta.data)
skimr::skim(meta_data) %>% 
  tibble::as_tibble() %>%
  select("Qc_type" = "skim_variable","mean"="numeric.mean","median"="numeric.p50") %>%
  drop_na() -> meta_data
write.csv(meta_data,paste(QC_dir,"/","1.2_cell_QC_before_filter.csv",sep = ""))
#QC1 以基因表达量，线粒体基因数目过滤
#nFeature_RNA < 3000 筛选表达在3000个细胞以下的基因，防止双细胞或多细胞读数
#percent.mt < 10 筛选线粒体基因表达较低的细胞，低质量/濒临死亡的细胞经常表现出线粒体污染
stereo <- subset(stereo, subset =  nFeature_RNA < 7500 & percent.mt < 40)
#绘制过滤后QC图---------
VlnPlot(stereo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"),sort = "increasing",pt.size = 0,ncol = 5)*
  geom_boxplot(colour = "#002333",position=position_dodge(0.6),size=0.5,width=0.2)
ggplot2::ggsave(paste(QC_dir, "/","2.2_cell_QC-VlnPlot_after_filter.pdf",sep = ""),width = 18,height = 8,dpi = 300)
#  qc输出
meta_data = as.data.frame(stereo@meta.data)
skimr::skim(meta_data) %>% 
  tibble::as_tibble() %>%
  select("Qc_type" = "skim_variable","mean"="numeric.mean","median"="numeric.p50") %>%
  drop_na() -> meta_data
write.csv(meta_data,paste(QC_dir,"/","2.2_cell_QC_after_filter.csv",sep = ""))

# 数据标准化--------------------
stereo <- NormalizeData(stereo)
stereo <- FindVariableFeatures(stereo, selection.method = "vst", nfeatures = 2000)
stereo <- ScaleData(stereo)
stereo <- SCTransform(object = stereo, verbose = TRUE)
#聚类分析----------
DefaultAssay(stereo) <- "integrated"

# Run the standard workflow for visualization and clustering
stereo <- ScaleData(stereo, verbose = FALSE)
stereo <- RunPCA(stereo, npcs = 30, verbose = FALSE)
stereo <- RunUMAP(stereo, reduction = "pca", dims = 1:30)
stereo <- FindNeighbors(stereo, reduction = "pca", dims = 1:30)
stereo <- FindClusters(stereo, resolution = 0.5)
stereo <- RunTSNE(stereo, dims = 1:30)
# 绘制图形--------------------
stereo = readRDS("s1_next.rds")
# 选择颜色
clu_num = length(levels(stereo@active.ident))
colpal1 = colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(clu_num)
#
DimPlot(stereo, reduction = "umap",label = TRUE,pt.size = 0.5,cols = colpal1)
ggplot2::ggsave(paste(cell_dir,"/","1.1_cluster_umap_point.pdf",sep = ""),
                width = 8,height = 6,dpi = 300)
DimPlot(stereo, reduction = "tsne",label = TRUE,pt.size = 0.5,cols = colpal1)
ggplot2::ggsave(paste(cell_dir,"/","1.2_cluster_tsne_point.pdf",sep = ""),
                width = 8,height = 6,dpi = 300)
# 分细胞群看质控--------------
VlnPlot(stereo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3,cols = colpal1)
ggplot2::ggsave(paste(cell_dir, "/","2.2_cluster_QC-VlnPlot_vlnplot.pdf",sep = ""),width = 12,height = 8,dpi = 300)

DimPlot(stereo, reduction = "umap",label = TRUE,pt.size = 0.5,group.by = "orig.ident")
ggplot2::ggsave(paste(cell_dir,"/","1.3_sample_umap_point.pdf",sep = ""),
                width = 8,height = 6,dpi = 300)

# 细胞群之间的相似性-----------
cluster_exp_all <- AverageExpression(
  stereo,
  assays = "SCT",
  features = NULL,
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  slot = "data",
  verbose = TRUE
)

cluster_exp = cluster_exp_all[["SCT"]]
cor(cluster_exp) -> cor_matrix
redgreen = colorRampPalette(c("#00441b", "white", "#e31a1c"))(100)

pheatmap::pheatmap(cor_matrix,cluster_cols = T,cluster_rows = T,
                   show_rownames = T,color=redgreen,fontsize_col = 8,
                   cellheight = 30,scale = "row",
                   cellwidth = 30,fontsize_row  = 8,treeheight_col = 8,
                   angle_col = 45,display_numbers = T,
                   filename = paste(cell_dir,"/","3.1_cluster_exp_cor_heatmap.pdf",sep = ""))
dev.off()
DefaultAssay(stereo) = "SCT"
stereo.markers <- FindAllMarkers(stereo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(stereo.markers, file = paste(diff_dir,"/","cluster_DEG.rds",sep = ""))
write.csv(stereo.markers, paste(diff_dir,"/","1.diff_exp_gene_exp_in_every_cell_cluster.csv",sep = ""))

#保存结果--------------
saveRDS(stereo,file = paste(parameter_workdir,"/",object_name,".rds",sep = ""))
#-------------
require(Spat)
load_spat_env()

stereo =  readRDS(paste(parameter_workdir,"/",object_name,".rds",sep = ""))
DefaultAssay(stereo) = "RNA"

stereo.markers = readRDS(file = paste(diff_dir,"/","cluster_DEG.rds",sep = ""))
write.csv(stereo.markers, paste(diff_dir,"/","1.diff_exp_gene_exp_in_every_cell_cluster.csv",sep = ""))

stereo.markers %>%
  dplyr::count(gene) %>%
  filter(n == 1) %>%
  left_join(stereo.markers, by = "gene") %>%
  select(-c(2)) -> stereo.markers.unique_raw

write.csv(stereo.markers.unique_raw, paste(diff_dir,"/","2.unique_gene_diff_exp_in_every_cell_cluster.csv",sep = ""))

# 展示marker基因-------------------------
stereo.markers.unique_raw %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> stereo.markers.unique
# 基因表达量可视化-小提琴图
write.csv(stereo.markers.unique, paste(diff_dir,"/","3.top20_unique_diff_gene_exp_gene.csv",sep = ""))

maker_gene <- stereo.markers.unique %>%
  group_by(cluster) %>%
  select("gene")

for(i in 1:dim(maker_gene)[1]){
  VlnPlot(stereo, features = maker_gene$gene[i], slot = "data",cols = colpal1)
  ggplot2::ggsave(paste(diff_dir_vlnplot,"/",maker_gene$cluster[i],"_",maker_gene$gene[i],".pdf",sep = ""),
                  width = 8,height = 6,dpi = 300)
  cat(paste("[",i,"]","the gene", maker_gene$gene[i] ,"has finised !","\n"))
}
# slot = "count" 等同与用原始count值来画图,=scale.data等同于二次标准化的数据
# 基因表达细胞分群热图

for(i in 1:dim(maker_gene)[1]){
  FeaturePlot(stereo, features = maker_gene$gene[i], slot = "data",keep.scale = "feature",pt.size = 1,cols = c("lightgrey","#02786A","#EB4B17"))
  ggplot2::ggsave(paste(diff_dir_sheatmap,"/",maker_gene$cluster[i],"_",maker_gene$gene[i],"exp_in_umap.pdf",sep = ""),
                  width = 8,height = 6,dpi = 300)
  cat(paste("[",i,"]","the gene", maker_gene$gene[i] ,"has finised !","\n"))
}


# 基因表达热图

DoHeatmap(stereo, features = maker_gene$gene,assay = "RNA",slot = "data")
ggplot2::ggsave(paste(diff_dir,"/","4.top_diff_gene_heatmap.pdf",sep = ""),
                width = 8,height = 11,dpi = 300)

# 富集分析

library(org.Hs.eg.db)
library(clusterProfiler)

for(i in levels(stereo@active.ident)){
  stereo.markers.unique_raw %>%
    dplyr::filter(cluster == i) %>%
    dplyr::select("gene") -> geneid
  gene1 = bitr(geneid$gene,fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)
  tryCatch({ego <- enrichGO(gene          = gene1$ENTREZID,
                            OrgDb         = org.Hs.eg.db,
                            keyType = "ENTREZID",
                            pvalueCutoff = 0.5)
  dotplot(ego)
  ggplot2::ggsave(paste(diff_dir_enrich,"/cluster",i,"_","go.pdf",sep = ""),
                  width = 8,height = 10,dpi = 300)
  }, error = function(e){})
  cat(paste("[",i,"]","the cluster" ,"has finised !","\n"))
}

for(i in levels(stereo@active.ident)){
  stereo.markers.unique_raw %>%
    dplyr::filter(cluster == i) %>%
    dplyr::select("gene") -> geneid
  gene1 = bitr(geneid$gene,fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)
  tryCatch({
    ek <- enrichKEGG(gene = gene1$ENTREZID,organism = "hsa",pvalueCutoff = 0.5)
    dotplot(ek)
    ggplot2::ggsave(paste(diff_dir_enrich,"/cluster",i,"_","ko.pdf",sep = ""),width = 8,height = 10,dpi = 300)
  }, error=function(e){})
  cat(paste("[",i,"]","the cluster" ,"has finised !","\n"))
}

# 细胞注释--------------------------
# marker list
ab_marker = read_csv("info/markergene_0528.csv")
colnames(ab_marker) <- c("id","module")
ab_marker =  unique(ab_marker) %>% drop_na()
ab_marker$id = sub(" ","",ab_marker$id)
genename = rownames(stereo)

for(i in unique(ab_marker$module)){
  module_dir = paste(annot_dir,"/",i,sep = "")
  if(!dir.exists(module_dir)){dir.create(module_dir)}
  data = ab_marker %>% filter(module == i)
  for(j in data$id){
    if(any(grepl(j,genename,ignore.case = TRUE))){
      tryCatch({
        FeaturePlot(stereo,features = j,slot = "data")+
        VlnPlot(stereo,features = j,slot = "data",cols = colpal1)
        ggsave(paste(module_dir,"/","ve_",j,"_marker_gene_display_chip.pdf",sep = ""),width = 10,height = 8,dpi = 300)
      }, error=function(e){})
    } else {
      cat(paste("gene",j,"not find !","\n"))
    }
  }
  cat(paste("[",i,"]","cluster" ,"marker gene has finised !","\n"))
}


im_marker = read_csv("info/markergene_0528.csv")
#im_marker = im_marker[,c(4,2)]
colnames(im_marker) <- c("id","module")
im_marker =  unique(im_marker) %>% drop_na()
im_marker$id = sub(" ","",im_marker$id)

ident_rows = function(x){
  m = apply(x, 1, FUN  = function(x) quantile(x,probs = c(0.6),names = F))
  return(x > m)
}

cluster_exp_all <- AverageExpression(
  stereo,
  assays = "RNA",
  features = NULL,
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  slot = "data",
  verbose = TRUE
)
cluster_exp = cluster_exp_all[["RNA"]]

cluster_exp_ident = ident_rows(cluster_exp)

as.data.frame(cluster_exp_ident)  %>%
  rownames_to_column("id") %>%
  right_join(im_marker, by = "id") %>%
  select(-"id") %>%
  replace(is.na(.), 0) %>%
  group_by(module) %>%
  summarise(across(everything(), list(sum))) -> all_cluster_avg

all_cluster_avg %>% column_to_rownames("module") -> all_cluster_avg
redgreen = colorRampPalette(c("#00441b", "white", "#e31a1c"))(100)

pheatmap::pheatmap(all_cluster_avg,cluster_cols = F,cluster_rows = F,
                   show_rownames = T,color=redgreen,fontsize_col = 8,
                   cellheight = 15,scale = "row",
                   cellwidth = 45,fontsize_row  = 8,treeheight_row = 0,treeheight_col = 8,
                   angle_col = 45,filename = paste(annot_dir,"/","1.immune_gene_exp_ratio0.6_module_heatmap.pdf",sep = ""))
dev.off()


rat_dir = paste("./",parameter_workdir,"/","6.Ratio",sep = "")
dir.create(rat_dir)
setwd(rat_dir)


#细胞注释-----------
setwd("../")
ann_dir = paste("./",parameter_workdir,"/","7.Anotation",sep = "")
dir.create(ann_dir)
setwd(ann_dir)

stereo$orig_cluster = Idents(stereo)
# stereo$second_cluster = factor(stereo$sub_cluster,levels = unique(stereo$sub_cluster))
# stereo@active.ident = stereo$second_cluster

current.cluster.ids <- seq(0,22)

name1 <- c("Enterocytes", "Enterocytes", "Enterocytes", "Naive T Cells", "Gut Cells", "Endothelial Cells", "DCs", "Enterocytes", "M1 Macrophages", "CD14+ Monocytes", "CD8+ T cells", "Smooth Muscle Cells", "Inflammatory Associated Fibroblasts", "Endothelial Cells", "Gut Cells", "Inflammatory Associated Fibroblasts", "M2 Macrophages", "Fibroblasts", "NK Cells", "Blood Endothelial Cells", "Gut Cells", "B Cells", "Enteroendocrines")

stereo$cellname1 <- plyr::mapvalues(x = stereo$orig_cluster, from = current.cluster.ids, to = name1)

stereo@active.ident = stereo$cellname1

levels(stereo) = c("Gut Cells", "Enterocytes", "Enteroendocrines", "Inflammatory Associated Fibroblasts", "Fibroblasts", "Smooth Muscle Cells", "Endothelial Cells", "Blood Endothelial Cells", "M1 Macrophages", "M2 Macrophages", "CD14+ Monocytes", "DCs", "B Cells", "NK Cells", "CD8+ T cells", "Naive T Cells")

saveRDS(stereo ,"~/Data/rdata/chen/merge_report_2022-05-28/merge.rds")

colpal1 = scCancer::getDefaultColors(n = length(levels(stereo)),type = 2)
DimPlot(object = stereo, reduction = "tsne", pt.size = 0.5,group.by = "cellname1",cols = colpal1,label = T,repel = T)

ggplot2::ggsave("annotation_tsne_cluster.pdf",width = 10,height = 6 ,dpi = 300)

DimPlot(object = stereo, reduction = "umap", pt.size = 0.5,group.by = "cellname1",cols = colpal1,label = T,repel = T)

ggplot2::ggsave("annotation_UMAP_cluster.pdf",width = 10,height = 6 ,dpi = 300)
#-------
stereo =  readRDS("merge_data/merge_report_2022-07-20/merge.rds")
library(reshape2)
library(ggplot2)
df<-(t(table(stereo@meta.data[c('orig.ident', 'cellname2')])))
df<-as.data.frame(df)
df<-dcast(df, cellname2~orig.ident)
write.table(df, file='./clusters_count.csv', sep='\t', quote=F, row.names=F)
for(col in colnames(df)){
  if(col=='cellname2'){
    next
  }
  print(col)
  all_sum = sum(df[col])
  print(df[col])
  df[col]<-prop.table(df[col])
}
write.table(df, file='./clusters_ratio.xls', sep='\t', quote=F, row.names=F)

#细胞堆积图
cell.prop<-as.data.frame(prop.table(table(
  stereo@meta.data[c('cellname2','orig.ident')])),stringsAsFactors = FALSE)
colnames(cell.prop)<-c("cluster","patient","proportion")
cell.prop$cluster<- factor(cell.prop[,1],
                           levels = levels(stereo@meta.data$cellname2))

clu_num=length(levels(stereo$cellname2))
colpal1 = colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(clu_num)

ggplot(cell.prop,aes(patient,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("cell_prop")+
  theme_bw(base_size = 10)+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,vjust=1,size=14))+
  guides(fill=guide_legend(title="cluster"))+
  scale_fill_manual(values = colpal1)
ggplot2::ggsave("cell_prop.jpeg",
                width = 6,height = 6 ,dpi = 300)
require(tidyverse)
df %>% 
  pivot_longer(-cellname2,names_to = "sample",values_to = "cell_ratio") %>%
  ggplot(aes(x = cellname2, y = cell_ratio, fill = sample))+
  geom_col(position = "dodge")+
  #scale_fill_manual(values = c("#e82c45","#34569d"))+
  facet_wrap(cellname2~.,scales = "free")+
  theme(strip.text.x = element_text(size = 10),
        axis.title.x = element_blank())+
  ggdark::dark_theme_minimal()
ggplot2::ggsave("cell_prop_facet.jpeg",
                width = 12,height = 10 ,dpi = 300)
#样本堆积图
ggplot(cell.prop,aes(cluster,proportion,fill=patient))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("patient_prop")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,vjust=1,size=14))+
  guides(fill=guide_legend(title=NULL))
ggplot2::ggsave("patient_prop.pdf",
                width = 8,height = 6 ,dpi = 300)

stereo$group = factor(gsub("\\d*-","",stereo$orig.ident),levels = c("ZC","CY"))

df = as.data.frame(table(stereo$group,stereo$cellname2))
colnames(df) = c("group","cellname","freq")

ggplot(df,aes(x = group, y = freq, fill = group)) +
  geom_col()+
  facet_wrap(cellname~.,scales = "free")+
  scale_fill_manual(values = xbox::need_colors("two5"))+
  theme(strip.text.x = element_text(size = 10),
        axis.title.x = element_blank())+
  ggdark::dark_theme_minimal()
ggplot2::ggsave("cell_prop_facet.jpeg",
                width = 8,height = 10 ,dpi = 300)
  
# 比较两个rds之间的区别
rm(list = ls())
setwd("/Users/dalena/Data/rdata/chen/merge_data/")
st1 = readRDS("merge_report_2022-07-20/merge.rds")

current.cluster.ids <- levels(st1@meta.data[["seurat_clusters"]])
name1 <-c("CD8+ T cells", "Enterocytes", "Enterocytes", "M2 Macrophages", "Endothelial cells", "T memory cells", "Enterocytes", "CD4+ T cells", "Smooth muscle cells/Enteric glia cells", "M1 Macrophages", "B cells", "Enterocytes", "Inflammatory Associated Fibroblasts", "Enterocytes", "Tuft cells", "T memory cells", "Inflammatory Associated Fibroblasts", "Smooth muscle cells", "NK cells", "Enterocytes", "Endothelial cells", "Enterocytes", "Enterocytes", "Enterocytes", "Erythroid-like and erythroid precursor cells")
st1$cellname2 <- plyr::mapvalues(x = st1$seurat_clusters, from = current.cluster.ids, to = name1)
require(Spat)
load_spat_env()
DimPlot(st1,group.by = "cellname2",reduction = "umap")
ggsave("0720-anno-umap.pdf")
saveRDS(st1,"merge_report_2022-07-20/merge.rds")
rename_data1 = data.frame("id" = st1@assays$RNA@counts@Dimnames[[2]],"name" = st1@meta.data$cellname2)
colnames(rename_data1) = c("id","name")

levels(st2@meta.data[["orig.ident"]])

levels(st1@meta.data[["orig.ident"]])

rename_data1$id = sub("_3","_5",rename_data1$id)
rename_data1$id = sub("_4","_6",rename_data1$id)

st2 = readRDS("merge_report_2022-08-16/merge.rds")
rename_data2 = data.frame("id" = st2@assays$RNA@counts@Dimnames[[2]],"name" = st2@meta.data$seurat_clusters)
colnames(rename_data2) = c("id","name")
#rename_data2$id = sub("_\\d","",rename_data2$id)
data = left_join(rename_data2,rename_data1,by = "id")

all(st2@assays$RNA@counts@Dimnames[[2]] == data$id)

st2@meta.data$map_cluster = data$name.y

DimPlot(st2,group.by = "map_cluster",cols = colpal1)
ggsave("0720-maping-0816-anno-umap.pdf",width = 12,height = 6)

clu_num=length(levels(st2$seurat_clusters))
colpal1 = colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(clu_num)

as.data.frame(table(st2$map_cluster,st2$seurat_clusters))%>% 
  rename("celltype" = Var1,"cluster" = Var2) %>%
  ggplot(aes(x = celltype,y = Freq,fill = as.character(cluster)))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values = colpal1)+
  cowplot::theme_cowplot(20)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("0720-maping-0816-anno-ratio_bar.pdf",width = 8,height = 6)

write.csv(as.data.frame(table(st2$map_cluster,st2$seurat_clusters)),"0720-maping-0816-anno-ratio_bar.csv")

rm(list = ls())
setwd("/Users/dalena/Data/rdata/chen/merge_data/")
st1 = readRDS("merge_report_2022-08-16/merge.rds")

current.cluster.ids <- levels(st1@meta.data[["seurat_clusters"]])
name1 <-c("CD4+ T cells", "T memory cells", "Enterocytes", "B cells", "CD8+ T cells", "Endothelial cells", "Inflammatory Associated Fibroblasts", "M2 Macrophages", "Smooth muscle cells/Enteric glia cells", "Enterocytes", "Enterocytes", "M1 Macrophages", "T memory cells", "Tuft cells", "T memory cells", "Enterocytes", "T memory cells", "NK cells", "Smooth muscle cells", "T memory cells", "T memory cells", "Endothelial cells", "Enterocytes")

  
st1$cellname2 <- plyr::mapvalues(x = st1$seurat_clusters, from = current.cluster.ids, to = name1)
require(Spat)
load_spat_env()
clu_num=length(levels(st1$seurat_clusters))
colpal1 = colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(clu_num)

DimPlot(st1,group.by = "cellname2",reduction = "umap",cols = colpal1)
ggsave("0820-anno-umap.pdf",width = 12,height = 8)
saveRDS(st1,"merge_report_2022-08-16/merge.rds")


###----------
require(Seurat)
stereo = readRDS("merge_report_2022-07-20/merge.rds")
stereo = st1
rm(st1)
tcell = subset(stereo,idents = c(0L,1L,4L,12L, 14L, 16L,19L, 20L))

# 数据标准化--------------------
tcell <- NormalizeData(tcell)
tcell <- FindVariableFeatures(tcell, selection.method = "vst", nfeatures = 2000)
tcell <- ScaleData(tcell)
tcell <- SCTransform(object = tcell, verbose = TRUE)
#聚类分析----------
DefaultAssay(stereo) <- "integrated"
require(Seurat)
# Run the standard workflow for visualization and clustering
tcell <- ScaleData(tcell, verbose = FALSE)
tcell <- RunPCA(tcell, npcs = 30, verbose = FALSE)
tcell <- FindNeighbors(tcell, reduction = "pca", dims = 1:30)
tcell <- FindClusters(tcell, resolution = 0.5)
tcell <- RunUMAP(tcell, reduction = "pca", dims = 1:30)
tcell <- RunTSNE(tcell, dims = 1:30)
# 绘制图形--------------------
#stereo = readRDS("s1_next.rds")
# 选择颜色
clu_num = length(levels(tcell@active.ident))
colpal1 = colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(clu_num)
#
DimPlot(tcell, reduction = "umap",label = TRUE,pt.size = 0.5,cols = colpal1)
ggplot2::ggsave("tcell_umap_point.pdf",
                width = 8,height = 6,dpi = 300)

FeaturePlot(tcell,features = c("CD4","CD8A"),reduction = "umap")


DimPlot(tcell, reduction = "tsne",label = TRUE,pt.size = 0.5,cols = colpal1)
ggplot2::ggsave("tcell_tsne_point.pdf",
                width = 8,height = 6,dpi = 300)

saveRDS(tcell,"only_tcell_2_cluster.rds")
