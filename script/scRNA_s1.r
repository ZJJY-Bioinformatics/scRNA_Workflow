#.libPaths("/hwfssz5/ST_INFECTION/HPV/wangjiaxuan/biosoft/miniconda3/envs/R/lib/R/library/")

library("optparse")
option_list <- list(
  make_option(c("-i","--input"),type="character"),
  make_option(c("-n","--name"),type = "character")
)

args <- parse_args(OptionParser(option_list=option_list))

data_dir = args$input
object_name = args$name

if(is.null(data_dir)){
  stop("no input dir path")
}

if(is.null(object_name)){
  stop("the parameter -n is needed")
}

# data_dir = "cellranger_count_out/metastatic/metastatic_result/outs/filtered_feature_bc_matrix"
# name = "meta"
# 加载R包-----------------
require(Spat)
load_spat_env()
# 本地运行需要加载--------
file_path = data_dir
f_min = 200
f_max = 10000
m_max = 10
time = Sys.Date()
# 文件夹设置---------------
if(!dir.exists("cRNA_data/")){ dir.create("scRNA_data/")}
parameter_workdir = paste("scRNA_data/",object_name,"_report",time,sep = "")
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
stereo.data <- Read10X(data.dir = file_path)
stereo <- CreateSeuratObject(counts = stereo.data,project = object_name,min.cells = 3,min.features = 200)
rm(stereo.data)
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
VlnPlot(stereo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"),sort = "increasing",pt.size = 0,ncol = 5,cols = "#4d8657")*
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
stereo <- RunPCA(stereo)
stereo <- RunUMAP(stereo, dims = 1:10)

# 去除双胞----------------------
#QC2 去除双细胞
library(DoubletFinder)
#去除多重
## pK Identification (no ground-truth) 
sweep.res.list <- paramSweep_v3(stereo, PCs = 1:10, sct = F)
# gt.calls <- stereo@meta.data[rownames(sweep.res.list[[1]]), "GT"]
# sweep.stats <- summarizeSweep(sweep.res.list, GT = TRUE, GT.calls = gt.calls)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpk<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
## Homotypic Doublet Proportion Estimate ------------------------------------------------------------
annotations <- stereo@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: 
nExp_poi <- round(0.075*nrow(stereo@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
stereo <- doubletFinder_v3(stereo, PCs = 1:10, pN = 0.25, pK = mpk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
# stereo <- doubletFinder_v3(stereo, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
rm(sweep.res.list,sweep.stats,bcmvn,annotations
   ,homotypic.prop,nExp_poi,nExp_poi.adj)
doublet_name = grep("DF.classifications",colnames(stereo@meta.data),value = T)
stereo@meta.data[["doublet"]] = stereo@meta.data[[doublet_name]]
stereo@meta.data[[doublet_name]] = NULL

stereo_raw = stereo
stereo = subset(stereo_raw, subset = doublet == "Singlet")
DimPlot(stereo_raw,group.by = "doublet")
ggplot2::ggsave(paste(dd_dir,"/","1.1_doublet_siglet_umap.pdf",sep = ""),
                width = 8,height = 6,dpi = 300)
write_tsv(as.data.frame(table(stereo_raw$doublet)),paste(dd_dir,"/","1.1_doublet_siglet_group.xls",sep = ""))
# 标准化分群------------------
stereo <- SCTransform(object = stereo, verbose = TRUE)
stereo <- RunPCA(stereo, features = VariableFeatures(object = stereo))
ElbowPlot(stereo)
ggplot2::ggsave(paste(QC_dir,"/","3.1_PC_select.pdf",sep = ""),
                width = 8,height = 6,dpi = 300)
# 细胞分群--------------------
stereo <- FindNeighbors(stereo, dims = 1:15)
stereo <- FindClusters(stereo, resolution = 0.5)
stereo <- RunUMAP(stereo, dims = 1:15)
stereo <- RunTSNE(stereo, dims = 1:15)
# 绘制图形--------------------
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

# 不怎么需要分析的--------
# stereo.markers <- FindAllMarkers(stereo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(stereo.markers, file = paste(diff_dir,"/","cluster_DEG.rds",sep = ""))
# write.csv(stereo.markers, paste(diff_dir,"/","1.diff_exp_gene_exp_in_every_cell_cluster.csv",sep = ""))

#保存结果--------------
saveRDS(stereo,file = paste(parameter_workdir,"/",object_name,".rds",sep = ""))

