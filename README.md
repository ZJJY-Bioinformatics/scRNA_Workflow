
# 单细胞基本流程分析

## 原始数据通过Cellranger进行基因定量

如果你有公司给的，或者从网上download的单细胞下机测序的原始数据，是可以通过10x公司的Cellranger进行流程处理。

# 1. 第一步：整理测序数据和meta表

> 如果是一个建库样本，无论测多少次(多次加测)，都可以直接用`cellranger count`来做分析，只是需要把多次测序所对应的fq文件放在一个文件下
> 
> 注意：一个样本必须要**有且只有一个**对应的文件夹，里面存放这个样本的下机测序数据

将Meta表整理成**不含表头的两列tsv表格**，第一列是这个样本存放fq文件的目录路径，第二列是这个样本的分析名称

例如：
```shell
/home/wangjiaxuan/scRNA/sample_1  case1
/home/wangjiaxuan/scRNA/sample_2  case2
/home/wangjiaxuan/scRNA/sample_3  case3
/home/wangjiaxuan/scRNA/sample_1  control1
/home/wangjiaxuan/scRNA/sample_2  control2
/home/wangjiaxuan/scRNA/sample_3  control3
```
> 注意必选是个文件夹


如上表格所示，我们就拥有了一个标准的meta表

# 2. 第二步：运行cellranger

目前脚本只支持**细胞悬液制备的，人组织**的样本分析。如果是细胞核悬液制备or其他物种组织的分析，请联系脚本的维护人员。或者查看附录，了解cellranger的参数和命令，自己单独运行

将上一步的meta表用`-i`参数输入即可

```shell
./run_scRNA.py -i sample.tsv
```
脚本会自动qsub投递，等投递任务跑完，就可以坐等结果，期间如遇到任何报错请联系脚本维护人员.

结果会放在`cellranger_count_out/{样本名称的文件夹种}`，其中最有用的，下一步分析需要的是`cellranger_count_out/{样本名称的文件夹种}/outs/filtered_feature_bc_matrix`文件夹中的文件。
以下就是outs文件夹结果输出的模版，其中filtered_feature_bc_matrix就是我们所说下一步导入Seurat中所需的文件

```
drwxr-sr-x 2  4.0K Apr 12 21:50 filtered_feature_bc_matrix
-rw-r--r-- 1  8.6M Apr 12 21:50 filtered_feature_bc_matrix.h5
-rw-r--r-- 1   651 Apr 12 21:57 metrics_summary.csv
-rw-r--r-- 1   87M Apr 12 21:51 molecule_info.h5
drwxr-sr-x 2  4.0K Apr 12 21:47 raw_feature_bc_matrix
-rw-r--r-- 1   17M Apr 12 21:47 raw_feature_bc_matrix.h5
-rw-r--r-- 1  2.4M Apr 12 21:57 web_summary.html
```

# 3. 第三步: 输入到Seurat中分析

## 3.1 针对单个样本的分析

```
/data/wangjiaxuan/biosoft/miniconda3/bin/Rscript ./script/scRNA_s1.r -i cellranger_count_out/{样本名称}/{样本名称}_result/outs/filtered_feature_bc_matrix
```

> 注意有多少个样本，就要运行多少次

输出结果目录：

```
.
├── 1.QC # 质控文件夹
├── 2.cell cluster # 细胞聚类
├── 2.remove_doublet # 去重双细胞
├── 3.cluster filter # 没有的文件，除非又要过滤的cluster
├── 4.Diff expression # 差异表达文件
└── 5.Marker_gene_display # marker基因展示，一般没结果，除非指定marker基因
```

## 3.2 针对多个样本的整合分析

当然不可能只有一个样本的，分析中多个样本需要整合在一起，作细胞注释、差异分析等。运行整合时候，脚本需要输入一个文件夹，文件中必须放入所有需要整合样本的rds文件。

```
mkdir -p scRNA_data/merge_out
cp scRNA_data/{sample1}_report2023-04-21/{sample1}.rds scRNA_data/merge_out
cp scRNA_data/{sample2}_report2023-04-21/{sample2}.rds scRNA_data/merge_out

ls scRNA_data/merge_out
# 可以看到scRNA_data/merge_out下有两个rds，分别代表不同样本
# 最后
/data/wangjiaxuan/biosoft/miniconda3/bin/Rscript ./script/scRNA_s2.r -i scRNA_data/merge_out
```
> 注意：需要整合样本的rds放在目录下，不需要一个样本一个文件夹，这样会识别不到。

# 完结~

# 附录（下载）

### 下载基因组文件
human GRCH38 （大小11G:tea:，这个文件不仅有fa和gtf文件，而且star的指引文件都构建好了。而且10x公司针对自己的产品做了部分优化）<sup>[[3]](#ref03)</sup>

```
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
```

mouse mm10 大于10G大小:tea:

```
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
```

### 下载测序文件

接下来，从 10X 官网上的一个公开数据集下载 FASTQ 文件。本示例使用来自人类外周血单核细胞 （PBMC） 的1000 个 细胞 数据集，这些数据集由淋巴细胞（T 细胞、B 细胞和 NK 杀伤）和单核细胞组成:dart::dart::dart:。

```
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
tar -xvf pbmc_1k_v3_fastqs.tar
```

### 定量分析

我们需要用到`cellranger count`命令来进行定量分析，可以先看命令的帮助文档，帮助我们理解其参数的含义。

```
cellranger count --help
```
作为参考，我用的参数是
```
cellranger count --no-bam\ #不需要bam文件
--nosecondary\ # 不需要它分析聚类和细胞分群
--disable-ui\ # 不需要web展示
--id=run_count_1kpbmcs \ #分析名称
--fastqs=/pbmc_1k_v3_fastqs \ #fq文件目录
--sample=pbmc_1k_v3 \ #样本名称
--transcriptome=/refdata-gex-GRCh38-2020-A #刚刚下载的基因组文件
--include-introns # 一般不需要，除非是细胞核测序
```

>  明显的filtered_feature_bc_matrix文件下的是可以直接作为Seurat的输入文件路径的（就是这么方便！）。

# 维护人员

齐灿灿

王家轩

# 维护人员
罗蓉蓉

# 参考文献:

<p id = "ref01">[1] [Cell Ranger “Aggr” read-depth normalization vs. Seurat NormalizeData](https://github.com/satijalab/seurat/issues/672)

<p id = "ref02">[2] [Cellranger aggr versus Seurat](https://www.biostars.org/p/433937/)

<p id = "ref03">[3] [Build Notes for Reference Packages](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#mm10_#{files.refdata_mm10.version})

