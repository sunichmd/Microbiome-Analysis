#!/usr/bin/env Rscript

# Copyright 2016-2023 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Yong-Xin Liu, Lei Chen, Tengfei Ma, Xiaofang Li, Maosheng Zheng, Xin Zhou, Liang Chen, Xubo Qian, Jiao Xi, Hongye Lu, Huiluo Cao, Xiaoya Ma, Bian Bian, Pengfan Zhang, Jiqiu Wu, Ren-You Gan, Baolei Jia, Linyang Sun, Zhicheng Ju, Yunyun Gao, Tao Wen, Tong Chen. 2023. EasyAmplicon: An easy-to-use, open-source, reproducible, and community-based pipeline for amplicon data analysis in microbiome research. iMeta 2: e83. https://doi.org/10.1002/imt2.83

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

#----1. 参数 Parameters#----

#----1.1 功能描述 Function description#----

# 程序功能：物种组成弦图
# Functions: Taxonomy circlize

options(warn = -1) # Turn off warning


## 设置输入输出文件和参数

# 修改下面`default=`后面的文件和参数。
#
# 输入文件为特征表(otutab.txt)+分组信息(metadata.tsv)
#
# 输入文件"-i", "--input"，otutab.txt; 特征表
#
# 实验设计"-d", "--design"，默认`metadata.tsv`，可手动修改文件位置；
#
# 分组列名"-n", "--group"，默认将metadata.tsv中的Group列作为分组信息，可修改为任意列名；
#
# 分组列名"-c", "--compare_pair"，默认将比较metadata.tsv中的Group列的前两个值，建议手动设定；
#
# 分组列名"-t", "--threhold"，丰度筛选阈值，默认千分之1
#
# 图片宽"-w", "--width"，默认89 mm，根据图像布局可适当增大或缩小
#
# 图片高"-e", "--height"，默认59 mm，根据图像布局可适当增大或缩小


#----1.2 参数缺少值 Default values#----
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse")
  require("optparse",character.only=T)
}
# 基于Bioconductor安装R包
library(BiocManager)
if (!requireNamespace("edgeR", quietly = TRUE))
  BiocManager::install("edgeR")
# 基于github安装，检测没有则安装
library(devtools)
if(!requireNamespace("amplicon", quietly = TRUE))
  install_github("microbiota/amplicon")

# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="result/otutab.txt",
                help="Feature table [default %default]"),
    make_option(c("-t", "--threshold"), type="numeric", default=0.1,
                help="Threshold of relative abundance 0.1% [default %default]"),
    make_option(c("-d", "--design"), type="character", default="result/metadata.txt",
                help="Design file or metadata [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]"),
    make_option(c("-c", "--compare"), type="character", default="KO-WT",
                help="Groups comparison [default %default]"),
    make_option(c("-m", "--method"), type="character", default="wilcox",
                help="Compare method, default wilcox, alternative edgeR or t.test [default %default]"),
    make_option(c("-s", "--scale"), type="logical", default=TRUE,
                help="Normalize to 100 [default %default]"),
    make_option(c("-p", "--pvalue"), type="numeric", default=0.05,
                help="Threshold of P-value [default %default]"),
    make_option(c("-f", "--fdr"), type="numeric", default=0.1,
                help="Threshold of FDR [default %default]"),
    make_option(c("-l", "--log2Foldchange"), type="numeric", default=0,
    						help="Threshold of log2Foldchange [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/compare/",
                help="Output predix; name according to input [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}
# 设置输出文件缺省值，如果为空，则为输入+pcoa.pdf
# if(opts$output==""){opts$output=opts$input}
# suppressWarnings(dir.create(dirname(opts$output), showWarnings = F))


#----1.3. 加载包 Load packages#----

suppressWarnings(suppressMessages(library(amplicon)))


#----2. 读取文件 Read files#----

#----2.1 实验设计 Metadata#----
metadata = read.table(opts$design, header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)

#----2.1 特征表 Feature table#----
otutab = read.table(opts$input, header=T, row.names=1, sep="\t", comment.char="", quote="")



#----3. 统计保存 Stat and saving#----

#----3.1 统计 Stat#----
# 输入特征表、元数据和分组列，比较组名，方法，丰度、p和fdr值
output = compare(data = otutab, metadata = metadata,
                 group = opts$group, compare_pair = opts$compare,
                 method = opts$method, RA = opts$threshold,
                 pvalue = opts$pvalue, fdr = opts$fdr, normalize = opts$scale)

#----3.2 保存表格 Saving#----
filename = paste0(opts$output, opts$compare, ".txt")
print(filename)
write.table("ID\t", file = filename,
            append = F, quote = F, eol = "", row.names = F, col.names = F)
suppressWarnings(write.table(output, file = filename,
                             append = T, quote = F, sep = '\t', row.names = T, col.names = T))
