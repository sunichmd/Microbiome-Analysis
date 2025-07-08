#!/usr/bin/env Rscript
# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Jingying Zhang, Yong-Xin Liu, et. al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录



#----1.1 程序功能描述和主要步骤#----

# 程序功能：差异比较热图展示
# Functions: pheatmap for different compare
# Main steps:
# - Reads data matrix and design
# - Prepare annotation for rows and columns
# - Save heatmap in input folder

# 清空工作环境 Clean enviroment object
rm(list=ls())


#----1.2 安装CRAN来源常用包#----
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("ggplot2","pheatmap","dplyr")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="Cervico_9m-CON.txt",
                help="Feature table [default %default]"),
    make_option(c("-d", "--design"), type="character", default="metadata_total.txt",
                help="Design file or metadata [default %default]"),
    make_option(c("-A", "--group"), type="character", default="Group",
                help="Group name [default %default]"),
    make_option(c("-l", "--remove"), type="numeric", default=7,
                help="remove the first column [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="Output directory; name according to input [default %default]"),
    make_option(c("-s", "--fontsize"), type="numeric", default=7,
                help="Figure width [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=8,
                help="Figure width [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=8,
                help="Figure heidth [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  suppressWarnings(dir.create(opts$output))
}

#----2. 读取输入文件#----

# 读取比较列表
input = read.table(opts$input, header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors=F)
# 筛选显著差异部分
input = subset(input, level %in% c("Enriched","Depleted"))
res_filename = paste0(opts$output, ".heatmap.pdf", sep="")
plot_wid = opts$width
plot_hei = opts$height
plot_size = opts$fontsize
design = read.table(opts$design, header=T, row.names=1, sep="\t", comment.char="")
# 统一改实验列为group
design$group = design[[opts$group]]

norm = input[,-(1:opts$remove)]

if (dim(norm)[1]>0){

  idx = rownames(design) %in% colnames(norm)
  design = design[idx,]

  anno_row = data.frame(Level = input$level, row.names = rownames(input))
  anno_col = data.frame(Group = design$group, row.names = rownames(design))
  norm = norm[apply(norm,1,sum)>0,]

  ## 注释文件存在时，添加物种注释，不聚类分组
  if (file.exists("11hehesvvv.txt")){
    taxonomy = read.table("11hehesvvv.txt", sep = "\t", row.names=1, header=T)
    # per = read.delim("result/taxonomy.txt", sep = "\t", row.names=1, header=T)
    # mean = rowMeans(per)
    # per = as.data.frame(mean[order(mean, decreasing = T)])
    # top_tax=head(rownames(per), n=10)
    #
    # x = input
    #
    # # 将低丰度的门变为Low Abundance
    # x$tax = x$Phylum# factor(x$Phylum, levels=c(as.vector(unique(x$Phylum)),"Low Abundance"))
    # # 将门中 proteobacteria 替换为纲
    # x[x$tax %in% "Proteobacteria",]$tax =  x[x$tax %in% "Proteobacteria",]$Class # no level can get value
    #
    # # 判断是否有需要替换为低丰度的类，没有报错的解决
    # if (length(x[!(x$tax %in% top_tax),]$tax > 0)){
    #   x[!(x$tax %in% top_tax),]$tax = "Low Abundance" # no level can get value
    # }
    #
    # # 颜色还是不能保证一致，因为不同组门数量不同？？
    # x$tax = factor(x$tax, levels=sort(c(top_tax,"Low Abundance")))
    input$Phylum = taxonomy[rownames(input),"Phylum"]
    anno_row = data.frame(Level = input$level, Taxonomy = input$Phylum, row.names = rownames(input))

	  pheatmap(norm,
	           scale = "row",
	           cutree_rows=2,cutree_cols = 2,
	           cluster_cols = TRUE,
	           annotation_col = anno_col,
	           annotation_row = anno_row,
	           filename = res_filename,
	           width=plot_wid, height=plot_hei,
	           annotation_names_row= T,annotation_names_col=T,
	           show_rownames=T,show_colnames=T,
	           fontsize=plot_size,display_numbers=F)
   }else{
	  pheatmap(norm,
	           scale = "row",
	           cutree_rows=2,cutree_cols = 2,
	           cluster_cols = F, cluster_rows = T, 
	           annotation_col = anno_col,
	           annotation_row = anno_row,
	           filename = res_filename,
	           width=plot_wid, height=plot_hei,
	           annotation_names_row= T,annotation_names_col=T,
	           show_rownames=T,show_colnames=T,
	           fontsize=plot_size,display_numbers=F)
	 }
  # # 提示工作完成
  print(paste("Output in ",res_filename," finished.", sep = ""))
}
