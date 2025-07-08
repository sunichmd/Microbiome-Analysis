#!/usr/bin/env Rscript

# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Jingying Zhang, Yong-Xin Liu, et. al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

#----1. 参数 Parameters#----

#----1.1 功能描述 Function description#----

# 程序功能：物种组成弦图
# Functions: Taxonomy circlize

options(warn = -1) # Turn off warning


## 设置输入输出文件和参数

# 修改下面`default=`后面的文件和参数。
#
# 输入文件为物种相对丰度矩阵(sum_p.txt)+分组信息(metadata.tsv)
#
# 输入文件"-i", "--input"，sum_p.txt; 物种组成表
#
# 实验设计"-d", "--design"，默认`metadata.tsv`，可手动修改文件位置；
#
# 分组列名"-n", "--group"，默认将metadata.tsv中的Group列作为分组信息，可修改为任意列名；
#
# 图片宽"-w", "--width"，默认89 mm，根据图像布局可适当增大或缩小
#
# 图片高"-e", "--height"，默认59 mm，根据图像布局可适当增大或缩小
#
# 无法指定输出文件名，默认保存于输入目录circlize.pdf和circlize_legend.pdf

#----1.2 参数缺少值 Default values#----
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="result/tax/sum_p.txt",
                help="Taxonomy composition [default %default]"),
    make_option(c("-l", "--legend"), type="numeric", default=8,
                help="Legend number [default %default]"),
    make_option(c("-d", "--design"), type="character", default="result/metadata.txt",
                help="Design file or metadata [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="Output predix; name according to input [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=89,
                help="Figure width in mm [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=59,
                help="Figure heidth in mm [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  #suppressWarnings(dir.create(opts$output))
}
# 设置输出文件缺省值，如果为空，则为输入+pcoa.pdf
if(opts$output==""){opts$output=opts$input}


#----1.3. 加载包 Load packages#----

suppressWarnings(suppressMessages(library(amplicon)))


#----2. 读取文件 Read files#----

#----2.1 实验设计 Metadata#----
metadata = read.table(opts$design, header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)

#----2.1 物种组成Taxonomy table#----
taxonomy = read.table(opts$input, header=T, row.names=1, sep="\t", comment.char="")


#----3. 绘图保存 Plotting and saving#----

#----3.1 绘图样品 Plotting#----
# 输入物种组成、元数据和分组列，保存circlize.pdf和circlize_legend.pdf
# tax_sum=taxonomy
# groupID = "Group"
# topN = opts$legend
# groupID = opts$group
# output=opts$output
# width = opts$width
# height = opts$height
tax_circlize = function(tax_sum, metadata, topN = 5, groupID = "Group",output,width,height) {
	
	# 依赖关系检测与安装
	p_list = c("ggplot2", "reshape2", "circlize")
	for(p in p_list){
		if (!requireNamespace(p)){install.packages(p)}
		suppressPackageStartupMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
	}
	
	# 测试默认参数
	# library(amplicon)
	# tax_sum = tax_phylum
	# topN = 5
	# groupID = "Group"
	
	# 交叉筛选
	idx = rownames(metadata) %in% colnames(tax_sum)
	metadata = metadata[idx,,drop=F]
	tax_sum = tax_sum[, rownames(metadata)]
	
	# 提取样品组信息,默认为group可指定
	sampFile = as.data.frame(metadata[, groupID],row.names = row.names(metadata))
	colnames(sampFile)[1] = "group"
	
	#----按丰度降序排序#----
	mean_sort = as.data.frame(tax_sum[(order(-rowSums(tax_sum))), ])
	# 筛选前N类，其它归为Other，可设置不同组数
	other = colSums(mean_sort[topN:dim(mean_sort)[1], ])
	mean_sort = mean_sort[1:(topN - 1), ]
	mean_sort = rbind(mean_sort,other)
	rownames(mean_sort)[topN] = c("Other")
	# 保存变量备份，并输出至文件
	merge_tax=mean_sort
	
	#----按组合并求均值#----
	
	# 转置样品名添加组名，并去除多余的两个样品列
	mat_t = t(merge_tax)
	mat_t2 = merge(sampFile, mat_t, by="row.names")
	mat_t2 = mat_t2[,c(-1)]
	
	# 按组求均值，转置，再添加列名
	mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
	# df = do.call(rbind, mat_mean)[-1,]
	df = t(mat_mean[,-1])
	geno = mat_mean$group
	colnames(df) = mat_mean$group
	
	#----默认参数绘图-颜色随机#----
	# 颜色设定
	library(RColorBrewer)
	grid.col = NULL
	# 分类学颜色，最多12种
	qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
	col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
	set.seed(623);manual1 = sample(col_vector,topN)
	grid.col[rownames(df)] = manual1
	# 定义分组颜色
	grid.col[colnames(df)] = brewer.pal(dim(df)[2], "Accent")
	# 设置图片文件名、长宽和字体大小
	pdf(file=paste0(output,".circlize.pdf"), width=width/25.4, height=height/25.4, pointsize=8)
	# 上方绘图和图例代码
	chordDiagram(t(df), directional = T,diffHeight = 0.02, grid.col = grid.col, transparency = 0.5,
							 annotationTrack = c("grid"), preAllocateTracks = 1)
	circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
		xlim = get.cell.meta.data("xlim")
		ylim = get.cell.meta.data("ylim")
		sector.name = get.cell.meta.data("sector.index")
		circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
								niceFacing = T, adj = c(-0.1, 0.5),cex=0.6) # 字体大小调整cex # 字体展示格式调整facing
		circos.axis(h="bottom",direction="inside",labels.cex=0.4)
		}, bg.border = NA,
	)
	circos.clear()
	# 绘图结束后写入文件
	dev.off()
	
	#----指定颜色绘图+图例#----
	pdf(file=paste0(output,".circlize_legend.pdf"), width=width/25.4, height=height/25.4, pointsize=8)
	# 上方绘图和图例代码
	chordDiagram(t(df), directional = TRUE,diffHeight = 0.03, grid.col = grid.col, transparency = 0.5,
							 annotationTrack = c("grid"), preAllocateTracks = 1)
	# 添加行-物种图例
	legend("left",pch=20,legend=rownames(df),col=grid.col[rownames(df)],bty="n",cex=1,pt.cex=3,border="black")
	# 添加列-分组图例
	legend("right",pch=20,legend=colnames(df),col=grid.col[colnames(df)],bty="n",cex=1,pt.cex=3,border="black")
	# 绘图结束后写入文件
	circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
		xlim = get.cell.meta.data("xlim")
		ylim = get.cell.meta.data("ylim")
		sector.name = get.cell.meta.data("sector.index")
		circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
								niceFacing = TRUE, adj = c(-0.1, 0.5),cex=0.8) # 字体大小调整cex # 字体展示格式调整facing
		circos.axis(h="bottom",direction="inside",labels.cex=0.4)
	}, bg.border = NA)
	circos.clear()
	dev.off()
}

tax_circlize(taxonomy, metadata, topN = opts$legend, groupID = opts$group,output=opts$output,width = opts$width, height = opts$height)
