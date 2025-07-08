#!/usr/bin/env Rscript

# Copyright 2016-2021 Yong-Xin Liu <yxliu@genetics.ac.cn / metagenome@126.com>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, Yang Bai. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein Cell 2021(12) 5:315-330 doi: 10.1007/s13238-020-00724-8

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

#----1. 参数 Parameters#----

#----1.1 功能描述 Function description#----

# 程序功能：物种组成堆叠柱状图
# Functions: Taxonomy stackplot

options(warn = -1) # Turn off warning


## 设置输入输出文件和参数

# 修改下面`default=`后面的文件和参数。
#
# 输入文件为物种相对丰度矩阵(sum_p/c/o/f/g.txt)+分组信息(metadata.txt)
#
# 输入文件"-i", "--input"，tax/sum_p.txt; 物种组成表
#
# 实验设计"-d", "--design"，默认`metadata.txt`，可手动修改文件位置；
#
# 分组列名"-n", "--group"，默认将metadata.txt中的Group列作为分组信息，可修改为任意列名；
#
# 图片宽"-w", "--width"，默认89 mm，根据图像布局可适当增大或缩小
#
# 图片高"-e", "--height"，默认59 mm，根据图像布局可适当增大或缩小


#----1.2 参数缺省值 Default values#----
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
    make_option(c("-d", "--design"), type="character", default="result/metadata.txt",
                help="Design file or metadata [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]"),
    make_option(c("-s", "--subgroup"), type="character", default=NULL,
    						help="Sub Group name [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="Output directory; name according to input [default %default]"),
    make_option(c("-l", "--legend"), type="numeric", default=12,
                help="Legend number [default %default]"),
    make_option(c("-c", "--color"), type="character", default="Paired",
                help="color ggplot, manual1, Paired or Set3 [default %default]"),
    make_option(c("-g", "--groupOrder"), type="character", default=NULL,
    						help="Group value level [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=181,
                help="Figure width in mm [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=118,
                help="Figure heidth in mm [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
  # suppressWarnings(dir.create(opts$output))
}
# 设置输出文件缺省值，如果为空，则为输入+pcoa.pdf
if(opts$output==""){opts$output=opts$input}


#----1.3. 加载包 Load packages#----

suppressWarnings(suppressMessages(library(amplicon)))


#----2. 读取文件 Read files#----

#----2.1 实验设计 Metadata#----
metadata = read.table(opts$design, header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)

#----2.2 物种组成矩阵Taxonomy matrix#----
taxonomy = read.table(opts$input, header=T, row.names=1, sep="\t", comment.char="", quote = "")



#----3. 绘图保存 Plotting and saving#----
# 设置颜色22种
colorset1 = c('#98d66d','#45b1f9','#ffa6f0','#f76fc3','#85c1b8','#a584ff','#ffb444','#c45cc0','#7ebfe5','#cec0c9','#467584','#005ff9','#bc8c38','#bcba6d','#91749e','#b2798d','#fcef5d','#b23301','#235931',"#892e4f",'#fabc75','#f75c39')
# 调颜色包：Paired/Set3
library(RColorBrewer)

#----3.1 绘图样品 Plotting#----
# 输入矩阵矩阵、元数据和分组列，返回ggplot2对象
tax_stackplot1 = function(tax_sum,
				 metadata,
				 topN = 8,
				 groupID = "Group",
				 style = "group",
				 sorted = "abundance",
				 subgroup = NULL,
				 groupOrder = NULL) {
	# 依赖关系检测与安装
	p_list = c("ggplot2", "reshape2")
	for (p in p_list) {
		if (!requireNamespace(p)) {
			install.packages(p)
		}
		library(
			p,
			character.only = TRUE,
			quietly = TRUE,
			warn.conflicts = FALSE
		)
	}
	
	# 测试默认参数
	# library(amplicon)
	# tax_sum = tax_phylum
	# topN = 8
	# groupID = "Group"
	# style = "group"
	# sorted = "abundance"
	
	# 交叉筛选
	idx = rownames(metadata) %in% colnames(tax_sum)
	metadata = metadata[idx, , drop = F]
	tax_sum = tax_sum[, rownames(metadata)]
	
	# 提取样品组信息,默认为group可指定
	sampFile = as.data.frame(metadata[, c(groupID, subgroup)], row.names = row.names(metadata))
	colnames(sampFile)[1] = "group"
	
	# 按丰度降序排序
	mean_sort = as.data.frame(tax_sum[(order(-rowSums(tax_sum))),])
	# 把末分类调整到最后面
	idx = grepl("unassigned|unclassified|unknown",
							rownames(mean_sort),
							ignore.case = T)
	mean_sort = rbind(mean_sort[!idx, ], mean_sort[idx, ])
	# 筛选前N类，其它归为Other，可设置不同组数
	other = colSums(mean_sort[topN:dim(mean_sort)[1],])
	mean_sort = mean_sort[1:(topN - 1),]
	mean_sort = rbind(mean_sort, other)
	rownames(mean_sort)[topN] = c("Other")
	# 保存变量备份，并输出至文件
	merge_tax = mean_sort
	
	if (!is.null(subgroup)){
		style = "subgroup"
	}
	
	if (style == "sample") {
		# 添加分类学并设置排序方式，默认字母，abundancer按丰度
		mean_sort$Taxonomy = rownames(mean_sort)
		data_all = as.data.frame(melt(mean_sort, id.vars = c("Taxonomy")))
		if (sorted == "abundance") {
			data_all$Taxonomy  = factor(data_all$Taxonomy, levels = rownames(mean_sort))
		}
		# set group facet
		data_all = merge(data_all, sampFile, by.x = "variable", by.y = "row.names")
		
		# 按group分面，x轴范围自由否则位置异常，swith设置标签底部，并调置去除图例背景
		# 关闭x轴刻度和标签
		if (!is.null(groupOrder)){
			library(stringr)
			group_order = as.vector(str_split(groupOrder,",",simplify = T))
			data_all$group = factor(data_all$group,levels = group_order)
		}
		p = ggplot(data_all, aes(x = variable, y = value, fill = Taxonomy)) +
			geom_bar(stat = "identity",
							 position = "fill",
							 width = 1) +
			scale_y_continuous(labels = scales::percent) +
			facet_grid(~ group, scales = "free_x", switch = "x") +
			theme(strip.background = element_blank()) +
			theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
			xlab("Groups") + ylab("Percentage (%)") +
			theme_classic() + theme(axis.text.x = element_text(
				angle = 45,
				vjust = 1,
				hjust = 1
			)) +
			theme(text = element_text(family = "sans", size = 7))
		p
	} else{
		if (!is.null(subgroup)) {
			mat_t = t(merge_tax)
			mat_t2 = merge(sampFile, mat_t, by = "row.names")
			mat_t2 = mat_t2[, c(-1)]
			
			ncol_sampFile <- ncol(sampFile)
			# �������ֵ��ת�ã�����������
			mat_mean = aggregate(mat_t2[, -(1:ncol_sampFile)], by = mat_t2[c(1:ncol_sampFile)], FUN =
													 	mean) # mean
			#mat_mean_final = do.call(rbind, mat_mean)[-1,]
			#geno = mat_mean$group
			#colnames(mat_mean_final) = geno
			#mean_sort=as.data.frame(mat_mean_final)
			
			# ���ӷ���ѧ����������ʽ��Ĭ����ĸ��abundancer�����
			#mean_sort$Taxonomy = rownames(mean_sort)
			#data_all = as.data.frame(melt(mean_sort, id.vars=c("Taxonomy")))
			data_all = melt(mat_mean, variable.name = "Taxonomy")
			
			if (sorted == "abundance") {
				data_all$Taxonomy  = factor(data_all$Taxonomy, levels = rownames(mean_sort))
			}
			data_all$value = as.numeric(data_all$value)
			
			variable_en = sym(subgroup)
			if (!is.null(groupOrder)){
				library(stringr)
				group_order = as.vector(str_split(groupOrder,",",simplify = T))
				data_all$group = factor(data_all$group,levels = group_order)
			}
			p = ggplot(data_all, aes(
				x = !!variable_en,
				y = value,
				fill = Taxonomy
			)) +
				geom_bar(stat = "identity",
								 position = "fill",
								 width = 0.7) +
				scale_y_continuous(labels = scales::percent) +
				facet_grid(~ group, scales = "free_x", switch = "x") +
				xlab("Groups") + ylab("Percentage (%)") + theme_classic() +
				theme(text = element_text(family = "sans", size = 7))
		} else {
			# 按组合并求均值
			
			# 转置样品名添加组名，并去除多余的两个样品列
			mat_t = t(merge_tax)
			mat_t2 = merge(sampFile, mat_t, by = "row.names")
			mat_t2 = mat_t2[, c(-1)]
			
			# 按组求均值，转置，再添加列名
			mat_mean = aggregate(mat_t2[, -1], by = mat_t2[1], FUN = mean) # mean
			mat_mean_final = do.call(rbind, mat_mean)[-1, ]
			geno = mat_mean$group
			colnames(mat_mean_final) = geno
			mean_sort = as.data.frame(mat_mean_final)
			write.table(mean_sort,paste0(opts$output,".plot_Group_relative_abu.txt"),sep="\t",quote=F,col.names = NA)
			
			
			# 添加分类学并设置排序方式，默认字母，abundancer按丰度
			mean_sort$Taxonomy = rownames(mean_sort)
			#res_mean <<- mean_sort
			#print("变量res_mean已创建！")
			
			data_all = as.data.frame(melt(mean_sort, id.vars = c("Taxonomy")))
			if (sorted == "abundance") {
				data_all$Taxonomy  = factor(data_all$Taxonomy, levels = rownames(mean_sort))
			}
			data_all$value = as.numeric(data_all$value)
			p = ggplot(data_all, aes(x = variable, y = value, fill = Taxonomy)) +
				geom_bar(stat = "identity",
								 position = "fill",
								 width = 0.7) +
				scale_y_continuous(labels = scales::percent) +
				xlab("Groups") + ylab("Percentage (%)") + theme_classic() +
				theme(text = element_text(family = "sans", size = 7))
		}
		p
	}
}

if (!is.null(opts$subgroup)){
	p = tax_stackplot1(taxonomy, metadata, topN = opts$legend, groupID = opts$group, subgroup = opts$subgroup, sorted = "abundance")
	if (opts$color == "manual1"){
		p = p + scale_fill_manual(values = colorset1)
	} else if (opts$color == "Paired"){
		p = p + scale_fill_brewer(palette = "Paired")
	}else if (opts$color == "Set3"){
		p = p + scale_fill_brewer(palette = "Set3")
	}
	ggsave(paste0(opts$output,".subgroup.pdf"), p, width = opts$width, height = opts$height, units = "mm")
}else{
	p = tax_stackplot1(taxonomy, metadata, topN = opts$legend, groupID = opts$group, style = "sample", sorted = "abundance",groupOrder = opts$groupOrder)
	
	if (opts$color == "manual1"){
		p = p + scale_fill_manual(values = colorset1)
	} else if (opts$color == "Paired"){
		p = p + scale_fill_brewer(palette = "Paired")
	}else if (opts$color == "Set3"){
		p = p + scale_fill_brewer(palette = "Set3")
	}
	
	#---3.2 保存 Saving#----
	# 大家可以修改图片名称和位置，长宽单位为毫米
	ggsave(paste0(opts$output,".sample.pdf"), p, width = opts$width, height = opts$height, units = "mm")
	
	#----3.3 绘图多层次分组 Plotting#----
	
	
	
	#----3.3 绘图分组 Plotting#----
	# 输入矩阵矩阵、元数据和分组列，返回ggplot2对象
	p = tax_stackplot1(taxonomy, metadata, topN = opts$legend, groupID = opts$group, style = "group", sorted = "abundance")
	if (opts$color == "manual1"){
		p = p + scale_fill_manual(values = colorset1)
	} else if (opts$color == "Paired"){
		p = p + scale_fill_brewer(palette = "Paired")
	}else if (opts$color == "Set3"){
		p = p + scale_fill_brewer(palette = "Set3")
	}
	
	if (!is.null(opts$groupOrder)){
		library(stringr)
		group_order = as.vector(str_split(opts$groupOrder,",",simplify = T))
		p = p + scale_x_discrete(limits = group_order)
	}
	#---3.4 保存 Saving#----
	# 大家可以修改图片名称和位置，长宽单位为毫米
	ggsave(paste0(opts$output,".group.pdf"), p, width = opts$width, height = opts$height, units = "mm")
	
}
