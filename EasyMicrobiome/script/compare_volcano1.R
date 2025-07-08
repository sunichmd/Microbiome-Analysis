#!/usr/bin/env Rscript

# Copyright 2016-2020 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
# Jingying Zhang, Yong-Xin Liu, et. al. NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice. Nature Biotechnology 37, 676-684, doi:10.1038/s41587-019-0104-4 (2019).

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录

#----1. 参数 Parameters#----

#----1.1 功能描述 Function description#----

# 程序功能：差异比较火山图绘制
# Functions: Volcano plot

options(warn = -1) # Turn off warning


## 设置输入输出文件和参数

# 修改下面`default=`后面的文件和参数。
#
# 输入文件为特征表(A-B.txt)+分组信息(metadata.tsv)
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
# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages("optparse", repos=site)
  require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="result/compare/KO-WT.txt",
                help="Feature table [default %default]"),
    make_option(c("-t", "--threshold"), type="numeric", default=0.1,
                help="Threshold of relative abundance 0.1% [default %default]"),
    make_option(c("-d", "--design"), type="character", default="result/metadata.tsv",
                help="Design file or metadata [default %default]"),
    make_option(c("-n", "--group"), type="character", default="Group",
                help="Group name [default %default]"),
    make_option(c("-c", "--compare"), type="character", default="KO-WT",
                help="Groups comparison [default %default]"),
    make_option(c("-m", "--method"), type="character", default="wilcox",
                help="Compare method, alternative edgeR or t.test [default %default]"),
    make_option(c("-p", "--pvalue"), type="numeric", default=0.05,
                help="Threshold of P-value [default %default]"),
    make_option(c("-f", "--fdr"), type="numeric", default=0.1,
                help="Threshold of FDR [default %default]"),
    make_option(c("-l", "--log2foldChange"), type="numeric", default=0,
    						help="Threshold of log2 fold change [default %default]"),
    make_option(c("-o", "--output"), type="character", default="",
                help="Output directory; name according to input [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=89,
                help="Figure width [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=59,
                help="Figure heidth [default %default]")
    )
  opts = parse_args(OptionParser(option_list=option_list))
  # suppressWarnings(dir.create(opts$output))
}
# 设置输出文件缺省值，如果为空，则为输入+pcoa.pdf
if(opts$output == ""){
  opts$output = paste0(opts$input, ".volcano.pdf")}


#----1.3. 加载包 Load packages#----

suppressWarnings(suppressMessages(library(amplicon)))
library(stringr)
# 5.23修改 加载包
package_list <- c('cowplot')
# 判断R包加载是否成功来决定是否安装后再加载
for(pac in package_list){
  if(!suppressWarnings(suppressMessages(require(pac, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(pac, repos=site)
    suppressWarnings(suppressMessages(library(pac, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

#----2. 读取文件 Read files#----

#----2.1 读取差异比较结果#----
diff = read.table(opts$input, header=T, row.names=1, sep="\t", comment.char="")
# diffN = as.data.frame(table(diff$level))
# 5.23修改 对区间赋值
if(F){
	diff$log2FC[diff$log2FC> 4]=4
	diff$log2FC[diff$log2FC< -4]=-4
	breaks= c(-4,-3,-2,-1,0,1,2,3,4)
}



#----3. 统计绘图#----

main_theme = theme(panel.background = element_blank(),
                   panel.grid = element_blank(),
                   axis.line.x = element_line(size = .5, colour = "black"),
                   axis.line.y = element_line(size = .5, colour = "black"),
                   axis.ticks = element_line(size = .5, color = "black"),
                   axis.text = element_text(size = 7, color = "black"),
                   legend.position = "top",
                   legend.background = element_blank(),
                   legend.key = element_blank(),
                   legend.text = element_text(size = 7),
                   text = element_text(family = "sans", size = 7))
table(diff$level)
# diff$color = ifelse(diff$PValue < 0.05,
# 										ifelse(diff$log2FC > opts$log2foldChange, 
# 													 "red",
# 													 ifelse(diff$log2FC < -opts$log2foldChange, 
# 													 			 "blue",
# 													 			 "grey")),
# 										"grey")
diff$color = ifelse(diff$level == "Enriched","red",
										ifelse(diff$level == "Depleted",
													 "green",
													 "grey"))
diff$color = factor(diff$color,
											 levels = c('red','green','grey'))
table(diff$color)
# 创建火山图
sampleA = str_split(opts$compare,"-",simplify = T)[1]
sampleB = str_split(opts$compare,"-",simplify = T)[2]
sampleA_num = ifelse(is.na(table(diff$level)["Enriched"]),0,table(diff$level)["Enriched"])
sampleB_num = ifelse(is.na(table(diff$level)["Depleted"]),0,table(diff$level)["Depleted"])

if(length(unique(diff$color))==1){
	point_col = "grey"
	diff$gene_num = factor(diff$color,
												 levels = point_col,
												 labels = "NotSig")
}else if(length(unique(diff$color))==2){
	if("red" %in% unique(diff$color)){
		point_col = c("red","grey")
		diff$gene_num = factor(diff$color,
													 levels = point_col,
													 labels = c(paste(sampleA,sampleA_num),"NotSig"))
	}else{
		point_col = c("green","grey")
		diff$gene_num = factor(diff$color,
													 levels = point_col,
													 labels = c(paste(sampleB,sampleB_num),"NotSig"))
	}
}else{
	point_col = c('red',"green","grey")
	diff$gene_num = factor(diff$color,
												 levels = point_col,
												 labels = c(paste(sampleA,sampleA_num),paste(sampleB,sampleB_num),"NotSig"))
}

p = ggplot(diff, aes(x = log2FC, y = -log10(PValue),color=gene_num)) +
	geom_point(size = 1.5) +
	scale_color_manual(values = point_col) +
	guides(colour=guide_legend(title=NULL)) +
	geom_point(data = subset(diff, level == "Enriched"), color = "red",size = 1.5) +
	geom_point(data = subset(diff, level == "Depleted"), color = "green",size = 1.5) +
	geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
	labs(x = "log2(Fold Change)", y = "-log10(p-value)", title = paste(sampleA,"vs",sampleB))+
	main_theme
p
# p =ggplot(diff, aes(x=log2FC, y=PValue, color=level)) +
#   geom_point(color = level)+
#   theme_classic()+
#   #加抖动
#   geom_jitter(data = subset(diff, log2FC > opts$log2foldChange & PValue < 0.05), width = 0.2, height = 0.2, color = "red")+
#   geom_jitter(data = subset(diff, log2FC < -opts$log2foldChange & PValue < 0.05), width = 0.2, height = 0.2, color = "green")+
#   # 5.23修改 重设刻度间隔
#   #scale_x_continuous(breaks=seq(-4,4,1))+ 
#   theme_classic()+
#   labs(x="log2(fold change)", y="log2(count per million)",
#        title=gsub(".txt","",basename(opts$input)))+
#   annotate("text",x=-3,y=15,label=paste(names(table(diff$level))[1],table(diff$level)[1]))+
#   annotate("text",x=3,y=15,paste(names(table(diff$level))[2],table(diff$level)[2])) +
#   main_theme
# # 5.23修改 加图例
# legend= ggplot(diff, aes(x=log2FC, y=log2CPM, color=level)) +
#   geom_point() + theme_classic()+
#   scale_colour_manual(values=c("green","red","grey")) +
#   main_theme
# 
# p1=plot_grid(p,
#              plot_grid(get_legend(legend),
#                        ncol = 1,
#                        align = "hv"),
#              ncol = 2,align = "hv",rel_widths=c(9,1))
p
ggsave(opts$output, p, width = opts$width, height = opts$height, units = "mm")

