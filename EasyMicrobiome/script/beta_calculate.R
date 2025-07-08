#!/usr/bin/env Rscript

# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
	install.packages("optparse", repos=site)
	require("optparse",character.only=T) 
}

# 解析参数-h显示帮助信息
option_list <- list(
	make_option(c("-i", "--input"), type="character", default="result/otutab_rare.txt",
							help="Input reads count rare file; such as OTU table, kraken2 taxonomy counts table [default %default]"),
	make_option(c("-b", "--beta_index"), type="character", default="all",
							help="Input beta matric if input all will calculate bray_curtis jaccard euclidean and unifrac [default %default]"),
	make_option(c("-t", "--tree"), type="character", default=NULL,
							help="Input tree filename [default %default]"),
	make_option(c("-o", "--output"), type="character", default="result/beta",
							help="Output beta diversity prefix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

# 显示输入输出确认是否正确
print(paste("The input feature table is ", opts$input,  sep = ""))
print(paste("Output prefix beta diversity: ", opts$output, sep = ""))

# 1.3 安装CRAN来源常用包
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list <- c("tidyverse","plyr","ape")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install.packages(p, repos=site)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
	}
}

if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager", repos = site)

a = rownames(installed.packages())

install_bioc <- c("phyloseq")

for (i in install_bioc) {
	if (!i %in% a)
		BiocManager::install(i, update = F, site_repository=site)
}
suppressWarnings(suppressMessages(library(phyloseq)))

# 1.4 读取输入文件

# 默认的quote会跳过2/3的数据行减少，产生NA，改为空
otu <- read.table(opts$input, header=T, sep="\t", quote = "", row.names=1, comment.char="")
if(!is.null(opts$tree)){
	phy <- read.tree(opts$tree)
	GP <- suppressWarnings(phyloseq(otu_table(as.matrix(otu), taxa_are_rows=TRUE),
																	phy_tree(phy)))
}else{
	GP <- suppressWarnings(phyloseq(otu_table(as.matrix(otu),
																						taxa_are_rows=TRUE)))
}




#unifrac、wunifrac、jaccard、bray_curtis、
if (opts$beta_index == "all"){
	dist_beta_metric = c("bray_curtis","unifrac","jaccard","euclidean")
}else{
	dist_beta_metric = opts$beta_index
}

for(i in dist_beta_metric){
	if(i == "bray_curtis") {name = i;i="bray"}else{name = i}
	if(i == "unifrac"){
		suppressWarnings(beta_dist <- as.matrix(phyloseq::distance(GP,method = 'unifrac')))
	}else{
		suppressWarnings(beta_dist <- as.matrix(phyloseq::distance(GP,method = i,binary = T)))
	}
	filename = paste0(opts$output,"/",name,"_binary.txt")
	write.table(paste0(name,"_binary\t"), file=filename,
							append = F, quote = F, eol = "", row.names = F, col.names = F)
	# 保存统计结果，有waring正常
	suppressWarnings(write.table(beta_dist, file=filename, append = T, quote = F,
															 sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))
	
	if(i == "unifrac"){
		beta_dist <- suppressWarnings(as.matrix(phyloseq::distance(GP,method = 'wunifrac')))
	}else{
		beta_dist <- suppressWarnings(as.matrix(phyloseq::distance(GP,method = i,binary = F)))
	}
	filename = paste0(opts$output,"/",name,".txt")
	write.table(paste0(name,"\t"), file=filename,
							append = F, quote = F, eol = "", row.names = F, col.names = F)
	# 保存统计结果，有waring正常
	suppressWarnings(write.table(beta_dist, file=filename, append = T, quote = F,
															 sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))
	
}

# vegan包可以计算除unifrac的距离
# library(vegan)
# otu = read.table(opts$input, header=T, sep="\t", quote = "", row.names=1, comment.char="") 
# dis_beta <- as.matrix(vegdist(t(otu),method = "bray"))



