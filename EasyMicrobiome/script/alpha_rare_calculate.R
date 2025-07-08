#!/usr/bin/env Rscript

# 设置清华源加速下载
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 判断命令行解析是否安装，安装并加载
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
	install.packages("optparse", repos=site)
	require("optparse",character.only=T) 
}

package_list <- c("vegan","ape")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
	if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
		install.packages(p, repos=site)
		suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
	}
}
# 解析参数-h显示帮助信息
option_list <- list(
	make_option(c("-i", "--input"), type="character", default="result/otutab_rare.txt",
							help="Input alpha diversity file [default %default]"),
	make_option(c("-s", "--seed"), type="numeric", default=1,
							help="Random sample seed, can set any integer [default %default]"),
	make_option(c("-m", "--sample_method"), type="character", default="freq",
							help="Input the method of sample freq or seq [default %default]"),
	make_option(c("-a", "--alpha_index"), type="character", default="richness",
							help="Input alpha matric [default %default]"),
	make_option(c("-o", "--output"), type="character", default="result/alpha/alpha_rare.txt",
							help="Output file [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

# 显示输入输出确认是否正确
print(paste("The input file is ", opts$input,  sep = ""))
print(paste("Output file is: ", opts$output, sep = ""))

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

alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
	if (method == 'richness') result <- rowSums(x > 0)    #丰富度指数
	else if (method == 'chao1') result <- estimateR(x)[3, ]    #Chao1 指数
	else if (method == 'ace') result <- estimateR(x)[5, ]    #ACE 指数
	else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)    #Shannon 指数
	else if (method == 'simpson') result <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
	else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)    #Pielou 均匀度
	else if (method == 'gc') result <- 1 - rowSums(x == 1) / rowSums(x)    #goods_coverage
	else if (method == 'pd' & !is.null(tree)) {    #PD_whole_tree
		pd <- pd(x, tree, include.root = FALSE)
		result <- pd[ ,1]
		names(result) <- rownames(pd)
	}
	result
}

# 可以rare设置固定总序列数来抽样
alpha_curves <- function(x, step, method = 'richness',sample_method = "freq",
												 rare = NULL, tree = NULL, base = exp(1)) {
	x_nrow <- nrow(x)
	if (is.null(rare)) rare <- rowSums(x) else rare <- rep(rare, x_nrow)
	alpha_rare <- list()
	
	if(sample_method == "seq"){
		for (i in 1:x_nrow) {
			step_num <- seq(0, rare[i], step)
			if (max(step_num) < rare[i]) step_num <- c(step_num, rare[i])
			
			alpha_rare_i <- NULL
			set.seed(opts$seed);sub_data = rrarefy(x[i, ], step_num_n)
			for (step_num_n in step_num) alpha_rare_i <- c(alpha_rare_i, alpha_index(x = sub_data,
																																							 method = method, tree = tree, base = base))
			names(alpha_rare_i) <- step_num
			alpha_rare <- c(alpha_rare, list(alpha_rare_i))
		}
		names(alpha_rare) <- rownames(x)
		richness_curves <- alpha_rare
		plot_richness <- c()
		for (i in names(richness_curves)) {
			richness_curves_i <- (richness_curves[[i]])
			richness_curves_all <- data.frame(rare = names(richness_curves_i), 
																				alpha = richness_curves_i,
																				sample = i,
																				stringsAsFactors = FALSE)
			plot_richness <- rbind(plot_richness, richness_curves_all)
		}
	}else if(sample_method == "freq"){
		for (i in 1:x_nrow) {
			# 生成1%-100%要抽样的序列数，小数四舍五入处理
			step_num <- round(quantile(1:rare[i],probs = seq(0, 1, 0.01)))[-1]
			
			alpha_rare_i <- NULL
			# 计算不同抽样下的序列算出的alpha多样性
			for (step_num_n in step_num) {
				set.seed(opts$seed);sub_data = rrarefy(x[i, ], step_num_n)
				alpha_rare_i <- c(alpha_rare_i, alpha_index(sub_data,
																										method = method, tree = tree, base = base))
			}
			alpha_rare <- c(alpha_rare, list(alpha_rare_i))
		}
		names(alpha_rare) <- rownames(x)
		richness_curves <- alpha_rare

		plot_richness <- data.frame(row.names = 1:100)
		for (i in names(richness_curves)) {
			richness_curves_i <- data.frame(richness_curves[[i]])
			plot_richness <- cbind(plot_richness, richness_curves_i)
		}
		colnames(plot_richness) = rownames(x)
	}
	plot_richness
}
# 每次加2000reads抽
#richness_curves <- alpha_curves(t(otu), step = 2000, method = "richness")
# 按照序列1%-100%抽
rare_alpha <- alpha_curves(t(otu),method = opts$alpha_index,
																sample_method = opts$sample_method)

write.table(paste0(opts$alpha_index,"\t"), file=opts$output,
						append = F, quote = F, eol = "", row.names = F, col.names = F)
# 保存统计结果，有waring正常
suppressWarnings(write.table(rare_alpha, file=opts$output, append = T, quote = F,
														 sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))