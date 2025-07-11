#!/usr/bin/env Rscript

# Copyright 2016-2021 Yong-Xin Liu <yxliu@genetics.ac.cn / metagenome@126.com>

# If used this script, please cited:
# Yong-Xin Liu, Yuan Qin, Tong Chen, Meiping Lu, Xubo Qian, Xiaoxuan Guo, Yang Bai. A practical guide to amplicon and metagenomic analysis of microbiome data. Protein Cell 2021(12) 5:315-330 doi: 10.1007/s13238-020-00724-8
# 
# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录


# 1.1 程序功能描述和主要步骤

# 程序功能：根据sintax注释结果过滤16S扩增子，去除线粒体、叶绿体和非细菌序列
# Functions: Remove chloroplast, mitocondria, and non-bacteria

options(warn = -1) # Turn off warning
## 设置输入输出文件和参数
# 修改下面`default=`后面的文件和参数。
# 输入文件"-i", "--input"，输入文件OTU表
# 物种注释"-t", "--taxonomy"，sintax物种注释文件位置；
# 输出文件"-o", "--output"，输出文件筛选OTU表
# 分组列名"-s", "--stat"，筛选物种注释和统计文件；
# 过滤物种注释"-d", "--discard"，筛选去除的物种。

# 1.2 解析命令行
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
        make_option(c("-i", "--input"), type="character", default="result/raw/otutab_CQ.txt",
                    help="OTU table [default %default]"),
        make_option(c("-t", "--taxonomy"), type="character", default="result/raw/otus_CQ.sintax",
                    help="sintax taxonomy [default %default]"),
        make_option(c("-s", "--stat"), type="character", default="result/raw/otutab_nonType.stat",
                    help="Filter stat result [default %default]"),
        make_option(c("-d", "--discard"), type="character", default="result/raw/otus.sintax.discard",
                    help="Filter stat result [default %default]"),
        make_option(c("-o", "--output"), type="character", default="result/otutab.txt",
                    help="Filtered OTU table [default %default]")
    )
    opts = parse_args(OptionParser(option_list=option_list))
#    suppressWarnings(dir.create(opts$output))
}


# # 设置输入文件OTU表和物种注释文件位置
# opts$input = "result/raw/otutab.txt"
# opts$taxonomy = "result/raw/otus.sintax"
# 
# # 设置输出文件筛选OTU表、及筛选物种注释和统计文件位置
# opts$output = "result/otutab.txt"
# opts$stat = "result/raw/otutab_nonBac.txt"

# 读取OTU表和sintax物种注释
otutab = read.table(opts$input, header=T, row.names=1, sep="\t", comment.char="")
sintax = read.table(opts$taxonomy, header=F, row.names=1, sep="\t", fill = TRUE, comment.char="")

print(paste0("Input feature table is ", opts$input))
print(paste0("Input sintax taxonomy table is ", opts$taxonomy))

# OTU表汇总
total_reads = colSums(otutab)

# 筛选无注释的信息
nonspecific = sintax[sintax$V4=="",]
nonspecific_reads = colSums(otutab[rownames(nonspecific),,drop=F])
sintax = sintax[sintax$V4!="",]

# 筛选OTU并按丰度由大到小排序
# 如果直接使用sintax的行名，如果与outtab数量不对等，如过多将出现NA行错误，改为交叉筛选
# otutab = otutab[rownames(sintax),]
idx = rownames(otutab) %in% rownames(sintax)
otutab = otutab[idx,,drop=F]
idx = order(rowSums(otutab), decreasing = T)
otutab = otutab[idx,,drop=F]
filtered_reads = colSums(otutab)

# 输出异常OTU编号和总量统计
#write.table(nonspecific, file=opts$discard, append = F, sep="\t", quote=F, row.names=T, col.names=F, eol = "\n")
df = as.data.frame(cbind(total_reads, nonspecific_reads,filtered_reads))
df = rbind(df,df/total_reads*100)
print(paste0("Fileter porpotion: "))
print(df)
suppressWarnings(write.table(paste("SampleID\t",  sep=""), file=opts$stat, append = F, sep="\t", quote=F, row.names=F, col.names=F, eol = ""))
suppressWarnings(write.table(df, file=paste(opts$stat, sep=""), append = T, sep="\t", quote=F, row.names=T, col.names=T, eol = "\n"))

# 筛选排序后的OTU表
write.table(paste("#OTUID\t",  sep=""), file=paste(opts$output, sep=""), append = F, sep="\t", quote=F, row.names=F, col.names=F, eol = "")
suppressWarnings(write.table(otutab, file=paste(opts$output, sep=""), append = T, sep="\t", quote=F, row.names=T, col.names=T, eol = "\n"))

# 输出样本基本信息
#print(paste0("Summary of samples size in final feature table: "))
#summary(filtered_reads)

print(paste0("Onput feature table is ", opts$output))
print(paste0("Detail and statistics in ", opts$stat))

