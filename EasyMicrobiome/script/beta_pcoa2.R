if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
	install.packages(p, repos=site)
	require("optparse",character.only=T)
}
# 解析参数-h显示帮助信息
if (TRUE){
	option_list = list(
		make_option(c("-i", "--input"), type="character", default="result/beta/bray_curtis.txt",
								help="Beta distance matrix [default %default]"),
		make_option(c("-d", "--design"), type="character", default="result/metadata.txt",
								help="Design file or metadata [default %default]"),
		make_option(c("-n", "--group"), type="character", default="Group",
								help="Group name [default %default]"),
		make_option(c("-l", "--label"), type="logical", default=FALSE,
								help="Design file or metadata [default %default]"),
		make_option(c("-o", "--output"), type="character", default="",
								help="Output directory; name according to input [default %default]"),
		make_option(c("-w", "--width"), type="numeric", default=89,
								help="Figure width in mm [default %default]"),
		make_option(c("-e", "--height"), type="numeric", default=59,
								help="Figure heidth in mm [default %default]")
	)
	opts = parse_args(OptionParser(option_list=option_list))
	# suppressWarnings(dir.create(opts$output))
}
# 设置输出文件缺省值，如果为空，则为输入+pcoa.pdf
if(opts$output==""){opts$output=paste0(opts$input,".pcoa.boxplot.pdf")}
opts$output = gsub(".txt","",opts$output)
suppressWarnings(dir.create(dirname(opts$output), showWarnings = F))


#----1.3. 加载包 Load packages#----

suppressWarnings(suppressMessages(library(amplicon)))
library(ggplot2)
library(ggforce)
library(ggrepel)
library(dplyr)
#----2. 读取文件 Read files#----

#----2.1 实验设计 Metadata#----
metadata = read.table(opts$design, header=T, row.names=1, sep="\t", comment.char="",
											stringsAsFactors = F,check.names = F)
dis_mat = read.csv(opts$input,sep = "\t",header = T,
									 row.names = 1,check.names = F,comment.char = "")

dis_mat=dis_mat[rownames(metadata), rownames(metadata)]

# 提取样品组信息,默认为group可指定
groupID = opts$group
sampFile=as.data.frame(metadata[, groupID],row.names=row.names(metadata))

# PCoA
pcoa=cmdscale(dis_mat, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points=as.data.frame(pcoa$points) # get coordinate string, format to dataframme
eig=pcoa$eig
points=cbind(points, sampFile[rownames(points),])
colnames(points)=c("x", "y", "z","group")
points$group = as.factor(points$group)

library(vegan)
permanova_result <- adonis(dis_mat ~ oldGroup, data = metadata)
permanova_result$aov.tab$Df[1]

text1 = paste0("PERMANOVA:\ndf = ",permanova_result$aov.tab$Df[1],
							"\nR2 = ",round(permanova_result$aov.tab$R2[1],3),
							"\nP-value = ",permanova_result$aov.tab$`Pr(>F)`[1])
text1

empty_plot <- ggplot() + 
	annotate("text", x = 0, y = 0, label = text1,
					 size = 4, fontface = "bold") + 
	theme(
		# 移除网格线
		panel.grid = element_blank(),
		# 设置图形区域为白色
		panel.background = element_rect(fill = "white", color = "black"),
		# 移除坐标轴
		axis.line = element_blank(),
		axis.text = element_blank(),
		axis.ticks = element_blank(),
		axis.title = element_blank()
	) 


anova_result <- aov(y ~ group, data = points)
# 使用Tukey HSD进行多重比较
tukey_result <- TukeyHSD(anova_result)

# 提取比较结果
Tukey_HSD_table = as.data.frame(tukey_result$group)
Tukey.levels = tukey_result$group[,4]
library(multcompView)
Tukey.labels = data.frame(multcompLetters(Tukey.levels)['Letters'])

# 设置分组位置为各组y最大值+高的5%
y2 = as.data.frame(points %>% group_by(group) %>% summarise(Max=max(y)))
rownames(y2) = y2$group
y2$label = Tukey.labels[rownames(y2),"Letters"]
y2 = y2[levels(points$group),]

p3 = ggplot(points, aes(x = group, y = y)) +
	geom_boxplot(aes(color = group), size = 1, width = 0.5,alpha=.7) +
	labs(x=NULL,y=NULL) +
	theme(
		# 移除网格线
		panel.grid = element_blank(),
		# 设置图形区域为白色
		panel.background = element_rect(fill = "white", color = "black"),
		# 移除坐标轴
		axis.line = element_blank(),
		axis.text = element_blank(),
		axis.ticks = element_blank(),
		axis.title = element_blank()
	) +
	theme(legend.position = "none") +
	geom_text(data = y2,
						aes(x = group, 
								y = Max+0.05, label = label), size = 5)

p3

# 添加显著性结果
anova_result <- aov(x ~ group, data = points)
# 使用Tukey HSD进行多重比较
tukey_result <- TukeyHSD(anova_result)

# 提取比较结果
Tukey_HSD_table = as.data.frame(tukey_result$group)
Tukey.levels = tukey_result$group[,4]
Tukey.labels = data.frame(multcompLetters(Tukey.levels)['Letters'])

# 设置分组位置为各组y最大值+高的5%
y1 = as.data.frame(points %>% group_by(group) %>% summarise(Max=max(x)))
rownames(y1) = y1$group
y1$label = Tukey.labels[rownames(y1),"Letters"]
y1 = y1[levels(points$group),]

p2 = ggplot(points, aes(x = group, y = x)) +
	geom_boxplot(aes(color=group), size = 1, width = 0.5) +
	labs(x=NULL,y = NULL) +
	theme(
		# 移除网格线
		panel.grid = element_blank(),
		# 设置图形区域为白色
		panel.background = element_rect(fill = "white", color = "black"),
		# 移除坐标轴
		axis.line = element_blank(),
		axis.text = element_blank(),
		axis.ticks = element_blank(),
		axis.title = element_blank()
	) +
	coord_flip()+
	theme(legend.position = "none") +
	geom_text(data = y1,
						aes(x = group, angle = 270,
								y = Max+0.05, label = label), size = 5)
p2


p1 = ggplot(points, aes(x=x, y=y,group = group,color=group))  +
	labs(x=paste("PCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
			 y=paste("PCo 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
	geom_point(aes(color=group),alpha=.7) + 
	theme_bw() + 
	theme(legend.position = "bottom",
				panel.grid = element_blank()) +
	guides(color = guide_legend(order = 2),
				 shape = guide_legend(order = 3),
				 size = guide_legend(order = 1)
	)+ stat_ellipse(level=0.68) 
p1

library(patchwork)

## 出图-----------
p4 = p2 + empty_plot +p1 +p3 +plot_layout(ncol = 2,
																		 heights = c(1, 3),
																		 widths  = c(3, 1))
	
p4
ggsave(opts$output,p4,width = opts$width,height = opts$height,
			 units = "mm")




