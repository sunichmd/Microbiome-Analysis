# Microbiome-Analysis

本项目整理并优化了多个微生物组分析流程，主要包括：

## EasyAmplicon

本流程基于 [YongxinLiu/EasyAmplicon](https://github.com/YongxinLiu/EasyAmplicon) 仓库（版本：2025年4月）进行修改与优化。原项目使用 GPL v3 协议发布，本项目遵循相同协议，并在其基础上加入了更丰富的绘图与结果输出模块。

### Qiime2 流程修改
- qiime2数据库文件更新
 - 添加数据库：符合2024.10版qiime2的silva数据库
 - 添加数据库：greenggene2_24_09数据库
 - 添加数据库：口腔微生物数据库HOMD v15.23数据库
- pipeline_qiime2.sh 进行优化为pipeline_qiime2_4ana.sh
 - 添加设置：export UNIFRAC_MAX_CPU=basic  #CPU不支持 AVX2的解决方法
 - 更新代码：物种注释数据库训练代码
 - 添加代码：整理qiime2结果变为常见人类可读文件格式
 - 添加代码：基于R对qiime2结果进行可视化代码
 - 添加代码：整理数据为lefse输入数据格式；实现lefse分析
 - 添加代码：qiime2 longitudinal对alpha多样性、beta多样性和物种丰度的分析
 - 添加代码：alpha.txt beta.txt 转为qza格式

### Usearch+Vsearch 流程（window系统和Linux系统）
- pipeline.sh 进行优化为pipeline_4ana.sh
 - 增加可选的数据库：silva-138-99-seqs.usearch.fa greengene2_2022.10.seqs.usearch.fna
 - 使用${db}/script/alpha_boxplot1.R；调整alpha统计结果输出以及支持修改图中分组顺序
 - 使用${db}/script/beta_pcoa2.R：输出pcoa+boxplot 展示分组是否在不同pc上存在差异
 - 使用${db}/script/tax_stackplot1.R：输出画图所用数据，支持修改图中分组顺序以及支持输入第二分组对图形进行分面绘制
 - 使用${db}/script/tax_circlize_xlim90.R：图形X轴文字旋转90度另外可指定输出结果位置
 - 使用${db}/script/compare_volcano1.R: 解决颜色容易错误分配以及部分结果展示不全的问题
 - 使用${db}/script/compare_heatmap1.sh：设置输入行的注释信息非必须
 - 附录增加2.1：可以解决样本序列文件为fasta文件的情况

### Vsearch流程（Mac系统）
- 更新了pipeline.sh
 - 修改cat -A为cat -vet 兼容Mac系统
 - 添加R计算alpha稀释曲线
 - 添加R计算beta多样性结果

### 通用分析
- 添加advanced/BatchCorrect文件夹：
 - 评估批次影响
 - 4种批次矫正的方法：MMUPHin、removeBatchEffect、ComBat、Percentile Normalisation
- 修改advanced/SourceTrackerFeast文件夹：
 - 增加SourceTracker1.r修正SourceTracker.r报错（由于一些R包更新数据类型的报错）
- 修改advanced/SourceTrackerFeast文件夹：
 - 增加FEAST1.Rmd 时间个体样本一一对应溯源分析
- 添加advanced/Multi-factorAnalysis文件夹：
 - 分析多因素情况下alpha多样性变化
 - 分析多因素情况下beta多样性变化
 - 分析多因素情况下物种丰度变化
 - 对于时间+分组分类可考虑qiime的longitudinal分析
- 修改advanced/RandomForestClassification文件夹：
 - 增加RF_classification_onegroup.Rmd调整设置随机抽样数据作为测试集和验证集，增加AUC曲线；绘制所有样本热图
- 添加advanced/Net_SingleMulti_Vis文件夹：
 - 添加One_cmd_calculate_net_onefile_twofile.Rmd，实现一行代码计算单文件或两个文件相关性网络并可视化
  - 计算每组相关性网络
  - 可视化网络
  - 计算zipi
  - 计算网络属性+可视化网络边数
  - 可进行网络模块分析
  - 网络稳定性评估
  - 基于qgraph包实现类graph的布局
  - 整理结果可导入graphi可视化
  - 整理结果可导入cytoscape可视化
  - 可用R调用cytoscape可视化
  - 统计核心物种
 - sparcc_vis_as_cellParper.Rmd，基于sparcc计算网络并且在R进行可视化（复现cell文章）
 - cor_2file.Rmd 计算简单相关性网络；生成可交互结果
- 添加Phytopathogen_Bac_predict文件夹：实现细菌的致病性功能预测以及差异分析
- 添加advanced/calulate_core_otu_by_usearch_treeplot.Rmd：基于Usearch计算核心物种并绘制进化树
- 添加advanced/adonis_lost_variable_solution.Rmd：解决adnois分析中无法评估一些变量影响
- 添加advanced/calulate_each_tax_raw_count.Rmd：计算物种水平数据的alpha多样性和beta多样性
- 添加advanced/PCoA_density_compare_multi-factor.Rmd：绘制PCOA图+密度曲线以及比较不同分组在PCo上是否存在差异
- 添加advanced/Fungi_function_predict.Rmd：真菌功能注释
- 添加advanced/calculate_Bray-Curtis_dissimilarities.Rmd：计算Bray-Curtis dissimilarities
- 添加advanced/Other_ana.Rmd（不常用分析）
 - adonism分析
 - R实现lefse分析
 - 绘制三元图（只适合三组数据）
 - NMDS展示beta多样性差异
 - 整理greengene2数据库为usearch注释可用
 - 计算相对丰度的alpha多样性
 - 计算物种累积曲线评估样本个数是否达到饱和
 - 物种相对丰度冲击图
 - svd分析
## EasyMetagenome

本流程基于 [YongxinLiu/EasyMetagenome](https://github.com/YongxinLiu/EasyMetagenome) 仓库（版本：2024年12月）进行修改优化，主要针对文件结构与兼容性进行了整理。

## EasyMicrobiome

本流程基于 [YongxinLiu/EasyMicrobiome](https://github.com/YongxinLiu/EasyMicrobiome) 仓库（版本：2024年12月）修改，并重构部分参数与默认模块以适配批量运行需求。

更新内容：
- usearch文件夹增加 greengene2数据库（基于2022.10月更新的数据库整理）、和silva-138数据库（基于qiime2整理好的数据库进行整理）
- script文件夹
 - 添加alpha_calculate_pd.R；R计算7种alpha多样性【pd多样性计算需提供进化树】

> 本项目所有衍生流程保留原始作者署名与协议，并明确标注改动内容，旨在促进微生物组数据处理的复用性与模块化。

This project is a derivative of EasyAmplicon, EasyMetagenome, and EasyMicrobiome by Yongxin Liu, originally distributed under the GNU General Public License v3.0.
All modifications to the original code are also distributed under the same license.
=======
>>>>>>> 1114af2cc93995c9647546431d51a5ac3c660e58
