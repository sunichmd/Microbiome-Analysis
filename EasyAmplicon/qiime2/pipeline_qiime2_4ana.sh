[TOC]

# QIIME2 2024.10标准分析流程

## 0. 软件安装(附录1)

    # 在Linux、Mac、Windows内置Linux子系统(支持右键粘贴)下安装并运行
    # 详细教程参阅官网https://docs.qiime2.org/2024.2/
    # 安装Windows子系统，https://mp.weixin.qq.com/s/0PfA0bqdvrEbo62zPVq4kQ

## 1. 准备工作
    # ubuntu 黑白的界面 -- qiime2流程
    # git for windows 彩色的界面 -- usearch流程
    # 可以在tools 的global option 的termina选择对应的界面；新建的就会是对应的页面
    
    # 设置工作目录，如服务器为~/amplicon/qiime2，Win子系统如下：
    wd=/mnt/c/0amplicon/EasyAmplicon/qiime2/
    # 进入工作目录
    mkdir -p ${wd}
    cd ${wd}
    # 激活QIIME2工作环境，旧版conda使用source替换conda运行
    conda activate qiime2-amplicon-2024.10
    
    # 准备样本元数据metadata.txt、原始数据seq/*.fq.gz
    
    ## 此处代码基于metadata.txt从公共数据下载测序数据，按GSA的CRA(批次)和CRR(样品)编号下载数据
    mkdir -p seq
    # # 公网下载
    # cd seq
    # awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_f1.fq.gz -O "$1"_1.fq.gz")}' \
    #     <(tail -n+2 ../metadata.txt)
    #  awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_r2.fq.gz -O "$1"_2.fq.gz")}' \
    #     <(tail -n+2 ../metadata.txt)
    # cd .. && ls -lsh seq
    # 从其他地方链接(不额外占用空间)
    ln /mnt/c/0amplicon/EasyAmplicon/seq/* seq/
    ln /mnt/c/0amplicon/EasyAmplicon/result/metadata.txt ./
    # 根据metadata生成manifest文件,还需要根据自己的情况适当修改
    
    awk 'NR==1{print "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"} \
      NR>1{print $1"\t$PWD/seq/"$1"_1.fq.gz\t$PWD/seq/"$1"_2.fq.gz"}' \
      metadata.txt > manifest
    head -n3 manifest
    
    # 数据导入qiime2，格式为双端33格式
    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path manifest \
      --output-path demux.qza \
      --input-format PairedEndFastqManifestPhred33V2
    # fq文件1G用时7m，fq.gz压缩格式仅34s，测试数据9s

## 2. 生成特征表和代表序列

### 方法1. DADA2(慢，依赖R和Python包多容易报错)

    # 支持多线程加速，90万条PE250数据，0/96p, 34m；24p, 44m；8p, 77m；1p, 462m
    # 27万，8p, 9m; 4p, 11m;
    time qiime dada2 denoise-paired \
      --i-demultiplexed-seqs demux.qza \
      --p-n-threads 4 \
      --p-trim-left-f 29 --p-trim-left-r 18 \
      --p-trunc-len-f 0 --p-trunc-len-r 0 \
      --o-table dada2-table.qza \
      --o-representative-sequences dada2-rep-seqs.qza \
      --o-denoising-stats denoising-stats.qza
    
    # 确定使用dada2结果并导入主流程
    cp dada2-table.qza table.qza
    cp dada2-rep-seqs.qza rep-seqs.qza

### 方法2. 外部导入特征表和代表序列(常用)

    # 上传我们生成的OTU表otutab.txt和代表序列otus.fa
    # 转换文本为Biom1.0，注意biom --version 2.1.5/8可以，2.1.7报错
    biom convert -i otutab.txt -o otutab.biom \
      --table-type="OTU table" --to-json
    # 导入特征表，9s
    qiime tools import --input-path otutab.biom \
      --type 'FeatureTable[Frequency]' --input-format BIOMV100Format \
      --output-path table.qza
    # 导入代表序列，8s
    qiime tools import --input-path otus.fa \
      --type 'FeatureData[Sequence]' \
      --output-path rep-seqs.qza
      
    biom convert -i 7ec3072f-5ebc-4fa9-a857-4dc258e98ce0/data/feature-table.biom -o otutab_qiime2.txt \
     --to-tsv

### 特征表和代表序列统计
    qiime feature-table summarize \
      --i-table table.qza \
      --o-visualization table.qzv \
      --m-sample-metadata-file metadata.txt
    qiime feature-table tabulate-seqs \
      --i-data rep-seqs.qza \
      --o-visualization rep-seqs.qzv
     
    unzip rep-seqs.qzv
    # 下载qzv在线查看，dada2只有1千多个ASV


## 3. Alpha和beta多样性分析

### 构建进化树用于多样性分析 53s

    qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences rep-seqs.qza \
      --o-alignment aligned-rep-seqs.qza \
      --o-masked-alignment masked-aligned-rep-seqs.qza \
      --o-tree unrooted-tree.qza \
      --o-rooted-tree rooted-tree.qza

### 计算核心多样性 

    # 13s，采样深度通常选择最小值，来自table.qzv
    export UNIFRAC_MAX_CPU=basic  #CPU不支持 AVX2的解决方法
    rm -rf core-metrics-results
    qiime diversity core-metrics-phylogenetic \
      --i-phylogeny rooted-tree.qza \
      --i-table table.qza \
      --p-sampling-depth 7439 \
      --m-metadata-file metadata.txt \
      --output-dir core-metrics-results

### Alpha多样性组间显著性分析和可视化

    # 7s, 可选的alpha指数有 faith_pd、shannon、observed_features、evenness
    for index in observed_features faith_pd shannon evenness;do
    qiime diversity alpha-group-significance \
      --i-alpha-diversity core-metrics-results/${index}_vector.qza \
      --m-metadata-file metadata.txt \
      --o-visualization core-metrics-results/${index}-group-significance.qzv
    done
    
    # 解压qzv文件方便后续用r脚本进行可视化
    mkdir -p result/alpha
    for index in observed_features faith_pd shannon evenness;do
    mv core-metrics-results/${index}_vector.qza result/alpha
    unzip -d result/alpha/ result/alpha/${index}_vector.qza
    done
    
    # 合并所有 alpha多样性结果
    awk 'BEGIN{OFS=FS="\t"}ARGIND==1{if(FNR==1) b[ARGIND]="sample";else a[$1]=$1}\
    ARGIND>1{if(FNR==1) b[ARGIND]=$2;else a[$1]=a[$1]"\t"$2}\
    END{print b[1],b[2],b[3],b[4],b[5];for(i in a) print a[i]}' \
    metadata.txt `ls result/alpha/*/data/alpha-diversity.tsv` > result/alpha/alpha.txt
    
    db=/C/EasyMicrobiome
    PATH=$PATH:${db}/win
    for i in `head -n1 result/alpha/alpha.txt|cut -f 2-`;do
      Rscript ${db}/script/alpha_boxplot1.R --alpha_index ${i} \
        --input result/alpha/alpha.txt --design metadata.txt \
        --group Group --output result/alpha/ --groupOrder "AA.0m,AA.6m,AA.12m,FA.0m,FA.6m,FA.12m,FF.0m,FF.6m,FF.12m"\
        --width 150 --height 80
    done
    mv alpha_boxplot_TukeyHSD.txt result/alpha/
    
    # for i in result/alpha/*.pdf;do start $i;done
    
    # rm -r `find result/alpha/ -maxdepth 1 -mindepth 1 -type d`
    # 
    # for index in observed_features faith_pd shannon evenness;do
    # unzip -d result/alpha/${index}  result/alpha/${index}-group-significance.qzv
    # done
    
### Alpha多样性稀疏曲线

    # 25s, max-depth选最大值，来自table.qzv
    qiime diversity alpha-rarefaction \
      --i-table table.qza \
      --i-phylogeny rooted-tree.qza \
      --p-max-depth 10299 \
      --m-metadata-file metadata.txt \
      --o-visualization alpha-rarefaction.qzv
    # 结果有observed_otus, shannon, 和faith_pd三种指数可选
    unzip -d result/alpha/ result/alpha/alpha-rarefaction.qzv 
    
    
### Beta多样性组间显著性分析和可视化

    # 可选的beta指数有 unweighted_unifrac、bray_curtis、weighted_unifrac和jaccard
    # 7s, 指定分组是减少计算量，置换检验较耗时
    for distance in weighted_unifrac unweighted_unifrac bray_curtis jaccard;do
    column=Group
    qiime diversity beta-group-significance \
      --i-distance-matrix core-metrics-results/${distance}_distance_matrix.qza \
      --m-metadata-file metadata.txt \
      --m-metadata-column ${column} \
      --o-visualization core-metrics-results/${distance}-${column}-significance.qzv \
      --p-pairwise
    done
    
    mkdir -p result/beta
    
    # 整理已有结果
    for distance in weighted_unifrac unweighted_unifrac bray_curtis jaccard;do
    column=Group
    mv core-metrics-results/${distance}* result/beta/

    #distance=weighted_unifrac
    unzip -d result/beta/${distance} result/beta/${distance}-${column}-significance.qzv
    mv result/beta/${distance}/*/data/* result/beta/${distance}/
    
    unzip -d result/beta/${distance}/pcoa result/beta/${distance}_emperor.qzv
    mv result/beta/${distance}/pcoa/*/data/* result/beta/${distance}/pcoa/
    done
    
    # 转beta多样性距离矩阵为普通txt file
    for distance in weighted_unifrac unweighted_unifrac bray_curtis jaccard;do
    qiime tools export  \
       --input-path result/beta/${distance}_distance_matrix.qza \
       --output-path result/beta/${distance}
    mv result/beta/${distance}/distance-matrix.tsv result/beta/${distance}.txt
    done
    
    # 基于R脚本进行可视化
    cut -f 1-4 result/metadata.txt > group.txt
    head group.txt
    for i in weighted_unifrac unweighted_unifrac bray_curtis jaccard
    do
    # -P添加行注释文件，-Q添加列注释
    bash ${db}/script/sp_pheatmap.sh \
      -f result/beta/${i}.txt \
      -H 'TRUE' -u 18 -v 16 -d "none" \
      -P group.txt -Q group.txt
      
    Rscript ${db}/script/beta_pcoa.R \
      --input result/beta/${i}.txt --design result/metadata.txt \
      --group Group1 --label F --width 89 --height 59 --shapeID Time -s F \
      --output result/beta/${i}.pcoa.time_shape.pdf
    done
    
    for i in weighted_unifrac unweighted_unifrac bray_curtis jaccard
    do
    for type in AA FA FF;
    do
    Rscript ${db}/script/beta_pcoa.R \
      --input result/beta/${i}.txt --design result/metadata_${type}.txt \
      --group Group --label F --width 89 --height 59 \
      --output result/beta/${type}_${i}.pcoa.pdf
    Rscript ${db}/script/beta_cpcoa.R \
      --input result/beta/${i}.txt --design result/metadata_${type}.txt \
      --group Group --label F --width 89 --height 59 \
      --output result/beta/${type}_${i}.cpcoa.pdf
    done
    done
    
## 4. 物种组成分析

    # 物种注释，数据库见附录，可选silva-138-99-qiime24.10-nb-classifier.qza
    time qiime feature-classifier classify-sklearn \
      --i-classifier gg_2024.09.backbone.full-length.nb.qza \
      --i-reads rep-seqs.qza \
      --o-classification taxonomy.qza
    # 可视化物种注释
    qiime metadata tabulate \
      --m-input-file taxonomy.qza \
      --o-visualization taxonomy.qzv
    # 堆叠柱状图展示
    qiime taxa barplot \
      --i-table table.qza \
      --i-taxonomy taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization result/tax/taxa-bar-plots.qzv
    
    # 整理结果用R可视化
    unzip -d result/tax result/tax/taxa-bar-plots.qzv
    mv result/tax/*/data/* result/tax
    rm -r `find result/tax/ -maxdepth 1 -mindepth 1 -type d|grep -v -w dist|grep -v -w q2templateassets`
    rm result/tax/*.jsonp
    
    declare -A tax_array
    tax_array=([1]=k [2]=p [3]=c [4]=o [5]=f [6]=g [7]=s)
    
    for i in "${!tax_array[@]}";do
    echo $i ${tax_array[$i]}
    awk 'BEGIN{FS=",";OFS="\t"}{for(i=1;i<=NF-4;i++) {\
    if(FNR==1) {split($i,arr,"__");\
      if(arr[length(arr)]=="") $i=$i"Unassigned";name[arr[length(arr)]]++};
      a[NR,i]=$i};}\
    END{for(m=1;m<=NF-4;m++) {split(a[1,m],arr,"__");\
    if(name[arr[length(arr)]]==1) printf arr[length(arr)];else {gsub("__;","Unassigned;",a[1,m]);printf a[1,m]};\
    for(n=2;n<=NR;n++) printf "\t"a[n,m];print ""}}' result/tax/level-${i}.csv|\
    sed 's/[kpcofg]__//g'|sed 's/__Unassigned/Unassigned/'  >result/tax/sum_${tax_array[$i]}.txt
    done
    
    for i in p c o f g s; do
    Rscript ${db}/script/tax_stackplot1.R \
      --input result/tax/sum_${i}.txt --design result/metadata.txt \
      --group Group1 --color manual1 --output result/tax/sum_${i}.stackplot \
      --legend 20 --width 200 --height 150 --subgroup Time \
      --subgroupOrder "0m,6m,12m" --groupOrder "FF,FA,AA" 
      
    Rscript ${db}/script/tax_stackplot1.R \
      --input result/tax/sum_${i}.txt --design result/metadata.txt \
      --group Group --color manual1 --output result/tax/sum_${i}.stackplot \
      --legend 20 --width 200 --height 150 \
      --groupOrder "FF.0m,FF.6m,FF.12m,FA.0m,FA.6m,FA.12m,AA.0m,AA.6m,AA.12m," 
    done
  
    unzip -d result/tax result/tax/taxonomy.qzv
    mv result/tax/*/data/metadata.tsv result/tax/taxonomy_confidence.txt
    rm -r `find result/tax/ -maxdepth 1 -mindepth 1 -type d|grep -v -w dist|grep -v -w q2templateassets`
    
    
## 5. 差异分析ancom

    # 格式化特征表，添加伪计数，4s
    qiime composition add-pseudocount \
      --i-table table.qza \
      --o-composition-table comp-table.qza
    
    # 计算差异特征，指定分组类型比较，1m
    column=Group
    time qiime composition ancom \
      --i-table comp-table.qza \
      --m-metadata-file metadata.txt \
      --m-metadata-column ${column} \
      --o-visualization ancom-${column}.qzv
    
    # 按属水平合并，并统计
    ## 按属水平合并，6s
    qiime taxa collapse \
      --i-table table.qza \
      --i-taxonomy taxonomy.qza \
      --p-level 6 \
      --o-collapsed-table table-l6.qza
    # 格式化特征表，添加伪计数，6s
    qiime composition add-pseudocount \
      --i-table table-l6.qza \
      --o-composition-table comp-table-l6.qza
    # 计算差异属，指定分组类型比较，16s
    qiime composition ancom \
      --i-table comp-table-l6.qza \
      --m-metadata-file metadata.txt \
      --m-metadata-column ${column} \
      --o-visualization ancom-l6-${column}.qzv

## 6. lefse分析
    qiime tools export \
      --input-path result/table.qza \
      --output-path exported-feature-table
      
    cp exported-feature-table/feature-table.biom result/
    biom convert -i result/feature-table.biom \
      -o result/feature-table.txt --to-tsv
    sed -i '1d' result/feature-table.txt
    
    awk 'BEGIN{OFS=FS="\t"}FNR>2{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
      result/tax/taxonomy_confidence.txt > otus.tax
    
    head otus.tax
    sed 's/;/\t/g;s/.__//g;' otus.tax|cut -f 1-8 | \
      sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
      > result/taxonomy.txt
    head -n3 result/taxonomy.txt
    
    mkdir -p result/lefse
    # threshold控制丰度筛选以控制作图中的枝数量
    Rscript ${db}/script/format2lefse1.R --input result/feature-table.txt \
      --taxonomy result/taxonomy.txt --design result/metadata.txt \
      --group Group --threshold 0 \
      --output result/lefse/LEfSe
    
    source activate /anaconda3/envs/lefse
    lefse-format_input.py LEfSe3.txt input.in -c 1 -o 1000000
    #运行lefse
    run_lefse.py input.in input.res
    #绘制物种树注释差异
    ## dpi会增大图形圈大小 labeled_stop_lev指定legend展示的物种等级 abrv_stop/start_lev 指定图中差异物种背景色展示的等级
    # lefse-plot_cladogram.py c_lefse.txt c_cladogram.pdf --format pdf \
    # --dpi 350 --labeled_stop_lev 7 --abrv_stop_lev 7 --abrv_start_lev 1 \
    # --right_space_prop 0.2 --left_space_prop 0.05
    lefse-plot_cladogram.py input.res cladogram.pdf --format pdf
    #绘制所有差异features柱状图
    lefse-plot_res.py input.res res.pdf --format pdf
    #绘制单个features柱状图(同STAMP中barplot)
    head input.res #查看差异features列表
    #批量绘制所有差异features柱状图，慎用(几百张差异结果柱状图阅读也很困难)
    mkdir -p features
    lefse-plot_features.py -f diff --archive none --format pdf \
      input.in input.res features/
      
# 以下分析要求数据必须有两类分组 Group和Time分组（Time分组必须为数值）
# 另外如果有个体的数据也可以提供；可以矫正个体差异的影响
## 7. alpha_longitudinal分析

    ### 7.1 longitudinal linear-mixed-effects
    # 本质为绘制分组线图--散点图和回归线（基于线性混合效应模型）
    # 在不同分组中；数据是否随时间变化
    for i in observed_features faith_pd shannon_entropy pielou_evenness
    do
    
    mkdir -p result/alpha_longitudinal/$i/lme/
    # 波动性分析
    time qiime longitudinal volatility \
      --m-metadata-file result/metadata.txt \
      --m-metadata-file result/alpha_longitudinal/alpha_vector.qza \
      --p-default-metric $i \
      --p-state-column Stage1 \
      --p-individual-id-column Person \
      --o-visualization result/alpha_longitudinal/${i}/${i}-group-volatility.qzv
    unzip result/alpha_longitudinal/${i}/${i}-group-volatility.qzv -d result/alpha_longitudinal/${i}/volatility/
    mv -f result/alpha_longitudinal/${i}/volatility/*/data/* result/alpha_longitudinal/${i}/volatility/
    
    # 线性混合效应
    time qiime longitudinal linear-mixed-effects \
      --m-metadata-file result/metadata.txt \
      --m-metadata-file result/alpha_longitudinal/alpha_vector.qza \
      --p-metric $i \
      --p-state-column Stage1 \
      --p-individual-id-column Person \
      --o-visualization result/alpha_longitudinal/${i}/${i}-group-linear-mixed-effects.qzv
    unzip result/alpha_longitudinal/${i}/${i}-group-linear-mixed-effects.qzv -d result/alpha_longitudinal/${i}/lme/
    mv -f result/alpha_longitudinal/${i}/lme/*/data/* result/alpha_longitudinal/${i}/lme/
    
    # 成对比较 两个时间点比较
    time qiime longitudinal pairwise-differences \
      --m-metadata-file result/metadata.txt \
      --m-metadata-file result/alpha_longitudinal/alpha_vector.qza \
      --p-metric $i \
      --p-state-column Stage1 \
      --p-state-1 0 \
      --p-state-2 1 \
      --p-individual-id-column Person \
      --p-replicate-handling random \
      --o-visualization result/alpha_longitudinal/${i}/${i}-0-1-pairwise-differences.qzv
    unzip result/alpha_longitudinal/${i}/${i}-0-1-pairwise-differences.qzv -d result/alpha_longitudinal/${i}/pairwise0-1/
    
    mv -f result/alpha_longitudinal/${i}/pairwise0-1/*/data/* result/alpha_longitudinal/${i}/pairwise0-1/
    
    # 成对比较
    time qiime longitudinal pairwise-differences \
      --m-metadata-file result/metadata.txt \
      --m-metadata-file result/alpha_longitudinal/alpha_vector.qza \
      --p-metric $i \
      --p-state-column Stage1 \
      --p-state-1 1 \
      --p-state-2 7 \
      --p-individual-id-column Person \
      --p-replicate-handling random \
      --o-visualization result/alpha_longitudinal/${i}/${i}-1-7-pairwise-differences.qzv
    unzip result/alpha_longitudinal/${i}/${i}-1-7-pairwise-differences.qzv -d result/alpha_longitudinal/${i}/pairwise1-7/
    mv -f result/alpha_longitudinal/${i}/pairwise1-7/*/data/* result/alpha_longitudinal/${i}/pairwise1-7/
    
    ### 先分析连续时间点的变化率，再用linear-mixed-effects可视化
    # 比较其他时间点与第一个时间点的差异
    mkdir -p result/alpha_longitudinal/$i/first-differences/
    
    time qiime longitudinal first-differences \
      --m-metadata-file result/metadata.txt \
      --m-metadata-file result/alpha_longitudinal/alpha_vector.qza  \
      --p-metric ${i} \
      --p-state-column Stage1 \
      --p-individual-id-column Person \
      --p-replicate-handling random \
      --o-first-differences result/alpha_longitudinal/${i}/${i}-first-differences.qza
    
    time qiime longitudinal linear-mixed-effects \
      --m-metadata-file result/alpha_longitudinal/${i}/${i}-first-differences.qza \
      --m-metadata-file result/metadata.txt \
      --p-metric Difference \
      --p-state-column Stage1 \
      --p-individual-id-column Person \
      --o-visualization result/alpha_longitudinal/${i}/${i}-group-first-distances-LME.qzv
    unzip result/alpha_longitudinal/${i}/${i}-group-first-distances-LME.qzv -d result/alpha_longitudinal/${i}/first-differences
    mv -f result/alpha_longitudinal/${i}/first-differences/*/data/* result/alpha_longitudinal/${i}/first-differences
    
    done
    
    # for i in `ls -d */lme`;do start $i/index.html;done
    # for i in `ls -d */first-differences`;do start $i/index.html;done
    # for i in `ls -d */pairwise1-7`;do start $i/index.html;done
    # for i in `ls -d */pairwise0-1`;do start $i/index.html;done
    
    # 使用R代码实现基于分组，比较不同时间点alpha多样性是否存在显著差异
    ## 线图、箱线图、配对箱线图
    # 查看alpha_or_single_species_compare_in_time_with_group.Rmd


## 8. beta_longitudinal分析
    ### 先分析连续时间点的变化率，再用linear-mixed-effects可视化
    for i in weighted_unifrac unweighted_unifrac bray_curtis jaccard;do
    qiime longitudinal pairwise-distances \
    --i-distance-matrix result/beta/${i}_distance_matrix.qza \
    --m-metadata-file result/metadata.txt \
    --p-group-column id1 \
    --p-state-column Group1 \
    --p-state-1 Baseline \
    --p-state-2 Treatment_1month \
    --p-individual-id-column Person \
    --p-replicate-handling random \
    --o-visualization result/beta_longitudinal/${i}/${i}-0-1-pairwise-distances.qzv
    unzip result/beta_longitudinal/${i}/${i}-0-1-pairwise-distances.qzv -d result/beta_longitudinal/${i}/pairwise0-1/
    mv -f result/beta_longitudinal/${i}/pairwise0-1/*/data/* result/beta_longitudinal/${i}/pairwise0-1/
    
    qiime longitudinal pairwise-distances \
    --i-distance-matrix result/beta/${i}_distance_matrix.qza \
    --m-metadata-file result/metadata.txt \
    --p-group-column id1 \
    --p-state-column Group1 \
    --p-state-1 Treatment_1month \
    --p-state-2 Treatment_7month \
    --p-individual-id-column Person \
    --p-replicate-handling random \
    --o-visualization result/beta_longitudinal/${i}/${i}-1-7-pairwise-distances.qzv
    unzip result/beta_longitudinal/${i}/${i}-1-7-pairwise-distances.qzv -d result/beta_longitudinal/${i}/pairwise1-7/
    mv -f result/beta_longitudinal/${i}/pairwise1-7/*/data/* result/beta_longitudinal/${i}/pairwise1-7/
    done
    
    for i in weighted_unifrac unweighted_unifrac bray_curtis jaccard
    do
    mkdir -p result/beta_longitudinal/$i/
    qiime longitudinal first-distances \
    --i-distance-matrix result/beta/${i}_distance_matrix.qza \
    --m-metadata-file result/metadata.txt \
    --p-state-column Stage1 \
    --p-individual-id-column Person \
    --p-replicate-handling random \
    --o-first-distances result/beta_longitudinal/${i}/${i}-first-distances.qza
    
    qiime longitudinal linear-mixed-effects \
    --m-metadata-file result/beta_longitudinal/${i}/${i}-first-distances.qza \
    --m-metadata-file result/metadata.txt \
    --p-metric Distance  \
    --p-state-column Stage1 \
    --p-individual-id-column Person \
    --o-visualization result/beta_longitudinal/${i}/${i}-group-first-distances-LME.qzv
    unzip result/beta_longitudinal/${i}/${i}-group-first-distances-LME.qzv -d result/beta_longitudinal/${i}/
    mv -f result/beta_longitudinal/${i}/*/data/* result/beta_longitudinal/${i}/
    done
    
    # for i in `ls -d */index.html`;do start $i;done
    # for i in `ls -d */pairwise0-1/index.html`;do start $i;done
    # for i in `ls -d */pairwise1-7/index.html`;do start $i;done
    
    # 使用R代码实现基于分组，比较不同时间点beta多样性是否存在显著差异
    ## 箱线图、配对pcoa图
    # 查看beta_compare_in_time_with_group.Rmd
    
## 9. feature_longitudinal分析
    time qiime longitudinal feature-volatility \
      --i-table result/table.qza \
      --m-metadata-file result/metadata.txt \
      --p-state-column Stage1 \
      --p-individual-id-column Person \
      --output-dir result/compare_longitudinal
    #accuracy_results.qzv  feature_importance.qza  filtered_table.qza  sample_estimator.qza  volatility_plot.qzv
    unzip result/compare_longitudinal/accuracy_results.qzv -d result/compare_longitudinal/
    mv result/compare_longitudinal/*/data/predictions.pdf  result/compare_longitudinal/
    mv result/compare_longitudinal/*/data/predictive_accuracy.tsv  result/compare_longitudinal/
    rm -rf mv result/compare_longitudinal/*/data/
    unzip result/compare_longitudinal/volatility_plot.qzv -d result/compare_longitudinal/
    mv -f result/compare_longitudinal/*/data/*  result/compare_longitudinal/
    
    # 使用R代码实现基于分组，比较不同时间点物种是否存在显著差异；挑选显著物种分析
    ## 线图、箱线图、配对箱线图
    # 查看muliti_species_compare_in_time_with_group.Rmd

# 附录

## 1. qiime2 2024.10安装

### 安装Conda

    # 下载、安装和启动conda
    wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -f
    ~/miniconda3/condabin/conda init
    # 关闭终端重新打开

### 方法1. Conda在线安装QIIME

    # 附软件在线安装和打包代码
    n=qiime2-amplicon-2024.10
    # 下载软件列表
    wget -c https://data.qiime2.org/distro/amplicon/${n}-py38-linux-conda.yml
    # 备用链接
    wget -c http://www.imeta.science/db/conda/${n}-py38-linux-conda.yml
    # 新环境安装，可在不同电脑服务器上安装成功后打包分发
    conda env create -n ${n} --file ${n}-py38-linux-conda.yml
    # 环境打包(可选，1.2G)
    conda pack -n ${n} -o ${n}.tar.gz

### 方法2. 本地安装QIIME

    n=qiime2-amplicon-2024.10
    # 安装包下载链接 
    wget -c ftp://download.nmdc.cn/tools/conda/${n}.tar.gz
    # 新环境安装
    mkdir -p ~/miniconda3/envs/${n}
    tar -xzf ${n}.tar.gz -C ~/miniconda3/envs/${n}
    # 激活并初始化环境
    conda activate ${n}
    conda unpack

## 2. 物种注释数据训练集

### Silva 138 99% OTUs full-length sequences
    
    ## 需要重新训练
    ## 从qiime2官网下载silva的序列和注释文件
    wget -c https://data.qiime2.org/2024.10/common/silva-138-99-seqs.qza
    wget -c https://data.qiime2.org/2024.10/common/silva-138-99-tax.qza
    
    time qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads silva-138-99-seqs.qza \
      --i-reference-taxonomy silva-138-99-tax.qza \
      --o-classifier silva-138-99-qiime24.10-nb-classifier.qza
      
    ## 以下版本不适合 qiime2024-10版本使用
    # 官网下载
    wget -c https://data.qiime2.org/2024.2/common/silva-138-99-nb-classifier.qza
    # 备用链接
    wget -c ftp://download.nmdc.cn/tools/amplicon/silva/silva-138-99-nb-classifier.qza
    

### Greengenes2 2022.10 full length sequences

    ## 2024.9版本的greengene 为:https://ftp.microbio.me/greengenes_release/2024.09/
    # 官网下载
    ## 全长
    wget -c --no-check-certificate https://ftp.microbio.me/greengenes_release/2024.09/2024.09.backbone.full-length.nb.qza
    ## V4
    wget -c --no-check-certificate https://ftp.microbio.me/greengenes_release/2024.09/2024.09.backbone.v4.nb.qza
    
    ## 2022.10版本
    # 官网下载
    wget -c https://data.qiime2.org/classifiers/greengenes/gg_2022_10_backbone_full_length.nb.qza
    wget -c http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.nb.qza
    # 备用链接
    wget -c ftp://download.nmdc.cn/tools/amplicon/GreenGenes/2022.10.backbone.full-length.nb.qza
    mv 2022.10.backbone.full-length.nb.qza gg_2022_10_backbone_full_length.nb.qza 
    
## 3. 物种注释数据训练集

    wd=/mnt/c/amplicon/qiime2
    mkdir -p $wd
    cd $wd
    # 下载数据库文件(greengenes, 320M)
    # wget -c ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
    # 国内备用链接
    # wget -c ftp://download.nmdc.cn/tools/amplicon/GreenGenes/gg_13_8_otus.tar.gz
    # 国内备份核心99数据库(60M)
    wget -c ftp://download.nmdc.cn/tools/amplicon/GreenGenes/gg_13_8_otus_99.tar.gz
    mv gg_13_8_otus_99.tar.gz gg_13_8_otus.tar.gz
    # 解压
    tar -zxvf gg_13_8_otus.tar.gz
    
    # 使用rep_set文件中的99_otus.fasta数据和taxonomy中的99_OTU_taxonomy.txt数据作为参考物种注释
    # 导入参考序列，50s
    qiime tools import \
      --type 'FeatureData[Sequence]' \
      --input-path gg_13_8_otus/rep_set/99_otus.fasta \
      --output-path 99_otus.qza
    # Fontconfig error: Cannot load default config file 不影响结果
    
    # 导入物种分类信息，6s
    qiime tools import \
      --type 'FeatureData[Taxonomy]' \
      --input-format HeaderlessTSVTaxonomyFormat \
      --input-path gg_13_8_otus/taxonomy/99_otu_taxonomy.txt \
      --output-path ref-taxonomy.qza
    
    # Train the classifier（训练分类器）——全长，30m
    time qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads 99_otus.qza \
      --i-reference-taxonomy ref-taxonomy.qza \
      --o-classifier classifier_gg_13_8_99.qza
    # 备用下载链接	wget -c http://bailab.genetics.ac.cn/db/GreenGenes/classifier_gg_13_8_99.qza
      
    # 引物提取参考序列的扩增区段 Extract reference reads
    # 常用Greengenes 13_8 99% OTUs from 341F CCTACGGGNGGCWGCAG/805R GACTACHVGGGTATCTAATCC region of sequences（分类器描述），提供测序的引物序列，截取对应的区域进行比对，达到分类的目的。
    # 本次使用引物799F-1193R，请根据实际替换, 8m
    time qiime feature-classifier extract-reads \
      --i-sequences 99_otus.qza \
      --p-f-primer AACMGGATTAGATACCCKG \
      --p-r-primer ACGTCATCCCCACCTTCC \
      --o-reads ref-seqs.qza
    # Train the classifier（训练分类器）
    # 基于筛选的指定区段，生成实验特异的分类器，7 min
    time qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads ref-seqs.qza \
      --i-reference-taxonomy ref-taxonomy.qza \
      --o-classifier classifier_gg_13_8_99_V5-V7.qza
    
    # 常见问题1：scikit-learn版本不兼容，重新fit-classifier-naive-bayes构建训练集即可
    Plugin error from feature-classifier:
      The scikit-learn version (0.21.2) used to generate this artifact does not match the current version of scikit-learn installed (0.22.1). Please retrain your classifier for your current deployment to prevent data-corruption errors.
    Debug info has been saved to /tmp/qiime2-q2cli-err-5ngzk2hm.log
    
    # 常见问题2: dada2运行报错，环境配置不完全，运行conda unpack初始化
    Plugin error from dada2:
    An error was encountered while running DADA2 in R (return code 255), please inspect stdout and stderr to learn more.
    Debug info has been saved to /tmp/qiime2-q2cli-err-utwt1cmu.log

## 4. HOMD数据库训练
# 整理数据库
grep ">" HOMD_16S_rRNA_RefSeq_V15.23.fasta |awk 'BEGIN{OFS="\t";FS="|"}{split($1,a," ");split($3,b," ");sub(">","",a[1]);print a[1],b[1]}' >id-HMTid.txt
awk 'BEGIN{OFS=FS="\t";}ARGIND==1&&FNR>2{a[$9]="k__"$1";p__"$2";c__"$3";o__"$4";f__"$5";g__"$6";s__"$7}ARGIND==2{print $1,$2,a[$2]}' HOMD_taxon_table2023-12-20_1703099281.txt id-HMTid.txt >id_hmtid_ann.txt
awk 'BEGIN{OFS=FS="\t";}ARGIND==1&&FNR>2{a[$9]="k__"$1";p__"$2";c__"$3";o__"$4";f__"$5";g__"$6";s__"$7}\
ARGIND==2{tax=a[$2];gsub(" ","_",tax);print $1,tax}' HOMD_taxon_table2023-12-20_1703099281.txt id-HMTid.txt >id_ann.txt

# 训练数据库
 # 使用rep_set文件中的99_otus.fasta数据和taxonomy中的99_OTU_taxonomy.txt数据作为参考物种注释
    # 导入参考序列，50s
    qiime tools import \
      --type 'FeatureData[Sequence]' \
      --input-path HOMD_16S_rRNA_RefSeq_V15.23.fasta \
      --output-path HOMD_v15.23_otus.qza
    # Fontconfig error: Cannot load default config file 不影响结果
    
    # 导入物种分类信息，6s
    qiime tools import \
      --type 'FeatureData[Taxonomy]' \
      --input-format HeaderlessTSVTaxonomyFormat \
      --input-path id_ann.txt \
      --output-path ref-taxonomy.qza
    
    # Train the classifier（训练分类器）——全长，30m
    time qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads HOMD_v15.23_otus.qza \
      --i-reference-taxonomy ref-taxonomy.qza \
      --o-classifier classifier_HOMD_v15.23_otus.qza
      
    qiime tools import \
      --type 'FeatureData[Sequence]' \
      --input-path HOMD_16S_rRNA_RefSeq_V15.23.p9.fasta \
      --output-path HOMD_v15.23.p9_otus.qza
    time qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads HOMD_v15.23.p9_otus.qza \
      --i-reference-taxonomy ref-taxonomy.qza \
      --o-classifier classifier_HOMD_v15.23.p9_otus.qza

##  5. Qiime2文件转换
    # 输入文件修改格式
    ## fa:FeatureData[Sequence]或者FeatureData[AlignedSequence]（等长fasta文件含-）
    ## Feature table data:FeatureTable[Frequency]  + --input-format BIOMV100Format(biom v.1.0.0) BIOMV210Format(biom v.2.1.0)
    ## 系统发育树:Phylogeny[Unrooted]
    ## AlphaDiversityFormat:SampleData[AlphaDiversity]?
    ## :DistanceMatrix
    ## 查看qiime2支持的格式
    qiime tools import \
      --show-importable-formats
    
    qiime tools import \
      --show-importable-types
    
    # fasta - qza
    qiime tools import \
      --input-path sequences.fna \
      --output-path sequences.qza \
      --type 'FeatureData[Sequence]'
      
    # alpha.txt - qza
    qiime tools import \
      --input-path ../alpha/alpha.txt \
      --output-path alpha.qza \
      --input-format AlphaDiversityFormat \
      --type 'SampleData[AlphaDiversity]'
    
    # beta.txt - qza
    qiime tools import \
      --input-path ../beta/bray_curtis.txt \
      --output-path beta_bray_curtis.qza  \
      --type 'DistanceMatrix'
    
  
## ref
    https://library.qiime2.org/quickstart/amplicon#id-1-installing-miniconda
    https://docs.qiime2.org/2024.10/
    
## Citation引文

    If used this script, please cited:
    使用此脚本，请引用下文：
    
    Evan Bolyen, Jai Ram Rideout, Matthew R. Dillon, Nicholas A. Bokulich, Christian C. Abnet, 
    Gabriel A. Al-Ghalith, Harriet Alexander, Eric J. Alm, Manimozhiyan Arumugam, Francesco Asnicar, 
    Yang Bai, Jordan E. Bisanz, Kyle Bittinger, Asker Brejnrod, Colin J. Brislawn, C. Titus Brown, 
    Benjamin J. Callahan, Andrés Mauricio Caraballo-Rodríguez, John Chase, Emily K. Cope, 
    Ricardo Da Silva, Christian Diener, Pieter C. Dorrestein, Gavin M. Douglas, Daniel M. Durall, 
    Claire Duvallet, Christian F. Edwardson, Madeleine Ernst, Mehrbod Estaki, Jennifer Fouquier, 
    Julia M. Gauglitz, Sean M. Gibbons, Deanna L. Gibson, Antonio Gonzalez, Kestrel Gorlick, 
    Jiarong Guo, Benjamin Hillmann, Susan Holmes, Hannes Holste, Curtis Huttenhower, Gavin A. Huttley, 
    Stefan Janssen, Alan K. Jarmusch, Lingjing Jiang, Benjamin D. Kaehler, Kyo Bin Kang, 
    Christopher R. Keefe, Paul Keim, Scott T. Kelley, Dan Knights, Irina Koester, Tomasz Kosciolek, 
    Jorden Kreps, Morgan G. I. Langille, Joslynn Lee, Ruth Ley, **Yong-Xin Liu**, Erikka Loftfield, 
    Catherine Lozupone, Massoud Maher, Clarisse Marotz, Bryan D. Martin, Daniel McDonald, 
    Lauren J. McIver, Alexey V. Melnik, Jessica L. Metcalf, Sydney C. Morgan, Jamie T. Morton, 
    Ahmad Turan Naimey, Jose A. Navas-Molina, Louis Felix Nothias, Stephanie B. Orchanian, 
    Talima Pearson, Samuel L. Peoples, Daniel Petras, Mary Lai Preuss, Elmar Pruesse, 
    Lasse Buur Rasmussen, Adam Rivers, Michael S. Robeson, Patrick Rosenthal, Nicola Segata, 
    Michael Shaffer, Arron Shiffer, Rashmi Sinha, Se Jin Song, John R. Spear, Austin D. Swafford, 
    Luke R. Thompson, Pedro J. Torres, Pauline Trinh, Anupriya Tripathi, Peter J. Turnbaugh, 
    Sabah Ul-Hasan, Justin J. J. van der Hooft, Fernando Vargas, Yoshiki Vázquez-Baeza, 
    Emily Vogtmann, Max von Hippel, William Walters, Yunhu Wan, Mingxun Wang, Jonathan Warren, 
    Kyle C. Weber, Charles H. D. Williamson, Amy D. Willis, Zhenjiang Zech Xu, Jesse R. Zaneveld, 
    Yilong Zhang, Qiyun Zhu, Rob Knight, J. Gregory Caporaso. 
    2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. 
    **Nature Biotechnology** 37: 852-857. https://doi.org/10.1038/s41587-019-0209-9
    
    Copyright 2016-2024 Yong-Xin Liu <liuyongxin@caas.cn>