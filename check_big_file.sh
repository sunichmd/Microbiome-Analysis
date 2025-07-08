#!/bin/bash

# 输入任意参数即可检查大文件
echo "正在扫描大于 100MB 的文件..."
BIG_FILES=$(find . -type f -size +100M|grep -v "/.git/") # 单个文件不允许超过100M

if [[ -n "$BIG_FILES" ]]; then
    echo "检测到以下文件超过 100MB，已自动添加到 .gitignore："
    echo "$BIG_FILES"
    echo "$BIG_FILES" >> .gitignore
    sort -u .gitignore -o .gitignore  # 去重
	
	for i in `echo $BIG_FILES`;
	do
	  # 删除文件缓存
		# git rm --cached $i
		# 移动文件到上一层的big files里
		mv $i ../Big_file/
	done
	
fi