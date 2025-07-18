#!/bin/bash

# 登录Github
# ssh -T git@github.com

PROJECT_NAME="${PWD##*/}"
REMOTE_URL="https://github.com/sunichmd/${PROJECT_NAME}.git"  # 可按需改为动态参数
SSH_URL="git@github.com:sunichmd/${PROJECT_NAME}.git"

echo "项目名称（目录名）：$PROJECT_NAME"

if [ ! -d ".git" ]; then
    echo "第一次执行，初始化git仓库并推送"
    if [ ! -f "README.md" ]; then
        echo "# $PROJECT_NAME" > README.md
        echo "自动创建README.md，内容为项目名: $PROJECT_NAME"
    fi
	
    git init
    git add README.md .gitignore
    git commit -m "first commit: initialize $PROJECT_NAME repo"
	git branch -M main
    
    git remote add origin "$SSH_URL"
    git push -u origin main
else
    echo "已存在git仓库，执行关联远程仓库并推送"
	
	if git remote get-url origin > /dev/null 2>&1; then
        echo "远程 origin 已存在，跳过添加"
    else
		git remote add origin "$SSH_URL"
	fi
	# 增加缓冲区
	git config --global http.postBuffer 524288000
	# 使用更强的压缩减少体积
	git config --global core.compression 9
	
	echo -e "\n?? 添加未被 ignore 的文件..."
    git add -u  # 已追踪文件更新
    git add .   # 新文件（.gitignore 已过滤）

	git commit -m "Auto-commit: $(date '+%Y-%m-%d %H:%M:%S')"
    git branch -M main
    git push -u origin main
fi

