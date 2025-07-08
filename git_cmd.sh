#!/bin/bash

# ��¼Github
# ssh -T git@github.com

PROJECT_NAME="${PWD##*/}"
REMOTE_URL="https://github.com/sunichmd/${PROJECT_NAME}.git"  # �ɰ����Ϊ��̬����
SSH_URL="git@github.com:sunichmd/${PROJECT_NAME}.git"

echo "��Ŀ���ƣ�Ŀ¼������$PROJECT_NAME"

# ��������������ɼ����ļ�
if [[ "$1" != "" ]]; then
    echo "����ɨ����� 100MB ���ļ�..."
    BIG_FILES=$(find . -type f -size +100M) # �����ļ���������100M
fi

if [[ -n "$BIG_FILES" ]]; then
    echo "��⵽�����ļ����� 100MB�����Զ���ӵ� .gitignore��"
    echo "$BIG_FILES"
    echo "$BIG_FILES" >> .gitignore
    sort -u .gitignore -o .gitignore  # ȥ��
	# ɾ���ļ�����
	for i in `echo $BIG_FILES`;
	do
		git rm --cached $i
	done
fi

if [ ! -d ".git" ]; then
    echo "��һ��ִ�У���ʼ��git�ֿⲢ����"
    if [ ! -f "README.md" ]; then
        echo "# $PROJECT_NAME" > README.md
        echo "�Զ�����README.md������Ϊ��Ŀ��: $PROJECT_NAME"
    fi
	
    git init
    git add README.md .gitignore
    git commit -m "first commit: initialize $PROJECT_NAME repo"
	git branch -M main
    
    git remote add origin "$SSH_URL"
    git push -u origin main
else
    echo "�Ѵ���git�ֿ⣬ִ�й���Զ�ֿ̲Ⲣ����"
	
	if git remote get-url origin > /dev/null 2>&1; then
        echo "Զ�� origin �Ѵ��ڣ��������"
    else
		git remote add origin "$SSH_URL"
	fi
	# ���ӻ�����
	git config --global http.postBuffer 524288000
	# ʹ�ø�ǿ��ѹ���������
	git config --global core.compression 9
	
	echo -e "\n?? ���δ�� ignore ���ļ�..."
    git add -u  # ��׷���ļ�����
    git add .   # ���ļ���.gitignore �ѹ��ˣ�

	git commit -m "Auto-commit: $(date '+%Y-%m-%d %H:%M:%S')"
    git branch -M main
    git push -u origin main
fi

