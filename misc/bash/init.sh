#!/bin/bash

# 设置环境名称和要安装的包
ENV_NAME="r_env"
PACKAGES="r-essentials r-base=4.4"

# 函数：尝试创建环境并安装包
create_env() {
	echo "尝试创建 Conda 环境并安装包..."
	conda create -n $ENV_NAME $PACKAGES -y
}

# 无限循环，直到成功创建环境并安装包
while true; do
	create_env
	if [ $? -eq 0 ]; then
		echo "成功创建 Conda 环境并安装包。"
		break
	else
		echo "创建 Conda 环境失败，重试..."
	fi
done

