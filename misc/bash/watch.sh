#!/bin/bash

if [ -z "$1" ]; then
    echo "enter target:"
    read target
else
    target=$1
fi

# 日志文件路径
cd /hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer
if [ -f "log/${target}.log" ]; then
    logfile="log/${target}.log"
elif [ -f "./${target}" ]; then
    logfile="${target}"
else
    echo "Error: Neither log/${target}.log nor ./${target} exists." >&2
    exit 1
fi
last_mod_time=$(stat -c %Y "$logfile")

# 初始化显示文件内容
clear
cat "$logfile"

# 循环检测文件修改时间和进程状态
while true; do
    # 检查文件修改时间
    current_mod_time=$(stat -c %Y "$logfile")
    if [ "$current_mod_time" -ne "$last_mod_time" ]; then
        clear
        cat "$logfile"
        last_mod_time=$current_mod_time
    fi

    # 检查是否按下 q 键退出
    if read -t 0.1 -n 1 input && [ "$input" = "q" ]; then
        echo "Exiting watching program."
        break
    fi

    # 每次检查间隔100毫秒
    sleep 0.1
done
