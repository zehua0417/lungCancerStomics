#!/bin/bash

if [ -z "$1" ]; then
    echo "enter target:"
    read target
else
    target=$1
fi
if [ -z "$2" ]; then
    echo "enter pid:"
    read pid
else
    pid=$2
fi

# 日志文件路径
logfile=log/${target}.log
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

    # 检查进程是否结束
    if ! ps -p $pid > /dev/null; then
        echo "Process has finished."
        break
    fi

    # 检查是否按下 e 键退出
    if read -t 0.1 -n 1 input && [ "$input" = "e" ]; then
        echo "Exiting log monitor."
        echo "if you want to monitor the log again, please run: "
        echo "make reload_log then enter $pid"
        break
    fi

    # 检查是否按下 q 键退出
    if read -t 0.1 -n 1 input && [ "$input" = "q" ]; then
        echo "Exiting python program."
        kill -9 $pid
        break
    fi

    # 每次检查间隔100毫秒
    sleep 0.1
done
