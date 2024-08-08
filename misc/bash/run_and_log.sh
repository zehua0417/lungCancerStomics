#!/bin/bash

if [ -z "$1" ]; then
    echo "enter python path:"
    read python
else
    python=$1
fi
if [ -z "$2" ]; then
    echo "enter target:"
    read target
else
    target=$2
fi

if [ ! -d "log" ]; then
    mkdir log
fi

# 运行 make target
PYTHONUNBUFFERED=1 "${python}" src/${target}.py > log/${target}.log 2>&1 &

# 获取 target 的进程ID
target_pid=$!

# 创建日志文件 $target.log
logfile=log/${target}.log

# 初始化显示文件内容
clear
cat $logfile

# 获取文件修改时间
last_mod_time=$(stat -c %Y $logfile)

# 循环检测文件修改时间和进程状态
while true; do
    # 检查文件修改时间
    current_mod_time=$(stat -c %Y $logfile)
    if [ "$current_mod_time" -ne "$last_mod_time" ]; then
        clear
        cat "$logfile"
        last_mod_time=$current_mod_time
    fi

    # 检查进程是否结束
    if ! kill -0 $target_pid 2>/dev/null; then
        echo "Process $target has finished."
        break
    fi

    # 检查是否按下 e 键退出log monitor
    if read -t 0.1 -n 1 input && [ "$input" = "e" ]; then
        echo "exiting log monitor."
        echo "if you want to monitor the log again, please run: "
        echo "make reload_log $target_pid then enter $target_pid"
        break
    fi

    # 检查是否按下 q 键退出python程序
    if read -t 0.1 -n 1 input && [ "$input" = "q" ]; then
        echo "Exiting python program."
        kill -9 $target_pid
        break
    fi

    # 每次检查间隔100毫秒
    sleep 0.1
done
