# lung cancer stomics

> 本仓库用于处理肺癌时空组数据

```txt
.
|-- README.md
|-- LICENSE
|-- data
|-- log
|-- temp
|-- out
|-- figures
|-- makefile
|-- requirements.txt
|-- misc                        杂项
|   |-- R
|   |-- bash
|   `-- python
`-- src
    |-- NMF.py                  TODO
    |-- annotation.R            使用SPOTlight细胞注释
    |-- tissueAna.py            组织数据分析主函数
    `-- utils
        |-- config.py           各种参数设置和配置数据
        |-- loader.py           读取数据
        |-- filter.py           细胞过滤
        |-- annotation.py       手动注释, 主要是输出reference基因点图
        |-- converter.py        格式转换
        `-- imgcatcher.py       用于保存stereopy图片
```

