#!/bin/bash

# 遍历当前目录中的所有文件
for file in cds_data/*.fa; do
    python3 02.range.py "$file"
done

