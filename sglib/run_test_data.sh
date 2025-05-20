#!/bin/bash
# Test script for the fixed Robust Savitzky-Golay Filter library
# Author: Mapoet
# Date: 2025-05-20

echo "====== 测试修复后的Savitzky-Golay滤波器 ======"

# 确保当前在项目根目录
cd "$(dirname "$0")"

# 确保构建目录存在
mkdir -p build

# 清理旧的构建
make clean

# 编译库
echo "正在编译库..."
make release

if [ $? -ne 0 ]; then
    echo "错误: 编译失败"
    exit 1
fi

# 生成测试数据
echo "生成测试数据..."
python3 scripts/generate_test_data.py -n 300 --noise 0.1 --outlier-prob 0.07 --outlier-scale 5.0 -p

# 运行测试程序
echo "运行测试程序..."
cd build
../build/test_sg_filter
cd ..

# 绘制结果图表
echo "生成图表..."
python3 scripts/plot_results.py -i build/sg_filter_results.dat -o sg_filter_results.png

echo "测试完成！请检查生成的图表 sg_filter_results.png" 