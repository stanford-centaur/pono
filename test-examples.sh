#!/bin/bash

# Pono示例测试脚本

set -e

# 颜色输出
BLUE='\033[0;34m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_test() {
    echo -e "${YELLOW}[TEST]${NC} $1"
}

echo "🚀 Pono示例测试"
echo "================="

# 测试1: 简单的BTOR2文件
print_test "测试1: constarrfalse.btor2"
./run-pono.sh -v ./samples -- /workspace/constarrfalse.btor2
echo

# 测试2: 使用BMC引擎
print_test "测试2: 使用BMC引擎"
./run-pono.sh -v ./samples -- --engine bmc --bound 5 /workspace/uv_example.btor2
echo

# 测试3: 使用IC3引擎
print_test "测试3: 使用IC3引擎"
./run-pono.sh -v ./samples -- --engine ic3 /workspace/constarrfalse.btor2
echo

# 测试4: 详细输出
print_test "测试4: 详细输出模式"
./run-pono.sh -v ./samples -- --verbosity 2 /workspace/int_win.btor2
echo

# 测试5: 生成反例
print_test "测试5: 尝试生成反例"
./run-pono.sh -v ./samples -- --witness /workspace/uv_example.btor2
echo

print_success "所有测试完成！"

echo
echo "📚 可用的BTOR2文件："
ls samples/*.btor2 | sed 's|samples/|  - |g'

echo
echo "💡 使用提示："
echo "  - 使用 -v ./samples 挂载样例目录"
echo "  - 文件路径使用 /workspace/文件名"
echo "  - 尝试不同的引擎: bmc, ic3, kind"
echo "  - 调整边界: --bound N"
echo "  - 详细输出: --verbosity N"
