#!/bin/bash

# Pono远端编译和运行测试脚本

set -e

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_test() {
    echo -e "${YELLOW}[TEST]${NC} $1"
}

# 测试计数器
TESTS_PASSED=0
TESTS_FAILED=0

run_test() {
    local test_name="$1"
    local test_command="$2"
    
    print_test "运行测试: $test_name"
    
    if eval "$test_command" > /dev/null 2>&1; then
        print_success "✓ $test_name"
        ((TESTS_PASSED++))
    else
        print_error "✗ $test_name"
        ((TESTS_FAILED++))
    fi
}

print_info "开始Pono远端编译和运行测试..."

# 测试1: 检查脚本权限
run_test "检查脚本执行权限" "test -x ./docker-build-and-run.sh && test -x ./remote-build.sh && test -x ./run-pono.sh"

# 测试2: 检查Docker可用性
run_test "检查Docker可用性" "command -v docker && docker info"

# 测试3: 测试Docker脚本帮助
run_test "Docker脚本帮助" "./docker-build-and-run.sh --help"

# 测试4: 测试智能运行脚本帮助
run_test "智能运行脚本帮助" "./run-pono.sh --help"

# 测试5: 测试远端编译脚本帮助
run_test "远端编译脚本帮助" "./remote-build.sh --help"

# 测试6: 测试Pono版本（使用预构建镜像）
print_test "测试Pono版本信息"
if ./docker-build-and-run.sh -p -- --version 2>/dev/null | grep -q "v[0-9]"; then
    print_success "✓ Pono版本信息"
    ((TESTS_PASSED++))
else
    print_error "✗ Pono版本信息"
    ((TESTS_FAILED++))
fi

# 测试7: 测试Pono帮助信息
print_test "测试Pono帮助信息"
if ./docker-build-and-run.sh -p -- --help 2>/dev/null | grep -q "USAGE"; then
    print_success "✓ Pono帮助信息"
    ((TESTS_PASSED++))
else
    print_error "✗ Pono帮助信息"
    ((TESTS_FAILED++))
fi

# 测试8: 测试样例文件验证
if [[ -f "samples/constarrfalse.btor2" ]]; then
    print_test "测试BTOR2文件验证"
    if ./run-pono.sh -v ./samples -- /workspace/constarrfalse.btor2 2>/dev/null | grep -q "unknown\|sat\|unsat"; then
        print_success "✓ BTOR2文件验证"
        ((TESTS_PASSED++))
    else
        print_error "✗ BTOR2文件验证"
        ((TESTS_FAILED++))
    fi
else
    print_warning "跳过BTOR2文件验证测试（样例文件不存在）"
fi

# 测试9: 测试智能运行脚本自动检测
print_test "测试智能运行脚本自动检测"
if ./run-pono.sh -- --version 2>/dev/null | grep -q "v[0-9]"; then
    print_success "✓ 智能运行脚本自动检测"
    ((TESTS_PASSED++))
else
    print_error "✗ 智能运行脚本自动检测"
    ((TESTS_FAILED++))
fi

# 测试10: 检查文档文件
run_test "检查文档文件存在" "test -f REMOTE_BUILD_GUIDE.md"

# 输出测试结果
echo
print_info "测试完成!"
print_info "通过: $TESTS_PASSED"
print_info "失败: $TESTS_FAILED"
print_info "总计: $((TESTS_PASSED + TESTS_FAILED))"

if [[ $TESTS_FAILED -eq 0 ]]; then
    print_success "所有测试通过! 🎉"
    echo
    print_info "您现在可以使用以下方式运行Pono:"
    echo
    echo "1. 智能运行（推荐）:"
    echo "   ./run-pono.sh -- --help"
    echo "   ./run-pono.sh -v ./samples -- /workspace/constarrfalse.btor2"
    echo
    echo "2. 直接使用Docker:"
    echo "   ./docker-build-and-run.sh -p -- --help"
    echo "   ./docker-build-and-run.sh -v ./samples -- /workspace/uv_example.btor2"
    echo
    echo "3. 远端服务器编译:"
    echo "   ./remote-build.sh"
    echo "   ./build/pono --help"
    echo
    echo "更多信息请查看: REMOTE_BUILD_GUIDE.md"
    exit 0
else
    print_error "部分测试失败，请检查配置"
    exit 1
fi
