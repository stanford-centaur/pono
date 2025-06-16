#!/bin/bash

# Pono运行脚本 - 支持本地和远端执行
# 自动检测环境并选择最佳运行方式

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

show_usage() {
    cat << EOF
Pono智能运行脚本

此脚本会自动检测环境并选择最佳的运行方式：
1. 如果存在本地编译的pono，直接使用
2. 如果安装了Docker，使用Docker容器运行
3. 否则提示用户编译或安装

使用方法:
    $0 [选项] [-- pono参数]

选项:
    -h, --help          显示此帮助信息
    -f, --force-docker  强制使用Docker运行
    -b, --build         如果使用Docker，强制重新构建镜像
    -v, --volume DIR    挂载本地目录到Docker容器
    --examples          显示使用示例

Pono参数示例:
    --help              显示Pono帮助信息
    --engine bmc        使用BMC引擎
    --bound 10          设置展开深度为10
    --verbosity 2       设置详细输出级别

示例:
    $0 -- --help                           # 显示Pono帮助
    $0 -- samples/counter-false.btor2      # 验证BTOR2文件
    $0 -v ./samples -- counter-false.btor2 # 使用Docker并挂载samples目录
    $0 -- --engine bmc --bound 20 test.btor2  # 使用BMC引擎，深度20

EOF
}

show_examples() {
    cat << EOF
Pono使用示例

1. 基本模型检查:
   $0 -- samples/counter-false.btor2

2. 使用不同的验证引擎:
   $0 -- --engine bmc samples/counter-false.btor2      # BMC
   $0 -- --engine ic3 samples/counter-false.btor2      # IC3
   $0 -- --engine kind samples/counter-false.btor2     # K-induction

3. 设置参数:
   $0 -- --bound 50 --engine bmc samples/counter-false.btor2
   $0 -- --verbosity 3 samples/counter-false.btor2

4. 生成反例:
   $0 -- --witness samples/counter-false.btor2

5. 使用Docker运行:
   $0 -f -v ./samples -- counter-false.btor2

更多信息请运行: $0 -- --help

EOF
}

# 默认参数
FORCE_DOCKER=false
BUILD_DOCKER=false
VOLUME_DIR=""
PONO_ARGS=""

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_usage
            exit 0
            ;;
        --examples)
            show_examples
            exit 0
            ;;
        -f|--force-docker)
            FORCE_DOCKER=true
            shift
            ;;
        -b|--build)
            BUILD_DOCKER=true
            shift
            ;;
        -v|--volume)
            VOLUME_DIR="$2"
            shift 2
            ;;
        --)
            shift
            PONO_ARGS="$@"
            break
            ;;
        *)
            print_error "未知选项: $1"
            show_usage
            exit 1
            ;;
    esac
done

print_info "Pono智能运行脚本启动..."

# 检查本地编译版本
LOCAL_PONO=""
if [[ -f "./build/pono" ]]; then
    LOCAL_PONO="./build/pono"
elif [[ -f "./pono" ]]; then
    LOCAL_PONO="./pono"
elif command -v pono &> /dev/null; then
    LOCAL_PONO="pono"
fi

# 决定运行方式
if [[ "$FORCE_DOCKER" == true ]]; then
    print_info "强制使用Docker运行"
    RUN_METHOD="docker"
elif [[ -n "$LOCAL_PONO" ]] && [[ "$FORCE_DOCKER" == false ]]; then
    print_info "发现本地Pono: $LOCAL_PONO"
    RUN_METHOD="local"
elif command -v docker &> /dev/null && docker info &> /dev/null 2>&1; then
    print_info "未找到本地Pono，使用Docker运行"
    RUN_METHOD="docker"
else
    print_error "未找到可用的运行方式"
    cat << EOF

请选择以下方式之一：

1. 编译本地版本:
   ./remote-build.sh

2. 安装Docker并使用容器:
   # 安装Docker后运行:
   $0 -f $@

3. 使用预构建的Docker镜像:
   ./docker-build-and-run.sh -p -- $PONO_ARGS

EOF
    exit 1
fi

# 执行运行
case $RUN_METHOD in
    local)
        print_info "使用本地Pono: $LOCAL_PONO"
        print_info "执行命令: $LOCAL_PONO $PONO_ARGS"
        exec $LOCAL_PONO $PONO_ARGS
        ;;
    docker)
        print_info "使用Docker运行Pono"
        DOCKER_ARGS=""
        if [[ "$BUILD_DOCKER" == true ]]; then
            DOCKER_ARGS="$DOCKER_ARGS -b"
        fi
        if [[ -n "$VOLUME_DIR" ]]; then
            DOCKER_ARGS="$DOCKER_ARGS -v $VOLUME_DIR"
        fi
        exec ./docker-build-and-run.sh $DOCKER_ARGS -- $PONO_ARGS
        ;;
esac
