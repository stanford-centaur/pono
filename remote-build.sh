#!/bin/bash

# Pono远端服务器编译脚本
# 适用于Ubuntu/Debian系统

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

# 显示使用说明
show_usage() {
    cat << EOF
Pono远端服务器编译脚本

使用方法:
    $0 [选项]

选项:
    -h, --help          显示此帮助信息
    -j, --jobs N        并行编译任务数 (默认: CPU核心数)
    --with-msat         编译时包含MathSAT支持 (需要手动下载许可证)
    --with-btor         编译时包含Boolector支持
    --debug             编译调试版本
    --clean             清理之前的构建
    --install           编译后安装到系统

示例:
    $0                  # 基本编译
    $0 -j 8             # 使用8个并行任务编译
    $0 --debug          # 编译调试版本
    $0 --clean --install # 清理重新编译并安装

EOF
}

# 默认参数
JOBS=$(nproc)
WITH_MSAT=false
WITH_BTOR=false
DEBUG=false
CLEAN=false
INSTALL=false

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_usage
            exit 0
            ;;
        -j|--jobs)
            JOBS="$2"
            shift 2
            ;;
        --with-msat)
            WITH_MSAT=true
            shift
            ;;
        --with-btor)
            WITH_BTOR=true
            shift
            ;;
        --debug)
            DEBUG=true
            shift
            ;;
        --clean)
            CLEAN=true
            shift
            ;;
        --install)
            INSTALL=true
            shift
            ;;
        *)
            print_error "未知选项: $1"
            show_usage
            exit 1
            ;;
    esac
done

print_info "开始Pono远端编译流程..."
print_info "使用 $JOBS 个并行任务"

# 检查操作系统
if [[ ! -f /etc/os-release ]]; then
    print_error "无法检测操作系统版本"
    exit 1
fi

source /etc/os-release
print_info "检测到操作系统: $PRETTY_NAME"

# 安装依赖
print_info "安装系统依赖..."
case $ID in
    ubuntu|debian)
        sudo apt-get update
        sudo apt-get install -y \
            build-essential \
            cmake \
            curl \
            git \
            libgmp-dev \
            m4 \
            meson \
            ninja-build \
            pkg-config \
            python3 \
            python3-pip \
            python3-venv \
            wget
        ;;
    centos|rhel|fedora)
        if command -v dnf &> /dev/null; then
            sudo dnf install -y \
                gcc gcc-c++ \
                cmake \
                curl \
                git \
                gmp-devel \
                m4 \
                meson \
                ninja-build \
                pkgconfig \
                python3 \
                python3-pip \
                wget
        else
            sudo yum install -y \
                gcc gcc-c++ \
                cmake \
                curl \
                git \
                gmp-devel \
                m4 \
                pkgconfig \
                python3 \
                python3-pip \
                wget
        fi
        ;;
    *)
        print_warning "未知的Linux发行版，请手动安装依赖"
        ;;
esac

# 清理之前的构建
if [[ "$CLEAN" == true ]]; then
    print_info "清理之前的构建..."
    rm -rf build deps
fi

# 设置依赖
print_info "设置Bison和Flex..."
if ! command -v bison &> /dev/null || ! command -v flex &> /dev/null; then
    ./contrib/setup-bison.sh
    ./contrib/setup-flex.sh
fi

print_info "设置SMT-Switch..."
SETUP_SMT_ARGS=""
if [[ "$WITH_MSAT" == true ]]; then
    SETUP_SMT_ARGS="$SETUP_SMT_ARGS --with-msat"
    print_warning "MathSAT需要手动下载并放置在 ./deps/mathsat/ 目录"
fi
if [[ "$WITH_BTOR" == true ]]; then
    SETUP_SMT_ARGS="$SETUP_SMT_ARGS --with-btor"
fi

./contrib/setup-smt-switch.sh $SETUP_SMT_ARGS

print_info "设置BTOR2Tools..."
./contrib/setup-btor2tools.sh

# 配置构建
print_info "配置构建系统..."
CONFIGURE_ARGS=""
if [[ "$WITH_MSAT" == true ]]; then
    CONFIGURE_ARGS="$CONFIGURE_ARGS --with-msat"
fi
if [[ "$WITH_BTOR" == true ]]; then
    CONFIGURE_ARGS="$CONFIGURE_ARGS --with-btor"
fi
if [[ "$DEBUG" == true ]]; then
    CONFIGURE_ARGS="$CONFIGURE_ARGS --debug"
fi

./configure.sh $CONFIGURE_ARGS

# 编译
print_info "开始编译 (使用 $JOBS 个并行任务)..."
cd build
make -j$JOBS

print_success "编译完成!"

# 运行测试
print_info "运行测试..."
if make check; then
    print_success "所有测试通过!"
else
    print_warning "部分测试失败，但编译成功"
fi

# 安装
if [[ "$INSTALL" == true ]]; then
    print_info "安装Pono到系统..."
    sudo make install
    print_success "Pono已安装到系统"
fi

cd ..

print_success "Pono远端编译完成!"
print_info "可执行文件位置: ./build/pono"
print_info "使用示例:"
print_info "  ./build/pono --help"
print_info "  ./build/pono samples/counter-false.btor2"
