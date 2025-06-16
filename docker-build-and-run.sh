#!/bin/bash

# Pono远端Docker编译和运行脚本
# 使用方法: ./docker-build-and-run.sh [选项]

set -e

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 打印带颜色的消息
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
Pono远端Docker编译和运行脚本

使用方法:
    $0 [选项] [-- pono参数]

选项:
    -h, --help          显示此帮助信息
    -b, --build         强制重新构建Docker镜像
    -p, --pull          从GitHub Container Registry拉取预构建镜像
    -t, --tag TAG       指定Docker镜像标签 (默认: pono-local)
    -v, --volume DIR    挂载本地目录到容器的/workspace
    --no-cache          构建时不使用缓存

示例:
    $0                                    # 构建并运行Pono
    $0 -p                                # 使用预构建镜像
    $0 -v ./samples                      # 挂载samples目录
    $0 -- --help                         # 显示Pono帮助信息
    $0 -- samples/counter-false.btor2    # 运行特定的BTOR2文件

EOF
}

# 默认参数
BUILD_IMAGE=false
PULL_IMAGE=false
IMAGE_TAG="pono-local"
VOLUME_DIR=""
NO_CACHE=""
PONO_ARGS=""

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_usage
            exit 0
            ;;
        -b|--build)
            BUILD_IMAGE=true
            shift
            ;;
        -p|--pull)
            PULL_IMAGE=true
            shift
            ;;
        -t|--tag)
            IMAGE_TAG="$2"
            shift 2
            ;;
        -v|--volume)
            VOLUME_DIR="$2"
            shift 2
            ;;
        --no-cache)
            NO_CACHE="--no-cache"
            shift
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

# 检查Docker是否安装
if ! command -v docker &> /dev/null; then
    print_error "Docker未安装。请先安装Docker。"
    exit 1
fi

# 检查Docker是否运行
if ! docker info &> /dev/null; then
    print_error "Docker未运行。请启动Docker服务。"
    exit 1
fi

print_info "开始Pono远端编译和运行流程..."

# 拉取预构建镜像
if [[ "$PULL_IMAGE" == true ]]; then
    print_info "从GitHub Container Registry拉取预构建镜像..."
    if docker pull ghcr.io/stanford-centaur/pono:latest; then
        docker tag ghcr.io/stanford-centaur/pono:latest "$IMAGE_TAG"
        print_success "成功拉取并标记预构建镜像为 $IMAGE_TAG"
    else
        print_warning "拉取预构建镜像失败，将构建本地镜像"
        BUILD_IMAGE=true
    fi
fi

# 检查镜像是否存在
if [[ "$BUILD_IMAGE" == true ]] || ! docker image inspect "$IMAGE_TAG" &> /dev/null; then
    print_info "构建Docker镜像: $IMAGE_TAG"
    print_info "这可能需要几分钟时间，请耐心等待..."
    
    if docker build $NO_CACHE -t "$IMAGE_TAG" .; then
        print_success "Docker镜像构建成功: $IMAGE_TAG"
    else
        print_error "Docker镜像构建失败"
        exit 1
    fi
else
    print_info "使用现有Docker镜像: $IMAGE_TAG"
fi

# 准备Docker运行参数
DOCKER_RUN_ARGS="--rm"
# 只在交互式终端时添加 -it 参数
if [[ -t 0 ]]; then
    DOCKER_RUN_ARGS="$DOCKER_RUN_ARGS -it"
fi

# 添加卷挂载
if [[ -n "$VOLUME_DIR" ]]; then
    if [[ ! -d "$VOLUME_DIR" ]]; then
        print_error "指定的目录不存在: $VOLUME_DIR"
        exit 1
    fi
    VOLUME_DIR=$(realpath "$VOLUME_DIR")
    DOCKER_RUN_ARGS="$DOCKER_RUN_ARGS -v $VOLUME_DIR:/workspace"
    print_info "挂载目录: $VOLUME_DIR -> /workspace"
fi

# 运行容器
print_info "启动Pono容器..."
print_info "Docker命令: docker run $DOCKER_RUN_ARGS $IMAGE_TAG $PONO_ARGS"

# 运行Docker容器并获取退出码
docker run $DOCKER_RUN_ARGS "$IMAGE_TAG" $PONO_ARGS
EXIT_CODE=$?

# 对于--help和--version命令，退出码可能非零但仍然是成功的
if [[ $EXIT_CODE -eq 0 ]] || [[ "$PONO_ARGS" == *"--help"* ]] || [[ "$PONO_ARGS" == *"--version"* ]]; then
    print_success "Pono执行完成"
else
    print_error "Pono执行失败 (退出码: $EXIT_CODE)"
    exit $EXIT_CODE
fi
