# Pono远端编译和运行指南

本指南提供了多种在远端服务器上编译和运行Pono的方法。

## 🚀 快速开始

### 方法1：一键智能运行（推荐）
```bash
# 自动检测环境并选择最佳运行方式
./run-pono.sh -- --help
./run-pono.sh -- samples/counter-false.btor2
```

### 方法2：Docker容器运行（最简单）
```bash
# 使用预构建镜像
./docker-build-and-run.sh -p -- --help

# 本地构建镜像
./docker-build-and-run.sh -b -- samples/counter-false.btor2
```

### 方法3：远端服务器直接编译
```bash
# 完整编译流程
./remote-build.sh
./build/pono --help
```

## 📋 详细说明

### 🐳 Docker方式（推荐）

**优点：**
- 环境隔离，不污染系统
- 依赖管理简单
- 跨平台兼容性好
- 可使用预构建镜像

**使用方法：**
```bash
# 基本使用
./docker-build-and-run.sh -- --help

# 使用预构建镜像（更快）
./docker-build-and-run.sh -p -- samples/counter-false.btor2

# 挂载本地目录
./docker-build-and-run.sh -v ./samples -- counter-false.btor2

# 强制重新构建
./docker-build-and-run.sh -b --no-cache -- --help
```

**选项说明：**
- `-p, --pull`: 从GitHub拉取预构建镜像
- `-b, --build`: 强制重新构建本地镜像
- `-v, --volume DIR`: 挂载本地目录到容器
- `--no-cache`: 构建时不使用缓存

### 🖥️ 远端服务器直接编译

**优点：**
- 性能最佳
- 可自定义编译选项
- 适合开发和调试

**系统要求：**
- Ubuntu 18.04+ / Debian 10+ / CentOS 7+ / Fedora 30+
- 至少4GB内存
- 至少2GB磁盘空间

**使用方法：**
```bash
# 基本编译
./remote-build.sh

# 使用8个并行任务编译
./remote-build.sh -j 8

# 编译调试版本
./remote-build.sh --debug

# 包含额外求解器支持
./remote-build.sh --with-btor

# 清理重新编译并安装
./remote-build.sh --clean --install
```

**选项说明：**
- `-j, --jobs N`: 并行编译任务数
- `--with-msat`: 包含MathSAT支持（需要许可证）
- `--with-btor`: 包含Boolector支持
- `--debug`: 编译调试版本
- `--clean`: 清理之前的构建
- `--install`: 安装到系统

### 🤖 智能运行脚本

**特点：**
- 自动检测运行环境
- 优先使用本地编译版本
- 自动回退到Docker方式
- 统一的命令行接口

**使用方法：**
```bash
# 显示帮助
./run-pono.sh -- --help

# 基本模型检查
./run-pono.sh -- samples/counter-false.btor2

# 强制使用Docker
./run-pono.sh -f -- samples/counter-false.btor2

# 使用Docker并挂载目录
./run-pono.sh -v ./samples -- counter-false.btor2

# 显示使用示例
./run-pono.sh --examples
```

## 🔧 Pono使用示例

### 基本验证
```bash
# 验证BTOR2文件
./run-pono.sh -- samples/counter-false.btor2

# 验证SMV文件
./run-pono.sh -- samples/simple_counter.smv
```

### 使用不同引擎
```bash
# BMC (Bounded Model Checking)
./run-pono.sh -- --engine bmc --bound 20 samples/counter-false.btor2

# IC3 (Incremental Construction of Inductive Clauses)
./run-pono.sh -- --engine ic3 samples/counter-false.btor2

# K-induction
./run-pono.sh -- --engine kind samples/counter-false.btor2
```

### 高级选项
```bash
# 生成反例
./run-pono.sh -- --witness samples/counter-false.btor2

# 详细输出
./run-pono.sh -- --verbosity 3 samples/counter-false.btor2

# 设置超时
./run-pono.sh -- --timeout 300 samples/counter-false.btor2
```

## 🐛 故障排除

### Docker相关问题

**问题：Docker镜像构建失败**
```bash
# 清理Docker缓存
docker system prune -a

# 重新构建
./docker-build-and-run.sh -b --no-cache
```

**问题：权限错误**
```bash
# 将用户添加到docker组
sudo usermod -aG docker $USER
# 重新登录或运行
newgrp docker
```

### 编译相关问题

**问题：依赖安装失败**
```bash
# 手动安装依赖（Ubuntu/Debian）
sudo apt-get update
sudo apt-get install -y build-essential cmake git libgmp-dev

# 手动安装依赖（CentOS/RHEL）
sudo yum install -y gcc gcc-c++ cmake git gmp-devel
```

**问题：内存不足**
```bash
# 减少并行任务数
./remote-build.sh -j 2

# 或者增加交换空间
sudo fallocate -l 2G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

### 运行时问题

**问题：找不到输入文件**
```bash
# 使用绝对路径
./run-pono.sh -- /full/path/to/file.btor2

# 或者挂载目录（Docker方式）
./run-pono.sh -v /path/to/files -- file.btor2
```

## 📚 更多资源

- [Pono官方文档](https://github.com/stanford-centaur/pono)
- [BTOR2格式说明](http://fmv.jku.at/btor2/)
- [SMT-Switch文档](https://github.com/stanford-centaur/smt-switch)

## 🤝 贡献

如果您发现问题或有改进建议，请提交Issue或Pull Request。

## 📄 许可证

本脚本遵循与Pono相同的BSD 3-Clause许可证。
