# Pono远端编译和运行 - 快速开始

## 🚀 一分钟快速开始

```bash
# 1. 使用预构建Docker镜像运行Pono
./run-pono.sh -- --help

# 2. 验证一个BTOR2文件
./run-pono.sh -v ./samples -- /workspace/constarrfalse.btor2

# 3. 使用不同的验证引擎
./run-pono.sh -v ./samples -- --engine bmc --bound 20 /workspace/uv_example.btor2
```

## 📁 文件说明

我已经为您创建了以下脚本和文档：

### 🔧 主要脚本

1. **`run-pono.sh`** - 智能运行脚本（推荐使用）
   - 自动检测环境并选择最佳运行方式
   - 优先使用本地编译版本，自动回退到Docker
   - 统一的命令行接口

2. **`docker-build-and-run.sh`** - Docker容器运行脚本
   - 支持使用预构建镜像或本地构建
   - 支持目录挂载和各种Docker选项
   - 适合容器化部署

3. **`remote-build.sh`** - 远端服务器编译脚本
   - 自动安装依赖并编译Pono
   - 支持多种Linux发行版
   - 可自定义编译选项

4. **`test-remote-build.sh`** - 测试脚本
   - 验证所有脚本和功能是否正常工作
   - 自动化测试流程

### 📚 文档文件

1. **`REMOTE_BUILD_GUIDE.md`** - 详细使用指南
   - 完整的使用说明和示例
   - 故障排除指南
   - 高级配置选项

2. **`QUICK_START.md`** - 本文件，快速开始指南

## 🎯 使用场景

### 场景1：快速验证模型
```bash
# 使用智能脚本，自动选择最佳运行方式
./run-pono.sh -v ./samples -- /workspace/your-model.btor2
```

### 场景2：开发和调试
```bash
# 在远端服务器编译，获得最佳性能
./remote-build.sh --debug
./build/pono --verbosity 3 samples/your-model.btor2
```

### 场景3：CI/CD集成
```bash
# 使用Docker确保环境一致性
./docker-build-and-run.sh -p -v ./models -- /workspace/test-model.btor2
```

### 场景4：批量验证
```bash
# 编写脚本批量处理多个模型
for model in models/*.btor2; do
    echo "验证: $model"
    ./run-pono.sh -v ./models -- "/workspace/$(basename "$model")"
done
```

## 🔍 验证结果说明

Pono的输出结果含义：
- **`sat`** - 找到反例，属性不成立
- **`unsat`** - 未找到反例，属性在给定边界内成立
- **`unknown`** - 无法确定结果（可能需要增加边界或使用不同引擎）

## ⚡ 性能优化建议

1. **选择合适的引擎：**
   - `bmc` - 适合寻找浅层反例
   - `ic3` - 适合证明属性成立
   - `kind` - 适合归纳证明

2. **调整参数：**
   - `--bound N` - 设置展开深度
   - `--verbosity N` - 控制输出详细程度
   - `--timeout N` - 设置超时时间

3. **使用本地编译：**
   - 如果需要频繁使用，建议编译本地版本以获得最佳性能

## 🆘 常见问题

**Q: Docker镜像下载慢怎么办？**
A: 使用国内镜像源或本地构建：
```bash
./docker-build-and-run.sh -b -- --help
```

**Q: 如何处理大型模型文件？**
A: 使用目录挂载并增加Docker内存限制：
```bash
./docker-build-and-run.sh -v ./large-models --memory=8g -- /workspace/large-model.btor2
```

**Q: 编译失败怎么办？**
A: 检查系统依赖并查看详细错误信息：
```bash
./remote-build.sh --clean 2>&1 | tee build.log
```

## 📞 获取帮助

- 查看详细文档：`cat REMOTE_BUILD_GUIDE.md`
- 运行测试脚本：`./test-remote-build.sh`
- 查看Pono帮助：`./run-pono.sh -- --help`
- 查看脚本帮助：`./run-pono.sh --help`

## 🎉 开始使用

现在您已经准备好开始使用Pono进行模型检查了！

```bash
# 验证您的第一个模型
./run-pono.sh -v ./samples -- /workspace/constarrfalse.btor2
```

祝您使用愉快！ 🚀
