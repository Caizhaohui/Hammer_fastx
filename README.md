# Hammer_fastx

> 用Rust语言编写的一个处理FASTA和FASTQ文件的小工具

## 项目简介

Hammer_fastx 是一个用 Rust 开发的轻量级命令行工具，旨在高效处理生物信息学中的 FASTA 和 FASTQ 文件。适用于基因组数据的预处理、格式转换等常见场景。

## 功能特性

- 支持 FASTA 和 FASTQ 格式的读取与处理
- 高性能，适合大规模数据
- 跨平台支持（Linux/macOS/Windows）
- 简单易用的命令行界面

## 安装方法

1. **通过 cargo 构建安装**  
   请确保已安装 Rust 环境（推荐使用 [rustup](https://rustup.rs/)）：

   ```bash
   git clone https://github.com/Caizhaohui/Hammer_fastx.git
   cd Hammer_fastx
   cargo build --release
   ```

2. **直接运行**  
   或者在项目目录下直接运行：

   ```bash
   cargo run -- <参数>
   ```

## 使用示例

假设你已经编译完成，在 `target/release` 目录下有可执行文件：

```bash
./hammer_fastx --help
```

常见命令示例：

```bash
# 读取 fasta 文件并统计序列数
./hammer_fastx stat -i example.fasta

# fastq 转换为 fasta
./hammer_fastx fq2fa -i example.fastq -o output.fasta
```

> 具体支持的参数和子命令请使用 `--help` 查看。

## 输入输出格式

- **输入**：支持标准 FASTA (.fasta, .fa) 和 FASTQ (.fastq, .fq) 文件
- **输出**：可生成 FASTA/FASTQ 格式文件，或命令行统计信息

## 依赖环境

- Rust 1.60 及以上

## 贡献和反馈

欢迎 issue 和 PR！如有建议或 bug 可在 [GitHub Issues](https://github.com/Caizhaohui/Hammer_fastx/issues) 提出。

## 许可证

本项目采用 MIT License，详见 [LICENSE](./LICENSE)。

---

**作者**：[Caizhaohui](https://github.com/Caizhaohui)
