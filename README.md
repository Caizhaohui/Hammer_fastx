# Hammer_fastx

一个使用 Rust 编写的，用于处理 FASTA/FASTQ 文件的高性能、多功能命令行工具集。

---

## 简介

`Hammer_fastx` 旨在为生物信息学分析提供一个速度极快、内存占用低的序列文件（FASTA/FASTQ）处理方案。它利用了 Rust 的内存安全和并发特性，可以高效地处理大规模测序数据。

目前，该工具集成了序列拆分（demux）、统计（stats）和过滤（filter）三个核心功能。

## 主要特性

* **多功能集成**: 在一个可执行文件中提供了多个实用的子命令。
* **极致性能**: 基于 Rust 和多线程（针对 `demux` 功能）构建，处理速度远超脚本语言。
* **格式兼容**: 智能识别并处理 FASTA 和 FASTQ 格式，同时支持 `.gz` 压缩文件。
* **清晰的命令行接口**: 采用标准的子命令结构，用法清晰，易于上手。
* **跨平台**: 可在 Linux, macOS 和 Windows 上编译和运行。

## 安装

确保您已经安装了 [Rust 工具链](https://www.rust-lang.org/tools/install)。然后，您可以使用 `cargo` 直接从源码安装 `Hammer_fastx`：

```bash
# 首先，克隆本仓库
git clone [https://github.com/Caizhaohui/Hammer_fastx.git](https://github.com/Caizhaohui/Hammer_fastx.git)

# 进入项目目录
cd Hammer_fastx

# 使用 cargo install 进行安装
cargo install --path .
