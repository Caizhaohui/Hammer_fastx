
# Hammer_fastx

Hammer_fastx 是一个使用 Rust 编写的高性能、多功能命令行工具集，用于处理 FASTA/FASTQ 文件。

## 简介

Hammer_fastx 旨在为生物信息学分析提供一个速度极快、内存占用低的序列文件（FASTA/FASTQ）处理方案。它利用了 Rust 的内存安全和并发特性，可以高效地处理大规模测序数据。目前，该工具集成了序列拆分（demux）、统计（stats）和过滤（filter）三个核心功能。

## 主要特性

- **多功能集成**: 在一个可执行文件中提供了多个实用的子命令。
- **极致性能**: 基于 Rust 和多线程（针对 demux 功能）构建，处理速度远超脚本语言。
- **格式兼容**: 智能识别并处理 FASTA 和 FASTQ 格式，同时支持 `.gz` 压缩文件。
- **清晰的命令行接口**: 采用标准的子命令结构，用法清晰，易于上手。
- **跨平台**: 可在 Linux, macOS 和 Windows 上编译和运行。

## 安装

确保您已经安装了 Rust 工具链。然后，您可以使用 `cargo` 直接从源码安装 Hammer_fastx：

```bash
# 克隆本仓库
git clone https://github.com/Caizhaohui/Hammer_fastx.git

# 进入项目目录
cd Hammer_fastx

# 使用 cargo install 进行安装
cargo install --path .
markdown

安装完成后，Hammer_fastx 命令将在您的终端中全局可用。
使用方法

Hammer_fastx 的使用遵循 Hammer_fastx <子命令> [参数] 的标准格式。
1. demux - 序列拆分

根据双端 barcode 将一个 FASTQ 文件拆分成多个样本文件。
用法:

bash
Hammer_fastx demux --inputfile <输入文件> --output <输出目录> --tags <标签文件> [其他选项]

参数:

--inputfile <路径>: 必需。输入的 FASTQ 文件，可为 .gz 格式。
--output <目录>: 必需。存放拆分后文件的输出目录。
--tags <路径>: 必需。包含样本和 barcode 的 CSV/TSV 文件。文件需包含 SampleID, F_tag, R_tag 这三列。
--threads <数字>: (可选) 指定使用的线程数，默认为机器的物理核心数。
--tag-len <数字>: (可选) 指定 barcode 的长度，默认为 8。
--trim: (可选) 激活此选项以在拆分后裁剪序列两端的 barcode。
--out-fasta: (可选) 将输出文件格式转换为 FASTA。

示例:

bash
Hammer_fastx demux \
    --inputfile ./raw_data/all_reads.fastq.gz \
    --output ./demux_results \
    --tags ./metadata/barcodes.csv \
    --threads 16 \
    --trim

2. stats - 序列统计

快速计算一个 FASTA 或 FASTQ 文件的基本统计信息。
用法:

bash
Hammer_fastx stats --inputfile <输入文件>

参数:

--inputfile <路径>: 必需。输入的 FASTA 或 FASTQ 文件，可为 .gz 格式。

示例:

bash
Hammer_fastx stats --inputfile ./sequences.fasta.gz

输出示例:

==================== 统计结果 ====================
序列总数 (Total sequences): 150000
总碱基数 (Total bases):    75000000
最大序列长度 (Max length):   501
最小序列长度 (Min length):   498
平均序列长度 (Avg length):   500.00
================================================

3. filter - 序列过滤

根据序列的最小和/或最大长度来过滤文件。
用法:

bash
Hammer_fastx filter --inputfile <输入文件> [过滤选项] [输出选项]

参数:

--inputfile <路径>: 必需。输入的 FASTA 或 FASTQ 文件，可为 .gz 格式。
-m, --min-len <长度>: (可选) 保留长度大于或等于该值的序列。
-M, --max-len <长度>: (可选) 保留长度小于或等于该值的序列。
--outfile <路径>: (可选) 将过滤结果输出到指定文件。如果省略，结果将直接输出到屏幕（标准输出）。

示例:

保留长度在 100bp 到 500bp 之间的序列，并输出到新文件:

bash
Hammer_fastx filter --inputfile reads.fastq -m 100 -M 500 --outfile filtered_reads.fastq

只保留长度大于等于 300bp 的序列，并输出到屏幕:

bash
Hammer_fastx filter --inputfile my.fasta -m 300

使用管道符进行高级操作（过滤后直接压缩）:

bash
Hammer_fastx filter --inputfile my.fasta -m 100 | gzip > filtered.fasta.gz

从源码构建

bash
# 克隆仓库
git clone https://github.com/Caizhaohui/Hammer_fastx.git
cd Hammer_fastx

# 以发布模式编译
cargo build --release

# 编译好的可执行文件位于 ./target/release/Hammer_fastx

许可证

本项目采用 MIT License。

这次我确保了后续部分保持了完整的 Markdown 格式，方便您进行粘贴。
