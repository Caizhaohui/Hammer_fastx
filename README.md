```markdown
# Hammer_fastx

[![Version](https://img.shields.io/badge/version-v0.7.1-blue)](https://github.com/your-repo/hammer_fastx)

`Hammer_fastx` 是一个专为处理高通量测序数据（FASTA/FASTQ）而设计的多功能命令行工具集。它集成了质控、合并、拆分、统计和分析等多种功能，旨在简化从原始数据到下游分析的完整工作流程。

---

## ������ 功能概览

`Hammer_fastx` 提供了多个子命令，覆盖了常见的生物信息学处理步骤：

### ������ 主要流程
- **`demux_all`**: 一键式完整流程，自动执行 `fastp` 质控 → `flash2` 合并 → `demux_only` 拆分。
- **`mergePE`**: 将双端测序数据（Paired-End）进行质控后合并为单条序列。

### ������ 单步工具
- **`fastp`**: 包装 `fastp` 工具，对双端 FASTQ 文件进行快速质控。
- **`flash2`**: 包装 `flash2` 工具，将双端 reads 合并为单条序列。
- **`demux_only`**: 根据用户提供的标签文件，将已合并的 FASTQ 文件拆分为多个样本。
- **`stats`**: 统计一个或多个 FASTA/FASTQ 文件的序列信息（数量、长度分布等）。
- **`filter`**: 根据序列长度过滤 FASTA/FASTQ 文件。
- **`Ns_count`**: 将 reads 比对到含有 N 区域的参考序列上，提取并统计 N 区域的组合（Haplotype）。

---

## ������ 快速开始

### 1. 安装

`Hammer_fastx` 是一个 Rust 程序，您需要先安装 [Rust 工具链](https://www.rust-lang.org/tools/install)。

```bash
# 克隆仓库
git clone https://github.com/your-repo/hammer_fastx.git
cd hammer_fastx

# 编译并安装
cargo install --path .
```

或者，您可以直接下载预编译的二进制文件（如果可用）。

### 2. 依赖

本工具依赖以下外部软件，请确保它们已安装并位于您的 `PATH` 环境变量中：
- [`fastp`](https://github.com/OpenGene/fastp)
- [`flash2`](https://github.com/jaekss/flash2) (或原版 `FLASH`)

### 3. 查看帮助

```bash
# 查看所有子命令
hammer_fastx --help

# 查看特定子命令的帮助
hammer_fastx demux_all --help
```

---

## ������️ 使用示例

### 示例 1: 运行完整拆分流程 (`demux_all`)

此命令将自动完成质控、合并和样本拆分。

```bash
hammer_fastx demux_all \
  -i raw_reads/R1.fastq.gz \
  -I raw_reads/R2.fastq.gz \
  --tags sample_tags.csv \
  -o ./analysis_results \
  --cleanup
```

### 示例 2: 合并双端数据 (`mergePE`)

将质控后的双端数据合并为单条序列。

```bash
hammer_fastx mergePE \
  -i cleaned_R1.fastq.gz \
  -I cleaned_R2.fastq.gz \
  -o merged.fasta \
  --out_fasta \
  --cleanup
```

### 示例 3: 拆分已合并的文件 (`demux_only`)

根据标签文件拆分一个已合并的 FASTQ 文件。

```bash
hammer_fastx demux_only \
  --inputfile merged_reads.fastq \
  --output ./demux_results \
  --tags sample_tags.csv \
  --tag-len 8 \
  --trim \
  --out_fasta
```

### 示例 4: 统计文件信息 (`stats`)

统计多个文件的序列信息。

```bash
hammer_fastx stats --inputfile *.fastq.gz
```

### 示例 5: 提取 N 区域组合 (`Ns_count`)

比对 reads 到带有 N 的参考序列，并统计 N 区域的组合频率。

```bash
hammer_fastx Ns_count \
  --reads input_reads.fasta \
  --refSEQ reference_with_Ns.fasta \
  --output ./ns_results \
  --group "SampleGroup" \
  --extract_matches
```

---

## ������ 输入文件格式

### 样本标签文件 (`--tags`)

用于 `demux_only` 和 `demux_all` 的 CSV 文件，必须包含以下三列：

| SampleID | F_tag | R_tag |
| :--- | :--- | :--- |
| Sample1 | ATGCATGC | TCGATCGA |
| Sample2 | GCGCGCGC | ATATATAT |

- `F_tag`: 正向引物/标签序列。
- `R_tag`: 反向引物/标签序列。
- 程序会自动处理反向互补。

---

## ������ 输出说明

- **`demux_all`**: 输出目录包含 `01_fastp_out`, `02_flash2_out`, `03_demux_out` 三个子目录。
- **`demux_only`**: 在指定输出目录中为每个样本（包括 `unmatched`）生成单独的 FASTQ/FASTA 文件。
- **`Ns_count`**: 为每个含有 N 的参考序列生成一个 CSV 文件，报告所有观察到的组合及其频率。

---

## ⚙️ 开发

### 构建

```bash
cargo build --release
```

构建后的可执行文件位于 `target/release/hammer_fastx`。

### 版本更新

当前版本: `v0.7.1`
- **主要修复**: 彻底修复了 `Ns_count` 子命令因序列长度判断不当导致的 panic bug。

---

## ������ 致谢

- 本项目由 CZH 开发，并借助了 Google Gemini 的辅助。
- 感谢 `fastp`, `flash2`, `bio` crate 等开源项目。

---

## ������ 许可

本项目采用 MIT 许可证。
```
