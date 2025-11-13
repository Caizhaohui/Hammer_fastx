# Hammer_fastx

> 用 Rust 语言编写的高性能 FASTA/FASTQ 文件处理工具套件

## 项目概述

Hammer_fastx 是一个用 Rust 开发的轻量级命令行工具套件，专为高效处理生物信息学中的 FASTA 和 FASTQ 文件而设计。该工具集成了从原始测序数据的质控、合并、拆分到序列统计和特殊分析的完整工作流程，适用于基因组测序数据的预处理、格式转换、质控、合并及拆分等多种场景。

### 主要特点

- **高性能**：利用 Rust 语言的内存安全性和高性能特性，结合多线程并行处理，大幅提升数据处理速度
- **完整工作流**：集成质控、合并、拆分等多个步骤，支持一键式全流程处理
- **多线程优化**：使用 Rayon 库实现高效并行计算，充分利用多核 CPU 资源
- **内存友好**：采用分块处理策略，有效处理超大数据文件而不占用过多内存
- **功能丰富**：提供多种子命令，涵盖序列统计、长度过滤、翻译和突变分析等功能
- **灵活配置**：支持丰富的命令行参数，满足不同的处理需求

### 典型应用场景

- 高通量测序数据的质控和预处理
- 双端测序数据的合并与拆分
- 基于 Barcode 的样本多路分解
- 大规模序列文件的统计分析
- 序列长度过滤和筛选
- DNA 到氨基酸序列的翻译
- 蛋白质序列突变分析

## 安装指南

### 系统要求

- 操作系统：Linux, macOS, Windows (WSL 推荐)
- Rust 编译环境：1.60 或更高版本
- 外部依赖工具：部分子命令需要 `fastp` 和 `flash2`（详见下方说明）

### 步骤 1: 安装 Rust

如果您的系统上尚未安装 Rust，可以使用以下命令安装：

```bash
# Linux/macOS
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# 或者手动下载并安装
# 安装完成后，确保将 Rust 添加到环境变量
```

### 步骤 2: 克隆项目

```bash
git clone https://github.com/Caizhaohui/Hammer_fastx.git
cd Hammer_fastx
```

### 步骤 3: 编译安装

```bash
# 编译项目
cargo build --release

# 可执行文件将位于 target/release/ 目录下
# 可以将其复制到系统路径中方便使用
cp target/release/hammer_fastx /usr/local/bin/
```

### 步骤 4: 安装外部依赖

部分子命令依赖于外部工具：

#### 安装 fastp

```bash
# 使用 conda 安装
easy_install pip
pip install fastp

# 或者从源码编译
git clone https://github.com/OpenGene/fastp.git
cd fastp
make
cp fastp /usr/local/bin/
```

#### 安装 flash2

```bash
# 使用 conda 安装
conda install -c bioconda flash2

# 或者从源码编译
git clone https://github.com/dstreett/FLASH2.git
cd FLASH2
make
cp flash /usr/local/bin/
```

### 步骤 5: 验证安装

```bash
# 验证 hammer_fastx 安装成功
hammer_fastx --version

# 验证外部依赖安装成功
fastp --version
flash --version
```

## 子命令与功能详细说明

### 1. demux_all

**[主流程] 一键运行质控(fastp)、合并(flash2)、barcode拆分**

- **功能**：执行完整的数据处理流程，包括质控、合并和样本拆分
- **输入**：成对的原始双端 fastq 文件、barcode 文件等
- **输出**：根据 barcode 拆分后的样本文件
- **常用参数**：
  - `--in1`/`--in2`：输入的 R1/R2 fastq 文件
  - `--barcode`：barcode 文件路径
  - `--outdir`：输出目录
  - `--fastp-threads`：fastp 线程数
  - `--flash2-threads`：flash2 线程数
  - `--out-fasta`：输出为 FASTA 格式（默认 FASTQ）
  - `--min-length`：合并后序列的最小长度阈值

### 2. mergePE

**[流程] 先用 fastp 质控，再用 flash2 合并双端数据**

- **功能**：对双端测序数据进行质控后合并
- **输入**：成对的原始双端 fastq 文件
- **输出**：质控并合并后的序列文件
- **常用参数**：
  - `--in1`/`--in2`：输入的 R1/R2 fastq 文件
  - `--out1`/`--out2`：输出的质控后 R1/R2 文件
  - `--merged`：输出的合并后文件
  - `--threads`：总线程数
  - `--out-fasta`：输出为 FASTA 格式（默认 FASTQ）
  - `--html`/`--json`：fastp 报告输出路径

### 3. demux_only

**[单步骤] 仅基于 barcode 对已合并 fastq 文件进行拆分**

- **功能**：快速根据 barcode 对合并后的序列进行样本拆分
- **输入**：合并后的 fastq 文件和 barcode 文件
- **输出**：根据 barcode 拆分的样本文件
- **常用参数**：
  - `--infile`：合并后的 fastq 文件
  - `--barcode`：barcode 文件
  - `--outdir`：结果输出目录
  - `--threads`：线程数
  - `--allow-mismatch`：允许的 mismatch 数量
  - `--min-length`：输出序列的最小长度

### 4. Fastp

**调用 fastp 对双端 fastq 文件进行质控**

- **功能**：包装外部 fastp 工具，提供统一的命令行接口
- **输入**：双端 fastq 文件
- **输出**：质控后的 fastq 文件和质控报告
- **常用参数**：
  - `-i/--in1`：输入文件1 (Read1)
  - `-I/--in2`：输入文件2 (Read2)
  - `-o/--out1`：输出文件1 (Read1)
  - `-O/--out2`：输出文件2 (Read2)
  - `-h/--html`：指定 HTML 报告输出路径
  - `-j/--json`：指定 JSON 报告输出路径
  - `-R/--report-title`：报告标题
  - `-t/--threads`：线程数
  - `--cut_front`：前端裁剪长度
  - `--cut_tail`：末端裁剪长度

### 5. Flash2

**调用 flash2 合并双端 reads**

- **功能**：包装外部 flash2 工具，用于合并双端测序数据
- **输入**：成对的 fastq 文件
- **输出**：合并后的序列文件
- **常用参数**：
  - `-i/--in1`：输入文件1 (Read1)
  - `-I/--in2`：输入文件2 (Read2)
  - `-o/--out`：输出文件前缀
  - `-t/--threads`：线程数
  - `-M/--max-overlap`：最大重叠长度
  - `-m/--min-overlap`：最小重叠长度
  - `--max-mismatch-density`：最大错配密度

### 6. Stats

**统计FASTA/FASTQ文件中的序列信息**

- **功能**：计算序列文件的基本统计信息，包括序列数、总碱基数、平均长度等
- **输入**：一个或多个 FASTA/FASTQ 文件
- **输出**：统计信息表格
- **常用参数**：
  - `-i/--input`：输入文件或目录
  - `-o/--output`：输出文件（可选）

### 7. Filter

**根据序列长度过滤FASTA/FASTQ文件**

- **功能**：根据指定的长度范围过滤序列
- **输入**：一个或多个 FASTA/FASTQ 文件
- **输出**：过滤后的序列文件
- **常用参数**：
  - `-i/--input`：输入文件或目录
  - `-o/--output`：输出文件或目录
  - `--min-length`：最小序列长度（默认：0）
  - `--max-length`：最大序列长度（默认：无限制）
  - `--threads`：线程数

### 8. Ns_count

**将 reads 比对到含 N 的参考序列并提取组合**

- **功能**：将序列比对到包含 N 碱基的参考序列，分析 N 区域的碱基组成
- **输入**：reads 文件和参考序列文件
- **输出**：比对结果统计
- **常用参数**：
  - `--reads`：输入 reads 文件
  - `--ref-seq`：参考序列文件
  - `--output`：输出文件
  - `--threads`：线程数

### 9. dna2aa（新增）

**DNA 序列转氨基酸序列翻译工具**

- **功能**：将 FASTA 格式的 DNA 序列翻译为氨基酸序列
- **输入**：FASTA 格式的 DNA 序列文件或目录
- **输出**：翻译后的氨基酸序列文件
- **常用参数**：
  - `--input`：输入文件或目录
  - `--output`：输出目录
  - `--aa-length`：目标氨基酸序列长度（可选）

### 10. count_aa（新增）

**比较序列与参考蛋白序列的氨基酸突变统计**

- **功能**：统计序列与参考蛋白序列的氨基酸突变情况
- **输入**：参考蛋白序列文件和样本 FASTA 文件目录
- **输出**：突变统计 CSV 文件
- **常用参数**：
  - `-r/--reference`：参考蛋白 FASTA 序列
  - `-i/--input-dir`：包含多个 FASTA 文件的目录
  - `-o/--output-dir`：输出 CSV 文件的目录
  - `-A/--aa-offset`：位置偏移量（默认：0）
  - `-c/--config`：CSV 配置文件，包含 protected_sites 列
  - `--match_len`：匹配参考起始位置的氨基酸数量（默认：6）
  - `--threads`：每个文件的并行线程数（默认：8）
  - `--chunk_size`：每块 reads 数量（默认：100000）

## 使用示例

### 1. 完整工作流示例

```bash
# 使用 demux_all 执行完整流程
hammer_fastx demux_all \
    --in1 raw_data/sample_R1.fastq.gz \
    --in2 raw_data/sample_R2.fastq.gz \
    --barcode barcodes.txt \
    --outdir results \
    --fastp-threads 8 \
    --flash2-threads 8 \
    --out-fasta

# 参数说明：
# --in1/--in2: 输入的双端测序文件
# --barcode: 包含样本barcode信息的文件
# --outdir: 输出目录
# --fastp-threads/--flash2-threads: 指定工具使用的线程数
# --out-fasta: 输出FASTA格式而不是默认的FASTQ
```

### 2. 分步处理示例

```bash
# 步骤1: 使用Fastp进行质控
hammer_fastx Fastp \
    --in1 raw_data/sample_R1.fastq.gz \
    --in2 raw_data/sample_R2.fastq.gz \
    --out1 clean_data/sample_R1.clean.fq.gz \
    --out2 clean_data/sample_R2.clean.fq.gz \
    --html reports/fastp_report.html \
    --json reports/fastp_report.json \
    --threads 8

# 步骤2: 使用Flash2合并质控后的双端数据
hammer_fastx Flash2 \
    --in1 clean_data/sample_R1.clean.fq.gz \
    --in2 clean_data/sample_R2.clean.fq.gz \
    --out merged/sample \
    --threads 8

# 步骤3: 使用demux_only拆分合并的数据
hammer_fastx demux_only \
    --infile merged/sample.extendedFrags.fastq \
    --barcode barcodes.txt \
    --outdir demultiplexed \
    --threads 12
```

### 3. 统计分析示例

```bash
# 统计序列文件信息
hammer_fastx Stats \
    --input demultiplexed/*.fasta \
    --output stats.csv

# 根据序列长度过滤文件
hammer_fastx Filter \
    --input demultiplexed \
    --output filtered \
    --min-length 200 \
    --max-length 1000 \
    --threads 8
```

### 4. 新增功能使用示例

```bash
# DNA 转氨基酸序列
hammer_fastx dna2aa \
    --input dna_sequences/ \
    --output protein_sequences/ \
    --aa-length 100

# 氨基酸突变统计
hammer_fastx count_aa \
    --reference reference_protein.fasta \
    --input-dir protein_sequences/ \
    --output-dir mutation_stats/ \
    --aa-offset 1 \
    --config protected_sites.csv \
    --threads 12 \
    --chunk_size 500000
```

### 5. 常见任务组合示例

```bash
# 批量处理多个样本的完整流程
for sample in sampleA sampleB sampleC; do
    hammer_fastx demux_all \
        --in1 raw_data/${sample}_R1.fastq.gz \
        --in2 raw_data/${sample}_R2.fastq.gz \
        --barcode barcodes.txt \
        --outdir results/${sample} \
        --threads 8
        
    echo "处理完成: $sample"
done

# 对处理后的样本进行统计分析
hammer_fastx Stats \
    --input results/*/*.fasta \
    > all_samples_stats.txt
```

## 输入输出格式

### 支持的输入格式

- **FASTA 格式**：`.fasta`, `.fa`, `.fna`
- **FASTQ 格式**：`.fastq`, `.fq`
- **压缩格式**：支持 `.gz` 压缩的 FASTA/FASTQ 文件
- **Barcode 文件**：制表符分隔的文本文件，格式：`barcode\t样本名`

### 输出格式

- **序列文件**：FASTA 或 FASTQ 格式（根据命令行参数选择）
- **统计报告**：文本表格或 CSV 格式
- **质控报告**：HTML 和 JSON 格式（由 fastp 生成）
- **突变分析**：CSV 格式，包含突变位点和计数

## 依赖说明

### Rust 内部依赖

项目使用以下 Rust 库：

- **bio**: 生物序列处理库
- **clap**: 命令行参数解析
- **csv**: CSV 文件处理
- **rayon**: 并行处理
- **crossbeam-channel**: 线程通信
- **dashmap**: 并发哈希映射
- **indicatif**: 命令行进度显示
- **flate2**: 压缩文件处理
- **anyhow**: 错误处理
- **glob**: 文件匹配
- **num_cpus**: CPU 核心数检测

### 外部工具依赖

部分子命令需要以下外部工具：

- **fastp**: 用于 FASTQ 文件质控
- **flash2**: 用于合并双端 reads

请确保这些工具已正确安装并添加到系统路径中。

## 性能优化建议

1. **线程数设置**：根据系统 CPU 核心数调整 `--threads` 参数，通常设置为 CPU 核心数的 1/2 到全部核心数

2. **内存管理**：对于超大文件，可调整 `--chunk_size` 参数（特别是 count_aa 子命令），减少内存占用

3. **并行策略**：
   - 对于多个文件，使用多线程并行处理每个文件
   - 对于单个大文件，内部会自动分块并行处理

4. **I/O 优化**：将输入输出文件放在不同的物理磁盘上，以减少 I/O 竞争

## 贡献指南

欢迎对项目进行贡献！如果您想参与开发，请遵循以下步骤：

1. **Fork 仓库**：在 GitHub 上 Fork 项目仓库到您自己的账户

2. **克隆仓库**：
   ```bash
   git clone https://github.com/您的用户名/Hammer_fastx.git
   cd Hammer_fastx
   ```

3. **创建分支**：
   ```bash
   git checkout -b feature/您的功能名称
   ```

4. **进行开发**：实现新功能或修复 bug

5. **测试**：确保您的代码可以正常编译和运行
   ```bash
   cargo build
   cargo test
   ```

6. **提交更改**：
   ```bash
   git add .
   git commit -m "添加描述性的提交信息"
   ```

7. **推送到远程**：
   ```bash
   git push origin feature/您的功能名称
   ```

8. **创建 Pull Request**：在 GitHub 上从您的分支创建 Pull Request 到主仓库

### 代码规范

- 遵循 Rust 的代码风格（使用 `rustfmt` 格式化代码）
- 为新功能添加适当的文档和注释
- 确保新代码有良好的错误处理

## 问题反馈

如有任何问题、建议或发现 bug，请在 [GitHub Issues](https://github.com/Caizhaohui/Hammer_fastx/issues) 页面提交。提交问题时，请包含：

- 问题描述
- 复现步骤
- 错误信息（如果有）
- 您的系统环境信息

## 许可证

本项目采用 MIT License。详见 [LICENSE](./LICENSE) 文件。

## 作者

**Caizhaohui**

- GitHub: [https://github.com/Caizhaohui](https://github.com/Caizhaohui)

---

*Hammer_fastx - 让生物信息学数据分析更加高效！*
