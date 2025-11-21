# Hammer_fastx

> 用 Rust 编写的高性能 FASTA/FASTQ 处理工具套件

## 概述

Hammer_fastx 是一个面向 FASTA/FASTQ 的命令行工具集，提供从原始测序数据的质控、合并、样本拆分到统计分析、长度过滤、DNA→AA 翻译、蛋白突变统计，以及基序查找与上下游片段提取等功能。所有命令均支持大文件、`.gz` 压缩，以及合理的错误提示。

### 特点

- 高性能与内存安全，适合大规模数据处理
- 完整工作流（质控→合并→拆分）与丰富单步工具
- 并行优化与流式读取，支持 `.gz` 压缩
- 输出格式清晰统一，便于下游分析

## 安装

### 环境要求

- Linux/macOS/WSL
- Rust 1.60+（建议使用最新版）
- 外部工具：部分命令需要 `fastp` 与 `flash2`

### 安装步骤

```bash
# 安装 Rust（Linux/macOS）
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# 获取并编译
git clone https://github.com/Caizhaohui/Hammer_fastx.git
cd Hammer_fastx
cargo build --release
cp target/release/hammer_fastx /usr/local/bin/
```

### 外部依赖

### 1. 完整工作流示例

```bash
# fastp（建议使用 conda）
conda install -c bioconda fastp

# flash2（建议使用 conda）
conda install -c bioconda flash2
```

## 命令总览

```bash
hammer_fastx --help
```

- `demux_all`：完整流程（fastp→flash2→demux）
- `mergePE`：质控并合并双端数据
- `demux_only`：对已合并 FASTQ 进行样本拆分
- `fastp`：包装 fastp 进行质控
- `flash2`：包装 flash2 进行合并
- `stats`：统计 FASTA/FASTQ 基本信息
- `filter`：按长度过滤（支持批量或拼接）
- `Ns_count`：N 区域锚定比对与组合统计
- `DNA2AA`：DNA FASTA 批量翻译到 AA FASTA
- `count_AA`：参考蛋白突变统计（并行）
- `find_seq`：查找基序并提取上下游片段（支持反向互补）

---

## 子命令详解

### demux_all（完整流程）

- 功能：一键执行 `fastp` 质控、`flash2` 合并、`demux_only` 样本拆分
- 输入：`--in1`/`--in2` 原始双端 FASTQ、`--tags` 样本标签 CSV
- 输出：`01_fastp_out/`、`02_flash2_out/`、`03_demux_out/`
- 常用参数：
  - `--in1`、`--in2`：原始 R1/R2 FASTQ
  - `--tags`：样本标签 CSV（格式见下）
  - `-o/--output-dir`：主输出目录
  - `--cleanup`：流程成功后删除中间文件
  - `--fastp-threads`、`--flash-threads`：fastp/flash2 线程数
  - `--min-overlap`、`--max-overlap`：flash2 合并重叠范围
  - `--demux-threads`、`-l/--tag-len`、`--trim`、`--out-fasta`：拆分阶段参数
- 使用示例：
```bash
hammer_fastx demux_all \
  --in1 raw/sample_R1.fastq.gz \
  --in2 raw/sample_R2.fastq.gz \
  --tags metadata/tags.csv \
  --output-dir results/sampleA \
  --fastp-threads 8 --flash-threads 8 \
  --min-overlap 10 --max-overlap 300 \
  --demux-threads 12 -l 8 --trim --out-fasta
```
- 样本标签 CSV 示例：
```csv
SampleID,F_tag,R_tag
S1,ACGTACGT,TGCATGCA
S2,AAAACCCC,GGGGTTTT
```
- 输出示例结构：
  - `results/sampleA/01_fastp_out/filtered_R1.fastq.gz`、`filtered_R2.fastq.gz`
  - `results/sampleA/02_flash2_out/merged.extendedFrags.fastq`
  - `results/sampleA/03_demux_out/<SampleID>.fasta`

### mergePE（质控并合并）

- 功能：对双端测序数据先质控后合并，得到最终输出（FASTA/FASTQ）
- 参数：`-i/--in1`、`-I/--in2`、`-o/--outfile`、`--out-fasta`、`--cleanup`、`--temp-dir`、`--fastp-threads`、`--flash-threads`、`--min-overlap`、`--max-overlap`
- 使用示例：
```bash
hammer_fastx mergePE \
  --in1 raw/R1.fastq.gz --in2 raw/R2.fastq.gz \
  --outfile merged/merged.fastq \
  --fastp-threads 8 --flash-threads 8 \
  --min-overlap 10 --max-overlap 300
```
- 输出：`merged.fastq` 或 `merged.fasta`（取决于 `--out-fasta`）
- FASTA 输出示例：
```bash
hammer_fastx mergePE \
  --in1 raw/R1.fastq.gz --in2 raw/R2.fastq.gz \
  --outfile merged/merged.fasta \
  --out-fasta
```

### demux_only（样本拆分）

- 功能：对已合并的 FASTQ 根据样本标签进行拆分
- 参数：`--inputfile`、`--output`、`--threads`、`--tags`、`-l/--tag-len`、`--trim`、`--out-fasta`
- 使用示例：
```bash
hammer_fastx demux_only \
  --inputfile merged/extendedFrags.fastq \
  --tags metadata/tags.csv \
  --output demux_out \
  --threads 12 -l 8 --trim --out-fasta
```
- 输出：`demux_out/SampleID.(fastq|fasta)`
- Barcode 文件格式（制表符分隔）也支持：
```text
ACGTACGT\tS1
AAAACCCC\tS2
```

### fastp（质控包装）

- 功能：调用 `fastp` 对双端 FASTQ 进行质控并输出报告
- 参数：`-i/--in1`、`-I/--in2`、`-o/--out1`、`-O/--out2`、`-h/--html`、`-j/--json`、`-R/--report-title`、`-t/--threads`
- 使用示例：
```bash
hammer_fastx fastp \
  --in1 raw/R1.fastq.gz --in2 raw/R2.fastq.gz \
  --out1 clean/R1.clean.fastq.gz --out2 clean/R2.clean.fastq.gz \
  --html reports/fastp.html --json reports/fastp.json \
  --threads 8
```
- 输出：质控后的 FASTQ 与 HTML/JSON 报告

### flash2（合并包装）

- 功能：调用 `flash2` 合并双端 reads
- 参数：`-i/--in1`、`-I/--in2`、`-o/--out-prefix`、`-d/--out-dir`、`-m/--min-overlap`、`-M/--max-overlap`、`-t/--threads`
- 使用示例：
```bash
hammer_fastx flash2 \
  --in1 clean/R1.clean.fastq.gz --in2 clean/R2.clean.fastq.gz \
  --out-prefix merged --out-dir flash_out \
  --min-overlap 10 --max-overlap 300 --threads 8
```
- 输出：`flash_out/merged.extendedFrags.fastq` 等标准 flash2 文件

### stats（文件统计）

- 功能：统计 FASTA/FASTQ 基本信息（序列数、总碱基数、最短/最长、平均长度），并可将“序列种类与数量”按降序导出到 CSV
- 参数：
  - `--inputfile <files...>`：一个或多个输入文件（支持 `.gz`）
  - `--outfile <path>`：将每个唯一序列的计数导出为 CSV，列为 `filename,sequence,count`
  - 规范化：导出时会对序列做大小写归一（转大写）与首尾空白去除；结果按 `count` 降序排列，计数相同时按 `sequence` 升序
- 使用示例：
```bash
# 仅打印总体统计到标准输出
hammer_fastx stats --inputfile samples/a.fasta samples/b.fastq.gz

# 导出序列计数到 CSV
hammer_fastx stats --inputfile samples/a.fasta --outfile seq_counts.csv
```
- CSV 示例：
```csv
filename,sequence,count
a,ACGT,2
a,TTTT,1
```
- 说明：导出时会归一化大小写（转大写）并去除首尾空白；计数按降序排列

### filter（长度过滤）

- 功能：按长度过滤（两种模式：批量目录、或拼接多个文件后过滤）
- 参数：
  - 批量模式：`--input-dir`、`--output-dir`、`--min-len`、`--max-len`
  - 拼接模式：`--input-files <files...>`、`--outfile`、`--min-len`、`--max-len`
- 使用示例（批量）：
```bash
hammer_fastx filter --input-dir demux_out --output-dir filtered --min-len 200 --max-len 1000
```
- 使用示例（拼接）：
```bash
hammer_fastx filter --input-files a.fasta b.fasta --outfile ab.filtered.fasta --min-len 200
```
- 输出：对应目录或单文件结果

### Ns_count（N 区域组合统计）

- 功能：将 reads 比对到含 `N` 的参考序列，基于锚定区域统计每个 `N` 区块的组合，支持导出匹配 reads
- 参数：`--reads`、`--refSEQ`、`--output`、`--threads`、`--group`、`--dig`、`--mismatches`、`--anchor-len`、`--extract-matches`
- 使用示例：
```bash
hammer_fastx Ns_count \
  --reads reads.fasta \
  --refSEQ ref_with_N.fasta \
  --output ns_out \
  --threads 12 --group T0 --dig 2 --mismatches 2 --anchor-len 15 --extract-matches
```
- 输出：
  - `ns_out/<ref_id>_combo_counts.csv`
  - 选项开启时：`ns_out/<ref_id>_matched_reads.fasta`
- 参考 FASTA 示例（含 N 块）：
```fasta
>ref1
AAAANNNCCCTTTNNNGGG
```

### DNA2AA（DNA→蛋白）

- 功能：批量将目录中的 DNA FASTA 翻译为 AA FASTA
- 参数：`-i/--input`、`-o/--output`、`--aa-length`（默认 50）
- 使用示例：
```bash
hammer_fastx DNA2AA --input dna_dir --output aa_dir --aa-length 80
```
- 输出：`aa_dir/<input_basename>.fasta`
- 批量目录示例：
```bash
hammer_fastx DNA2AA --input data/dna --output data/aa
```

### count_AA（蛋白突变统计）

- 功能：比对样本蛋白 FASTA 与参考蛋白，统计突变与受保护位点；并行处理多文件
- 参数：`-r/--reference`、`-i/--input-dir`、`-o/--output-dir`、`-A/--aa-offset`、`-c/--config`、`--match_len`、`--threads`、`--chunk_size`
- 使用示例：
```bash
hammer_fastx count_AA \
  --reference ref_protein.fasta \
  --input-dir aa_dir \
  --output-dir aa_stats \
  --aa-offset 1 --config protected_sites.csv --threads 12 --chunk_size 500000
```
- 输出：每个输入 FASTA 对应一个 CSV（突变计数与频率）
- 保护位点配置 CSV 示例：
```csv
protected_sites
12
55
```

### find_seq（基序查找与片段提取）

- 功能：在 FASTA/FASTQ（支持 `.gz`）中查找指定 `motif` 及其反向互补，截取上/下游片段并按每读段唯一窗口计数，生成 CSV
- 参数：`--inputfile`、`--output`、`--motif`、`--up-flank`、`--down-flank`
  - 当两者都未给出时默认 `40`；若只给一侧，另一侧取同值
- 使用示例：
```bash
hammer_fastx find_seq \
  --inputfile merged.fasta \
  --output motif_windows.csv \
  --motif ATGC \
  --up-flank 40 --down-flank 40
```
- 输出 CSV 示例：
```csv
Sequence,UpFlank,DownFlank,ReadsCount
NAAATGCCTT,NAA,CTT,1
NAAATGCTTT,NAA,TTT,1
```
- 反向互补示例：当读段中出现 `GCAT`（`ATGC` 的反向互补）时，会归一到正向窗口并计数

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

## 输入/输出格式示例

### FASTA 示例
```fasta
>read1
ACGTNNNACGT
>read2
TTTTAAAAGGGG
```

### FASTQ 示例
```fastq
@read1
ACGTNNNACGT
+
FFFFFFFFFFFF
@read2
TTTTAAAAGGGG
+
FFFFFFFFFFFF
```

### demux 标签 CSV 示例
```csv
SampleID,F_tag,R_tag
S1,ACGTACGT,TGCATGCA
S2,AAAACCCC,GGGGTTTT
```

### Ns_count 组合统计 CSV（示意）
```csv
Combo,Count,Frequency
AAAA_TTTT,123,0.45
AAAC_TTGA,15,0.06
```

### count_AA 结果 CSV（示意）
```csv
Mutation,Count,Frequency
S12A,34,0.12
P55L,5,0.02
```

- **FASTA 格式**：`.fasta`, `.fa`, `.fna`
- **FASTQ 格式**：`.fastq`, `.fq`
- **压缩格式**：支持 `.gz` 压缩的 FASTA/FASTQ 文件
- **Barcode 文件**：制表符分隔的文本文件，格式：`barcode\t样本名`

## 依赖与优化建议

- Rust crates：bio、clap、csv、rayon、crossbeam-channel、dashmap、indicatif、flate2、anyhow、glob、num_cpus
- 外部工具：fastp、flash2（需在 PATH 中）
- 建议：
  - 根据 CPU 核心数调整线程参数
  - 大文件使用目录批量模式与分块参数（如 `count_AA --chunk_size`）
  - 输入输出放不同磁盘以优化 I/O

## 贡献与许可证

- 欢迎贡献，建议使用分支与 Pull Request 流程
- 运行 `cargo build` 进行编译，必要时添加测试后再提交
- 许可证：MIT（见 `LICENSE`）

## 作者

- Caizhaohui（GitHub: https://github.com/Caizhaohui）

## 依赖说明

Hammer_fastx - 让生物信息学数据分析更加高效！
