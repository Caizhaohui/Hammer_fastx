# Hammer_fastx

> 用 Rust 语言编写的多功能 FASTA/FASTQ 文件处理工具

## 项目简介

Hammer_fastx 是一个用 Rust 开发的轻量级命令行工具，旨在高效处理生物信息学中的 FASTA 和 FASTQ 文件。适用于基因组测序数据的预处理、格式转换、质控、合并及拆分等场景。

---

## 命令行参数与 `--help` 说明

运行以下命令可查看完整帮助信息：

```bash
./hammer_fastx --help
```

主要输出格式示例（不同版本略有差异）：

```
Hammer_fastx v0.7.1
CZH with the help of Gemini
一个用于处理FASTX文件、集成质控和合并功能的多功能工具集

USAGE:
    hammer_fastx <SUBCOMMAND>

OPTIONS:
    -h, --help       打印帮助信息
    -V, --version    打印版本号

SUBCOMMANDS:
    demux_all    [主流程] 运行从质控、合并到拆分的完整流程
    mergePE      [流程] 先质控后合并双端数据，可选择输出格式
    demux_only   [单步骤] 仅根据barcode拆分已合并的FASTQ文件
    Fastp        (包装器) 使用 fastp 对双端fastq文件进行质控
    Flash2       (包装器) 使用 flash2 合并双端 reads
    Stats        统计一个或多个FASTA/FASTQ文件中的序列信息
    Filter       根据序列长度过滤一个或多个FASTA/FASTQ文件
    Ns_count     将reads比对到含有N的参考序列上并提取组合
```

---

## 子命令与功能详细说明

### 1. demux_all
**[主流程] 一键运行质控(fastp)、合并(flash2)、barcode拆分**

- 输入：成对的原始双端 fastq 文件、barcode 文件等
- 步骤：依次调用 fastp 进行质控、flash2 合并 reads、barcode 拆分
- 常用参数：
    - `--in1`/`--in2`：输入的 R1/R2 fastq 文件
    - `--barcode`：barcode 文件
    - `--outdir`：输出目录
    - `--fastp-threads`：fastp 线程数
    - `--out-fasta`：输出为 FASTA 格式（默认 FASTQ）

### 2. mergePE
**[流程] 先用 fastp 质控，再用 flash2 合并双端数据，适于 PE 数据的融合处理**

- 主要参数同上（具体参数可通过 `mergePE --help` 查看）

### 3. demux_only
**[单步骤] 仅基于 barcode 对已合并 fastq 文件进行拆分**

- 用于合并后数据的快速 barcode 拆分
- 主要参数：
    - `--infile`：合并后的 fastq 文件
    - `--barcode`：barcode 文件
    - `--outdir`：结果输出目录

### 4. Fastp
**调用 fastp 对双端 fastq 文件进行质控**

- 实际为外部 fastp 工具的包装，适用于高效过滤与质量统计
- 常用参数：
    - `-i/--in1`：输入文件1 (Read1)
    - `-I/--in2`：输入文件2 (Read2)
    - `-o/--out1`：输出文件1 (Read1)
    - `-O/--out2`：输出文件2 (Read2)
    - `-h/--html`：指定 HTML 报告输出路径
    - `-j/--json`：指定 JSON 报告输出路径
    - `-R/--report-title`：报告标题
    - `-t/--threads`：线程数

### 5. Flash2
**调用 flash2 合并双端 reads**

- 合并 PE reads，提升拼接效率
- 参数详见 `Flash2 --help`

### 6. Stats
**统计FASTA/FASTQ文件中的序列信息**

- 包括序列条数、总碱基数、平均长度等
- 输入支持多个文件

### 7. Filter
**根据序列长度过滤FASTA/FASTQ文件**

- 可自定义最小/最大长度阈值

### 8. Ns_count
**将 reads 比对到含 N 的参考序列并提取组合**

- 适合特殊定制需求

---

## 使用示例

```bash
# 查看帮助
./hammer_fastx --help

# 统计序列文件信息
./hammer_fastx Stats -i example.fasta

# fastq 转换为 fasta（如有支持）
./hammer_fastx fq2fa -i example.fastq -o output.fasta

# 一键处理流程
./hammer_fastx demux_all --in1 R1.fq.gz --in2 R2.fq.gz --barcode barcodes.txt --outdir result_dir
```

---

## 输入输出格式

- **输入**：支持标准 FASTA (.fasta, .fa) 与 FASTQ (.fastq, .fq) 文件
- **输出**：支持 FASTA/FASTQ，或统计、质控、合并、拆分等中间结果

---

## 依赖环境

- Rust 1.60 及以上
- 需依赖外部工具 fastp、flash2（部分子命令）

---

## 贡献与反馈

欢迎 issue 与 PR！如有建议或 bug 可在 [GitHub Issues](https://github.com/Caizhaohui/Hammer_fastx/issues) 提出。

## 许可证

MIT License，详见 [LICENSE](./LICENSE)。

---

**作者**：[Caizhaohui](https://github.com/Caizhaohui)
