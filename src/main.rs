use anyhow::Result;
use clap::Parser;
use std::process::{Command, Stdio}; // 导入用于执行外部命令的模块

// ==================================================================================
// 主程序入口和子命令定义 (CLI Structure)
// ==================================================================================

#[derive(Parser, Debug)]
#[command(
    name = "Hammer_fastx",
    version = "v0.7.1", // 版本号提升, 最终修复 Ns_count panic bug
    author = "CZH with the help of Gemini",
    about = "一个用于处理FASTX文件、集成质控和合并功能的多功能工具集"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Parser, Debug)]
enum Commands {
    /// [主流程] 运行从质控、合并到拆分的完整流程
    #[command(name = "demux_all")]
    DemuxAll(pipeline::Args),

    /// [流程] 先质控后合并双端数据，可选择输出格式
    #[command(name = "mergePE")]
    MergePE(merge_pe::Args),

    /// [单步骤] 仅根据barcode拆分已合并的FASTQ文件
    #[command(name = "demux_only")]
    DemuxOnly(demux::Args),

    /// (包装器) 使用 fastp 对双端fastq文件进行质控
    Fastp(fastp::Args),

    /// (包装器) 使用 flash2 合并双端 reads
    Flash2(flash2::Args),

    /// 统计一个或多个FASTA/FASTQ文件中的序列信息
    Stats(stats::Args),

    /// 根据序列长度过滤一个或多个FASTA/FASTQ文件
    Filter(filter::Args),

    /// 将reads比对到含有N的参考序列上并提取组合
    #[command(name = "Ns_count")]
    NsCount(ns_count::Args),
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::DemuxAll(args) => pipeline::run(args),
        Commands::MergePE(args) => merge_pe::run(args),
        Commands::DemuxOnly(args) => demux::run(args),
        Commands::Fastp(args) => fastp::run(args),
        Commands::Flash2(args) => flash2::run(args),
        Commands::Stats(args) => stats::run(args),
        Commands::Filter(args) => filter::run(args),
        Commands::NsCount(args) => ns_count::run(args),
    }
}

// ==================================================================================
// `pipeline` 子命令模块 (对应 demux_all)
// ==================================================================================
mod pipeline {
    use super::{demux, fastp, flash2};
    use anyhow::{Context, Result};
    use clap::Parser;
    use std::fs;
    use std::path::PathBuf;
    use std::time::Instant;

    #[derive(Parser, Debug)]
    #[command(name = "demux_all", about = "[主流程] 运行从质控(fastp)、合并(flash2)到拆分(demux)的完整流程")]
    pub struct Args {
        // === 主要输入 ===
        #[arg(short = 'i', long, help = "原始输入文件1 (Raw Read1)")]
        pub in1: PathBuf,

        #[arg(short = 'I', long, help = "原始输入文件2 (Raw Read2)")]
        pub in2: PathBuf,

        #[arg(long, help = "样本标签文件 (用于 demux)")]
        pub tags: PathBuf,

        // === 主要输出 ===
        #[arg(short = 'o', long, help = "总输出目录，所有结果和中间文件将存放于此")]
        pub output_dir: PathBuf,

        // === 流程控制 ===
        #[arg(long, help = "流程成功结束后删除 fastp 和 flash2 的中间文件")]
        pub cleanup: bool,

        // === fastp 参数 ===
        #[arg(long, help = "fastp 使用的线程数", default_value_t = 4)]
        pub fastp_threads: usize,

        // === flash2 参数 ===
        #[arg(long, help = "flash2 使用的线程数", default_value_t = 4)]
        pub flash_threads: usize,
        #[arg(long, help = "flash2 最小重叠长度", default_value_t = 10)]
        pub min_overlap: usize,
        #[arg(long, help = "flash2 最大重叠长度", default_value_t = 300)]
        pub max_overlap: usize,

        // === demux_only 参数 ===
        #[arg(long, help = "demux_only 使用的线程数", default_value_t = num_cpus::get_physical())]
        pub demux_threads: usize,
        #[arg(short = 'l', long, help = "demux_only 使用的标签长度", default_value_t = 8)]
        pub tag_len: usize,
        #[arg(long, help = "demux_only 时激活修剪标签的功能")]
        pub trim: bool,
        #[arg(long, help = "demux_only 后输出为 FASTA 格式 (默认为 FASTQ)")]
        pub out_fasta: bool,
    }

    pub fn run(args: Args) -> Result<()> {
        let total_start_time = Instant::now();
        println!("🚀 [主流程] 开始执行 hammer_fastx demux_all 完整流程...");

        // --- 1. 设置目录结构 ---
        let fastp_dir = args.output_dir.join("01_fastp_out");
        let flash_dir = args.output_dir.join("02_flash2_out");
        let demux_dir = args.output_dir.join("03_demux_out");

        fs::create_dir_all(&args.output_dir)
            .with_context(|| format!("无法创建总输出目录: {:?}", args.output_dir))?;
        fs::create_dir_all(&fastp_dir)
            .with_context(|| format!("无法创建 fastp 输出目录: {:?}", fastp_dir))?;
        fs::create_dir_all(&flash_dir)
            .with_context(|| format!("无法创建 flash2 输出目录: {:?}", flash_dir))?;
        // demux::run 会自己创建目录，这里无需预创建

        // --- 2. 执行 fastp ---
        println!("\n[步骤 1/3] ➡️  执行 fastp 质控...");
        let fastp_out1 = fastp_dir.join("filtered_R1.fastq.gz");
        let fastp_out2 = fastp_dir.join("filtered_R2.fastq.gz");
        let fastp_args = fastp::Args {
            in1: args.in1.clone(),
            in2: args.in2.clone(),
            out1: fastp_out1.clone(),
            out2: fastp_out2.clone(),
            html: Some(fastp_dir.join("fastp_report.html")),
            json: Some(fastp_dir.join("fastp_report.json")),
            report_title: "Hammer_fastx demux_all pipeline: fastp report".to_string(),
            threads: Some(args.fastp_threads),
        };
        fastp::run(fastp_args)?;

        // --- 3. 执行 flash2 ---
        println!("\n[步骤 2/3] ➡️  执行 flash2 合并...");
        let flash_prefix = "merged";
        let flash_args = flash2::Args {
            in1: fastp_out1.clone(),
            in2: fastp_out2.clone(),
            out_prefix: flash_prefix.to_string(),
            out_dir: flash_dir.clone(),
            min_overlap: args.min_overlap,
            max_overlap: args.max_overlap,
            threads: args.flash_threads,
        };
        flash2::run(flash_args)?;

        // --- 4. 执行 demux_only ---
        println!("\n[步骤 3/3] ➡️  执行 demux_only 拆分...");
        let demux_input = flash_dir.join(format!("{}.extendedFrags.fastq", flash_prefix));
        let demux_args = demux::Args {
            inputfile: demux_input,
            output: demux_dir.clone(),
            threads: args.demux_threads,
            tags: args.tags.clone(),
            tag_len: args.tag_len,
            trim: args.trim,
            out_fasta: args.out_fasta,
        };
        demux::run(demux_args)?;

        // --- 5. 清理中间文件 ---
        if args.cleanup {
            println!("\n[清理] 正在删除中间文件...");
            fs::remove_dir_all(&fastp_dir)
                .with_context(|| format!("清理 fastp 目录失败: {:?}", fastp_dir))?;
            fs::remove_dir_all(&flash_dir)
                .with_context(|| format!("清理 flash2 目录失败: {:?}", flash_dir))?;
            println!("✔ 清理完成。");
        }

        println!("\n🎉 [主流程] 所有步骤成功完成！总耗时: {:.2?}", total_start_time.elapsed());
        println!("最终拆分结果位于: {}", demux_dir.display());

        Ok(())
    }
}

// ==================================================================================
// `merge_pe` 子命令模块 (无变化)
// ==================================================================================
mod merge_pe {
    use super::{fastp, flash2};
    use anyhow::{anyhow, Context, Result};
    use bio::io::{fasta, fastq};
    use clap::Parser;
    use std::fs;
    use std::io::BufReader;
    use std::path::PathBuf;
    use std::time::Instant;

    #[derive(Parser, Debug)]
    #[command(name = "mergePE", about = "[流程] 先用 fastp 质控，再用 flash2 合并双端数据")]
    pub struct Args {
        // === 主要输入 ===
        #[arg(short = 'i', long, help = "原始输入文件1 (Raw Read1)")]
        pub in1: PathBuf,
        #[arg(short = 'I', long, help = "原始输入文件2 (Raw Read2)")]
        pub in2: PathBuf,

        // === 主要输出 ===
        #[arg(short = 'o', long, help = "最终合并文件的输出路径")]
        pub outfile: PathBuf,

        // === 流程控制 ===
        #[arg(long, help = "将最终输出格式转换为FASTA (默认: FASTQ)")]
        pub out_fasta: bool,
        #[arg(long, help = "流程成功结束后删除中间文件")]
        pub cleanup: bool,
        #[arg(long, help = "存放中间文件的目录 (默认在输出文件所在目录下创建 'intermediates')")]
        pub temp_dir: Option<PathBuf>,

        // === fastp 参数 ===
        #[arg(long, help = "fastp 使用的线程数", default_value_t = 4)]
        pub fastp_threads: usize,

        // === flash2 参数 ===
        #[arg(long, help = "flash2 使用的线程数", default_value_t = 4)]
        pub flash_threads: usize,
        #[arg(long, help = "flash2 最小重叠长度", default_value_t = 10)]
        pub min_overlap: usize,
        #[arg(long, help = "flash2 最大重叠长度", default_value_t = 300)]
        pub max_overlap: usize,
    }

    pub fn run(args: Args) -> Result<()> {
        let total_start_time = Instant::now();
        println!("🚀 [流程] 开始执行 hammer_fastx mergePE 流程...");

        // --- 1. 设置目录结构 ---
        let output_parent_dir = args.outfile.parent().ok_or_else(|| anyhow!("无法获取输出文件的父目录"))?;
        fs::create_dir_all(output_parent_dir)
            .with_context(|| format!("无法创建输出目录: {:?}", output_parent_dir))?;

        let temp_dir = args.temp_dir.clone().unwrap_or_else(|| output_parent_dir.join("intermediates"));
        fs::create_dir_all(&temp_dir)
            .with_context(|| format!("无法创建中间目录: {:?}", temp_dir))?;

        // --- 2. 执行 fastp ---
        println!("\n[步骤 1/3] ➡️  执行 fastp 质控...");
        let fastp_out1 = temp_dir.join("filtered_R1.fastq.gz");
        let fastp_out2 = temp_dir.join("filtered_R2.fastq.gz");
        let fastp_args = fastp::Args {
            in1: args.in1.clone(),
            in2: args.in2.clone(),
            out1: fastp_out1.clone(),
            out2: fastp_out2.clone(),
            html: Some(temp_dir.join("fastp_report.html")),
            json: Some(temp_dir.join("fastp_report.json")),
            report_title: "Hammer_fastx mergePE: fastp report".to_string(),
            threads: Some(args.fastp_threads),
        };
        fastp::run(fastp_args)?;

        // --- 3. 执行 flash2 ---
        println!("\n[步骤 2/3] ➡️  执行 flash2 合并...");
        let flash_prefix = "merged";
        let flash_args = flash2::Args {
            in1: fastp_out1.clone(),
            in2: fastp_out2.clone(),
            out_prefix: flash_prefix.to_string(),
            out_dir: temp_dir.clone(),
            min_overlap: args.min_overlap,
            max_overlap: args.max_overlap,
            threads: args.flash_threads,
        };
        flash2::run(flash_args)?;

        // --- 4. 写入最终输出文件 (并转换格式) ---
        println!("\n[步骤 3/3] ➡️  写入最终输出文件...");
        let merged_fastq_path = temp_dir.join(format!("{}.extendedFrags.fastq", flash_prefix));
        
        let in_file = fs::File::open(&merged_fastq_path)
            .with_context(|| format!("无法打开合并后的文件: {:?}", merged_fastq_path))?;
        let in_reader = BufReader::new(in_file);
        let fastq_reader = fastq::Reader::new(in_reader);

        let out_file = fs::File::create(&args.outfile)
            .with_context(|| format!("无法创建最终输出文件: {:?}", args.outfile))?;

        let mut records_written = 0;
        if args.out_fasta {
            let mut fasta_writer = fasta::Writer::new(out_file);
            for result in fastq_reader.records() {
                let record = result?;
                let fasta_record = fasta::Record::with_attrs(record.id(), record.desc(), record.seq());
                fasta_writer.write_record(&fasta_record)?;
                records_written += 1;
            }
        } else {
            let mut fastq_writer = fastq::Writer::new(out_file);
            for result in fastq_reader.records() {
                let record = result?;
                fastq_writer.write_record(&record)?;
                records_written += 1;
            }
        }
        println!("✔ 成功写入 {} 条记录到 {}", records_written, args.outfile.display());

        // --- 5. 清理中间文件 ---
        if args.cleanup {
            println!("\n[清理] 正在删除中间文件...");
            fs::remove_dir_all(&temp_dir)
                .with_context(|| format!("清理中间目录失败: {:?}", temp_dir))?;
            println!("✔ 清理完成。");
        }

        println!("\n🎉 [流程] mergePE 流程成功完成！总耗时: {:.2?}", total_start_time.elapsed());
        Ok(())
    }
}


// ==================================================================================
// `fastp` 子命令模块 (无变化)
// ==================================================================================
mod fastp {
    use super::{Command, Stdio}; // 引用顶层导入
    use anyhow::{anyhow, Context, Result};
    use clap::Parser;
    use std::path::PathBuf;

    #[derive(Parser, Debug)]
    #[command(
        name = "fastp",
        about = "包装器: 调用 fastp 对双端fastq文件进行质控"
    )]
    pub struct Args {
        #[arg(short = 'i', long, help = "输入文件1 (Read1)")]
        pub in1: PathBuf,

        #[arg(short = 'I', long, help = "输入文件2 (Read2)")]
        pub in2: PathBuf,

        #[arg(short = 'o', long, help = "输出文件1 (Read1)")]
        pub out1: PathBuf,

        #[arg(short = 'O', long, help = "输出文件2 (Read2)")]
        pub out2: PathBuf,

        #[arg(short = 'h', long, help = "指定HTML报告的输出路径")]
        pub html: Option<PathBuf>,

        #[arg(short = 'j', long, help = "指定JSON报告的输出路径")]
        pub json: Option<PathBuf>,

        #[arg(short = 'R', long, help = "报告的标题", default_value = "fastp report")]
        pub report_title: String,

        #[arg(short = 't', long, help = "线程数 (默认: 自动检测)")]
        pub threads: Option<usize>,
    }

    /// 检查指定的外部命令是否存在于系统的PATH中
    fn command_exists(cmd: &str) -> bool {
        Command::new(cmd)
            .arg("--version")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok()
    }

    pub fn run(args: Args) -> Result<()> {
        println!("---> 启动 fastp 质控流程...");

        if !command_exists("fastp") {
            return Err(anyhow!(
                "错误: 未找到 'fastp' 可执行文件。\n请确保 fastp 已安装并存在于您的系统 PATH 环境变量中。"
            ));
        }

        let mut cmd = Command::new("fastp");
        cmd.arg("-i").arg(&args.in1);
        cmd.arg("-I").arg(&args.in2);
        cmd.arg("-o").arg(&args.out1);
        cmd.arg("-O").arg(&args.out2);
        cmd.arg("-R").arg(&args.report_title);

        if let Some(html_path) = &args.html {
            cmd.arg("-h").arg(html_path);
        }
        if let Some(json_path) = &args.json {
            cmd.arg("-j").arg(json_path);
        }
        if let Some(threads) = args.threads {
            cmd.arg("-t").arg(threads.to_string());
        }

        println!("🔧 执行命令: {:?}", cmd);

        let status = cmd
            .status()
            .with_context(|| "执行 fastp 命令失败。请检查 fastp 是否已正确安装。")?;

        if status.success() {
            println!("\n✔ fastp 质控成功完成！");
            println!("   - 清理后的 R1: {}", args.out1.display());
            println!("   - 清理后的 R2: {}", args.out2.display());
            if let Some(html_path) = &args.html {
                println!("   - HTML 报告: {}", html_path.display());
            }
            Ok(())
        } else {
            Err(anyhow!(
                "fastp 执行失败，退出码: {:?}\n请检查 fastp 的输出日志以获取详细错误信息。",
                status.code()
            ))
        }
    }
}

// ==================================================================================
// `flash2` 子命令模块 (无变化)
// ==================================================================================
mod flash2 {
    use super::{Command, Stdio}; // 引用顶层导入
    use anyhow::{anyhow, Context, Result};
    use clap::Parser;
    use std::path::PathBuf;

    #[derive(Parser, Debug)]
    #[command(
        name = "flash2",
        about = "包装器: 调用 flash2 合并双端 reads"
    )]
    pub struct Args {
        #[arg(help = "输入的 Read1 文件")]
        pub in1: PathBuf,
        
        #[arg(help = "输入的 Read2 文件")]
        pub in2: PathBuf,

        #[arg(short = 'o', long, help = "输出文件的前缀", default_value = "out")]
        pub out_prefix: String,
        
        #[arg(short = 'd', long, help = "输出文件的目录", default_value = ".")]
        pub out_dir: PathBuf,

        #[arg(short = 'm', long, help = "最小重叠长度 (min overlap)", default_value_t = 10)]
        pub min_overlap: usize,

        #[arg(short = 'M', long, help = "最大重叠长度 (max overlap)", default_value_t = 300)]
        pub max_overlap: usize,

        #[arg(short = 't', long, help = "线程数 (默认: 1)", default_value_t = 1)]
        pub threads: usize,
    }
    
    /// 检查指定的外部命令是否存在于系统的PATH中
    fn command_exists(cmd: &str) -> bool {
        // flash2 没有 --version, 我们尝试直接运行它，但这有风险
        // 一个更安全的方法是使用 which 命令 (Unix-like)
        if cfg!(unix) {
            Command::new("which")
                .arg(cmd)
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .status()
                .map_or(false, |s| s.success())
        } else {
            // 对于 Windows, 可以使用 where
            Command::new("where")
                .arg(cmd)
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .status()
                .map_or(false, |s| s.success())
        }
    }

    pub fn run(args: Args) -> Result<()> {
        println!("---> 启动 flash2 合并流程...");
        
        if !command_exists("flash2") {
            return Err(anyhow!(
                "错误: 未找到 'flash2' 可执行文件。\n请确保 flash2 已安装并存在于您的系统 PATH 环境变量中。"
            ));
        }

        let mut cmd = Command::new("flash2");
        cmd.arg(&args.in1);
        cmd.arg(&args.in2);
        cmd.arg("-o").arg(&args.out_prefix);
        cmd.arg("-d").arg(&args.out_dir);
        cmd.arg("-m").arg(args.min_overlap.to_string());
        cmd.arg("-M").arg(args.max_overlap.to_string());
        cmd.arg("-t").arg(args.threads.to_string());

        println!("🔧 执行命令: {:?}", cmd);

        let status = cmd
            .status()
            .with_context(|| "执行 flash2 命令失败。请检查 flash2 是否已正确安装。")?;

        if status.success() {
            println!("\n✔ flash2 合并成功完成！");
            println!("   - 输出目录: {}", args.out_dir.display());
            println!("   - 输出文件前缀: {}", args.out_prefix);
            println!("   - 合并后文件: {}", args.out_dir.join(format!("{}.extendedFrags.fastq", args.out_prefix)).display());
            Ok(())
        } else {
            Err(anyhow!(
                "flash2 执行失败，退出码: {:?}\n请检查 flash2 的输出日志以获取详细错误信息。",
                status.code()
            ))
        }
    }
}


// ==================================================================================
// `common` 模块: 存放共享的工具函数 (无变化)
// ==================================================================================
mod common {
    use anyhow::{anyhow, Result};
    use flate2::bufread::MultiGzDecoder;
    use std::fs::File;
    use std::io::{BufRead, BufReader, Read}; // 引入 Read trait
    use std::path::Path;

    pub enum Format {
        Fasta,
        Fastq,
    }

    pub fn detect_format(path: &Path) -> Result<Format> {
        let file = File::open(path)?;
        let buf_reader = BufReader::new(file);
        let mut first_char_reader: Box<dyn BufRead> =
            if path.extension().map_or(false, |ext| ext == "gz") {
                Box::new(BufReader::new(MultiGzDecoder::new(buf_reader)))
            } else {
                Box::new(buf_reader)
            };
        let mut buf = [0; 1];
        // 使用 read_exact 确保我们准确读取一个字节
        match first_char_reader.read_exact(&mut buf) {
            Ok(_) => match buf[0] {
                b'>' => Ok(Format::Fasta),
                b'@' => Ok(Format::Fastq),
                _ => Err(anyhow!(
                    "无法识别的文件格式: {:?}，请确保文件以 '>' (FASTA) 或 '@' (FASTQ) 开头。",
                    path
                )),
            },
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
                Err(anyhow!("文件为空或无法读取: {:?}", path))
            }
            Err(e) => Err(e.into()),
        }
    }
}

// ==================================================================================
// `demux` 子命令模块 (对应 demux_only)
// ==================================================================================
mod demux {
    use anyhow::{anyhow, Context, Result};
    use bio::io::{
        fasta,
        fastq::{self, Record},
    };
    use clap::Parser;
    use csv::ReaderBuilder;
    use flate2::bufread::MultiGzDecoder;
    use indicatif::{ProgressBar, ProgressStyle};
    use rayon::prelude::*;
    use std::collections::{HashMap, HashSet};
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::{Path, PathBuf};
    use std::sync::Arc;
    use std::thread;
    use std::time::Instant;

    const CHUNK_SIZE: usize = 8192;

    #[derive(Parser, Debug)]
    pub struct Args {
        #[arg(long, help = "输入的FASTQ文件 (可为.gz压缩)")]
        pub inputfile: PathBuf,

        #[arg(long, help = "输出目录")]
        pub output: PathBuf,

        #[arg(long, help = "线程数", default_value_t = num_cpus::get_physical())]
        pub threads: usize,
        
        #[arg(short, long, help = "样本标签文件 (CSV格式, 列: SampleID,F_tag,R_tag)")]
        pub tags: PathBuf,
        
        #[arg(short = 'l', long, default_value_t = 8, help = "标签的长度")]
        pub tag_len: usize,
        
        #[arg(long, help = "激活此选项以修剪拆分后序列两端的标签")]
        pub trim: bool,
        
        #[arg(long, help = "将输出格式转换为FASTA (默认: FASTQ)")]
        pub out_fasta: bool,
    }

    #[derive(Debug, Clone)]
    struct MatchInfo {
        sample_id: String,
        orientation: Orientation,
    }
    #[derive(Debug, Clone, PartialEq)]
    enum Orientation {
        Forward,
        Reverse,
    }
    type RawChunk = Vec<Record>;
    type ProcessedChunk = HashMap<String, Vec<Record>>;
    enum GenericWriter {
        Fastq(fastq::Writer<File>),
        Fasta(fasta::Writer<File>),
    }
    impl GenericWriter {
        fn write_record(&mut self, record: &Record) -> Result<()> {
            match self {
                GenericWriter::Fastq(writer) => writer.write_record(record)?,
                GenericWriter::Fasta(writer) => {
                    let fasta_record =
                        fasta::Record::with_attrs(record.id(), record.desc(), record.seq());
                    writer.write_record(&fasta_record)?;
                }
            }
            Ok(())
        }
    }
    fn load_tags(
        tag_file: &Path,
        tag_len: usize,
    ) -> Result<(HashMap<(Vec<u8>, Vec<u8>), MatchInfo>, HashSet<String>)> {
        let mut lookup_map = HashMap::new();
        let mut all_samples = HashSet::new();
        let file = File::open(tag_file)
            .with_context(|| format!("无法打开标签文件: {:?}", tag_file))?;
        let mut rdr = ReaderBuilder::new()
            .has_headers(true)
            .flexible(true)
            .delimiter(b',')
            .from_reader(file);
        let headers = rdr.headers()?.clone();
        if !headers.iter().any(|h| h == "SampleID")
            || !headers.iter().any(|h| h == "F_tag")
            || !headers.iter().any(|h| h == "R_tag")
        {
            return Err(anyhow!(
                "标签文件必须包含 'SampleID', 'F_tag', 'R_tag' 这几列。"
            ));
        }
        for result in rdr.records() {
            let record = result?;
            let sample_id = record.get(0).ok_or_else(|| anyhow!("缺少 SampleID"))?.to_string();
            let f_tag = record.get(1).ok_or_else(|| anyhow!("缺少 F_tag"))?.as_bytes().to_ascii_uppercase();
            let r_tag = record.get(2).ok_or_else(|| anyhow!("缺少 R_tag"))?.as_bytes().to_ascii_uppercase();
            if f_tag.len() != tag_len || r_tag.len() != tag_len {
                return Err(anyhow!("样本 {} 的标签长度与指定的 --tag-len {} 不匹配", sample_id, tag_len));
            }
            all_samples.insert(sample_id.clone());
            let r_tag_rc = bio::alphabets::dna::revcomp(&r_tag);
            let fwd_key = (f_tag.clone(), r_tag_rc.clone());
            lookup_map.insert(fwd_key, MatchInfo { sample_id: sample_id.clone(), orientation: Orientation::Forward });
            
            let f_tag_rc = bio::alphabets::dna::revcomp(&f_tag);
            let rev_key = (r_tag_rc, f_tag_rc);
            lookup_map.insert(rev_key, MatchInfo { sample_id, orientation: Orientation::Reverse });
        }
        Ok((lookup_map, all_samples))
    }
    fn reader_thread(
        input_path: PathBuf,
        tx: crossbeam_channel::Sender<RawChunk>,
        pb: ProgressBar,
    ) -> Result<()> {
        let file = File::open(&input_path)?;
        let buf_reader = BufReader::new(file);
        let boxed_buf_reader: Box<dyn BufRead> =
            if input_path.extension().map_or(false, |ext| ext == "gz") {
                Box::new(BufReader::new(MultiGzDecoder::new(buf_reader)))
            } else {
                Box::new(buf_reader)
            };
        let mut records_iter = fastq::Reader::new(boxed_buf_reader).records();
        loop {
            let mut chunk = Vec::with_capacity(CHUNK_SIZE);
            for _ in 0..CHUNK_SIZE {
                match records_iter.next() {
                    Some(Ok(record)) => chunk.push(record),
                    Some(Err(e)) => return Err(e.into()),
                    None => break,
                }
            }
            if chunk.is_empty() {
                break;
            }
            pb.inc(chunk.len() as u64);
            if tx.send(chunk).is_err() {
                break;
            }
        }
        pb.finish_with_message("✔ 文件读取完成");
        Ok(())
    }
    fn worker_thread(
        rx_raw: crossbeam_channel::Receiver<RawChunk>,
        tx_processed: crossbeam_channel::Sender<ProcessedChunk>,
        lookup_map: Arc<HashMap<(Vec<u8>, Vec<u8>), MatchInfo>>,
        args: Arc<Args>,
    ) -> Result<()> {
        rx_raw.into_iter().par_bridge().for_each(|chunk| {
            let processed_results: Vec<(String, Record)> = chunk
                .into_par_iter()
                .filter_map(|record| process_record(&record, &lookup_map, &args))
                .collect();
            let mut processed_chunk: ProcessedChunk = HashMap::new();
            for (sample_id, record) in processed_results {
                processed_chunk.entry(sample_id).or_default().push(record);
            }
            if !processed_chunk.is_empty() {
                let _ = tx_processed.send(processed_chunk);
            }
        });
        Ok(())
    }
    fn process_record<'a>(
        record: &'a Record,
        lookup_map: &'a HashMap<(Vec<u8>, Vec<u8>), MatchInfo>,
        args: &'a Args,
    ) -> Option<(String, Record)> {
        let seq = record.seq();
        if seq.len() < args.tag_len * 2 {
            return Some(("unmatched".to_string(), record.clone()));
        }
        let read_start = seq[..args.tag_len].to_ascii_uppercase();
        let read_end = seq[seq.len() - args.tag_len..].to_ascii_uppercase();
        let lookup_key = (read_start, read_end);
        match lookup_map.get(&lookup_key) {
            Some(match_info) => {
                let final_record = if args.trim {
                    let trimmed_seq = &seq[args.tag_len..seq.len() - args.tag_len];
                    let trimmed_qual = &record.qual()[args.tag_len..record.qual().len() - args.tag_len];
                    if match_info.orientation == Orientation::Reverse {
                        let rc_seq = bio::alphabets::dna::revcomp(trimmed_seq);
                        let mut rc_qual = trimmed_qual.to_vec();
                        rc_qual.reverse();
                        Record::with_attrs(record.id(), record.desc(), &rc_seq, &rc_qual)
                    } else {
                        Record::with_attrs(record.id(), record.desc(), trimmed_seq, trimmed_qual)
                    }
                } else {
                    record.clone()
                };
                Some((match_info.sample_id.clone(), final_record))
            }
            None => Some(("unmatched".to_string(), record.clone())),
        }
    }
    fn writer_thread(
        rx_processed: crossbeam_channel::Receiver<ProcessedChunk>,
        output_dir: PathBuf,
        mut all_samples: HashSet<String>,
        out_fasta: bool,
    ) -> Result<HashMap<String, u64>> {
        let mut writers: HashMap<String, GenericWriter> = HashMap::new();
        let extension = if out_fasta { "fasta" } else { "fastq" };
        
        // 确保 unmatched 样本也被包含
        all_samples.insert("unmatched".to_string());

        for sample_id in &all_samples {
            let path = output_dir.join(format!("{}.{}", sample_id, extension));
            let file = File::create(&path)?;
            let writer = if out_fasta {
                GenericWriter::Fasta(fasta::Writer::new(file))
            } else {
                GenericWriter::Fastq(fastq::Writer::new(file))
            };
            writers.insert(sample_id.clone(), writer);
        }

        let mut counts: HashMap<String, u64> = HashMap::new();
        for chunk in rx_processed {
            for (sample_id, records) in chunk {
                *counts.entry(sample_id.clone()).or_insert(0) += records.len() as u64;
                let writer = writers.get_mut(&sample_id).expect("Writer for sample not found!");
                for record in records {
                    writer.write_record(&record)?;
                }
            }
        }
        Ok(counts)
    }
    fn print_summary(counts: HashMap<String, u64>, start_time: Instant, output_dir: &PathBuf) {
        let duration = start_time.elapsed();
        let total_reads = counts.values().sum::<u64>();
        let matched_reads = total_reads - *counts.get("unmatched").unwrap_or(&0);
        println!("\n\n==================== 拆分总结 (多线程) ====================");
        println!("处理耗时: {:.2?}", duration);
        println!("总处理 Reads: {}", total_reads);
        if total_reads > 0 {
            let matched_percent = matched_reads as f64 * 100.0 / total_reads as f64;
            let unmatched_percent = *counts.get("unmatched").unwrap_or(&0) as f64 * 100.0 / total_reads as f64;
            println!("  - 匹配上的 Reads: {:>10} ({:.2}%)", matched_reads, matched_percent);
            println!("  - 未匹配的 Reads: {:>10} ({:.2}%)", counts.get("unmatched").unwrap_or(&0), unmatched_percent);
            println!("--------------------------------------------------");
            let mut sorted_samples: Vec<_> = counts.into_iter().collect();
            sorted_samples.sort_by(|a, b| b.1.cmp(&a.1));
            for (sample, count) in sorted_samples {
                if sample != "unmatched" {
                    let sample_percent = count as f64 * 100.0 / total_reads as f64;
                    println!("  - 样本 {}: {:>10} reads ({:.2}%)", sample, count, sample_percent);
                }
            }
        }
        println!("============================================================");
        println!("✔ 程序执行完毕! 结果已输出至: {}", output_dir.display());
    }
    pub fn run(args: Args) -> Result<()> {
        let start_time = Instant::now();
        let output_dir = args.output.clone();
        std::fs::create_dir_all(&output_dir)
            .with_context(|| format!("无法创建输出目录: {:?}", output_dir))?;
        println!("---> 正在加载标签...");
        let (lookup_map, all_samples) = load_tags(&args.tags, args.tag_len)?;
        let lookup_map = Arc::new(lookup_map);
        let args_arc = Arc::new(args);
        let channel_capacity = args_arc.threads * 2;
        let (raw_tx, raw_rx) = crossbeam_channel::bounded::<RawChunk>(channel_capacity);
        let (processed_tx, processed_rx) = crossbeam_channel::bounded::<ProcessedChunk>(channel_capacity);
        let pb = ProgressBar::new_spinner();
        pb.enable_steady_tick(std::time::Duration::from_millis(120));
        pb.set_style(
            ProgressStyle::default_spinner()
                .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
                .template("{spinner:.blue} [{elapsed_precise}] {msg} {pos:>10} reads")?,
        );
        pb.set_message("正在处理...");
        thread::scope(|s| -> Result<()> {
            let out_fasta_flag = args_arc.out_fasta;
            let output_dir_for_writer = output_dir.clone();
            let writer_handle = s.spawn(move || {
                writer_thread(processed_rx, output_dir_for_writer, all_samples, out_fasta_flag)
            });
            for _ in 0..args_arc.threads {
                let (rx_raw, tx_processed) = (raw_rx.clone(), processed_tx.clone());
                let (lookup_map_clone, args_clone) = (lookup_map.clone(), Arc::clone(&args_arc));
                s.spawn(move || worker_thread(rx_raw, tx_processed, lookup_map_clone, args_clone));
            }
            drop(processed_tx);
            reader_thread(args_arc.inputfile.clone(), raw_tx, pb)?;
            match writer_handle.join() {
                Ok(Ok(counts)) => print_summary(counts, start_time, &output_dir),
                Ok(Err(e)) => eprintln!("写入线程出错: {:?}", e),
                Err(e) => eprintln!("写入线程发生Panic: {:?}", e),
            }
            Ok(())
        })?;
        Ok(())
    }
}

// ==================================================================================
// `stats` 子命令模块 (无变化)
// ==================================================================================
mod stats {
    use super::common::{detect_format, Format};
    use anyhow::Result;
    use bio::io::{fasta, fastq};
    use clap::Parser;
    use flate2::bufread::MultiGzDecoder;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::{Path, PathBuf};

    #[derive(Parser, Debug)]
    pub struct Args {
        #[arg(long, help = "一个或多个输入文件 (支持通配符, e.g., '*.fasta')", required = true, num_args = 1..)]
        inputfile: Vec<PathBuf>,
    }
    
    struct FileStats {
        filename: String,
        count: u64,
        total_len: u64,
        min_len: usize,
        max_len: usize,
    }

    fn get_sample_name(path: &Path) -> String {
        path.file_name()
            .unwrap_or_default()
            .to_str()
            .unwrap_or_default()
            .split('.')
            .next()
            .unwrap_or("")
            .to_string()
    }
    
    fn print_stats_table(stats: &[FileStats]) {
        if stats.is_empty() {
            println!("未处理任何文件或未找到任何序列。");
            return;
        }

        println!("\n====================================== 序列统计汇总 ======================================");
        println!("{:<25} {:>15} {:>18} {:>10} {:>10} {:>12}",
                 "样品名 (Sample)", "序列总数", "总碱基数", "最短长度", "最长长度", "平均长度");
        println!("{:-<25} {:-<15} {:-<18} {:-<10} {:-<10} {:-<12}",
                 "", "", "", "", "", "");

        for s in stats {
            let avg_len = if s.count > 0 {
                s.total_len as f64 / s.count as f64
            } else {
                0.0
            };
            println!("{:<25} {:>15} {:>18} {:>10} {:>10} {:<12.2}",
                       s.filename, s.count, s.total_len, s.min_len, s.max_len, avg_len);
        }
        println!("==========================================================================================");
    }

    pub fn run(args: Args) -> Result<()> {
        let mut all_stats: Vec<FileStats> = Vec::new();

        for input_path in &args.inputfile {
            println!("---> 正在处理: {}", input_path.display());
            let format = detect_format(input_path)?;

            let file = File::open(input_path)?;
            let buf_reader = BufReader::new(file);
            let input_reader: Box<dyn BufRead> =
                if input_path.extension().map_or(false, |ext| ext == "gz") {
                    Box::new(BufReader::new(MultiGzDecoder::new(buf_reader)))
                } else {
                    Box::new(buf_reader)
                };

            let mut count = 0;
            let mut total_len = 0;
            let mut min_len = usize::MAX;
            let mut max_len = 0;

            match format {
                Format::Fasta => {
                    let reader = fasta::Reader::new(input_reader);
                    for result in reader.records() {
                        let record = result?;
                        count += 1;
                        let len = record.seq().len();
                        total_len += len as u64;
                        if len < min_len { min_len = len; }
                        if len > max_len { max_len = len; }
                    }
                }
                Format::Fastq => {
                    let reader = fastq::Reader::new(input_reader);
                    for result in reader.records() {
                        let record = result?;
                        count += 1;
                        let len = record.seq().len();
                        total_len += len as u64;
                        if len < min_len { min_len = len; }
                        if len > max_len { max_len = len; }
                    }
                }
            };
            
            all_stats.push(FileStats {
                filename: get_sample_name(input_path),
                count,
                total_len,
                min_len: if count > 0 { min_len } else { 0 },
                max_len,
            });
        }
        
        print_stats_table(&all_stats);
        
        Ok(())
    }
}

// ==================================================================================
// `filter` 子命令模块 (无变化)
// ==================================================================================
mod filter {
    use super::common::{detect_format, Format};
    use anyhow::Result;
    use bio::io::{fasta, fastq};
    use clap::Parser;
    use flate2::bufread::MultiGzDecoder;
    use std::fs::File;
    use std::io::{self, BufRead, BufReader, BufWriter, Write};
    use std::path::PathBuf;

    #[derive(Parser, Debug)]
    pub struct Args {
        #[arg(long, help = "一个或多个输入文件 (支持通配符, e.g., '*.fasta')", required = true, num_args = 1..)]
        inputfile: Vec<PathBuf>,

        #[arg(long, help = "输出文件 (默认: 标准输出)")]
        outfile: Option<PathBuf>,

        #[arg(short = 'm', long, help = "过滤掉小于该长度的序列")]
        min_len: Option<usize>,
        
        #[arg(short = 'M', long, help = "过滤掉大于该长度的序列")]
        max_len: Option<usize>,
    }

    pub fn run(args: Args) -> Result<()> {
        let min_len = args.min_len.unwrap_or(0);
        let max_len = args.max_len.unwrap_or(usize::MAX);

        let mut writer: Box<dyn Write> = if let Some(path) = args.outfile {
            Box::new(BufWriter::new(File::create(path)?))
        } else {
            Box::new(BufWriter::new(io::stdout().lock()))
        };

        for input_path in &args.inputfile {
            let format = detect_format(input_path)?;
            
            let file = File::open(input_path)?;
            let buf_reader = BufReader::new(file);
            let input_reader: Box<dyn BufRead> =
                if input_path.extension().map_or(false, |ext| ext == "gz") {
                    Box::new(BufReader::new(MultiGzDecoder::new(buf_reader)))
                } else {
                    Box::new(buf_reader)
                };

            match format {
                Format::Fasta => {
                    let reader = fasta::Reader::new(input_reader);
                    let mut fasta_writer = fasta::Writer::new(&mut writer);
                    for result in reader.records() {
                        let record = result?;
                        let len = record.seq().len();
                        if len >= min_len && len <= max_len {
                            fasta_writer.write_record(&record)?;
                        }
                    }
                }
                Format::Fastq => {
                    let reader = fastq::Reader::new(input_reader);
                    let mut fastq_writer = fastq::Writer::new(&mut writer);
                    for result in reader.records() {
                        let record = result?;
                        let len = record.seq().len();
                        if len >= min_len && len <= max_len {
                            fastq_writer.write_record(&record)?;
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

// ==================================================================================
// `ns_count` 子命令模块 (已重写并简化)
// ==================================================================================
mod ns_count {
    use anyhow::{Context, Result};
    use bio::io::fasta::{self, Record};
    use clap::Parser;
    use flate2::bufread::MultiGzDecoder;
    use indicatif::{ProgressBar, ProgressStyle};
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::PathBuf;
    use std::sync::Arc;
    use std::thread;

    const CHUNK_SIZE: usize = 4096;

    #[derive(Parser, Debug)]
    pub struct Args {
        #[arg(long, help = "包含待比对reads的FASTA文件 (可为.gz)")]
        reads: PathBuf,
        #[arg(long = "refSEQ", help = "包含带N区域参考序列的FASTA文件")]
        ref_seq: PathBuf,
        #[arg(long, help = "存放输出CSV文件的目录")]
        output: PathBuf,
        #[arg(long, help = "线程数", default_value_t = num_cpus::get_physical())]
        threads: usize,
        #[arg(long, help = "输出CSV文件中的列名前缀", default_value = "T0")]
        group: String,
        #[arg(long, help = "频率的小数位数", default_value_t = 2)]
        dig: u8,
        #[arg(long, help = "允许的最大错配数", default_value_t = 2)]
        mismatches: usize,
        #[arg(long, help = "Read需达到的相对于Ref的最小长度比例", default_value_t = 0.9)]
        min_len_ratio: f64,
        #[arg(long, help = "提取所有比对上的reads并输出到单独的FASTA文件")]
        extract_matches: bool,
    }

    struct MatchResult {
        ref_id: String,
        combo: Vec<u8>,
        read_record: Record,
    }

    struct RefData {
        id: String,
        seq: Vec<u8>,
        len: usize,
        n_blocks: Vec<(usize, usize)>,
    }

    fn find_n_blocks(seq: &[u8]) -> Vec<(usize, usize)> {
        let mut blocks = Vec::new();
        let mut in_block = false;
        let mut start = 0;
        for (i, &base) in seq.iter().enumerate() {
            if base == b'N' {
                if !in_block {
                    start = i;
                    in_block = true;
                }
            } else if in_block {
                blocks.push((start, i - start));
                in_block = false;
            }
        }
        if in_block {
            blocks.push((start, seq.len() - start));
        }
        blocks
    }

    // This is the completely rewritten, safer alignment function.
    fn find_alignment(
        read_seq: &[u8], 
        ref_data: &RefData, 
        max_mismatches: usize, 
        min_len_ratio: f64,
        is_rc_read: bool
    ) -> Option<Vec<u8>> {
        let read_len = read_seq.len();
        let ref_len = ref_data.len;

        // Condition 1: Check if the read is long enough compared to the reference.
        if (read_len as f64) < (ref_len as f64 * min_len_ratio) {
            return None;
        }

        // Condition 2: A read can't be longer than the reference for a valid alignment.
        // This is the CRITICAL fix that prevents underflow and subsequent panics.
        if read_len > ref_len {
            return None;
        }

        let mut best_mismatches = usize::MAX;
        let mut best_ref_start: Option<usize> = None;

        // Find the best alignment position by sliding the read over the reference.
        for ref_start in 0..=(ref_len - read_len) {
            let mut current_mismatches = 0;
            for i in 0..read_len {
                if read_seq[i] != ref_data.seq[ref_start + i] {
                    current_mismatches += 1;
                }
            }
            if current_mismatches < best_mismatches {
                best_mismatches = current_mismatches;
                best_ref_start = Some(ref_start);
            }
        }

        // If a potential best alignment was found and it's within the threshold...
        if let Some(start_pos) = best_ref_start {
            if best_mismatches <= max_mismatches {
                // Condition 3: Check if this best alignment covers ALL N-blocks.
                let mut all_n_blocks_covered = true;
                for &(n_start, n_len) in &ref_data.n_blocks {
                    if !(n_start >= start_pos && (n_start + n_len) <= (start_pos + read_len)) {
                        all_n_blocks_covered = false;
                        break;
                    }
                }

                if all_n_blocks_covered {
                    let mut combo_parts = Vec::new();
                    for &(n_start, n_len) in &ref_data.n_blocks {
                        let read_idx_start = n_start - start_pos;
                        let segment = &read_seq[read_idx_start .. read_idx_start + n_len];
                        if is_rc_read {
                            combo_parts.push(bio::alphabets::dna::revcomp(segment));
                        } else {
                            combo_parts.push(segment.to_vec());
                        }
                    }
                    return Some(combo_parts.join(&b'-'));
                }
            }
        }
        
        None
    }

    fn collector_thread(
        rx: crossbeam_channel::Receiver<MatchResult>,
        output_dir: PathBuf,
        group: String,
        dig: u8,
        extract_matches: bool,
        ref_data_map: HashMap<String, Vec<(usize, usize)>>,
    ) -> Result<()> {
        let mut counters: HashMap<String, HashMap<Vec<u8>, u64>> = HashMap::new();
        let mut writers: HashMap<String, fasta::Writer<File>> = HashMap::new();

        for result in rx {
            let counter = counters.entry(result.ref_id.clone()).or_default();
            *counter.entry(result.combo).or_insert(0) += 1;

            if extract_matches {
                let writer = writers.entry(result.ref_id.clone()).or_insert_with(|| {
                    let out_path = output_dir.join(format!("{}_matched_reads.fasta", result.ref_id));
                    fasta::Writer::to_file(out_path).expect("Failed to create writer")
                });
                writer.write_record(&result.read_record)?;
            }
        }

        for (ref_id, counter) in counters {
            let total: u64 = counter.values().sum();
            if total > 0 {
                let n_blocks = ref_data_map.get(&ref_id).unwrap();
                let n_label = (1..=n_blocks.len()).map(|i| format!("N{}", i)).collect::<Vec<_>>().join("_");
                let out_csv_path = output_dir.join(format!("{}_combo_counts.csv", ref_id));
                let mut csv_writer = csv::Writer::from_path(out_csv_path)?;
                csv_writer.write_record(&[format!("{}_{}_combo", group, n_label), "Count".to_string(), "Frequency (%)".to_string()])?;
                
                let mut sorted_combos: Vec<_> = counter.iter().collect();
                sorted_combos.sort_by(|a, b| b.1.cmp(a.1));

                for (combo, count) in sorted_combos {
                    let freq = (*count as f64 / total as f64) * 100.0;
                    csv_writer.write_record(&[String::from_utf8_lossy(combo).to_string(), count.to_string(), format!("{:.1$}", freq, dig as usize)])?;
                }
                println!("[完成] {}: 共找到 {} 个匹配, {} 种独特组合。", ref_id, total, counter.len());
            }
        }

        for (_, mut writer) in writers {
            writer.flush()?;
        }

        Ok(())
    }

    pub fn run(args: Args) -> Result<()> {
        std::fs::create_dir_all(&args.output)
            .with_context(|| format!("无法创建输出目录: {:?}", args.output))?;
        
        let ref_file = File::open(&args.ref_seq)?;
        let ref_reader = BufReader::new(ref_file);
        let ref_records: Vec<_> = fasta::Reader::new(ref_reader).records().collect::<Result<_,_>>()?;
        
        let ref_data_vec: Vec<RefData> = ref_records.into_iter().filter_map(|rec| {
            let seq = rec.seq().to_ascii_uppercase();
            let n_blocks = find_n_blocks(&seq);
            if n_blocks.is_empty() {
                println!("[跳过] {}: 参考序列中未找到 'N' 区块。", rec.id());
                return None;
            }
            Some(RefData {
                id: rec.id().to_string(),
                len: seq.len(),
                seq,
                n_blocks,
            })
        }).collect();
        
        println!("---> 开始并行比对 {} 条有效参考序列...", ref_data_vec.len());
        
        rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global()?;
        
        let pb = ProgressBar::new_spinner();
        pb.enable_steady_tick(std::time::Duration::from_millis(120));
        pb.set_style(
            ProgressStyle::default_spinner()
                .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
                .template("{spinner:.blue} [{elapsed_precise}] {msg} {pos:>10} reads")?,
        );
        pb.set_message("正在读取reads...");

        let args_arc = Arc::new(args);
        let ref_data_arc = Arc::new(ref_data_vec);

        thread::scope(|s| -> Result<()> {
            let (reads_tx, reads_rx) = crossbeam_channel::bounded::<Vec<Record>>(args_arc.threads * 2);
            let (results_tx, results_rx) = crossbeam_channel::bounded::<MatchResult>(1024);

            let ref_data_for_collector: HashMap<_, _> = ref_data_arc.iter().map(|d| (d.id.clone(), d.n_blocks.clone())).collect();
            
            let collector_args = Arc::clone(&args_arc);
            let collector_handle = s.spawn(move || {
                collector_thread(results_rx, collector_args.output.clone(), collector_args.group.clone(), collector_args.dig, collector_args.extract_matches, ref_data_for_collector)
            });

            for _ in 0..args_arc.threads {
                let rx = reads_rx.clone();
                let tx = results_tx.clone();
                let refs = Arc::clone(&ref_data_arc);
                let args_clone = Arc::clone(&args_arc);

                s.spawn(move || {
                    for read_chunk in rx {
                        for read_record in read_chunk {
                            let read_seq = read_record.seq().to_ascii_uppercase();
                            if read_seq.contains(&b'N') { continue; }

                            for ref_data in refs.iter() {
                                if let Some(combo) = find_alignment(&read_seq, ref_data, args_clone.mismatches, args_clone.min_len_ratio, false) {
                                    if tx.send(MatchResult { ref_id: ref_data.id.clone(), combo, read_record: read_record.clone() }).is_ok() {
                                        break;
                                    }
                                }
                                let rc_read = bio::alphabets::dna::revcomp(&read_seq);
                                if let Some(combo) = find_alignment(&rc_read, ref_data, args_clone.mismatches, args_clone.min_len_ratio, true) {
                                     if tx.send(MatchResult { ref_id: ref_data.id.clone(), combo, read_record: read_record.clone() }).is_ok() {
                                        break;
                                    }
                                }
                            }
                        }
                    }
                });
            }
            drop(results_tx);

            let reads_file = File::open(&args_arc.reads)?;
            let reads_reader = BufReader::new(reads_file);
            let boxed_reads_reader: Box<dyn BufRead> = if args_arc.reads.extension().map_or(false, |ext| ext == "gz") {
                Box::new(BufReader::new(MultiGzDecoder::new(reads_reader)))
            } else {
                Box::new(reads_reader)
            };
            let mut records_iter = fasta::Reader::new(boxed_reads_reader).records();
            
            loop {
                let mut chunk = Vec::with_capacity(CHUNK_SIZE);
                for _ in 0..CHUNK_SIZE {
                    match records_iter.next() {
                        Some(Ok(record)) => chunk.push(record),
                        Some(Err(e)) => return Err(e.into()),
                        None => break,
                    }
                }
                if chunk.is_empty() { break; }
                pb.inc(chunk.len() as u64);
                if reads_tx.send(chunk).is_err() { break; }
            }
            drop(reads_tx);
            pb.finish_with_message("✔ reads读取完毕，等待比对完成...");

            collector_handle.join().unwrap()?;
            Ok(())
        })?;

        println!("\n✔ 所有比对任务已完成。");
        Ok(())
    }
}
