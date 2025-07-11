// --- Cargo.toml ---
// [package]
// name = "Hammer_fastx"
// version = "3.2.1"
// edition = "2021"
//
// [dependencies]
// bio = "1.2.0"
// clap = { version = "4.5.7", features = ["derive"] }
// csv = "1.3.0"
// flate2 = "1.0.30"
// indicatif = "0.17.8"
// anyhow = "1.0.86"
// rayon = "1.10.0"
// crossbeam-channel = "0.5.13"
// num_cpus = "1.16.0"

// --- src/main.rs ---
use anyhow::Result; 
use clap::Parser;
use std::path::PathBuf;

// ==================================================================================
// 主程序入口和子命令定义 (CLI Structure)
// ==================================================================================

#[derive(Parser, Debug)]
#[command(
    name = "Hammer_fastx",
    version = "v0.1.1",
    author = "CZH with the help of Gemini pro 2.5",
    about = "处理FASTX文件的工具"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Parser, Debug)]
enum Commands {
    /// 根据barcode拆分FASTQ文件
    Demux(demux::Args),
    /// 统计FASTA/FASTQ文件中的序列信息
    Stats(stats::Args),
    /// 根据序列长度过滤FASTA/FASTQ文件
    Filter(filter::Args),
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Demux(args) => {
            demux::run(args)?;
        }
        Commands::Stats(args) => {
            stats::run(args)?;
        }
        Commands::Filter(args) => {
            filter::run(args)?;
        }
    }
    Ok(())
}

// ==================================================================================
// `common` 模块: 存放共享的结构和函数
// ==================================================================================
mod common {
    use super::*;
    use anyhow::{anyhow, Result};
    use flate2::bufread::MultiGzDecoder;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::Path;

    #[derive(Parser, Debug)]
    pub struct GlobalArgs {
        #[arg(long, help = "输入文件 (FASTA/FASTQ, 可为.gz压缩)")]
        pub inputfile: PathBuf,
    }
    
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
        first_char_reader.read_exact(&mut buf)?;
        match buf[0] {
            b'>' => Ok(Format::Fasta),
            b'@' => Ok(Format::Fastq),
            _ => Err(anyhow!(
                "无法识别的文件格式，请确保文件以 '>' (FASTA) 或 '@' (FASTQ) 开头。"
            )),
        }
    }
}

// ==================================================================================
// `demux` 子命令模块
// ==================================================================================
mod demux {
    use super::common::GlobalArgs;
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
        #[command(flatten)]
        global: GlobalArgs,

        #[arg(long, help = "输出目录")]
        output: PathBuf,

        #[arg(long, help = "线程数", default_value_t = num_cpus::get_physical())]
        threads: usize,
        
        #[arg(short, long, help = "样本标签文件 (CSV或TSV格式)")]
        tags: PathBuf,
        
        #[arg(short = 'l', long, default_value_t = 8, help = "标签的长度")]
        tag_len: usize,
        
        #[arg(long, help = "激活此选项以修剪拆分后序列两端的标签")]
        trim: bool,
        
        #[arg(long, help = "将输出格式转换为FASTA (默认: FASTQ)")]
        out_fasta: bool,
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
            let rev_key = (r_tag_rc, f_tag);
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
        all_samples: HashSet<String>,
        out_fasta: bool,
    ) -> Result<HashMap<String, u64>> {
        let mut writers: HashMap<String, GenericWriter> = HashMap::new();
        let extension = if out_fasta { "fasta" } else { "fastq" };
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
                let writer = writers.get_mut(&sample_id).unwrap();
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
        let (lookup_map, mut all_samples) = load_tags(&args.tags, args.tag_len)?;
        all_samples.insert("unmatched".to_string());
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
            reader_thread(args_arc.global.inputfile.clone(), raw_tx, pb)?;
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
// `stats` 子命令模块
// ==================================================================================
mod stats {
    use super::common::{detect_format, Format, GlobalArgs};
    use anyhow::Result;
    use bio::io::{fasta, fastq};
    use clap::Parser;
    use flate2::bufread::MultiGzDecoder;
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    #[derive(Parser, Debug)]
    pub struct Args {
        #[command(flatten)]
        global: GlobalArgs,
    }

    pub fn run(args: Args) -> Result<()> {
        let input_path = &args.global.inputfile;
        println!("---> 正在统计文件: {}", input_path.display());
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
                    total_len += len;
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
                    total_len += len;
                    if len < min_len { min_len = len; }
                    if len > max_len { max_len = len; }
                }
            }
        };

        println!("\n==================== 统计结果 ====================");
        if count > 0 {
            println!("序列总数 (Total sequences): {}", count);
            println!("总碱基数 (Total bases):    {}", total_len);
            println!("最大序列长度 (Max length):   {}", max_len);
            println!("最小序列长度 (Min length):   {}", min_len);
            println!("平均序列长度 (Avg length):   {:.2}", total_len as f64 / count as f64);
        } else {
            println!("文件中未找到任何序列。");
        }
        println!("================================================");

        Ok(())
    }
}

// ==================================================================================
// `filter` 子命令模块
// ==================================================================================
mod filter {
    use super::common::{detect_format, Format, GlobalArgs};
    use anyhow::Result;
    use bio::io::{fasta, fastq};
    use clap::Parser;
    use flate2::bufread::MultiGzDecoder;
    use std::fs::File;
    use std::io::{self, BufRead, BufReader, BufWriter, Write};
    use std::path::PathBuf;

    #[derive(Parser, Debug)]
    pub struct Args {
        #[command(flatten)]
        global: GlobalArgs,

        #[arg(long, help = "输出文件 (默认: 标准输出)")]
        outfile: Option<PathBuf>,

        #[arg(short = 'm', long, help = "过滤掉小于该长度的序列")]
        min_len: Option<usize>,
        
        #[arg(short = 'M', long, help = "过滤掉大于该长度的序列")]
        max_len: Option<usize>,
    }

    pub fn run(args: Args) -> Result<()> {
        let input_path = &args.global.inputfile;
        let format = detect_format(input_path)?;
        let min_len = args.min_len.unwrap_or(0);
        let max_len = args.max_len.unwrap_or(usize::MAX);

        let file = File::open(input_path)?;
        let buf_reader = BufReader::new(file);
        let input_reader: Box<dyn BufRead> =
            if input_path.extension().map_or(false, |ext| ext == "gz") {
                Box::new(BufReader::new(MultiGzDecoder::new(buf_reader)))
            } else {
                Box::new(buf_reader)
            };

        let mut writer: Box<dyn Write> = if let Some(path) = args.outfile {
            Box::new(BufWriter::new(File::create(path)?))
        } else {
            Box::new(BufWriter::new(io::stdout().lock()))
        };

        match format {
            Format::Fasta => {
                let reader = fasta::Reader::new(input_reader);
                let mut writer = fasta::Writer::new(&mut writer);
                for result in reader.records() {
                    let record = result?;
                    let len = record.seq().len();
                    if len >= min_len && len <= max_len {
                        writer.write_record(&record)?;
                    }
                }
            }
            Format::Fastq => {
                let reader = fastq::Reader::new(input_reader);
                let mut writer = fastq::Writer::new(&mut writer);
                for result in reader.records() {
                    let record = result?;
                    let len = record.seq().len();
                    if len >= min_len && len <= max_len {
                        writer.write_record(&record)?;
                    }
                }
            }
        }
        Ok(())
    }
}

