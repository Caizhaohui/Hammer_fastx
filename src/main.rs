use anyhow::Result;
use clap::Parser;
use std::process::{Command, Stdio}; // For executing external commands

// ==================================================================================
// 模块声明 (Module declarations) - 已被移除
//
// (原第 8-17 行的 "mod pipeline;", "mod merge_pe;" 等已被删除，
// 因为所有模块都已在下方本文件中定义。)
//
// ==================================================================================

// ==================================================================================
// Main entry point and CLI command definitions
// ==================================================================================

#[derive(Parser, Debug)]
#[command(
    name = "Hammer_fastx",
    version,
    author = "Caizhaohui",
    about = "A versatile toolkit for FASTX file processing, including QC, merging, and demultiplexing.",
    help_template = "{name} v{version}\n{about}\n\n{usage-heading} {usage}\n\n{all-args}\n"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Parser, Debug)]
enum Commands{
    /// [Workflow] Run the complete pipeline from QC and merging to demultiplexing
    #[command(name = "demux_all")]
    DemuxAll(pipeline::Args),

    /// [Workflow] Quality control and merge paired-end data, with optional output formats
    #[command(name = "mergePE")]
    MergePE(merge_pe::Args),

    /// [Single Step] Demultiplex a merged FASTQ file based on barcodes
    #[command(name = "demux_only")]
    DemuxOnly(demux::Args),

    /// (Wrapper) Quality control paired-end FASTQ files using fastp
    Fastp(fastp::Args),

    /// (Wrapper) Merge paired-end reads using flash2
    Flash2(flash2::Args),

    /// Get sequence statistics from one or more FASTA/FASTQ files
    Stats(stats::Args),

    /// Filter one or more FASTA/FASTQ files based on sequence length
    Filter(filter::Args),

    #[command(name = "merge_file")]
    MergeFile(merge_file::Args),

    /// [Anchor Logic] Align reads to a reference with Ns using a strict anchor-based method
    #[command(name = "Ns_count")]
    NsCount(ns_count::Args),

    /// Translate DNA FASTA files to Amino Acid FASTA files
    #[command(name = "DNA2AA")]
    DNA2AA(dna2aa::Args),

    /// [NEW] Count Amino Acid mutations against a reference protein sequence
    #[command(name = "count_AA")]
    CountAA(count_aa::Args), // <-- 新添加的命令

    /// Find motif occurrences and extract flanks; counts unique per-read windows; supports reverse complement
    #[command(name = "find_seq", about = "Find motif occurrences and extract flanks; counts unique per-read windows; supports reverse complement")]
    FindSeq(find_seq::Args),
}

fn main() -> Result<()> {
    // FIX: Changed Cli.parse() to Cli::parse()
    let cli = Cli::parse();

    match cli.command {
        Commands::DemuxAll(args) => pipeline::run(args),
        Commands::MergePE(args) => merge_pe::run(args),
        Commands::DemuxOnly(args) => demux::run(args),
        Commands::Fastp(args) => fastp::run(args),
        Commands::Flash2(args) => flash2::run(args),
        Commands::Stats(args) => stats::run(args),
        Commands::Filter(args) => filter::run(args),
        Commands::MergeFile(args) => merge_file::run(args),
        Commands::NsCount(args) => ns_count::run(args),
        Commands::DNA2AA(args) => dna2aa::run(args),
        Commands::CountAA(args) => count_aa::run(args),
        Commands::FindSeq(args) => find_seq::run(args), // <-- 新添加的分支
    }
}

// ==================================================================================
// `pipeline` subcommand module (for `demux_all`)
// ==================================================================================
mod pipeline {
    use super::{demux, fastp, flash2};
    use anyhow::{Context, Result};
    use clap::Parser;
    use std::fs;
    use std::path::PathBuf;
    use std::time::Instant;

    #[derive(Parser, Debug)]
    #[command(name = "demux_all", about = "[Workflow] Run the complete pipeline: QC (fastp), merge (flash2), and demultiplex (demux)")]
    pub struct Args {
        #[arg(short = 'i', long, help = "Input file 1 (Raw Read1)")]
        pub in1: PathBuf,

        #[arg(short = 'I', long, help = "Input file 2 (Raw Read2)")]
        pub in2: PathBuf,

        #[arg(long, help = "Sample tags file for demultiplexing (CSV format)")]
        pub tags: PathBuf,

        #[arg(short = 'o', long, help = "Main output directory for all results and intermediate files")]
        pub output_dir: PathBuf,

        #[arg(long, help = "Delete intermediate files from fastp and flash2 upon successful completion")]
        pub cleanup: bool,

        #[arg(long, help = "Number of threads for fastp", default_value_t = 4)]
        pub fastp_threads: usize,

        #[arg(long, help = "Number of threads for flash2", default_value_t = 4)]
        pub flash_threads: usize,
        #[arg(long, help = "Minimum overlap length for flash2", default_value_t = 10)]
        pub min_overlap: usize,
        #[arg(long, help = "Maximum overlap length for flash2", default_value_t = 300)]
        pub max_overlap: usize,

        #[arg(long, help = "Number of threads for demux_only", default_value_t = num_cpus::get_physical())]
        pub demux_threads: usize,
        #[arg(short = 'l', long, help = "Tag length for demux_only", default_value_t = 8)]
        pub tag_len: usize,
        #[arg(long, help = "Activate tag trimming for demux_only")]
        pub trim: bool,
        #[arg(long, help = "Output in FASTA format after demux_only (default: FASTQ)")]
        pub out_fasta: bool,
    }

    pub fn run(args: Args) -> Result<()> {
        let total_start_time = Instant::now();
        println!("🚀 [Workflow] Starting hammer_fastx demux_all pipeline...");

        let fastp_dir = args.output_dir.join("01_fastp_out");
        let flash_dir = args.output_dir.join("02_flash2_out");
        let demux_dir = args.output_dir.join("03_demux_out");

        fs::create_dir_all(&args.output_dir)
            .with_context(|| format!("Failed to create main output directory: {:?}", args.output_dir))?;
        fs::create_dir_all(&fastp_dir)
            .with_context(|| format!("Failed to create fastp output directory: {:?}", fastp_dir))?;
        fs::create_dir_all(&flash_dir)
            .with_context(|| format!("Failed to create flash2 output directory: {:?}", flash_dir))?;
        
        println!("\n[Step 1/3] ➡️  Running fastp for quality control...");
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

        println!("\n[Step 2/3] ➡️  Running flash2 to merge reads...");
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

        println!("\n[Step 3/3] ➡️  Running demux_only to demultiplex...");
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

        if args.cleanup {
            println!("\n[Cleanup] Removing intermediate files...");
            fs::remove_dir_all(&fastp_dir)
                .with_context(|| format!("Failed to clean up fastp directory: {:?}", fastp_dir))?;
            fs::remove_dir_all(&flash_dir)
                .with_context(|| format!("Failed to clean up flash2 directory: {:?}", flash_dir))?;
            println!("✔ Cleanup complete.");
        }

        println!("\n🎉 [Workflow] All steps completed successfully! Total time: {:.2?}", total_start_time.elapsed());
        println!("Final demultiplexed results are in: {}", demux_dir.display());

        Ok(())
    }
}

// ==================================================================================
// `merge_pe` subcommand module
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
    #[command(name = "mergePE", about = "[Workflow] QC with fastp, then merge paired-end data with flash2")]
    pub struct Args {
        #[arg(short = 'i', long, help = "Input file 1 (Raw Read1)")]
        pub in1: PathBuf,
        #[arg(short = 'I', long, help = "Input file 2 (Raw Read2)")]
        pub in2: PathBuf,

        #[arg(short = 'o', long, help = "Output path for the final merged file")]
        pub outfile: PathBuf,

        #[arg(long, help = "Convert final output to FASTA format (default: FASTQ)")]
        pub out_fasta: bool,
        #[arg(long, help = "Delete intermediate files upon successful completion")]
        pub cleanup: bool,
        #[arg(long, help = "Directory for intermediate files (default: 'intermediates' in the output file's directory)")]
        pub temp_dir: Option<PathBuf>,

        #[arg(long, help = "Number of threads for fastp", default_value_t = 4)]
        pub fastp_threads: usize,

        #[arg(long, help = "Number of threads for flash2", default_value_t = 4)]
        pub flash_threads: usize,
        #[arg(long, help = "Minimum overlap length for flash2", default_value_t = 10)]
        pub min_overlap: usize,
        #[arg(long, help = "Maximum overlap length for flash2", default_value_t = 300)]
        pub max_overlap: usize,
    }

    pub fn run(args: Args) -> Result<()> {
        let total_start_time = Instant::now();
        println!("🚀 [Workflow] Starting hammer_fastx mergePE workflow...");

        let output_parent_dir = args.outfile.parent().ok_or_else(|| anyhow!("Could not get parent directory of output file"))?;
        fs::create_dir_all(output_parent_dir)
            .with_context(|| format!("Failed to create output directory: {:?}", output_parent_dir))?;

        let temp_dir = args.temp_dir.clone().unwrap_or_else(|| output_parent_dir.join("intermediates"));
        fs::create_dir_all(&temp_dir)
            .with_context(|| format!("Failed to create temporary directory: {:?}", temp_dir))?;

        println!("\n[Step 1/3] ➡️  Running fastp for quality control...");
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

        println!("\n[Step 2/3] ➡️  Running flash2 to merge reads...");
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

        println!("\n[Step 3/3] ➡️  Writing final output file...");
        let merged_fastq_path = temp_dir.join(format!("{}.extendedFrags.fastq", flash_prefix));
        
        let in_file = fs::File::open(&merged_fastq_path)
            .with_context(|| format!("Failed to open merged file: {:?}", merged_fastq_path))?;
        let in_reader = BufReader::new(in_file);
        let fastq_reader = fastq::Reader::new(in_reader);

        let out_file = fs::File::create(&args.outfile)
            .with_context(|| format!("Failed to create final output file: {:?}", args.outfile))?;

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
        println!("✔ Successfully wrote {} records to {}", records_written, args.outfile.display());

        if args.cleanup {
            println!("\n[Cleanup] Removing intermediate files...");
            fs::remove_dir_all(&temp_dir)
                .with_context(|| format!("Failed to clean up temporary directory: {:?}", temp_dir))?;
            println!("✔ Cleanup complete.");
        }

        println!("\n🎉 [Workflow] mergePE workflow completed successfully! Total time: {:.2?}", total_start_time.elapsed());
        Ok(())
    }
}

// ==================================================================================
// `fastp` subcommand module
// ==================================================================================
mod fastp {
    use super::{Command, Stdio};
    use anyhow::{anyhow, Context, Result};
    use clap::Parser;
    use std::path::PathBuf;

    #[derive(Parser, Debug)]
    #[command(
        name = "fastp",
        about = "(Wrapper) Quality control paired-end FASTQ files using fastp"
    )]
    pub struct Args {
        #[arg(short = 'i', long, help = "Input file 1 (Read1)")]
        pub in1: PathBuf,

        #[arg(short = 'I', long, help = "Input file 2 (Read2)")]
        pub in2: PathBuf,

        #[arg(short = 'o', long, help = "Output file 1 (Read1)")]
        pub out1: PathBuf,

        #[arg(short = 'O', long, help = "Output file 2 (Read2)")]
        pub out2: PathBuf,

        #[arg(short = 'h', long, help = "Specify path for HTML report")]
        pub html: Option<PathBuf>,

        #[arg(short = 'j', long, help = "Specify path for JSON report")]
        pub json: Option<PathBuf>,

        #[arg(short = 'R', long, help = "Report title", default_value = "fastp report")]
        pub report_title: String,

        #[arg(short = 't', long, help = "Number of threads (default: auto-detect)")]
        pub threads: Option<usize>,
    }

    fn command_exists(cmd: &str) -> bool {
        Command::new(cmd)
            .arg("--version")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok()
    }

    pub fn run(args: Args) -> Result<()> {
        println!("---> Starting fastp quality control...");

        if !command_exists("fastp") {
            return Err(anyhow!(
                "Error: 'fastp' executable not found.\nPlease ensure fastp is installed and in your system's PATH environment variable."
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

        println!("🔧 Executing command: {:?}", cmd);

        let status = cmd
            .status()
            .with_context(|| "Failed to execute fastp command. Please check if fastp is installed correctly.")?;

        if status.success() {
            println!("\n✔ fastp quality control completed successfully!");
            println!("   - Cleaned R1: {}", args.out1.display());
            println!("   - Cleaned R2: {}", args.out2.display());
            if let Some(html_path) = &args.html {
                println!("   - HTML Report: {}", html_path.display());
            }
            Ok(())
        } else {
            Err(anyhow!(
                "fastp execution failed with exit code: {:?}\nPlease check the fastp logs for detailed error information.",
                status.code()
            ))
        }
    }
}

// ==================================================================================
// `flash2` subcommand module
// ==================================================================================
mod flash2 {
    use super::{Command, Stdio};
    use anyhow::{anyhow, Context, Result};
    use clap::Parser;
    use std::path::PathBuf;

    #[derive(Parser, Debug)]
    #[command(
        name = "flash2",
        about = "(Wrapper) Merge paired-end reads using flash2"
    )]
    pub struct Args {
        #[arg(help = "Input Read1 file")]
        pub in1: PathBuf,
        
        #[arg(help = "Input Read2 file")]
        pub in2: PathBuf,

        #[arg(short = 'o', long, help = "Output file prefix", default_value = "out")]
        pub out_prefix: String,
        
        #[arg(short = 'd', long, help = "Output directory", default_value = ".")]
        pub out_dir: PathBuf,

        #[arg(short = 'm', long, help = "Minimum overlap length", default_value_t = 10)]
        pub min_overlap: usize,

        #[arg(short = 'M', long, help = "Maximum overlap length", default_value_t = 300)]
        pub max_overlap: usize,

        #[arg(short = 't', long, help = "Number of threads (default: 1)", default_value_t = 1)]
        pub threads: usize,
    }
    
    fn command_exists(cmd: &str) -> bool {
        // Use the same robust check as the 'fastp' module
        Command::new(cmd)
            .arg("--version")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok()
    }

    pub fn run(args: Args) -> Result<()> {
        println!("---> Starting flash2 read merging...");
        
        if !command_exists("flash2") {
            return Err(anyhow!(
                "Error: 'flash2' executable not found.\nPlease ensure flash2 is installed and in your system's PATH environment variable."
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

        println!("🔧 Executing command: {:?}", cmd);

        let status = cmd
            .status()
            .with_context(|| "Failed to execute flash2 command. Please check if flash2 is installed correctly.")?;

        if status.success() {
            println!("\n✔ flash2 merging completed successfully!");
            println!("   - Output directory: {}", args.out_dir.display());
            println!("   - Output prefix: {}", args.out_prefix);
            println!("   - Merged file: {}", args.out_dir.join(format!("{}.extendedFrags.fastq", args.out_prefix)).display());
            Ok(())
        } else {
            Err(anyhow!(
                "flash2 execution failed with exit code: {:?}\nPlease check the flash2 logs for detailed error information.",
                status.code()
            ))
        }
    }
}

// ==================================================================================
// `common` module: Shared utility functions
// ==================================================================================
mod common {
    use anyhow::{anyhow, Result};
    use flate2::bufread::MultiGzDecoder;
    use std::fs::File;
    use std::io::{BufRead, BufReader, Read};
    use std::path::Path;

    #[derive(Debug, PartialEq, Eq, Clone, Copy)]
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
        match first_char_reader.read_exact(&mut buf) {
            Ok(_) => match buf[0] {
                b'>' => Ok(Format::Fasta),
                b'@' => Ok(Format::Fastq),
                _ => Err(anyhow!(
                    "Cannot identify file format for {:?}. Please ensure it starts with '>' (FASTA) or '@' (FASTQ).",
                    path
                )),
            },
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
                Err(anyhow!("File is empty or unreadable: {:?}", path))
            }
            Err(e) => Err(e.into()),
        }
    }
}

// ==================================================================================
// `demux` subcommand module (for `demux_only`)
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
        #[arg(long, help = "Input FASTQ file (can be gzipped)")]
        pub inputfile: PathBuf,

        #[arg(long, help = "Output directory")]
        pub output: PathBuf,

        #[arg(long, help = "Number of threads", default_value_t = num_cpus::get_physical())]
        pub threads: usize,
        
        #[arg(short, long, help = "Sample tags file (CSV format: SampleID,F_tag,R_tag)")]
        pub tags: PathBuf,
        
        #[arg(short = 'l', long, default_value_t = 8, help = "Length of the tags")]
        pub tag_len: usize,
        
        #[arg(long, help = "Activate this flag to trim tags from both ends of the sequence")]
        pub trim: bool,
        
        #[arg(long, help = "Convert output to FASTA format (default: FASTQ)")]
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
            .with_context(|| format!("Failed to open tag file: {:?}", tag_file))?;
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
                "Tag file must contain the columns 'SampleID', 'F_tag', and 'R_tag'."
            ));
        }
        for result in rdr.records() {
            let record = result?;
            let sample_id = record.get(0).ok_or_else(|| anyhow!("Missing SampleID"))?.to_string();
            let f_tag = record.get(1).ok_or_else(|| anyhow!("Missing F_tag"))?.as_bytes().to_ascii_uppercase();
            let r_tag = record.get(2).ok_or_else(|| anyhow!("Missing R_tag"))?.as_bytes().to_ascii_uppercase();
            if f_tag.len() != tag_len || r_tag.len() != tag_len {
                return Err(anyhow!("Tag length for sample {} does not match the specified --tag-len {}", sample_id, tag_len));
            }
            all_samples.insert(sample_id.clone());
            let r_tag_rc = bio::alphabets::dna::revcomp(&r_tag);

            // Forward key: 5'-[F_tag]...[R_tag_rc]-3'
            let fwd_key = (f_tag.clone(), r_tag_rc.clone());
            lookup_map.insert(fwd_key, MatchInfo { sample_id: sample_id.clone(), orientation: Orientation::Forward });
            
            // Reverse key: 5'-[R_tag_rc]...[F_tag]-3'
            // FIX: The original code used f_tag_rc here, which was incorrect.
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
        pb.finish_with_message("✔ File reading complete");
        Ok(())
    }

    // This worker function processes a record
    // MODIFIED: Takes ownership of Record to avoid clones
    fn process_record(
        record: Record, // Takes ownership
        lookup_map: &HashMap<(Vec<u8>, Vec<u8>), MatchInfo>,
        args: &Args,
    ) -> (String, Record) { // Returns tuple, not Option
        let seq = record.seq();
        if seq.len() < args.tag_len * 2 {
            return ("unmatched".to_string(), record); // Move record
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
                    record // Move record
                };
                (match_info.sample_id.clone(), final_record)
            }
            None => ("unmatched".to_string(), record), // Move record
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
        println!("\n\n==================== Demultiplexing Summary (Multi-threaded) ====================");
        println!("Processing Time: {:.2?}", duration);
        println!("Total Reads Processed: {}", total_reads);
        if total_reads > 0 {
            let matched_percent = matched_reads as f64 * 100.0 / total_reads as f64;
            let unmatched_percent = *counts.get("unmatched").unwrap_or(&0) as f64 * 100.0 / total_reads as f64;
            println!("  - Matched Reads:       {:>10} ({:.2}%)", matched_reads, matched_percent);
            println!("  - Unmatched Reads: {:>10} ({:.2}%)", counts.get("unmatched").unwrap_or(&0), unmatched_percent);
            println!("--------------------------------------------------");
            let mut sorted_samples: Vec<_> = counts.into_iter().collect();
            sorted_samples.sort_by(|a, b| b.1.cmp(&a.1));
            for (sample, count) in sorted_samples {
                if sample != "unmatched" {
                    let sample_percent = count as f64 * 100.0 / total_reads as f64;
                    println!("  - Sample {}: {:>10} reads ({:.2}%)", sample, count, sample_percent);
                }
            }
        }
        println!("===================================================================================");
        println!("✔ Done! Results written to: {}", output_dir.display());
    }

    // Optimization: This function combines the original worker_thread and the rayon::par_bridge logic
    fn parallel_processing(
        rx_raw: crossbeam_channel::Receiver<RawChunk>,
        tx_processed: crossbeam_channel::Sender<ProcessedChunk>,
        lookup_map: Arc<HashMap<(Vec<u8>, Vec<u8>), MatchInfo>>,
        args: Arc<Args>,
    ) {
        // Use rayon's par_bridge to consume chunks from the channel in parallel
        rx_raw.into_iter().par_bridge().for_each(|chunk| {
            let processed_results: Vec<(String, Record)> = chunk
                .into_par_iter() // Process records within the chunk in parallel (moves records)
                .map(|record| process_record(record, &lookup_map, &args)) // Use map
                .collect();
            
            let mut processed_chunk: ProcessedChunk = HashMap::new();
            for (sample_id, record) in processed_results {
                processed_chunk.entry(sample_id).or_default().push(record);
            }

            if !processed_chunk.is_empty() {
                let _ = tx_processed.send(processed_chunk);
            }
        });
    }

    pub fn run(args: Args) -> Result<()> {
        let start_time = Instant::now();
        let output_dir = args.output.clone();
        std::fs::create_dir_all(&output_dir)
            .with_context(|| format!("Failed to create output directory: {:?}", output_dir))?;
        
        println!("---> Loading tags...");
        let (lookup_map, all_samples) = load_tags(&args.tags, args.tag_len)?;
        let lookup_map = Arc::new(lookup_map);
        let args_arc = Arc::new(args);
        
        // Configure rayon thread pool
        rayon::ThreadPoolBuilder::new().num_threads(args_arc.threads).build_global()?;

        let channel_capacity = args_arc.threads * 2;
        let (raw_tx, raw_rx) = crossbeam_channel::bounded::<RawChunk>(channel_capacity);
        let (processed_tx, processed_rx) = crossbeam_channel::bounded::<ProcessedChunk>(channel_capacity);
        
        let pb = ProgressBar::new_spinner();
        pb.enable_steady_tick(std::time::Duration::from_millis(120));
        pb.set_style(
            ProgressStyle::default_spinner()
                .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "a"])
                .template("{spinner:.blue} [{elapsed_precise}] {msg} {pos:>10} reads")?,
        );
        pb.set_message("Processing...");

        thread::scope(|s| -> Result<()> {
            let out_fasta_flag = args_arc.out_fasta;
            let output_dir_for_writer = output_dir.clone();

            // 1. Writer Thread
            let writer_handle = s.spawn(move || {
                writer_thread(processed_rx, output_dir_for_writer, all_samples, out_fasta_flag)
            });

            // 2. Parallel Processing (consuming from raw_rx, sending to processed_tx)
            let (lookup_clone, args_clone) = (lookup_map.clone(), args_arc.clone());
            let processing_handle = s.spawn(move || {
                parallel_processing(raw_rx, processed_tx, lookup_clone, args_clone);
            });

            // 3. Reader Thread (Main thread role, feeds raw_tx)
            // This will block until reading is done, then drop raw_tx
            let reader_res = reader_thread(args_arc.inputfile.clone(), raw_tx, pb);
            if let Err(e) = reader_res {
                eprintln!("Error in reader thread: {:?}", e);
            }

            // Wait for processing to finish
            processing_handle.join().unwrap(); 

            // Wait for writer to finish
            match writer_handle.join().unwrap() {
                Ok(counts) => print_summary(counts, start_time, &output_dir),
                Err(e) => eprintln!("Writer thread error: {:?}", e),
            }
            Ok(())
        })?;
        Ok(())
    }
}

// ==================================================================================
// `stats` subcommand module
// ==================================================================================
mod stats {
    use super::common::{detect_format, Format};
    use anyhow::Result;
    use bio::io::{fasta, fastq};
    use clap::Parser;
    use csv::Writer;
    use flate2::bufread::MultiGzDecoder;
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::{Path, PathBuf};

    #[derive(Parser, Debug)]
    pub struct Args {
        #[arg(long, help = "One or more input files (wildcards supported, e.g., '*.fasta')", required = true, num_args = 1..)]
        inputfile: Vec<PathBuf>,
        #[arg(long, help = "Output CSV file for per-sequence counts")]
        outfile: Option<PathBuf>,
    }
    
    struct FileStats {
        filename: String,
        count: u64,
        total_len: u64,
        min_len: usize,
        max_len: usize,
    }

    fn get_sample_name(path: &Path) -> String {
        let filename = path.file_name()
            .unwrap_or_default()
            .to_str()
            .unwrap_or_default();

        // Try to handle .fastq.gz, .fa.gz etc.
        let stem1 = Path::new(filename)
            .file_stem()
            .unwrap_or_default()
            .to_str()
            .unwrap_or_default();
        
        if stem1.ends_with(".fastq") || stem1.ends_with(".fa") || stem1.ends_with(".fasta") || stem1.ends_with(".fq") {
             Path::new(stem1)
                .file_stem()
                .unwrap_or_default()
                .to_str()
                .unwrap_or_default()
                .to_string()
        } else {
            stem1.to_string()
        }
    }
    
    fn print_stats_table(stats: &[FileStats]) {
        if stats.is_empty() {
            println!("No files processed or no sequences found.");
            return;
        }

        println!("\n====================================== Sequence Statistics Summary ======================================");
        println!("{:<30} {:>15} {:>18} {:>10} {:>10} {:>12}",
                 "Sample Name", "Total Seqs", "Total Bases", "Min Length", "Max Length", "Avg Length");
        println!("{:-<30} {:-<15} {:-<18} {:-<10} {:-<10} {:-<12}",
                 "", "", "", "", "", "");

        for s in stats {
            let avg_len = if s.count > 0 {
                s.total_len as f64 / s.count as f64
            } else {
                0.0
            };
            println!("{:<30} {:>15} {:>18} {:>10} {:>10} {:<12.2}",
                     s.filename, s.count, s.total_len, s.min_len, s.max_len, avg_len);
        }
        println!("===================================================================================================");
    }

    pub fn run(args: Args) -> Result<()> {
        let mut all_stats: Vec<FileStats> = Vec::new();
        let mut wtr_opt: Option<Writer<File>> = if let Some(path) = args.outfile.clone() {
            let mut w = Writer::from_path(path)?;
            w.write_record(["filename", "sequence", "count"])?;
            Some(w)
        } else { None };

        for input_path in &args.inputfile {
            println!("---> Processing: {}", input_path.display());
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
            let mut seq_counts: HashMap<String, u64> = HashMap::new();

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
                        let seq = String::from_utf8(record.seq().to_vec()).unwrap().trim().to_uppercase();
                        *seq_counts.entry(seq).or_insert(0) += 1;
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
                        let seq = String::from_utf8(record.seq().to_vec()).unwrap().trim().to_uppercase();
                        *seq_counts.entry(seq).or_insert(0) += 1;
                    }
                }
            };

            if let Some(wtr) = wtr_opt.as_mut() {
                let fname = get_sample_name(input_path);
                let mut entries: Vec<(String, u64)> = seq_counts.into_iter().collect();
                entries.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
                for (seq, c) in entries {
                    wtr.write_record([fname.clone(), seq, c.to_string()])?;
                }
            }
            
            all_stats.push(FileStats {
                filename: get_sample_name(input_path),
                count,
                total_len,
                min_len: if count > 0 { min_len } else { 0 },
                max_len,
            });
        }

        if let Some(wtr) = wtr_opt.as_mut() { wtr.flush()?; }
        print_stats_table(&all_stats);
        Ok(())
    }
}

// ==================================================================================
// `filter` subcommand module (MODIFIED FOR BATCH PROCESSING)
// ==================================================================================
mod filter {
    use super::common::{detect_format, Format};
    use anyhow::{anyhow, Context, Result};
    use bio::io::{fasta, fastq};
    use clap::Parser;
    use flate2::bufread::MultiGzDecoder;
    use std::fs::{self, File};
    use std::io::{self, BufRead, BufReader, BufWriter, Write};
    use std::path::{Path, PathBuf};

    #[derive(Parser, Debug)]
    #[command(name = "filter", about = "Filter FASTA/FASTQ files by length, either individually or in batches.")]
    #[clap(group(
        clap::ArgGroup::new("input_mode")
            .required(true)
            .args(["input_files", "input_dir"]),
    ))]
    pub struct Args {
        #[arg(long, help = "One or more input files to concatenate and filter", num_args = 1..)]
        input_files: Vec<PathBuf>,

        #[arg(long, help = "Input directory to batch process files")]
        input_dir: Option<PathBuf>,

        #[arg(long, help = "Output file (default: stdout, used with --input-files)")]
        outfile: Option<PathBuf>,

        #[arg(long, help = "Output directory (required with --input-dir)")]
        output_dir: Option<PathBuf>,

        #[arg(short = 'm', long, help = "Filter out sequences shorter than this length")]
        min_len: Option<usize>,
        
        #[arg(short = 'M', long, help = "Filter out sequences longer than this length")]
        max_len: Option<usize>,
    }

    /// Helper function to process a single stream (file)
    fn process_file_stream(
        input_reader: Box<dyn BufRead>,
        writer: &mut Box<dyn Write>,
        format: &Format,
        min_len: usize,
        max_len: usize,
    ) -> Result<u64> { // Returns count of records written
        let mut records_written = 0;
        match format {
            Format::Fasta => {
                let reader = fasta::Reader::new(input_reader);
                let mut fasta_writer = fasta::Writer::new(writer);
                for result in reader.records() {
                    let record = result?;
                    let len = record.seq().len();
                    if len >= min_len && len <= max_len {
                        fasta_writer.write_record(&record)?;
                        records_written += 1;
                    }
                }
            }
            Format::Fastq => {
                let reader = fastq::Reader::new(input_reader);
                let mut fastq_writer = fastq::Writer::new(writer);
                for result in reader.records() {
                    let record = result?;
                    let len = record.seq().len();
                    if len >= min_len && len <= max_len {
                        fastq_writer.write_record(&record)?;
                        records_written += 1;
                    }
                }
            }
        }
        Ok(records_written)
    }

    /// Generates the output filename with `_filtered` suffix
    fn get_output_filename(input_path: &Path) -> Result<(String, bool)> {
        let file_name = input_path.file_name()
            .ok_or_else(|| anyhow!("Failed to get file name from path: {:?}", input_path))?
            .to_str()
            .ok_or_else(|| anyhow!("File name contains invalid UTF-8: {:?}", input_path))?;

        let (stem, extension) = match file_name.rfind('.') {
            Some(dot_pos) => (&file_name[..dot_pos], &file_name[dot_pos..]),
            None => (file_name, "")
        };
        
        let (real_stem, real_ext, is_gz) = if extension == ".gz" {
             if let Some(dot_pos) = stem.rfind('.') {
                 // e.g., ("input", ".fastq", true) from "input.fastq.gz"
                 (&stem[..dot_pos], &stem[dot_pos..], true)
             } else {
                 // e.g., ("archive", "", true) from "archive.gz"
                 (stem, "", true)
             }
        } else {
            // e.g., ("input", ".fasta", false) from "input.fasta"
            (stem, extension, false)
        };
        
        // Check if this is a file type we want to process
        let extensions_to_process = [".fasta", ".fa", ".fastq", ".fq", ".fna"];
        let should_process = extensions_to_process.contains(&real_ext);
        
        let new_file_name = if is_gz {
            format!("{}_filtered{}{}", real_stem, real_ext, extension) // input_filtered.fastq.gz
        } else {
            format!("{}_filtered{}", real_stem, real_ext) // input_filtered.fasta
        };

        Ok((new_file_name, should_process))
    }


    pub fn run(args: Args) -> Result<()> {
        let min_len = args.min_len.unwrap_or(0);
        let max_len = args.max_len.unwrap_or(usize::MAX);

        // --- BRANCH 1: Batch processing from a directory ---
        if let Some(input_dir) = args.input_dir {
            let output_dir = args.output_dir.ok_or_else(|| {
                anyhow!("--output-dir is required when using --input-dir")
            })?;
            
            if args.outfile.is_some() {
                return Err(anyhow!("--outfile cannot be used with --input-dir. Use --output-dir instead."));
            }

            fs::create_dir_all(&output_dir)
                .with_context(|| format!("Failed to create output directory: {:?}", output_dir))?;

            println!("---> Starting batch filter in directory: {}", input_dir.display());

            for entry in fs::read_dir(input_dir)? {
                let entry = entry?;
                let input_path = entry.path();
                
                if input_path.is_file() {
                    let (new_file_name, should_process) = match get_output_filename(&input_path) {
                        Ok((name, process)) => (name, process),
                        Err(e) => {
                             println!("---> Skipping file {}: {}", input_path.display(), e); // <-- 修复：将 input_PANTS 改为 input_path
                             continue;
                        }
                    };

                    if !should_process {
                         println!("---> Skipping unsupported file type: {}", input_path.display());
                         continue;
                    }
                    
                    let output_path = output_dir.join(new_file_name);

                    // 2. Open reader
                    let format = match detect_format(&input_path) {
                         Ok(f) => f,
                         Err(e) => {
                             println!("---> Skipping file {}: {}", input_path.display(), e);
                             continue;
                         }
                    };
                    
                    let file = File::open(&input_path)?;
                    let buf_reader = BufReader::new(file);
                    let input_reader: Box<dyn BufRead> =
                        if input_path.extension().map_or(false, |ext| ext == "gz") {
                            Box::new(BufReader::new(MultiGzDecoder::new(buf_reader)))
                        } else {
                            Box::new(buf_reader)
                        };
                    
                    // 3. Open writer
                    let mut writer: Box<dyn Write> = Box::new(BufWriter::new(File::create(&output_path)?));

                    // 4. Process
                    println!("---> Filtering {} -> {}", input_path.display(), output_path.display());
                    let count = process_file_stream(input_reader, &mut writer, &format, min_len, max_len)
                        .with_context(|| format!("Failed to process file: {:?}", input_path))?;
                    println!("✔ Wrote {} records to {}", count, output_path.display());
                }
            }
            println!("🎉 Batch filtering complete.");

        // --- BRANCH 2: Original logic (concatenate and filter) ---
        } else if !args.input_files.is_empty() {
            if args.output_dir.is_some() {
                 return Err(anyhow!("--output-dir can only be used with --input-dir."));
            }
            
            let mut writer: Box<dyn Write> = if let Some(path) = args.outfile {
                Box::new(BufWriter::new(File::create(path)?))
            } else {
                Box::new(BufWriter::new(io::stdout().lock()))
            };

            let mut first_format: Option<Format> = None;
            let mut total_records = 0;

            for input_path in &args.input_files {
                eprintln!("---> Processing (and appending): {}", input_path.display());
                let format = detect_format(input_path)?;
                
                // Ensure all files are the same format when concatenating
                if let Some(ref first) = first_format {
                    if *first != format {
                         return Err(anyhow!(
                            "Mismatched formats: Cannot concatenate FASTA and FASTQ files into one output."
                         ));
                    }
                } else {
                    first_format = Some(format);
                }
                
                let file = File::open(input_path)?;
                let buf_reader = BufReader::new(file);
                let input_reader: Box<dyn BufRead> =
                    if input_path.extension().map_or(false, |ext| ext == "gz") {
                        Box::new(BufReader::new(MultiGzDecoder::new(buf_reader)))
                    } else {
                        Box::new(buf_reader)
                    };

                total_records += process_file_stream(input_reader, &mut writer, first_format.as_ref().unwrap(), min_len, max_len)
                    .with_context(|| format!("Failed to process file: {:?}", input_path))?;
            }
            eprintln!("✔ Total records written: {}", total_records);
        }
        // No 'else' needed, as clap's 'input_mode' group ensures one branch is taken
        
        Ok(())
    }
}

// ==================================================================================
// `merge_file` subcommand module
// ==================================================================================
mod merge_file {
    use super::common::{detect_format, Format};
    use anyhow::{anyhow, Context, Result};
    use bio::io::{fasta, fastq};
    use clap::Parser;
    use flate2::bufread::MultiGzDecoder;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::fs::File;
    use std::io::{BufRead, BufReader, BufWriter, Write};
    use std::path::PathBuf;
    use indicatif::{ProgressBar, ProgressStyle};
    use rand::seq::SliceRandom;
    use rand::thread_rng;

    #[derive(Parser, Debug)]
    #[command(name = "merge_file", about = "Merge multiple FASTA/FASTQ files with optional shuffle, concurrency, progress, and fastq→fasta conversion.")]
    pub struct Args {
        #[arg(long, num_args = 1.., help = "Input FASTA/FASTQ files (gz supported)")]
        pub input_files: Vec<PathBuf>,

        #[arg(long, help = "Output file (.fasta/.fastq or .gz)")]
        pub outfile: Option<PathBuf>,

        #[arg(long, help = "Keep input file order (default)")]
        pub keep_order: bool,

        #[arg(long, help = "Shuffle record order before writing")]
        pub shuffle: bool,

        #[arg(long, default_value_t = num_cpus::get_physical(), help = "Parallel read workers")]
        pub threads: usize,

        #[arg(long, default_value_t = 10000, help = "Chunk size per read batch")]
        pub chunk_size: usize,

        #[arg(long, help = "Convert FASTQ to FASTA before merging (if inputs are FASTQ)")]
        pub fastq_to_fasta: bool,

        #[arg(long, help = "Only perform FASTQ→FASTA conversion and write output (no merge)")]
        pub convert_only: bool,
    }

    pub fn run(args: Args) -> Result<()> {
        if args.input_files.is_empty() {
            return Err(anyhow!("No input files provided"));
        }

        // Determine output path (required unless convert_only prints to stdout in future)
        let outfile = args.outfile
            .clone()
            .ok_or_else(|| anyhow!("--outfile is required unless future support for stdout is added"))?;

        let first_format = detect_format(&args.input_files[0])?;
        for p in &args.input_files[1..] {
            let f = detect_format(p)?;
            if f != first_format {
                return Err(anyhow!("Input files must have the same format"));
            }
        }

        // If fastq_to_fasta is set, we will treat output as FASTA even if inputs are FASTQ
        let target_format = if args.convert_only || args.fastq_to_fasta { Format::Fasta } else { first_format };

        if args.convert_only && args.input_files.len() != 1 {
            return Err(anyhow!("--convert-only 仅支持单输入文件。如需合并请不要使用该选项"));
        }

        // Combine and optionally shuffle the list of files respecting keep_order/shuffle
        let mut files = args.input_files.clone();
        if args.shuffle && !args.keep_order {
            files.shuffle(&mut thread_rng());
        }
        let out_file = File::create(&outfile)
            .with_context(|| format!("Failed to create output file: {:?}", outfile))?;
        let out_writer: Box<dyn Write> = if outfile.extension().map_or(false, |ext| ext == "gz") {
            Box::new(GzEncoder::new(BufWriter::new(out_file), Compression::default()))
        } else {
            Box::new(BufWriter::new(out_file))
        };
        let mut out_writer = out_writer;

        let pb = ProgressBar::new(0);
        pb.set_style(
            ProgressStyle::default_spinner()
                .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
                .template("{spinner:.blue} {msg} {pos} records")?
        );
        pb.set_message("Merging records...");

        let mut total = 0u64;

        match (first_format, target_format) {
            (Format::Fasta, Format::Fasta) => {
                let mut out = fasta::Writer::new(&mut out_writer);
                for input_path in files {
                    let in_file = File::open(&input_path)
                        .with_context(|| format!("Failed to open input file: {:?}", input_path))?;
                    let buf_reader = BufReader::new(in_file);
                    let input_reader: Box<dyn BufRead> = if input_path.extension().map_or(false, |ext| ext == "gz") {
                        Box::new(BufReader::new(MultiGzDecoder::new(buf_reader)))
                    } else { Box::new(buf_reader) };
                    let reader = fasta::Reader::new(input_reader);
                    // Optionally parallelize by collecting chunks; here sequential writing keeps order
                    for result in reader.records() { let record = result?; out.write_record(&record)?; total += 1; pb.inc(1); }
                }
            }
            (Format::Fastq, Format::Fastq) => {
                let mut out = fastq::Writer::new(&mut out_writer);
                for input_path in files {
                    let in_file = File::open(&input_path)
                        .with_context(|| format!("Failed to open input file: {:?}", input_path))?;
                    let buf_reader = BufReader::new(in_file);
                    let input_reader: Box<dyn BufRead> = if input_path.extension().map_or(false, |ext| ext == "gz") {
                        Box::new(BufReader::new(MultiGzDecoder::new(buf_reader)))
                    } else { Box::new(buf_reader) };
                    let reader = fastq::Reader::new(input_reader);
                    let chunk_size = args.chunk_size;
                    let mut records_iter = reader.records();
                    loop {
                        let mut chunk = Vec::with_capacity(chunk_size);
                        for _ in 0..chunk_size {
                            match records_iter.next() { Some(Ok(r)) => chunk.push(r), Some(Err(e)) => return Err(e.into()), None => break }
                        }
                        if chunk.is_empty() { break; }
                        if args.shuffle { chunk.shuffle(&mut thread_rng()); }
                        // Parallel write is unsafe due to single writer; we parallel map then write sequentially
                        for rec in chunk { out.write_record(&rec)?; total += 1; pb.inc(1); }
                    }
                }
            }
            (Format::Fastq, Format::Fasta) => {
                let mut out = fasta::Writer::new(&mut out_writer);
                for input_path in files {
                    let in_file = File::open(&input_path)
                        .with_context(|| format!("Failed to open input file: {:?}", input_path))?;
                    let buf_reader = BufReader::new(in_file);
                    let input_reader: Box<dyn BufRead> = if input_path.extension().map_or(false, |ext| ext == "gz") {
                        Box::new(BufReader::new(MultiGzDecoder::new(buf_reader)))
                    } else { Box::new(buf_reader) };
                    let reader = fastq::Reader::new(input_reader);
                    let chunk_size = args.chunk_size;
                    let mut records_iter = reader.records();
                    loop {
                        let mut chunk = Vec::with_capacity(chunk_size);
                        for _ in 0..chunk_size {
                            match records_iter.next() { Some(Ok(r)) => chunk.push(r), Some(Err(e)) => return Err(e.into()), None => break }
                        }
                        if chunk.is_empty() { break; }
                        if args.shuffle { chunk.shuffle(&mut thread_rng()); }
                        for rec in chunk {
                            let fasta_rec = fasta::Record::with_attrs(rec.id(), rec.desc(), rec.seq());
                            out.write_record(&fasta_rec)?; total += 1; pb.inc(1);
                        }
                    }
                }
            }
            (Format::Fasta, Format::Fastq) => {
                return Err(anyhow!("Cannot convert FASTA to FASTQ because quality scores are unavailable"));
            }
        }

        pb.finish_with_message("✔ Merging complete");
        println!("✔ Processed {} records into {}", total, outfile.display());
        Ok(())
    }
}

// ==================================================================================
// `ns_count` subcommand module (Restored v0.5.1 anchor-based logic with syntax fix)
// ==================================================================================
mod ns_count {
    use anyhow::{Context, Result};
    use bio::io::fasta::{self, Record};
    use clap::Parser;
    use flate2::bufread::MultiGzDecoder;
    use indicatif::{ProgressBar, ProgressStyle};
    use std::collections::{HashMap, HashSet};
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::PathBuf;
    use std::sync::Arc;
    use std::thread;

    const CHUNK_SIZE: usize = 4096;

    #[derive(Parser, Debug)]
    pub struct Args {
        #[arg(long, help = "FASTA file containing reads to be aligned (can be gzipped)")]
        reads: PathBuf,
        #[arg(long = "refSEQ", help = "FASTA file containing the reference sequence with N-regions")]
        ref_seq: PathBuf,
        #[arg(long, help = "Output directory for CSV files")]
        output: PathBuf,
        #[arg(long, help = "Number of threads", default_value_t = num_cpus::get_physical())]
        threads: usize,
        #[arg(long, help = "Prefix for column headers in the output CSV", default_value = "T0")]
        group: String,
        #[arg(long, help = "Number of decimal places for frequency", default_value_t = 2)]
        dig: u8,
        #[arg(long, help = "Maximum mismatches allowed in non-anchor regions", default_value_t = 2)]
        mismatches: usize,
        #[arg(long, help = "Length of the anchor region on each side of an N-block", default_value_t = 15)]
        anchor_len: usize,
        #[arg(long, help = "Extract all matching reads into a separate FASTA file")]
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
        anchor_indices: HashSet<usize>,
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

    fn calculate_anchor_indices(n_blocks: &[(usize, usize)], ref_len: usize, anchor_len: usize) -> HashSet<usize> {
        let mut indices = HashSet::new();
        for &(start, len) in n_blocks {
            let upstream_start = start.saturating_sub(anchor_len);
            indices.extend(upstream_start..start);

            let downstream_start = start + len;
            let downstream_end = (downstream_start + anchor_len).min(ref_len);
            indices.extend(downstream_start..downstream_end);
        }
        indices
    }

    fn find_alignment(read_seq: &[u8], ref_data: &RefData, args: &Arc<Args>, is_rc_read: bool) -> Option<Vec<u8>> {
        let read_len = read_seq.len();
        let ref_len = ref_data.len;

        for ref_start in 0..=ref_len.saturating_sub(read_len) {
            let overlap_len = read_len;

            if !ref_data.n_blocks.iter().all(|(n_start, n_len)| 
                *n_start >= ref_start && (*n_start + *n_len) <= (ref_start + overlap_len)
            ) {
                continue;
            }

            let mut anchor_mismatch = false;
            for &anchor_idx in &ref_data.anchor_indices {
                if anchor_idx >= ref_start && anchor_idx < (ref_start + overlap_len) {
                    let read_idx = anchor_idx - ref_start;
                    if read_seq[read_idx] != ref_data.seq[anchor_idx] {
                        anchor_mismatch = true;
                        break;
                    }
                }
            }
            if anchor_mismatch { continue; }

            let mut mismatches = 0;
            for i in 0..overlap_len {
                let ref_idx = ref_start + i;
                if ref_data.anchor_indices.contains(&ref_idx) || ref_data.seq[ref_idx] == b'N' {
                    continue;
                }
                let read_idx = i;
                if read_seq[read_idx] != ref_data.seq[ref_idx] {
                    mismatches += 1;
                }
            }

            if mismatches <= args.mismatches {
                let mut combo_parts = Vec::new();
                for &(n_start, n_len) in &ref_data.n_blocks {
                    let read_idx_start = n_start - ref_start;
                    let segment = &read_seq[read_idx_start..read_idx_start + n_len];
                    if is_rc_read {
                        combo_parts.push(bio::alphabets::dna::revcomp(segment));
                    } else {
                        combo_parts.push(segment.to_vec());
                    }
                }
                return Some(combo_parts.join(&b'-'));
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
                sorted_combos.sort_by(|a, b| b.1.cmp(&a.1));

                for (combo, count) in sorted_combos {
                    let freq = (*count as f64 / total as f64) * 100.0;
                    csv_writer.write_record(&[String::from_utf8_lossy(combo).to_string(), count.to_string(), format!("{:.1$}", freq, dig as usize)])?;
                }
                println!("[Done] {}: Found {} matches with {} unique combinations.", ref_id, total, counter.len());
            }
        }

        for (_, mut writer) in writers {
            writer.flush()?;
        }

        Ok(())
    }

    pub fn run(args: Args) -> Result<()> {
        std::fs::create_dir_all(&args.output)
            .with_context(|| format!("Failed to create output directory: {:?}", args.output))?;
        
        let ref_file = File::open(&args.ref_seq)?;
        let ref_reader = BufReader::new(ref_file);
        let ref_records: Vec<_> = fasta::Reader::new(ref_reader).records().collect::<Result<_,_>>()?;
        
        let args_arc = Arc::new(args);

        let ref_data_vec: Vec<RefData> = ref_records.into_iter().filter_map(|rec| {
            let seq = rec.seq().to_ascii_uppercase();
            let n_blocks = find_n_blocks(&seq);
            if n_blocks.is_empty() {
                println!("[Skipping] {}: No 'N' blocks found in reference sequence.", rec.id());
                return None;
            }
            let anchor_indices = calculate_anchor_indices(&n_blocks, seq.len(), args_arc.anchor_len);
            Some(RefData {
                id: rec.id().to_string(),
                len: seq.len(),
                seq,
                n_blocks,
                anchor_indices,
            })
        }).collect();
        
        println!("---> Starting parallel alignment against {} valid reference(s)...", ref_data_vec.len());
        
        rayon::ThreadPoolBuilder::new().num_threads(args_arc.threads).build_global()?;
        
        let pb = ProgressBar::new_spinner();
        pb.enable_steady_tick(std::time::Duration::from_millis(120));
        pb.set_style(
            ProgressStyle::default_spinner()
                .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
                .template("{spinner:.blue} [{elapsed_precise}] {msg} {pos:>10} reads")?,
        );
        pb.set_message("Reading reads...");

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

                            'ref_loop: for ref_data in refs.iter() {
                                if let Some(combo) = find_alignment(&read_seq, ref_data, &args_clone, false) {
                                    if tx.send(MatchResult { ref_id: ref_data.id.clone(), combo, read_record: read_record.clone() }).is_ok() {
                                        break 'ref_loop;
                                    }
                                }
                                let rc_read = bio::alphabets::dna::revcomp(&read_seq);
                                if let Some(combo) = find_alignment(&rc_read, ref_data, &args_clone, true) {
                                    if tx.send(MatchResult { ref_id: ref_data.id.clone(), combo, read_record: read_record.clone() }).is_ok() {
                                        break 'ref_loop;
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
            pb.finish_with_message("✔ Reads loaded, waiting for alignment to finish...");

            collector_handle.join().unwrap()?;
            Ok(())
        })?;

        println!("\n✔ All alignment tasks are complete.");
        Ok(())
    }
}

// ==================================================================================
// `dna2aa` subcommand module (NEW)
// ==================================================================================
mod dna2aa {
    use anyhow::{anyhow, Context, Result};
    use bio::io::fasta; // 只导入 FASTA 读写器
    use clap::Parser;
    use rayon::prelude::*;
    use std::collections::HashMap;
    use std::fs::{self, File};
    use std::io::BufReader;
    use std::path::{Path, PathBuf};
    use std::sync::Arc;

    #[derive(Parser, Debug)]
    #[command(name = "DNA2AA", about = "Translate DNA FASTA files in a directory to Amino Acid FASTA files")]
    pub struct Args {
        #[arg(long, short = 'i', help = "Input directory containing DNA FASTA files")]
        pub input: PathBuf,

        #[arg(long, short = 'o', help = "Output directory for translated protein FASTA files")]
        pub output: PathBuf,

        #[arg(long, default_value_t = 50, help = "Minimum amino acid length to keep")]
        pub aa_length: usize,
    }

    // --------------------------------------------------------------------------------
    // 修复：手动实现标准密码子表，移除对 `bio` 库翻译功能的依赖
    // --------------------------------------------------------------------------------
    type CodonTable = HashMap<[u8; 3], u8>;

    /// 构建一个标准的DNA密码子表
    fn build_codon_table() -> CodonTable {
        let mut table = HashMap::new();
        // 终止密码子 (Stop Codons)
        table.insert(*b"TAA", b'*');
        table.insert(*b"TAG", b'*');
        table.insert(*b"TGA", b'*');
        // 苯丙氨酸 (F)
        table.insert(*b"TTT", b'F');
        table.insert(*b"TTC", b'F');
        // 亮氨酸 (L)
        table.insert(*b"TTA", b'L');
        table.insert(*b"TTG", b'L');
        table.insert(*b"CTT", b'L');
        table.insert(*b"CTC", b'L');
        table.insert(*b"CTA", b'L');
        table.insert(*b"CTG", b'L');
        // 异亮氨酸 (I)
        table.insert(*b"ATT", b'I');
        table.insert(*b"ATC", b'I');
        table.insert(*b"ATA", b'I');
        // 甲硫氨酸 (M) - 起始
        table.insert(*b"ATG", b'M');
        // 缬氨酸 (V)
        table.insert(*b"GTT", b'V');
        table.insert(*b"GTC", b'V');
        table.insert(*b"GTA", b'V');
        table.insert(*b"GTG", b'V');
        // 丝氨酸 (S)
        table.insert(*b"TCT", b'S');
        table.insert(*b"TCC", b'S');
        table.insert(*b"TCA", b'S');
        table.insert(*b"TCG", b'S');
        table.insert(*b"AGT", b'S');
        table.insert(*b"AGC", b'S');
        // 脯氨酸 (P)
        table.insert(*b"CCT", b'P');
        table.insert(*b"CCC", b'P');
        table.insert(*b"CCA", b'P');
        table.insert(*b"CCG", b'P');
        // 苏氨酸 (T)
        table.insert(*b"ACT", b'T');
        table.insert(*b"ACC", b'T');
        table.insert(*b"ACA", b'T');
        table.insert(*b"ACG", b'T');
        // 丙氨酸 (A)
        table.insert(*b"GCT", b'A');
        table.insert(*b"GCC", b'A');
        table.insert(*b"GCA", b'A');
        table.insert(*b"GCG", b'A');
        // 酪氨酸 (Y)
        table.insert(*b"TAT", b'Y');
        table.insert(*b"TAC", b'Y');
        // 组氨酸 (H)
        table.insert(*b"CAT", b'H');
        table.insert(*b"CAC", b'H');
        // 谷氨酰胺 (Q)
        table.insert(*b"CAA", b'Q');
        table.insert(*b"CAG", b'Q');
        // 天冬酰胺 (N)
        table.insert(*b"AAT", b'N');
        table.insert(*b"AAC", b'N');
        // 赖氨酸 (K)
        table.insert(*b"AAA", b'K');
        table.insert(*b"AAG", b'K');
        // 天冬氨酸 (D)
        table.insert(*b"GAT", b'D');
        table.insert(*b"GAC", b'D');
        // 谷氨酸 (E)
        table.insert(*b"GAA", b'E');
        table.insert(*b"GAG", b'E');
        // 半胱氨酸 (C)
        table.insert(*b"TGT", b'C');
        table.insert(*b"TGC", b'C');
        // 色氨酸 (W)
        table.insert(*b"TGG", b'W');
        // 精氨酸 (R)
        table.insert(*b"CGT", b'R');
        table.insert(*b"CGC", b'R');
        table.insert(*b"CGA", b'R');
        table.insert(*b"CGG", b'R');
        table.insert(*b"AGA", b'R');
        table.insert(*b"AGG", b'R');
        // 甘氨酸 (G)
        table.insert(*b"GGT", b'G');
        table.insert(*b"GGC", b'G');
        table.insert(*b"GGA", b'G');
        table.insert(*b"GGG", b'G');
        table
    }

    /// Translates a DNA sequence until the first stop codon (which is not included).
    /// Mimics Biopython's `seq.translate(to_stop=True)`
    fn translate_to_stop(dna_seq: &[u8], table: &CodonTable) -> Vec<u8> {
        let mut protein = Vec::new();

        // 遍历3碱基的密码子
        for codon_bytes in dna_seq.chunks_exact(3) {
            // 将 &[u8] 转换为 [u8; 3]
            let codon: [u8; 3] = [
                codon_bytes[0].to_ascii_uppercase(), 
                codon_bytes[1].to_ascii_uppercase(), 
                codon_bytes[2].to_ascii_uppercase()
            ];

            match table.get(&codon) {
                Some(aa) => {
                    if *aa == b'*' {
                        // 找到终止密码子，停止翻译
                        break;
                    }
                    // 添加氨基酸
                    protein.push(*aa);
                }
                None => {
                    // 无法识别的密码子 (例如 'N')，翻译为 'X'
                    protein.push(b'X');
                }
            }
        }
        protein
    }
    // --------------------------------------------------------------------------------
    // 修复结束
    // --------------------------------------------------------------------------------

    /// Processes a single FASTA file: translates it and saves the result.
    fn process_single_file(
        input_path: &Path,
        output_dir: &Path,
        min_aa_length: usize,
        table: &CodonTable, // <-- 接收密码子表
    ) -> Result<()> {
        // 1. Determine output path
        let file_stem = input_path
            .file_stem()
            .ok_or_else(|| anyhow!("Could not get file stem for {:?}", input_path))?;
        let output_filename = format!("{}_protein.fasta", file_stem.to_string_lossy());
        let output_path = output_dir.join(output_filename);

        // 2. Setup reader and writer
        let file = File::open(input_path)
            .with_context(|| format!("Failed to open input file: {:?}", input_path))?;
        let reader = fasta::Reader::new(BufReader::new(file));
        let mut writer = fasta::Writer::to_file(&output_path)
            .with_context(|| format!("Failed to create output file: {:?}", output_path))?;

        let mut records_written = 0;

        // 3. Translation logic
        for result in reader.records() {
            let record = result?;
            
            // Translate the DNA sequence, stopping at the first STOP codon
            let protein = translate_to_stop(record.seq(), table); // <-- 传入密码子表

            if protein.len() >= min_aa_length {
                // Create a new FASTA record for the protein
                let aa_record =
                    fasta::Record::with_attrs(record.id(), None, &protein);
                writer.write_record(&aa_record)?;
                records_written += 1;
            }
        }

        if records_written > 0 {
             println!(
                "Processed {:?} -> {:?} (Wrote {} records)",
                input_path.file_name().unwrap_or_default(),
                output_path.file_name().unwrap_or_default(),
                records_written
            );
        }

        Ok(())
    }

    /// Main run function for the DNA2AA subcommand
    pub fn run(args: Args) -> Result<()> {
        let start_time = std::time::Instant::now();
        
        // 1. 创建密码子表并用 Arc 包装，以便安全地跨线程共享
        let codon_table = Arc::new(build_codon_table());

        // 2. Create output directory
        fs::create_dir_all(&args.output)
            .with_context(|| format!("Failed to create output directory: {:?}", args.output))?;

        // 3. Find all input files
        let input_files: Vec<PathBuf> = fs::read_dir(&args.input)
            .with_context(|| format!("Failed to read input directory: {:?}", args.input))?
            .filter_map(|entry_result| {
                let entry = entry_result.ok()?;
                let path = entry.path();
                if path.is_file() {
                    if let Some(ext_str) = path.extension().and_then(|s| s.to_str()) {
                        // Match common FASTA extensions
                        if ext_str == "fasta" || ext_str == "fa" || ext_str == "fna" {
                            return Some(path);
                        }
                    }
                }
                None
            })
            .collect();
        
        if input_files.is_empty() {
             println!("Warning: No FASTA files (.fasta, .fa, .fna) found in {:?}.", args.input);
             return Ok(());
        }

        println!(
            "---> Found {} FASTA files to process in parallel...",
            input_files.len()
        );

        // 4. Process files in parallel (similar to Python's ProcessPoolExecutor)
        input_files.par_iter().for_each(|input_path| {
            // 为每个线程克隆 Arc 引用（开销很小）
            let table_clone = Arc::clone(&codon_table);
            if let Err(e) = process_single_file(input_path, &args.output, args.aa_length, &table_clone) {
                // Print errors from within the parallel loop
                eprintln!("\n[Error] Failed to process file {:?}: {}\n", input_path.display(), e);
            }
        });

        println!("\n🎉 All files processed successfully! Total time: {:.2?}", start_time.elapsed());
        println!("Results are in: {}", args.output.display());
        Ok(())
    }
}

// ==================================================================================
// `count_AA` subcommand module (NEWLY ADDED)
// ==================================================================================
mod count_aa {
    use anyhow::{anyhow, Context, Result};
    use bio::io::fasta::{self, Record};
    use clap::Parser;
    use crossbeam_channel::bounded;
    use dashmap::DashMap; // For concurrent counting
    use glob::glob; // For file matching
    use rayon::prelude::*; // For parallel iteration
    use std::collections::HashSet; // <-- 修复：移除未使用的 HashMap
    use std::fs::{self, File};
    use std::io::BufReader;
    use std::path::{Path, PathBuf};
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::sync::Arc;
    use std::thread;
    use std::time::Instant;

    #[derive(Parser, Debug)]
    #[command(name = "count_AA", about = "[NEW] Count AA mutations against a reference protein sequence, replicating the logic from Count_AAmutants.py")]
    pub struct Args {
        #[arg(short = 'r', long, help = "参考蛋白FASTA序列 (Reference protein FASTA sequence)")]
        pub reference: PathBuf,

        #[arg(short = 'i', long, help = "包含多个FASTA文件的目录 (Directory containing multiple FASTA files)")]
        pub input_dir: PathBuf,

        #[arg(short = 'o', long, help = "输出CSV文件的目录 (Output directory for CSV files)")]
        pub output_dir: PathBuf,

        #[arg(short = 'A', long, help = "位置偏移量 (Position offset)", default_value_t = 0)]
        pub aa_offset: i32,

        #[arg(short = 'c', long, help = "CSV配置文件，包含protected_sites列 (CSV config file with 'protected_sites' column)")]
        pub config: Option<PathBuf>,

        #[arg(long, help = "匹配参考起始位置的氨基酸数量 (Number of AAs to match ref start)", default_value_t = 6)]
        pub match_len: usize,

        #[arg(long, help = "每个文件的并行线程数 (Parallel threads per file)", default_value_t = 8)]
        pub threads: usize,

        #[arg(long, help = "每块reads数量 (Number of reads per chunk)", default_value_t = 100000)]
        pub chunk_size: usize,
    }

    /// (Helper) Loads the first sequence from a FASTA file.
    fn load_reference_sequence(path: &Path) -> Result<Vec<u8>> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open reference file: {:?}", path))?;
        let reader = fasta::Reader::new(BufReader::new(file)); // <-- 修复：移除 mut
        let record = reader
            .records()
            .next()
            .ok_or_else(|| anyhow!("Reference FASTA file is empty: {:?}", path))??;
        
        Ok(record.seq().to_ascii_uppercase())
    }

    /// (Helper) Loads protected sites from the config CSV.
    fn load_config(path: &Option<PathBuf>) -> Result<HashSet<usize>> {
        let mut protected = HashSet::new();
        let config_path = match path {
            Some(p) => p,
            None => return Ok(protected), // No config, return empty set
        };

        println!("---> Loading config: {}", config_path.display());
        let file = File::open(config_path)
            .with_context(|| format!("Failed to open config file: {:?}", config_path))?;
        let mut rdr = csv::Reader::from_reader(file);

        // Find the 'protected_sites' column index
        let headers = rdr.headers()?.clone();
        let site_col_idx = headers.iter().position(|h| h == "protected_sites");

        let site_col_idx = match site_col_idx {
            Some(idx) => idx,
            None => {
                println!("Warning: Config file provided, but 'protected_sites' column not found.");
                return Ok(protected);
            }
        };

        // Iterate through records and parse sites
        for result in rdr.records() {
            let record = result?;
            if let Some(val_str) = record.get(site_col_idx) {
                if val_str.is_empty() { continue; }
                
                match val_str.trim().parse::<f64>() { // Parse as f64 to match Python's behavior
                    Ok(val) => {
                        let site_1_based = val as usize;
                        if site_1_based > 0 {
                            protected.insert(site_1_based - 1); // Convert to 0-based index
                        }
                    }
                    Err(_) => {
                         println!("Warning: Could not parse protected site value '{}'", val_str);
                    }
                }
            }
        }
        
        println!("Loaded {} protected sites from config.", protected.len());
        Ok(protected)
    }

    /// (Helper) This is the core logic from the Python `analyze_chunk` function.
    /// It processes a chunk of reads and updates the global concurrent counters.
    fn analyze_chunk(
        reference_seq: &Vec<u8>,
        reads: Vec<fasta::Record>,
        protected_sites: &HashSet<usize>,
        match_len: usize,
        aa_counts: &[DashMap<u8, AtomicU64>], // A slice of concurrent maps
        total_reads: &AtomicU64,
        total_valid: &AtomicU64,
    ) {
        let seq_len = reference_seq.len();
        total_reads.fetch_add(reads.len() as u64, Ordering::Relaxed);
        let mut local_valid_reads = 0;

        for record in reads {
            let read = record.seq().to_ascii_uppercase();
            if read.is_empty() { continue; }

            // Prevent panic on reads shorter than match_len
            let match_segment_len = match_len.min(read.len());
            let read_start_segment = &read[..match_segment_len];
            
            if read_start_segment.is_empty() { continue; }

            // Find start position (Rust equivalent of Python's `str.find()`)
            let ref_start_pos = reference_seq
                .windows(read_start_segment.len())
                .position(|window| window == read_start_segment);

            if ref_start_pos.is_none() {
                continue; // Not found
            }
            
            let ref_start = ref_start_pos.unwrap();
            let mut violate = false;
            let mut mutation_count = 0;
            
            // We store (pos, aa) pairs and only update counts *after* checking for violations.
            // This ensures we don't count AAs from invalid (violated) reads.
            let mut pos_counts: Vec<(usize, u8)> = Vec::with_capacity(read.len());

            for (i, &aa) in read.iter().enumerate() {
                let pos = ref_start + i;
                if pos >= seq_len {
                    break; // Read is longer than remaining ref
                }
                
                let ref_aa = reference_seq[pos];

                if aa != ref_aa {
                    if protected_sites.contains(&pos) {
                        violate = true;
                        break; // Stop checking this read immediately
                    }
                    mutation_count += 1;
                }
                pos_counts.push((pos, aa));
            }

            // Only if the read is valid do we add its counts to the global map
            if !violate {
                for (pos, aa) in pos_counts {
                    // Find the concurrent map for this position
                    // Get or create an AtomicU64 counter for this AA
                    // Increment the counter atomically
                    aa_counts[pos]
                        .entry(aa)
                        .or_insert_with(|| AtomicU64::new(0))
                        .fetch_add(1, Ordering::Relaxed);
                }
                
                // Check the *other* condition for a "valid read"
                if mutation_count <= 1 {
                    local_valid_reads += 1;
                }
            }
        }
        
        // Atomically update the global "valid" counter
        total_valid.fetch_add(local_valid_reads, Ordering::Relaxed);
    }


    /// Main run function for the count_AA subcommand
    pub fn run(args: Args) -> Result<()> {
        let main_start_time = Instant::now();
        
        fs::create_dir_all(&args.output_dir)
            .with_context(|| format!("Failed to create output directory: {:?}", args.output_dir))?;

        // 1. Load Reference and Config
        let reference_seq = Arc::new(load_reference_sequence(&args.reference)?);
        let protected_sites = Arc::new(load_config(&args.config)?);
        println!("Reference sequence loaded ({} AAs).", reference_seq.len());

        // 2. Find input FASTA files (using `glob` crate)
        let pattern1 = args.input_dir.join("*.fasta").to_string_lossy().to_string();
        let pattern2 = args.input_dir.join("*.fa").to_string_lossy().to_string();
        
        let fasta_files: Vec<PathBuf> = glob(&pattern1)?
            .filter_map(Result::ok)
            .chain(glob(&pattern2)?.filter_map(Result::ok))
            .collect();

        if fasta_files.is_empty() {
            println!("No FASTA files (.fasta, .fa) found in {:?}.", args.input_dir);
            return Ok(());
        }

        println!("Processing {} FASTA files in parallel ({} threads per file)...", fasta_files.len(), args.threads);

        // 3. Configure Rayon global thread pool
        // This sets the *total* number of threads Rayon will use.
        rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global()?;

        // 4. Process each file (sequentially, as in Python)
        // The parallelism is *within* each file's chunk processing.
        for fasta_file in fasta_files {
            let file_start_time = Instant::now();
            let file_stem = fasta_file.file_stem().unwrap_or_default().to_string_lossy();
            println!("\n---> Processing file: {}", fasta_file.display());

            // --- Setup concurrent data structures for this file ---
            let seq_len = reference_seq.len();
            // Create a Vec of DashMaps, one for each position in the reference
            // Each DashMap stores: AA (u8) -> AtomicU64 (count)
            let global_counts: Arc<Vec<DashMap<u8, AtomicU64>>> = 
                Arc::new((0..seq_len).map(|_| DashMap::new()).collect());
            
            let total_reads = Arc::new(AtomicU64::new(0));
            let total_valid = Arc::new(AtomicU64::new(0));

            // Create Arcs for data to be shared across threads
            let reference_seq_clone = Arc::clone(&reference_seq);
            let protected_sites_clone = Arc::clone(&protected_sites);
            let global_counts_clone = Arc::clone(&global_counts);
            let total_reads_clone = Arc::clone(&total_reads);
            let total_valid_clone = Arc::clone(&total_valid);
            
            let (tx, rx) = bounded::<Vec<fasta::Record>>(args.threads * 2); // Channel for chunks of records

            // --- Use thread::scope for structured concurrency ---
            let res: Result<()> = thread::scope(|s| {
                // --- 1. Reader Thread ---
                // This thread reads the FASTA file and sends chunks of records to the channel
                let fasta_file_clone = fasta_file.clone();
                let chunk_size = args.chunk_size;
                s.spawn(move || {
                    let file = match File::open(&fasta_file_clone) {
                        Ok(f) => f,
                        Err(e) => {
                            eprintln!("Error opening {}: {:?}", fasta_file_clone.display(), e);
                            return; // Exit thread
                        }
                    };
                    let reader = BufReader::new(file);
                    let fasta_reader = fasta::Reader::new(reader); // <-- 修复：移除 mut
                    let mut records_iter = fasta_reader.records();

                    loop {
                        let mut chunk = Vec::with_capacity(chunk_size);
                        for _ in 0..chunk_size {
                            match records_iter.next() {
                                Some(Ok(record)) => chunk.push(record),
                                Some(Err(e)) => {
                                    eprintln!("Error reading record from {}: {}", fasta_file_clone.display(), e);
                                    // Continue to next record
                                },
                                None => break, // End of file
                            }
                        }
                        
                        let is_empty = chunk.is_empty();
                        if tx.send(chunk).is_err() {
                            break; // Receiver hung up, stop reading
                        }
                        if is_empty {
                            break; // End of file
                        }
                    }
                });

                // --- 2. Worker Pool (using Rayon) ---
                // This consumes chunks from the channel (rx)
                // `par_bridge` turns the channel iterator into a parallel iterator
                // `for_each` processes each chunk in parallel using the Rayon thread pool
                rx.into_iter().par_bridge().for_each(|chunk: Vec<Record>| {
                    analyze_chunk(
                        &reference_seq_clone,
                        chunk,
                        &protected_sites_clone,
                        args.match_len,
                        &global_counts_clone,
                        &total_reads_clone,
                        &total_valid_clone,
                    );
                });

                Ok(())
            }); // --- End of thread::scope ---

            if let Err(e) = res {
                eprintln!("Error during processing {}: {:?}", fasta_file.display(), e);
                continue; // Skip to next file
            }

            // --- 3. Collate and Write Results for this file ---
            let total_r = total_reads.load(Ordering::Relaxed);
            let total_v = total_valid.load(Ordering::Relaxed);
            println!("{} - Valid reads: {} / {}", file_stem, total_v, total_r);

            let mut mutation_stats = Vec::new();
            for (i, counter_map) in global_counts.iter().enumerate() {
                let ref_aa = reference_seq[i]; // Get the reference AA at this position
                let adj_pos = (i as i32) + 1 + args.aa_offset; // Calculate the adjusted position
                
                for item in counter_map.iter() {
                    let aa = *item.key();
                    let count = item.value().load(Ordering::Relaxed);
                    if count > 0 {
                        // Format: e.g., "A123C"
                        let mutation_str = format!("{}{}{}", ref_aa as char, adj_pos, aa as char);
                        mutation_stats.push((mutation_str, count));
                    }
                }
            }
            
            // Sort by mutation string (e.g., "A10C" before "A11G")
            mutation_stats.sort_by(|a, b| a.0.cmp(&b.0));

            // Write to CSV
            let output_file_name = format!("{}_mutation.csv", file_stem);
            let output_path = args.output_dir.join(output_file_name);
            
            let mut wtr = csv::Writer::from_path(&output_path)
                .with_context(|| format!("Failed to create output CSV: {:?}", output_path))?;
            
            wtr.write_record(&["Mutation", "Count"])?;
            for (mutation, count) in mutation_stats {
                wtr.write_record(&[mutation, count.to_string()])?;
            }
            
            wtr.flush()?;
            println!("Results saved to: {}", output_path.display());
            println!("Time taken for {}: {:.2?}", file_stem, file_start_time.elapsed());
        }

        println!("\n🎉 All files have been processed. Total time: {:.2?}", main_start_time.elapsed());
        Ok(())
    }
}mod find_seq {
    use super::common::{detect_format, Format};
    use anyhow::Result;
    use bio::io::{fasta, fastq};
    use clap::Parser;
    use csv::Writer;
    use flate2::bufread::MultiGzDecoder;
    use std::collections::{HashMap, HashSet};
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::PathBuf;

    #[derive(Parser, Debug)]
    #[command(
        name = "find_seq",
        about = "Find a motif in FASTA/FASTQ (gz supported), consider reverse complement, extract upstream/downstream flanks, and count unique windows per read. Outputs CSV with Sequence, UpFlank, DownFlank, ReadsCount."
    )]
    pub struct Args {
        #[arg(long, help = "Input FASTA/FASTQ file (optionally .gz)")]
        pub inputfile: PathBuf,
        #[arg(long, help = "Output CSV file path")]
        pub output: PathBuf,
        #[arg(long, help = "Target motif sequence")]
        pub motif: String,
        #[arg(long, help = "Upstream flank length", default_value_t = 0)]
        pub up_flank: usize,
        #[arg(long, help = "Downstream flank length", default_value_t = 0)]
        pub down_flank: usize,
    }

    fn revcomp(s: &str) -> String {
        let mut out = String::with_capacity(s.len());
        for &b in s.as_bytes().iter().rev() {
            let c = match b {
                b'A' => 'T',
                b'T' => 'A',
                b'C' => 'G',
                b'G' => 'C',
                b'N' => 'N',
                _ => 'N',
            };
            out.push(c);
        }
        out
    }

    fn find_all(hay: &str, needle: &str) -> Vec<usize> {
        let mut res = Vec::new();
        let mut start = 0usize;
        while let Some(pos) = hay[start..].find(needle) {
            res.push(start + pos);
            start = start + pos + 1;
        }
        res
    }

    pub fn run(mut args: Args) -> Result<()> {
        let default_flank = 40usize;
        if args.up_flank == 0 && args.down_flank == 0 { args.up_flank = default_flank; args.down_flank = default_flank; }
        else if args.up_flank == 0 { args.up_flank = args.down_flank; }
        else if args.down_flank == 0 { args.down_flank = args.up_flank; }
        let up = args.up_flank; let down = args.down_flank;
        let motif = args.motif.to_uppercase();
        let motif_rc = revcomp(&motif);

        let format = detect_format(&args.inputfile)?;
        let file = File::open(&args.inputfile)?;
        let buf_reader = BufReader::new(file);
        let input_reader: Box<dyn BufRead> = if args.inputfile.extension().map_or(false, |ext| ext == "gz") { Box::new(BufReader::new(MultiGzDecoder::new(buf_reader))) } else { Box::new(buf_reader) };

        let mut counts: HashMap<String, usize> = HashMap::new();
        match format {
            Format::Fasta => {
                let reader = fasta::Reader::new(input_reader);
                for result in reader.records() {
                    let record = result?;
                    let seq = String::from_utf8(record.seq().to_vec()).unwrap().to_uppercase();
                    process_seq(&seq, &motif, &motif_rc, up, down, &mut counts);
                }
            }
            Format::Fastq => {
                let reader = fastq::Reader::new(input_reader);
                for result in reader.records() {
                    let record = result?;
                    let seq = String::from_utf8(record.seq().to_vec()).unwrap().to_uppercase();
                    process_seq(&seq, &motif, &motif_rc, up, down, &mut counts);
                }
            }
        }

        let mut wtr = Writer::from_path(&args.output)?;
        wtr.write_record(["Sequence", "UpFlank", "DownFlank", "ReadsCount"])?;
        for (seq, c) in counts.into_iter() {
            let up_seq = if up > 0 { seq[..up].to_string() } else { String::new() };
            let down_seq = if down > 0 { seq[seq.len() - down..].to_string() } else { String::new() };
            wtr.write_record([seq, up_seq, down_seq, c.to_string()])?;
        }
        wtr.flush()?;
        Ok(())
    }

    fn process_seq(seq: &str, motif: &str, motif_rc: &str, up: usize, down: usize, counts: &mut HashMap<String, usize>) {
        let mut per_read: HashSet<String> = HashSet::new();
        for idx in find_all(seq, motif) {
            let left = idx as isize - up as isize;
            let right = idx + motif.len() + down;
            if left < 0 || right > seq.len() { continue; }
            let w = &seq[left as usize..right];
            per_read.insert(w.to_string());
        }
        for idx in find_all(seq, motif_rc) {
            let left = idx as isize - down as isize;
            let right = idx + motif.len() + up;
            if left < 0 || right > seq.len() { continue; }
            let w = &seq[left as usize..right];
            let w_rc = revcomp(w);
            per_read.insert(w_rc);
        }
        for w in per_read { *counts.entry(w).or_insert(0) += 1; }
    }
}
