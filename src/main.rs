use anyhow::Result;
use clap::Parser;
use std::process::{Command, Stdio}; // For executing external commands

// ==================================================================================
// Main entry point and CLI command definitions
// ==================================================================================

#[derive(Parser, Debug)]
#[command(
    name = "Hammer_fastx",
    version = "v1.2.0", // Hybrid: Modern codebase with classic anchor-based Ns_count logic
    author = "CZH with the help of Gemini",
    about = "A versatile toolkit for FASTX file processing, including QC, merging, and demultiplexing."
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Parser, Debug)]
enum Commands {
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

    /// [Anchor Logic] Align reads to a reference with Ns using a strict anchor-based method
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
            println!("    - Cleaned R1: {}", args.out1.display());
            println!("    - Cleaned R2: {}", args.out2.display());
            if let Some(html_path) = &args.html {
                println!("    - HTML Report: {}", html_path.display());
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
        if cfg!(unix) {
            Command::new("which")
                .arg(cmd)
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .status()
                .map_or(false, |s| s.success())
        } else {
            Command::new("where")
                .arg(cmd)
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .status()
                .map_or(false, |s| s.success())
        }
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
            println!("    - Output directory: {}", args.out_dir.display());
            println!("    - Output prefix: {}", args.out_prefix);
            println!("    - Merged file: {}", args.out_dir.join(format!("{}.extendedFrags.fastq", args.out_prefix)).display());
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
        pb.finish_with_message("✔ File reading complete");
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
            println!("  - Matched Reads:   {:>10} ({:.2}%)", matched_reads, matched_percent);
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
    pub fn run(args: Args) -> Result<()> {
        let start_time = Instant::now();
        let output_dir = args.output.clone();
        std::fs::create_dir_all(&output_dir)
            .with_context(|| format!("Failed to create output directory: {:?}", output_dir))?;
        println!("---> Loading tags...");
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
        pb.set_message("Processing...");
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
                Ok(Err(e)) => eprintln!("Writer thread error: {:?}", e),
                Err(e) => eprintln!("Writer thread panicked: {:?}", e),
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
    use flate2::bufread::MultiGzDecoder;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::{Path, PathBuf};

    #[derive(Parser, Debug)]
    pub struct Args {
        #[arg(long, help = "One or more input files (wildcards supported, e.g., '*.fasta')", required = true, num_args = 1..)]
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
            println!("No files processed or no sequences found.");
            return;
        }

        println!("\n====================================== Sequence Statistics Summary ======================================");
        println!("{:<25} {:>15} {:>18} {:>10} {:>10} {:>12}",
                 "Sample Name", "Total Seqs", "Total Bases", "Min Length", "Max Length", "Avg Length");
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
        println!("===================================================================================================");
    }

    pub fn run(args: Args) -> Result<()> {
        let mut all_stats: Vec<FileStats> = Vec::new();

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
// `filter` subcommand module
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
        #[arg(long, help = "One or more input files (wildcards supported, e.g., '*.fasta')", required = true, num_args = 1..)]
        inputfile: Vec<PathBuf>,

        #[arg(long, help = "Output file (default: standard output)")]
        outfile: Option<PathBuf>,

        #[arg(short = 'm', long, help = "Filter out sequences shorter than this length")]
        min_len: Option<usize>,
        
        #[arg(short = 'M', long, help = "Filter out sequences longer than this length")]
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
                sorted_combos.sort_by(|a, b| b.1.cmp(a.1));

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
