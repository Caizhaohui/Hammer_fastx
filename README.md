[# Hammer_fastx

A versatile toolkit for FASTX file processing, including quality control, merging, demultiplexing, and sequence analysis.

## Overview

`Hammer_fastx` is a Rust-based command-line tool designed for processing FASTA and FASTQ files. It provides a comprehensive set of subcommands for quality control, paired-end read merging, demultiplexing, sequence statistics, filtering, and alignment to references with ambiguous regions (N-regions). The tool leverages external dependencies like `fastp` and `flash2` for specific tasks and includes parallelized workflows for efficient processing.

- **Version**: v1.2.0
- **Author**: CZH with the help of Gemini 2.5 pro
- **License**: MIT (or specify your preferred license)

## Features

- **Workflows**:
  - `demux_all`: Complete pipeline for quality control, merging, and demultiplexing.
  - `mergePE`: Quality control and merging of paired-end reads with customizable output formats.
- **Single-Step Commands**:
  - `demux_only`: Demultiplex merged FASTQ files based on barcodes.
  - `fastp`: Quality control for paired-end FASTQ files using `fastp`.
  - `flash2`: Merge paired-end reads using `flash2`.
  - `stats`: Generate sequence statistics for FASTA/FASTQ files.
  - `filter`: Filter FASTA/FASTQ files based on sequence length.
  - `Ns_count`: Align reads to a reference with N-regions using anchor-based logic.

## Installation

### Prerequisites

- **Rust**: Ensure you have Rust installed (`cargo` is used for building).
- **External Tools**:
  - `fastp`: Required for quality control (`fastp` subcommand and workflows).
  - `flash2`: Required for paired-end read merging (`flash2` subcommand and workflows).
- Ensure both `fastp` and `flash2` are in your system's PATH.

### Build from Source

1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/hammer_fastx.git
   cd hammer_fastx
   ```

2. Build the project:
   ```bash
   cargo build --release
   ```

3. Install the binary:
   ```bash
   cargo install --path .
   ```

The executable will be available as `Hammer_fastx`.

## Usage

Run `Hammer_fastx --help` to see all available subcommands and options.

```bash
Hammer_fastx v1.2.0
CZH with the help of Gemini
A versatile toolkit for FASTX file processing, including QC, merging, and demultiplexing.

USAGE:
    Hammer_fastx [SUBCOMMAND]

SUBCOMMANDS:
    demux_all    Run the complete pipeline from QC and merging to demultiplexing
    mergePE      Quality control and merge paired-end data, with optional output formats
    demux_only   Demultiplex a merged FASTQ file based on barcodes
    fastp        Quality control paired-end FASTQ files using fastp
    flash2       Merge paired-end reads using flash2
    stats        Get sequence statistics from one or more FASTA/FASTQ files
    filter       Filter one or more FASTA/FASTQ files based on sequence length
    Ns_count     Align reads to a reference with Ns using a strict anchor-based method
```

### Example Commands

1. **Run the full demultiplexing pipeline**:
   ```bash
   Hammer_fastx demux_all -i read1.fastq.gz -I read2.fastq.gz --tags samples.csv -o output_dir --fastp-threads 8 --flash-threads 8 --demux-threads 8
   ```

2. **Merge paired-end reads**:
   ```bash
   Hammer_fastx mergePE -i read1.fastq.gz -I read2.fastq.gz -o merged.fastq --out-fasta --fastp-threads 4
   ```

3. **Demultiplex a merged FASTQ file**:
   ```bash
   Hammer_fastx demux_only --inputfile merged.fastq --output demux_out --tags samples.csv --threads 16 --trim
   ```

4. **Generate sequence statistics**:
   ```bash
   Hammer_fastx stats --inputfile sample1.fasta sample2.fastq
   ```

5. **Filter sequences by length**:
   ```bash
   Hammer_fastx filter --inputfile input.fasta --outfile filtered.fasta --min-len 100 --max-len 1000
   ```

6. **Align reads to a reference with N-regions**:
   ```bash
   Hammer_fastx Ns_count --reads reads.fasta --refSEQ reference.fasta --output results --threads 8 --mismatches 2
   ```

## Input/Output Formats

- **Input**: Supports FASTA and FASTQ files, including gzipped files.
- **Output**:
  - Workflows and demultiplexing: FASTQ (default) or FASTA (optional).
  - Statistics: Tabular summary printed to the console.
  - Filtering: FASTA/FASTQ output to file or stdout.
  - Ns_count: CSV files with combination counts and optional FASTA files for matching reads.

### Tag File Format (for `demux_only` and `demux_all`)

The tag file for demultiplexing must be a CSV with the following columns:
- `SampleID`: Unique identifier for the sample.
- `F_tag`: Forward tag sequence.
- `R_tag`: Reverse tag sequence.

Example (`samples.csv`):
```csv
SampleID,F_tag,R_tag
sample1,ACGTACGT,TGCACTGC
sample2,GGTTCCAA,CCATGGTT
```

## Dependencies

- **Rust Crates**:
  - `anyhow`: Error handling.
  - `clap`: Command-line argument parsing.
  - `bio`: FASTA/FASTQ parsing and DNA sequence utilities.
  - `flate2`: Gzip file handling.
  - `indicatif`: Progress bars.
  - `rayon`: Parallel processing.
  - `crossbeam-channel`: Thread-safe communication.
  - `csv`: CSV file handling.
- **External Tools**:
  - `fastp`: For quality control.
  - `flash2`: For paired-end read merging.

## Contributing

Contributions are welcome! Please submit issues or pull requests to the GitHub repository. Ensure your code follows the project's coding style and includes appropriate tests.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.
](https://zread.ai/Caizhaohui/Hammer_fastx)
