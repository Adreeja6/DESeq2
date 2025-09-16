# RNA-Seq Analysis Pipeline

## Overview

This repository provides a complete RNA-Seq analysis workflow, combining a Bash pipeline for processing raw FASTQ files and a Shiny web application for differential expression analysis using DESeq2. The pipeline handles quality control, trimming, alignment, and gene quantification, producing a gene count matrix from paired-end RNA-Seq data. The Shiny app offers an interactive interface to perform DESeq2 analysis, generate visualizations (PCA, MA, volcano plots, and heatmaps), and download results. The workflow is designed for analyzing RNA-Seq data, such as the Myelodysplastic Syndromes (MDS) dataset provided in `MDS.txt`, which includes FASTQ file URLs from the European Nucleotide Archive (ENA).

## Repository Structure

- `pipeline.sh`: Bash script for the RNA-Seq pipeline, handling preprocessing, alignment, and gene quantification.
- `app.R`: R Shiny application for differential expression analysis and visualization using DESeq2.
- `MDS.txt`: Metadata file with URLs for downloading FASTQ files.
- `README.md`: This file, providing setup and usage instructions.

## Prerequisites

- **Software**:
  - `fastp` (v1.0.1)
  - `FastQC` (v0.11.9+)
  - `MultiQC` (v1.31)
  - `HISAT2` (v2.2.1+)
  - `samtools` (v1.13+)
  - `featureCounts` (Subread, v2.0.3+)
  - `Qualimap` (v2.2.2+)
  - `wget`
  - R (v4.0+) with packages: `shiny`, `DESeq2`, `SummarizedExperiment`, `biomaRt`, `ggplot2`, `RColorBrewer`, `pheatmap`, `circlize`, `ggrepel`, `DT`, `apeglm` (optional)
- **Reference Files**:
  - Human reference genome (e.g., hg38)
  - GTF annotation file (e.g., Gencode v44)
- **Hardware**:
  - Multi-core CPU (10+ threads recommended)
  - 16GB+ RAM for Qualimap
  - 100GB+ storage for FASTQ and BAM files
- **Internet**: Required for downloading reference files, FASTQ files, and R packages, and for `biomaRt` annotation (optional).

## Setup and Usage

Follow these steps to set up the pipeline, download FASTQ files, process them, and analyze results with the Shiny app. Run each command in a Linux terminal with sufficient permissions.

### 1. Clone the Repository
Clone the repository and navigate to the directory:
```bash
git clone https://github.com/Adreeja6/DESeq2.git
```
```bash
cd DESeq2
```

### 2. Create Required Directories
Set up the directory structure for the pipeline:
```bash
mkdir -p ref gtf reads aligned_reads bam_files sam_files fastqc_data trimmed_data results multiqc_output
```

### 3. Download Reference Genome and GTF File
Download the hg38 reference genome and Gencode v44 GTF file (update URLs if newer versions are available):
```bash
wget -P ref https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
```
```bash
gunzip ref/GRCh38.primary_assembly.genome.fa.gz
```
```bash
wget -P gtf https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
```
```bash
gunzip gtf/gencode.v44.annotation.gtf.gz
```

### 4. Download FASTQ Files
Download paired-end FASTQ files from `MDS.txt` into the `reads/` directory:
```bash
mkdir -p reads
```
```bash
tail -n +2 MDS.txt | cut -f8 | tr ';' '\n' | sed '/^$/d' | xargs -n1 -I {} wget -c -P reads {}
```

### 5. Update Pipeline Script Paths
Update `pipeline.sh` to use the correct paths for the reference genome, GTF, and directories:
```bash
sed -i 's|ref="/home/pls_trainee/adreeja/ref/hg38"|ref="$(pwd)/ref/GRCh38.primary_assembly.genome.fa"|' pipeline.sh
```
```bash
sed -i 's|gtf="/home/pls_trainee/adreeja/gtf/gencode.v44.annotation.gtf"|gtf="$(pwd)/gtf/gencode.v44.annotation.gtf"|' pipeline.sh
```
```bash
sed -i 's|/home/pls_trainee/adreeja/RNA/|$(pwd)/|' pipeline.sh
```

### 6. Make Your Script Executable
Set execute permissions for `script.sh`:
```bash
chmod +x script.sh
```

### 7. Run the RNA-Seq Pipeline
Execute the pipeline to process FASTQ files, align reads, and generate gene counts:
```bash
./script.sh
```

### 8. Install R and Required Packages
Install R and the necessary R packages for the Shiny app:
```bash
sudo apt-get update
```
```bash
sudo apt-get install -y r-base
```
```bash
R -e 'install.packages(c("shiny", "DESeq2", "SummarizedExperiment", "biomaRt", "ggplot2", "RColorBrewer", "pheatmap", "circlize", "ggrepel", "DT"))'
```
```bash
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("apeglm")'
```

### 9. Create Sample Metadata File
Generate `coldata.csv` for the Shiny app based on `MDS.txt` experiment titles:
```R
# Read MDS.txt
mds <- read.table("MDS.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Create data frame for coldata
coldata <- data.frame(
  Run = mds$run_accession,
  cell_type = ifelse(grepl("MDS-mesenchymal stem cell- conditioned monocytes", mds$experiment_title), 
                     "Monocytes", "Mesenchymal stem cells"),
  disease_state = ifelse(grepl("MDS", mds$experiment_title), "MDS", "healthy donor")
)

# Set row names and remove Run column from data
rownames(coldata) <- coldata$Run
coldata$Run <- NULL

# Write to coldata.csv
write.csv(coldata, "coldata.csv", row.names = TRUE)
```

### 10. Run the Shiny App
Start the Shiny app to analyze the gene counts:
```bash
Rscript app.R
```

### 11. Access the Shiny App
Open the URL displayed in the terminal (e.g., `http://127.0.0.1:XXXX`) in a web browser, then:
- Upload `gene_counts.txt` as the count matrix.
- Upload `coldata.csv` as the sample metadata.
- Click "Run DESeq2 Analysis" to view results.

## Pipeline Details

The `pipeline.sh` script automates:
1. **Logging**: Creates a timestamped log file (`pipeline_log_YYYYMMDD_HHMMSS.log`).
2. **FASTQ Trimming**: Uses `fastp` for adapter trimming and quality filtering.
3. **Quality Control**: Runs `FastQC` and aggregates results with `MultiQC`.
4. **Alignment and Sorting**: Aligns reads with `HISAT2`, converts SAM to BAM, sorts, and indexes with `samtools`.
5. **Gene Quantification**: Counts genes with `featureCounts` using the GTF file.
6. **Quality Assessment**: Generates RNA-Seq metrics with `Qualimap`.

### Outputs
- **Trimmed FASTQ**: `trimmed_data/*.trimmed.fastq.gz`
- **FastQC Reports**: `fastqc_data/*.html`
- **MultiQC Report**: `multiqc_output/multiqc_report.html`
- **BAM Files**: `bam_files/*.sorted.bam`
- **Gene Counts**: `bam_files/gene_counts.txt`
- **Qualimap Reports**: `bam_files/qualimap_<sample>/`

## Shiny App Details

The `app.R` script provides an interactive interface for DESeq2 analysis, offering:
- **Input**: Upload a count matrix (e.g., `gene_counts.txt`) and sample metadata (e.g., `coldata.csv`).
- **Analysis**: Performs DESeq2 differential expression analysis, variance stabilization, and gene annotation via `biomaRt`.
- **Visualizations**: PCA, MA, volcano plots, sample distance heatmap, and top 20 gene heatmap.
- **Downloads**: CSV files for all/significant genes and PNG files for plots.

### Input Files
- **Count Matrix** (`*.csv` or `*.txt`):
  ```csv
  gene_id,sample1,sample2,...
  ENSG00000000003,100,120,...
  ```
- **Sample Metadata** (`*.csv` or `*.txt`):
  ```csv
  sample,disease_state,cell_type
  SRR10416997,MDS,Mesenchymal stem cells
  SRR10417001,healthy donor,Mesenchymal stem cells
  ```

## Notes
- **Reference Files**: Ensure the hg38 genome is indexed for HISAT2:
  ```bash
  hisat2-build ref/GRCh38.primary_assembly.genome.fa ref/hg38
  ```
  Update `pipeline.sh` to point to `ref/hg38` if indexed.
- **FASTQ Downloads**: Verify files in `reads/` with `ls reads/*.fastq.gz`. Use `--limit-rate=500k` with `wget` for large datasets if needed.
- **Pipeline Paths**: Check `script.sh` paths manually with `nano pipeline.sh` if `sed` commands fail.
- **Metadata**: Verify `coldata.csv` sample names match `gene_counts.txt` columns:
  ```bash
  cat coldata.csv
  ```
- **Pipeline Errors**: Check the log file for issues:
  ```bash
  cat pipeline_log_*.log
  ```
- **Shiny App**: Increase `shiny.maxRequestSize` in `app.R` for large datasets. Internet is required for `biomaRt` (falls back if unavailable).
- **Dependencies**: Install tools with:
  ```bash
  sudo apt-get install -y fastp fastqc multiqc hisat2 samtools subread qualimap wget
  ```
- **Storage**: Ensure 100GB+ free space for FASTQ and BAM files.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Submit a pull request or open an issue for bugs, feature requests, or improvements.

## Contact

For questions or support, open an issue on this repository.
