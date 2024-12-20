# Snakemake Workflow for Cut&Run Data Analysis

This README provides instructions on how to run the Snakemake workflow for Cut&Run data analysis. The workflow includes steps for quality control, trimming, alignment, and filtering of paired-end sequencing data.

## Requirements

1. **Software**:
   - [Snakemake](https://snakemake.readthedocs.io/)
   - [Conda](https://docs.conda.io/en/latest/)
   - Bioinformatics tools (installed via Conda):
     - `fastqc`
     - `multiqc`
     - `fastp`
     - `bowtie2`
     - `samtools`

2. **Reference Genome**:
   - Indexed Bowtie2 genome for gibbon: `/home/groups/hoolock2/u0/genomes/ucsc/nomLeu3/indexes/bowtie2/nomLeu3`

3. **Input Data**:
   - Paired-end FASTQ files located in `/home/groups/hoolock2/u0/jimi/raw_data/241219_cutandrun_cora/LIB241018LC`

## Setting Up the Conda Environment

1. **Install Conda** (if not already installed):
   Download and install Miniconda or Anaconda from [https://docs.conda.io/en/latest/](https://docs.conda.io/en/latest/).

2. **Create the Environment**:
   Use the provided `cutrun_env.yml` file to create the environment:
   ```bash
   conda env create -f cutrun_env.yml
   ```

3. **Activate the Environment**:
   Activate the new environment:
   ```bash
   conda activate cutrun_env
   ```

4. **Verify Installation**:
   Ensure all tools are properly installed:
   ```bash
   snakemake --version
   fastqc --version
   multiqc --version
   fastp --version
   bowtie2 --version
   samtools --version
   ```

## Running the Workflow

1. **Dry Run**:
   Test the workflow without executing any commands to ensure the rules and dependencies are correctly defined:
   ```bash
   snakemake -n
   ```

2. **Execute the Workflow**:
   Run the workflow:
   ```bash
   snakemake --use-conda
   ```

3. **Generate Reports**:
   View the MultiQC reports:
   - Raw FastQC report: `multiqc/raw_fastqc_report.html`
   - Trimmed FastQC report: `multiqc/trimmed_fastqc_report.html`


## Workflow Steps

1. **FastQC on Raw Reads**:
   Performs quality control on the raw FASTQ files.

2. **MultiQC on Raw Reads**:
   Aggregates the FastQC results into a summary report.

3. **Trimming**:
   Trims adapters and low-quality regions from the raw reads using `fastp`.

4. **FastQC on Trimmed Reads**:
   Quality control on the trimmed FASTQ files.

5. **MultiQC on Trimmed Reads**:
   Summarizes the FastQC results of trimmed reads.

6. **Bowtie2 Alignment**:
   Aligns the trimmed reads to the gibbon genome and generates BAM files.

7. **Filtering BAM Files**:
   Filters low-quality alignments (MAPQ < 30) to produce high-quality BAM files.

## Outputs

- **FastQC Reports**:
  - Raw reads: `fastqc_raw/`
  - Trimmed reads: `fastqc_trimmed/`

- **MultiQC Reports**:
  - Raw reads summary: `multiqc/raw_fastqc_report.html`
  - Trimmed reads summary: `multiqc/trimmed_fastqc_report.html`

- **Alignment Outputs**:
  - Sorted BAM files: `aligned/{sample}.bam`
  - Filtered BAM files: `filtered/{sample}.bam`

## Troubleshooting

1. **Missing Files**:
   If Snakemake fails due to missing files, increase the latency wait time:
   ```bash
   snakemake -j 4 --latency-wait 60
   ```

2. **Environment Issues**:
   If dependencies are missing, recreate the environment using:
   ```bash
   conda env create -f environment.yml
   ```

3. **Debugging**:
   Run Snakemake in debug mode for detailed error messages:
   ```bash
   snakemake -j 4 --debug
   ```

