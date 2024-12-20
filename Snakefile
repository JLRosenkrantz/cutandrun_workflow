# Snakemake workflow for Cut&Run data analysis

# Define the location of raw FASTQ files
# FASTQ_DIR = "/home/groups/hoolock2/u0/jimi/raw_data/241219_cutandrun_cora/test"
FASTQ_DIR = "/home/groups/hoolock2/u0/jimi/raw_data/241219_cutandrun_cora/LIB241018LC"

# Define the list of samples
SAMPLES, = glob_wildcards(f"{FASTQ_DIR}/{{sample}}_R1_001.fastq.gz")

# Define the reference genome (for Bowtie2 alignment)
GENOME_INDEX = "/home/groups/hoolock2/u0/genomes/ucsc/nomLeu3/indexes/bowtie2/nomLeu3"

# Raw FastQ files (input)
RAW_FASTQ_R1 = expand(f"{FASTQ_DIR}/{{sample}}_R1_001.fastq.gz", sample=SAMPLES)
RAW_FASTQ_R2 = expand(f"{FASTQ_DIR}/{{sample}}_R2_001.fastq.gz", sample=SAMPLES)

# Trimmed FastQ files (output of fastp)
TRIMMED_FASTQ_R1 = expand("trimmed/{sample}_R1_trimmed.fastq.gz", sample=SAMPLES)
TRIMMED_FASTQ_R2 = expand("trimmed/{sample}_R2_trimmed.fastq.gz", sample=SAMPLES)

# Define FastQC reports
FASTQC_RAW = expand("fastqc_raw/{sample}_R1_001_fastqc.html", sample=SAMPLES) + expand("fastqc_raw/{sample}_R2_001_fastqc.html", sample=SAMPLES)
FASTQC_TRIMMED = expand("fastqc_trimmed/{sample}_R1_trimmed_fastqc.html", sample=SAMPLES) + expand("fastqc_trimmed/{sample}_R2_trimmed_fastqc.html", sample=SAMPLES)

# Define alignment outputs
ALIGNED_BAM = expand("aligned/{sample}.bam", sample=SAMPLES)
FILTERED_BAM = expand("filtered/{sample}.bam", sample=SAMPLES)

# Define the final rule
rule all:
    input:
        FASTQC_RAW + FASTQC_TRIMMED + ALIGNED_BAM + FILTERED_BAM + [
            "multiqc/raw_fastqc_report.html",
            "multiqc/trimmed_fastqc_report.html"
        ]

# FastQC analysis on raw fastq files
rule fastqc_raw:
    input:
        raw_fastq_r1=f"{FASTQ_DIR}/{{sample}}_R1_001.fastq.gz",
        raw_fastq_r2=f"{FASTQ_DIR}/{{sample}}_R2_001.fastq.gz"
    output:
        html_r1="fastqc_raw/{sample}_R1_001_fastqc.html",
        zip_r1="fastqc_raw/{sample}_R1_001_fastqc.zip",
        html_r2="fastqc_raw/{sample}_R2_001_fastqc.html",
        zip_r2="fastqc_raw/{sample}_R2_001_fastqc.zip"
    threads: 4
    log:
        "logs/fastqc_raw_{sample}.log"
    shell:
        "fastqc -t {threads} {input.raw_fastq_r1} {input.raw_fastq_r2} -o fastqc_raw/ &> {log}"
        

# MultiQC report on raw FastQC
rule multiqc_raw:
    input:
        expand("fastqc_raw/{sample}_R1_001_fastqc.zip", sample=SAMPLES) + expand("fastqc_raw/{sample}_R2_001_fastqc.zip", sample=SAMPLES)
    output:
        "multiqc/raw_fastqc_report.html"
    threads: 4
    log:
        "logs/multiqc_raw.log"
    shell:
        "multiqc fastqc_raw/ -o multiqc/ -n raw_fastqc_report &> {log}"

# Trimming with fastp
rule trim:
    input:
        priority="multiqc/raw_fastqc_report.html",
        raw_fastq_r1=f"{FASTQ_DIR}/{{sample}}_R1_001.fastq.gz",
        raw_fastq_r2=f"{FASTQ_DIR}/{{sample}}_R2_001.fastq.gz"
    output:
        trimmed_r1="trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_r2="trimmed/{sample}_R2_trimmed.fastq.gz"
    params:
        fastp_output="trimmed/{sample}_fastp.json",
        fastp_html="trimmed/{sample}_fastp.html"
    threads: 4
    log:
        "logs/trim_{sample}.log"
    shell:
        "fastp -i {input.raw_fastq_r1} -I {input.raw_fastq_r2} -o {output.trimmed_r1} -O {output.trimmed_r2} \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 15 \
        --length_required 36 \
        --json {params.fastp_output} \
        --html {params.fastp_html} &> {log}"

# FastQC analysis on trimmed fastq files
rule fastqc_trimmed:
    input:
        priority1=TRIMMED_FASTQ_R1,
        priority2=TRIMMED_FASTQ_R2,
        trimmed_r1="trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_r2="trimmed/{sample}_R2_trimmed.fastq.gz"
    output:
        html_r1="fastqc_trimmed/{sample}_R1_trimmed_fastqc.html",
        zip_r1="fastqc_trimmed/{sample}_R1_trimmed_fastqc.zip",
        html_r2="fastqc_trimmed/{sample}_R2_trimmed_fastqc.html",
        zip_r2="fastqc_trimmed/{sample}_R2_trimmed_fastqc.zip"
    threads: 4
    log:
        "logs/fastqc_trimmed_{sample}.log"
    shell:
        "fastqc -t {threads} {input.trimmed_r1} {input.trimmed_r2} -o fastqc_trimmed/ &> {log}"

# MultiQC report on trimmed FastQC
rule multiqc_trimmed:
    input:
        expand("fastqc_trimmed/{sample}_R1_trimmed_fastqc.zip", sample=SAMPLES) + expand("fastqc_trimmed/{sample}_R2_trimmed_fastqc.zip", sample=SAMPLES)
    output:
        "multiqc/trimmed_fastqc_report.html"
    threads: 4
    log:
        "logs/multiqc_trimmed.log"
    shell:
        "multiqc fastqc_trimmed/ -o multiqc/ -n trimmed_fastqc_report &> {log}"

# Bowtie2 alignment to gibbon genome
rule align_bowtie2:
    input:
        priority="multiqc/trimmed_fastqc_report.html",
        trimmed_r1="trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_r2="trimmed/{sample}_R2_trimmed.fastq.gz"
    output:
        sam="aligned/{sample}.sam",
        bam="aligned/{sample}.bam"
    params:
        bowtie2_index=GENOME_INDEX,
        min_insert=10,
        max_insert=700,
        threads=4
    threads: 4
    log:
        "logs/align_bowtie2_{sample}.log"
    shell:
        "bowtie2 -x {params.bowtie2_index} -1 {input.trimmed_r1} -2 {input.trimmed_r2} \
        --very-sensitive \
        -I {params.min_insert} -X {params.max_insert} \
        -p {threads} \
        -S {output.sam} &> {log} && \
        samtools view -bS {output.sam} | samtools sort -o {output.bam} && \
        samtools index {output.bam}"

# Filter low quality reads
rule filter_bam:
    input:
        bam="aligned/{sample}.bam"
    output:
        filtered_bam="filtered/{sample}.bam"
    params:
        quality=30
    threads: 4
    log:
        "logs/filter_bam_{sample}.log"
    shell:
        "samtools view -b -q {params.quality} {input.bam} > {output.filtered_bam} &> {log}"
