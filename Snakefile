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
TRIMMED_FASTQ_R1 = expand("02_trimmed/{sample}_R1_trimmed.fastq.gz", sample=SAMPLES)
TRIMMED_FASTQ_R2 = expand("02_trimmed/{sample}_R2_trimmed.fastq.gz", sample=SAMPLES)

# Define FastQC reports
FASTQC_RAW = expand("01_fastqc_raw/{sample}_R1_001_fastqc.html", sample=SAMPLES) + expand("01_fastqc_raw/{sample}_R2_001_fastqc.html", sample=SAMPLES)
FASTQC_TRIMMED = expand("03_fastqc_trimmed/{sample}_R1_trimmed_fastqc.html", sample=SAMPLES) + expand("03_fastqc_trimmed/{sample}_R2_trimmed_fastqc.html", sample=SAMPLES)

# Define alignment outputs
ALIGNED_BAM = expand("04_aligned/{sample}.bam", sample=SAMPLES)
FILTERED_BAM = expand("05_filtered/{sample}.bam", sample=SAMPLES)
BED_FILES = expand("06_bed/{sample}_filtered.bed", sample=SAMPLES)

# Define the final rule
rule all:
    input:
        FASTQC_RAW + FASTQC_TRIMMED + ALIGNED_BAM + FILTERED_BAM + BED_FILES + [
            "01_multiqc/raw_fastqc_report.html",
            "03_multiqc/trimmed_fastqc_report.html"
        ]

# FastQC analysis on raw fastq files
rule fastqc_raw:
    input:
        raw_fastq_r1=f"{FASTQ_DIR}/{{sample}}_R1_001.fastq.gz",
        raw_fastq_r2=f"{FASTQ_DIR}/{{sample}}_R2_001.fastq.gz"
    output:
        html_r1="01_fastqc_raw/{sample}_R1_001_fastqc.html",
        zip_r1="01_fastqc_raw/{sample}_R1_001_fastqc.zip",
        html_r2="01_fastqc_raw/{sample}_R2_001_fastqc.html",
        zip_r2="01_fastqc_raw/{sample}_R2_001_fastqc.zip"
    threads: 4
    log:
        "logs/fastqc_raw_{sample}.log"
    shell:
        "fastqc -t {threads} {input.raw_fastq_r1} {input.raw_fastq_r2} -o 01_fastqc_raw/ &> {log}"

# MultiQC report on raw FastQC
rule multiqc_raw:
    input:
        expand("01_fastqc_raw/{sample}_R1_001_fastqc.zip", sample=SAMPLES) + expand("01_fastqc_raw/{sample}_R2_001_fastqc.zip", sample=SAMPLES)
    output:
        "01_multiqc/raw_fastqc_report.html"
    threads: 4
    log:
        "logs/multiqc_raw.log"
    shell:
        "multiqc 01_fastqc_raw/ -o 01_multiqc/ -n raw_fastqc_report &> {log}"

# Trimming with fastp
rule trim:
    input:
        raw_fastq_r1=f"{FASTQ_DIR}/{{sample}}_R1_001.fastq.gz",
        raw_fastq_r2=f"{FASTQ_DIR}/{{sample}}_R2_001.fastq.gz"
    output:
        trimmed_r1="02_trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_r2="02_trimmed/{sample}_R2_trimmed.fastq.gz"
    params:
        fastp_output="02_trimmed/{sample}_fastp.json",
        fastp_html="02_trimmed/{sample}_fastp.html"
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
        trimmed_r1="02_trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_r2="02_trimmed/{sample}_R2_trimmed.fastq.gz"
    output:
        html_r1="03_fastqc_trimmed/{sample}_R1_trimmed_fastqc.html",
        zip_r1="03_fastqc_trimmed/{sample}_R1_trimmed_fastqc.zip",
        html_r2="03_fastqc_trimmed/{sample}_R2_trimmed_fastqc.html",
        zip_r2="03_fastqc_trimmed/{sample}_R2_trimmed_fastqc.zip"
    threads: 4
    log:
        "logs/fastqc_trimmed_{sample}.log"
    shell:
        "fastqc -t {threads} {input.trimmed_r1} {input.trimmed_r2} -o 03_fastqc_trimmed/ &> {log}"

# MultiQC report on trimmed FastQC
rule multiqc_trimmed:
    input:
        expand("03_fastqc_trimmed/{sample}_R1_trimmed_fastqc.zip", sample=SAMPLES) + expand("03_fastqc_trimmed/{sample}_R2_trimmed_fastqc.zip", sample=SAMPLES)
    output:
        "03_multiqc/trimmed_fastqc_report.html"
    threads: 4
    log:
        "logs/multiqc_trimmed.log"
    shell:
        "multiqc 03_fastqc_trimmed/ -o 03_multiqc/ -n trimmed_fastqc_report &> {log}"

# Bowtie2 alignment to genome
rule align_bowtie2:
    input:
        trimmed_r1="02_trimmed/{sample}_R1_trimmed.fastq.gz",
        trimmed_r2="02_trimmed/{sample}_R2_trimmed.fastq.gz"
    output:
        sam=temp("04_aligned/{sample}.sam"),
        bam="04_aligned/{sample}.bam"
    params:
        bowtie2_index=GENOME_INDEX,
        min_insert=10,
        max_insert=700
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

# Filter low-quality reads
rule filter_bam:
    input:
        bam="04_aligned/{sample}.bam"
    output:
        filtered_bam="05_filtered/{sample}.bam"
    params:
        quality=20
    threads: 4
    log:
        "logs/filter_bam_{sample}.log"
    shell:
        "samtools view -b -q {params.quality} -o {output.filtered_bam} {input.bam} &> {log}"

# Convert filtered BAM to BED files with UCSC track header
rule bam_to_bed:
    input:
        filtered_bam="05_filtered/{sample}.bam"
    output:
        bed="06_bed/{sample}_filtered.bed"
    params:
        track_name="Filtered {sample} BAM",
        track_description="Filtered aligned reads in BED format for UCSC Genome Browser",
        visibility="4",
        color="0,0,255"
    threads: 4
    log:
        "logs/bam_to_bed_{sample}.log"
    shell:
        """
        echo 'track name="{params.track_name}" description="{params.track_name}" visibility={params.visibility} color={params.color}' > {output.bed}
        bedtools bamtobed -i {input.filtered_bam} >> {output.bed} 2> {log}
        """
        
         
# MACS2 peak calling rule
# rule call_peaks_macs2:
#     input:
#         filtered_bam="filtered/{sample}.bam"
#     output:
#         narrowpeak="peaks/{sample}_peaks.narrowPeak",
#         summits="peaks/{sample}_summits.bed"
#     params:
#         genome="2.7e9",  # Replace with "hs" (human), "mm" (mouse), or 2.7e9 for gibbon genome size
#         name="{sample}",
#         pvalue=1e-3  # Adjust as needed
#     log:
#         "logs/macs2_{sample}.log"
#     threads: 4
#     shell:
#         "macs2 callpeak -t {input.filtered_bam} -f BAM -g {params.genome} -n {params.name} --outdir peaks/ \
#         --keep-dup all --call-summits --pvalue {params.pvalue} &> {log}" 
        
# macs2 callpeak -t filtered_S7.bam -f BAMPE -q 0.1 -g 2.7e9 -s 50 --keep-dup all -n peak4_S7
        
# macs2 callpeak -t sample.bam -c control.bam -f BAM -g hs -n sample_with_control --outdir peaks/ --pvalue 1e-5
# -t: Specifies the treatment sample (e.g., the BAM file with your enriched regions, such as ChIP or Cut&Run data).
# -c: Specifies the control or input sample (e.g., a BAM file for a sample without enrichment, like whole-cell extract or input DNA).
# --pvalue or --qvalue: Sets the statistical threshold for peak calling.
# --outdir: Specifies the output directory for MACS2 results.       
        
        
        
