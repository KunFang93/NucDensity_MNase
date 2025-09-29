import os
import pandas as pd

###########################################
## paths for pipeline and/or reference data
work_dir = config["work_dir"]
species = config["species"]
binsize = 5000 if species != "yeast" else 1

## read in sample and corresponding fq files table
SAMPLES = (
    pd.read_csv(config["samples"], sep="\t")
    .set_index("sample_id", drop=False)
    .sort_index()
)

# Check for samples where R1 and R2 file names are identical
duplicate_samples = SAMPLES[SAMPLES["R1"] == SAMPLES["R2"]]
if not duplicate_samples.empty:
    raise ValueError(
        f"Error: The following samples have the same R1 and R2 file names: {duplicate_samples['sample_id'].tolist()}"
    )

# effectiveSize = 2862010428 #Read length 150 hg38
effectiveSize = 2494787038 #Read length 150 mm10

## read in reference files' info
REF = pd.read_csv(config["ref_files"], sep="\t", header=None, index_col=0)

threads = config["threads"]

##  get corresponding bwa_index
def get_bwa_index():
    return REF.loc["bwa_index"][1]

def get_bowtie2_index():
    return REF.loc["bowtie2_index"][1]

def get_raw_fastq_se(wildcards):
    if config["paired-end"] == False:
        return SAMPLES.loc[wildcards.sample]["R1"].split(",")
    else:
        return ""

def get_raw_fastq_pe_R1(wildcards):
    if config["paired-end"]:
        return SAMPLES.loc[wildcards.sample]["R1"].split(",")
    else:
        return ""

def get_raw_fastq_pe_R2(wildcards):
    if config["paired-end"]:
        print()
        return SAMPLES.loc[wildcards.sample]["R2"].split(",")
    else:
        return ""

def get_renamed_fastq(wildcards):
    if config["paired-end"]:
        R1 = "renamed_fq/{}_R1.fastq.gz".format(wildcards.sample)
        R2 = "renamed_fq/{}_R2.fastq.gz".format(wildcards.sample)
        return R1 + R2
    else:
        return "{}/renamed_fq/{}.fastq.gz".format(work_dir,wildcards.sample)

def get_fastq_4trim(wildcards):
    if config["paired-end"]:
        R1 = "renamed_fq/{}_R1.fastq.gz".format(wildcards.sample),
        R2 = "renamed_fq/{}_R2.fastq.gz".format(wildcards.sample),
        return R1 + R2
    else:
        return "renamed_fq/{}.fastq.gz".format(wildcards.sample)

##################################
## get trimmed fastq files for BWA
def get_trimmed_fastq(wildcards):
    if config["paired-end"]:
        R1 = "trimmed_fq/{}_R1_val_1.fq.gz".format(wildcards.sample),
        R2 = "trimmed_fq/{}_R2_val_2.fq.gz".format(wildcards.sample),
        return R1 + R2
    else:
        return "trimmed_fq/{}_trimmed.fq.gz".format(wildcards.sample)

########################################
## get dedup bam files for meth_qc_quant
def get_dedup_bam(wildcards):
    if config["paired-end"]:
        return "dedup_bam_pe/{}_dedup.bam".format(wildcards.sample)
    else:
        return "dedup_bam_se/{}_dedup.bam".format(wildcards.sample)

def get_fastqc_stats(wildcards):
    if config["paired-end"]:
        r1_raw  = expand("fastqc_pe/{samples}_R1_fastqc.zip", samples = wildcards.sample),
        r2_raw  = expand("fastqc_pe/{samples}_R2_fastqc.zip", samples = wildcards.sample),
        r1_trim = expand("fastqc_pe/{samples}_R1_val_1_fastqc.zip", samples = wildcards.sample),
        r2_trim = expand("fastqc_pe/{samples}_R2_val_2_fastqc.zip", samples = wildcards.sample),
        return r1_raw + r2_raw + r1_trim + r2_trim
    else:
         r1_raw  = expand("fastqc_se/{samples}_fastqc.zip", samples = wildcards.sample),
         r1_trim = expand("fastqc_se/{samples}_trimmed_fastqc.zip", samples = wildcards.sample),
         return r1_raw + r1_trim

def get_dedup_bam_stats(wildcards):
    if config["paired-end"] and config["add_umi"]:
        return expand("dedup_bam_umi_pe/{samples}_dedup.bam.stats.txt", samples = wildcards.sample)
    elif config["paired-end"] and config["add_umi"] == False:
        return expand("dedup_bam_pe/{samples}_dedup.bam.stats.txt", samples = wildcards.sample)
    elif config["paired-end"] == False and config["add_umi"]:
        return expand("dedup_bam_umi_se/{samples}_dedup.bam.stats.txt", samples = wildcards.sample)
    else:
        return expand("dedup_bam_se/{samples}_dedup.bam.stats.txt", samples = wildcards.sample)

###########################################
## Update rule all input to include filtered BAM and remove plus/minus BW outputs
def get_rule_all_input():
    bam_out    = expand("dedup_bam_pe/{sample}_dedup.bam", sample=SAMPLES["sample_id"].values.tolist())
    frag_out   = expand("fragment_size/{sample}_insert_size_metrics.txt", sample=SAMPLES["sample_id"].values.tolist())
    bw_out1    = expand("bigWig/{sample}_dedup_filt.CPM.bw", sample=SAMPLES["sample_id"].values.tolist())
    filt_bam   = expand("filtered_bam/{sample}_dedup_filt.bam", sample=SAMPLES["sample_id"].values.tolist())
    return bam_out + frag_out + bw_out1 + filt_bam

rule all:
    input:
        get_rule_all_input()

###########################################
## Pipeline Rules

# Merge and rename FASTQ files for paired-end reads
rule merge_and_rename_fq_pe:
    input:
        R1 = get_raw_fastq_pe_R1,
        R2 = get_raw_fastq_pe_R2,
    output:
        "renamed_fq/{sample}_R1.fastq.gz",
        "renamed_fq/{sample}_R2.fastq.gz",
    shell:
        "ln -s {input.R1} {output[0]} && ln -s {input.R2} {output[1]}"

# Trim Galore for paired-end reads
rule trim_galore_pe:
    input:
        get_fastq_4trim
    output:
        temp("trimmed_fq/{sample}_R1_val_1.fq.gz"),
        temp("trimmed_fq/{sample}_R2_val_2.fq.gz"),
        "trimmed_fq/{sample}_R1.fastq.gz_trimming_report.txt",
        "trimmed_fq/{sample}_R2.fastq.gz_trimming_report.txt",
    params:
        path = work_dir + "/trimmed_fq"
    threads: threads
    log:
        "logs/{sample}_trim_galore_pe.log"
    shell:
        "(trim_galore -q 20 --stringency 3 --length 20 --cores {threads} --paired -o {params.path} {input}) 2> {log}"

# BWA alignment
rule bwa_map:
    input:
        get_bwa_index(),
        get_trimmed_fastq
    output:
        temp_bam = temp("raw_bam/{sample}.bam"),
        stats = "raw_bam/{sample}_bam.stats.txt"
    threads: threads
    resources: mem_mb=80000
    log:
        log_bwa="logs/{sample}_bwa_map.log",
        log_stats="logs/{sample}_stats.log"
    shell:
        """
        (bwa mem -M -K 100000000 -t {threads} {input} | samtools view -Sb --threads {threads} - > {output.temp_bam}) 2> {log.log_bwa}
        samtools stats -@ {threads} {output.temp_bam} > {output.stats}
        """

# Sort, fixmate, index and collect stats for BAM
rule samtools_sort_index_stats:
    input:
        "raw_bam/{sample}.bam"
    output:
        "raw_bam/{sample}_sorted.bam",
        "raw_bam/{sample}_sorted.bam.stats.txt"
    threads: threads
    resources: mem_mb=80000
    shell:
        "(samtools fixmate -@ {threads} -m {input} - | samtools sort -@ {threads} -o {output[0]} && samtools index -@ {threads} {output[0]} && samtools stats -@ {threads} {output[0]} > {output[1]})"

# Deduplication using samtools markdup
rule samtools_markdup_stats_pe:
    input:
        "raw_bam/{sample}_sorted.bam"
    output:
        bam  = "dedup_bam_pe/{sample}_dedup.bam",
        stat = "dedup_bam_pe/{sample}_dedup.bam.stats.txt"
    threads: threads
    shell:
        "(samtools view -b -f 2 -F 1804 --threads {threads} {input} | samtools markdup -@ {threads} -r - {output.bam} && samtools index -@ {threads} {output.bam} && samtools stats -@ {threads} {output.bam} > {output.stat})"


# Filter deduplicated BAM using alignmentSieve for MNase-seq fragments
rule alignmentSieve_filter:
    input:
        bam = "dedup_bam_pe/{sample}_dedup.bam"
    output:
        filt_bam = "filtered_bam/{sample}_dedup_filt.bam"
    threads: 4
    shell:
        """
        alignmentSieve -p 4 -b {input.bam} -o {output.filt_bam} --minFragmentLength 135 --maxFragmentLength 175 &&
        samtools index {output.filt_bam}
        """

# Generate non–strand-specific BigWig file
rule generate_dedup_bw:
    input:
        "filtered_bam/{sample}_dedup_filt.bam"
    output:
        "bigWig/{sample}_dedup_filt.CPM.bw"
    threads: threads
    params:
        eSize = effectiveSize
    shell:
        "bamCoverage --bam {input} -o {output} --binSize 1 --normalizeUsing CPM --extendReads -p {threads}"

# Collect insert size metrics using Picard
rule insert_size:
    input:
        get_dedup_bam
    output:
        txt = "fragment_size/{sample}_insert_size_metrics.txt",
        hist = "fragment_size/{sample}_insert_size_histogram.pdf"
    log:
        "logs/{sample}_picard_insert_size.log"
    shell:
        "(picard CollectInsertSizeMetrics M=0.05 I={input} O={output.txt} H={output.hist}) 2> {log}"