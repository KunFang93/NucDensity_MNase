import os
import pandas as pd

# ─────────────────────────────────────────────────────────────────────────────
# CONFIG & GLOBAL VARS
# ─────────────────────────────────────────────────────────────────────────────
work_dir  = config["work_dir"]
species   = config["species"]
threads   = config["threads"]
binsize   = 5000 if species != "yeast" else 1

# Effective genome size (e.g. mm10, read length 150)
effectiveSize = 2494787038 #Read length 150 mm10

# ─────────────────────────────────────────────────────────────────────────────
# READ SAMPLE & REFERENCE TABLES
# ─────────────────────────────────────────────────────────────────────────────
SAMPLES = (
    pd.read_csv(config["samples"], sep="\t")
      .set_index("sample_id", drop=False)
      .sort_index()
)
# Make sure no two samples point to the same FASTQ
dupes = SAMPLES[SAMPLES["R1"].duplicated(keep=False)]
if not dupes.empty:
    raise ValueError(f"Duplicate FASTQ entries found for: {dupes.index.tolist()}")

REF = pd.read_csv(config["ref_files"], sep="\t", header=None, index_col=0)

# ─────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────
def parse_ppqt_fragment_size(fn):
    """
    Parse phantompeakqualtools output to extract estimated fragment length.
    Expected in 3rd column of the first line.
    """
    with open(fn) as f:
        line = f.readline().strip()
        est = int(line.split("\t")[2].split(',')[0])
        return min(est, 250)  # Cap to avoid overextension

def get_bwa_index():
    return REF.loc["bwa_index"][1]

def get_raw_fastq(wildcards):
    # Return list for compatibility with Snakemake’s expand/str formatting
    return SAMPLES.loc[wildcards.sample]["R1"].split(",")

def get_renamed_fastq(wildcards):
    return f"renamed_fq/{wildcards.sample}.fastq.gz"

def get_trimmed_fastq(wildcards):
    return f"trimmed_fq/{wildcards.sample}_trimmed.fq.gz"

def get_dedup_bam(wildcards):
    # uses UMI‐dir if requested
    if config.get("add_umi", False):
        return f"dedup_bam_umi_se/{wildcards.sample}_dedup.bam"
    else:
        return f"dedup_bam_se/{wildcards.sample}_dedup.bam"

def get_rule_all_input():
    bams = expand("dedup_bam_se/{sample}_dedup.bam", sample=SAMPLES["sample_id"].tolist())
    frag = expand("fragment_size/{sample}_ppqt.stats.txt", sample=SAMPLES["sample_id"].tolist())
    bw   = expand("bigWig/{sample}_dedup.CPM.bw", sample=SAMPLES["sample_id"].tolist())
    return bams + frag + bw

# ─────────────────────────────────────────────────────────────────────────────
# MASTER TARGET
# ─────────────────────────────────────────────────────────────────────────────
rule all:
    input:
        get_rule_all_input()


# ─────────────────────────────────────────────────────────────────────────────
# RULES
# ─────────────────────────────────────────────────────────────────────────────

rule merge_and_rename_fq_se:
    """
    Symlink raw FASTQ → renamed_fq/{sample}.fastq.gz
    """
    input:
        raw = get_raw_fastq
    output:
        renamed = "renamed_fq/{sample}.fastq.gz"
    shell:
        "ln -s {input.raw} {output.renamed}"

rule trim_galore_se:
    """
    Trim Galore (single‐end)
    """
    input:
        fastq = get_renamed_fastq
    output:
        trimmed = "trimmed_fq/{sample}_trimmed.fq.gz",
        report  = "trimmed_fq/{sample}.fastq.gz_trimming_report.txt"
    params:
        outdir = work_dir + "/trimmed_fq"
    threads: threads
    log:
        "logs/{sample}_trim_galore_se.log"
    shell:
        "(trim_galore -q 20 --stringency 3 --length 20 --gzip -o {params.outdir} {input.fastq}) 2> {log}"

rule bwa_map_se:
    """
    BWA‐MEM alignment → raw_bam/{sample}.bam + stats
    """
    input:
        idx   = get_bwa_index(),
        fastq = get_trimmed_fastq
    output:
        temp_bam = temp("raw_bam/{sample}.bam"),
        stats    = "raw_bam/{sample}_bam.stats.txt"
    threads: threads
    log:
        bwa   = "logs/{sample}_bwa_map.log",
        stats = "logs/{sample}_stats.log"
    shell:
        """
        bwa mem -M -t {threads} {input.idx} {input.fastq} \
          | samtools view -Sb --threads {threads} - > {output.temp_bam} 2> {log.bwa}
        samtools stats -@ {threads} {output.temp_bam} > {output.stats}
        """

rule samtools_sort_index_stats:
    """
    Sort, index & stats on raw BAM
    """
    input:
        "raw_bam/{sample}.bam"
    output:
        sorted = "raw_bam/{sample}_sorted.bam",
        stats  = "raw_bam/{sample}_sorted.bam.stats.txt"
    threads: threads
    shell:
        """
        samtools sort -@ {threads} -o {output.sorted} {input}
        samtools index -@ {threads} {output.sorted}
        samtools stats -@ {threads} {output.sorted} > {output.stats}
        """

rule samtools_markdup_stats_se:
    """
    Mark duplicates (name‐sort → fixmate → pos‐sort → markdup → filter by MAPQ and flags)
    """
    input:
        "raw_bam/{sample}_sorted.bam"
    output:
        bam  = "dedup_bam_se/{sample}_dedup.bam",
        stat = "dedup_bam_se/{sample}_dedup.bam.stats.txt"
    threads: threads
    temp:
        temp_name_sorted = "raw_bam/{sample}_name_sorted.bam",
        temp_fixed       = "raw_bam/{sample}_fixed.bam",
        temp_pos_sorted  = "raw_bam/{sample}_pos_sorted.bam",
        temp_dedup_raw   = "dedup_bam_se/{sample}_dedup.raw.bam"
    shell:
        """
        samtools sort -n -@ {threads} -o {input}.name_sorted.bam {input}
        samtools fixmate -m {input}.name_sorted.bam {input}.fixed.bam
        samtools sort -@ {threads} -o {input}.pos_sorted.bam {input}.fixed.bam
        samtools markdup -r {input}.pos_sorted.bam {output.bam}.raw
        samtools view -@ {threads} -q 10 -F 3844 -b {output.bam}.raw > {output.bam}
        samtools index -@ {threads} {output.bam}
        samtools stats -@ {threads} {output.bam} > {output.stat}
        rm -f {output.bam}.raw
        """
        """

rule estimate_fragment_length_ppqt:
    """
    Estimate fragment length for single-end data using phantompeakqualtools (run_spp.R).
    """
    input:
        bam = "dedup_bam_se/{sample}_dedup.bam"
    output:
        stats = "fragment_size/{sample}_ppqt.stats.txt",
        plot  = "fragment_size/{sample}_crosscorr.pdf"
    log:
        "logs/{sample}_ppqt.log"
    shell:
        """
        mkdir -p fragment_size
        run_spp.R \
            -c={input.bam} \
            -savp={output.plot} \
            -out={output.stats} \
            > {log} 2>&1
        """

rule generate_dedup_bw_se:
    """
    Generate CPM-normalized bigWig from SE BAM, using estimated fragment length.
    """
    input:
        bam     = "dedup_bam_se/{sample}_dedup.bam",
        metrics = "fragment_size/{sample}_ppqt.stats.txt"
    output:
        bw = "bigWig/{sample}_dedup.CPM.bw"
    threads: threads
    params:
        extend_len = lambda wc: parse_ppqt_fragment_size(f"fragment_size/{wc.sample}_ppqt.stats.txt")
    shell:
        """
        bamCoverage \
          --bam {input.bam} \
          --outFileName {output.bw} \
          --binSize 1 \
          --normalizeUsing CPM \
          --extendReads {params.extend_len} \
          -p {threads}
        """