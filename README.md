# NucDensity_MNase
Calculate nucleosome density from MNase

## Required environment
Please use conda environment [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main)  
Make sure these tools are installed in your conda for __Step1__ and Step2
```
snakemake
deeptools
trim_galore
bwa
samtools
picard
bedtools
```
and these python packages need to be install when run Step3
```
pandas
numpy
seaborn
pybedtools
pyBigWig
click
```

## Tutorial
### Step1: Mapping (fastq to bam)
* **step1.1**: Create BWA index (if you don't have)
```
wget https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz
bwa index mm10_no_alt_analysis_set_ENCODE.fasta.gz
```
* **step1.2**: Prepare these three files: **samples.txt**, **ref.tsv** and **config.yaml**  

The repo contains the template for these three files, you can download it and revise their contents accordingly  

**Note**: if you have biological replicates, you could combine them at fastq-level (use cat) or bam-level (samtools merge); for simplicity, we combined at fastq-level; but theoritically, you should check the correlation of two replicates first (deeptools in bam-level) before combining.

* **step1.3**: Run snakemake pipeline
For paired-end data, use MNase_PE.smk (make --cores xx larger if you have more resources)
```
snakemake -s MNase_PE.smk --configfile config.yaml --cores 10 -p
```
For single-end data, contact me

### Step2: Nucleosome calling (bam to nucleosome.bed)
Once you finished Mapping step, you will see filtered_bam folder within the work_dir your entered in config.yaml
Then for each sample, do
* **step2.1**: convert bam into bed
```
bedtools bamtobed -i /work_dir/filtered_bam/mES_WT_ctr_dedup_filt.bam | grep -E '^chr([1-9]|1[0-9]|X|Y)[[:space:]]' - > mES_WT_ctr_dedup_filt.bed
# to make comparable bigwig later on, you might need to make the number of row in input.bed the same across different condition
shuf -n 55563772 mES_WT_ctr_dedup_filt.bed > mES_WT_ctr_dedup_filt_downsampled.bed
```
* **step2.2**: run iNPS
```
python iNPS_V1.2.3.py -i mES_WT_ctr_dedup_filt.bed -o wt_ctr/mES_WT_ctr --s_p s
```
* **step2.3**: convert like_wig into bigwig, make sure 1. the prefix is the same with the run iNPS step and 2. you are located in output folder of run iNPS step
```
bash likewig2bigwig_mm10.sh mES_WT_ctr mm10.chrom.size
```
* **step2.4**: converte Gathering.like_bed into tsv format
```
tail -n +23 mES_WT_ctr_Gathering.like_bed > mES_WT_ctr.tsv
```

### Step3: Nucleosome density measurement
* **step3.1**: measure the nucleosome density and occupancy (mm10 origins.bed and gene.txt are provided in the repo)
```
python NucDensity_cli_v1.py \
    -n Ctrl:/path/to/ctrl.tsv \
    -n IAA:/path/to/iaa.tsv \
    -n WT:/path/to/wt.tsv \
    -b Ctrl:/path/to/ctrl.bigWig \
    -b IAA:/path/to/iaa.bigWig \
    -b WT:/path/to/wt.bigWig \
    -o origins.bed \
    -g genes.txt \
    -w 1000,10000,50000,100000 \
    -d output_directory
```
If any bug pop up, please contact me :)

### Contact
Kun Fang: kf2799@cumc.columbia.edu
