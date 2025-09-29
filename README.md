# NucDensity_MNase
Calculate nucleosome density from MNase

## Required environment
Please use conda environment [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main)  
Make sure these tools are installed in your conda
```
snakemake
deeptools
trim_galore
bwa
samtools
picard
bedtools
```
and these python packages need to be install
```
pandas
numpy
seaborn
```

## Tutorial
### Mapping (fastq to bam)
* Create BWA index (if you don't have)
```
wget https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz
bwa index mm10_no_alt_analysis_set_ENCODE.fasta.gz
```
* Prepare these three files: **samples.txt**, **ref.tsv** and **config.yaml**  

The repo contains the template for these three files, you can download it and revise their contents accordingly  
**Note**: if you have biological replicates, you could combine them at fastq-level (use cat) or bam-level (samtools merge); for simplicity, we combined at fastq-level; but theoritically, you should check the correlation of two replicates first (deeptools in bam-level) before combining.
* Run snakemake pipeline
For paired-end data, use MNase_PE.smk
```
snakemake -s MNase_PE.smk --configfile config.yaml --cores 10 -p
```
For single-end data, contact me

### Nucleosome calling (bam to nucleosome.bed)
Once you finished Mapping step, you will see filtered_bam folder within the work_dir your entered in config.yaml
Then for each sample, do
* convert bam into bed
```
bedtools bamtobed -i /work_dir/filtered_bam/mES_WT_ctr_dedup_filt.bam | grep -E '^chr([1-9]|1[0-9]|X|Y)[[:space:]]' - > mES_WT_ctr_dedup_filt.bed
# to make comparable bigwig later on, you might need to make the number 
```
* run iNPS
```
python iNPS_V1.2.3.py -i mES_WT_ctr_dedup_filt.bed -o wt_ctr/mES_WT_ctr --s_p s
```
* convert like_wig into bigwig
```

```

### Nucleosome density measure


### Contact
Kun Fang: kf2799@cumc.columbia.edu
