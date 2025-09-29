# NucDensity_MNase
Calculate nucleosome density from MNase

# Required environment
Please use conda environment [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main)  
Make sure these tools are installed in your conda
```
snakemake
deeptools
trim_galore
bwa
samtools
picard
```
and these python packages need to be install
```
pandas
numpy
seaborn
```

# Tutorial
## Mapping (fastq to bam)
Before mapping, please prepare these three files: **mm10_samples.txt**, **mm10.tsv** and **config.yaml**  
The repo contains the template for these three files, you can download it and revise their contents accordingly  

For paired-end data, use MNase_PE.smk
```
snakemake -s MNase_PE.smk --configfile config.yaml --cores 10 -p
```


## Nucleosome calling (bam to nucleosome.bed)

## Nucleosome density measure


## Contact
Kun Fang: kf2799@cumc.columbia.edu
