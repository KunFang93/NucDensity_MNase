#!/usr/bin/env python3

import pandas as pd
import pybedtools as pybt
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pyBigWig
import click
from pathlib import Path


@click.command()
@click.option('--nuc-files', '-n', multiple=True, required=True, 
              help='Nucleosome TSV files (at least 2 required). Format: name:path')
@click.option('--bw-files', '-b', multiple=True, required=True,
              help='BigWig files corresponding to nucleosome files. Format: name:path')
@click.option('--origin-bed', '-o', required=True, type=click.Path(exists=True),
              help='Origins BED file')
@click.option('--gene-annotation', '-g', required=True, type=click.Path(exists=True),
              help='Gene annotation file')
@click.option('--origin-windows', '-w', default='1000,10000,50000,100000',
              help='Comma-separated list of window sizes for origin analysis (default: 10000,50000,100000)')
@click.option('--flank-size', '-f', default=100000, type=int,
              help='Flanking region size around origins (default: 100000)')
@click.option('--height-quantile', '-q', default=0.1, type=float,
              help='Height quantile cutoff for nucleosome filtering (default: 0.1)')
@click.option('--output-dir', '-d', default='.', type=click.Path(),
              help='Output directory for plots (default: current directory)')
def main(nuc_files, bw_files, origin_bed, gene_annotation, origin_windows, 
         flank_size, height_quantile, output_dir):
    """
    Nucleosome density analysis CLI tool.
    
    This tool analyzes nucleosome positioning around origins of replication
    and generates various plots including density distributions and profiles.
    
    Example usage:
    python NucDensity_cli_v1.py \\
        -n Ctrl:/path/to/ctrl.tsv \\
        -n IAA:/path/to/iaa.tsv \\
        -n WT:/path/to/wt.tsv \\
        -b Ctrl:/path/to/ctrl.bigWig \\
        -b IAA:/path/to/iaa.bigWig \\
        -b WT:/path/to/wt.bigWig \\
        -o origins.bed \\
        -g genes.txt \\
        -w 10000,50000,100000
    """
    
    # Validate inputs
    if len(nuc_files) < 2:
        raise click.BadParameter("At least 2 nucleosome files are required")
    
    if len(nuc_files) != len(bw_files):
        raise click.BadParameter("Number of nucleosome files must match number of BigWig files")
    
    # Parse file inputs
    nuc_data = {}
    bw_paths = {}
    
    for nuc_file in nuc_files:
        if ':' not in nuc_file:
            raise click.BadParameter(f"Nucleosome file must be in format 'name:path', got: {nuc_file}")
        name, path = nuc_file.split(':', 1)
        if not Path(path).exists():
            raise click.FileError(f"Nucleosome file not found: {path}")
        nuc_data[name] = path
    
    for bw_file in bw_files:
        if ':' not in bw_file:
            raise click.BadParameter(f"BigWig file must be in format 'name:path', got: {bw_file}")
        name, path = bw_file.split(':', 1)
        if not Path(path).exists():
            raise click.FileError(f"BigWig file not found: {path}")
        bw_paths[name] = path
    
    # Check that names match between nuc and bw files
    if set(nuc_data.keys()) != set(bw_paths.keys()):
        raise click.BadParameter("Sample names must match between nucleosome and BigWig files")
    
    # Parse origin windows
    origin_windows_list = [int(w.strip()) for w in origin_windows.split(',')]
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    click.echo(f"Loading data for {len(nuc_data)} samples: {', '.join(nuc_data.keys())}")
    
    # --- 1) Load nucleosome data ------------------------------
    nucs = {}
    for name, path in nuc_data.items():
        click.echo(f"Loading nucleosome data: {name}")
        nucs[name] = pd.read_table(path, header=None)
    
    # Load origins and gene annotation
    click.echo("Loading origins and gene annotation...")
    ori_bed = pd.read_csv(
        origin_bed, sep=r'\s+', usecols=[0,1,2],
        names=['chrom','start','end'], header=None
    )
    gene_df = pd.read_table(gene_annotation, header=None)
    
    # Filter unreliable nucleosomes using the first sample's quantile
    first_sample = list(nucs.keys())[0]
    height_cut = nucs[first_sample][5].quantile(height_quantile)
    click.echo(f"Filtering nucleosomes with height cutoff: {height_cut:.3f}")
    
    for name in nucs:
        before_count = len(nucs[name])
        nucs[name] = nucs[name][nucs[name][5] > height_cut]
        after_count = len(nucs[name])
        click.echo(f"  {name}: {before_count} -> {after_count} nucleosomes")
    
    # --- 2) Build flanking regions ----------------------
    click.echo(f"Building flanking regions (±{flank_size:,} bp)...")
    ori_bed['start_flank'] = (ori_bed['start'] - flank_size).clip(lower=0)
    ori_bed['end_flank'] = ori_bed['end'] + flank_size
    flank_bed = pybt.BedTool.from_dataframe(
        ori_bed[['chrom','start_flank','end_flank']].rename(
            columns={'start_flank':'start','end_flank':'end'}
        )
    )
    
    # --- 3) Make BedTools for nucleosome calls --------
    nuc_beds = {}
    for name, data in nucs.items():
        nuc_beds[name] = pybt.BedTool.from_dataframe(
            data[[0,1,2,5]].rename(columns={5:'score'})
        )
    
    # --- 4) Generate NumberNucs_bar.png --------------------------
    click.echo("Generating nucleosome count bar plot...")
    counts = {}
    for name, bed in nuc_beds.items():
        counts[name] = len(bed.intersect(flank_bed, u=True))
    
    plt.figure(figsize=(8,6))
    # Generate distinct colors for each sample
    colors = plt.cm.Set3(np.linspace(0, 1, len(counts)))
    bars = plt.bar(counts.keys(), counts.values(), color=colors)
    
    for bar in bars:
        y = bar.get_height()
        plt.text(bar.get_x()+bar.get_width()/2, y, f'{y:,}', ha='center', va='bottom')
    
    plt.ylabel('Nucleosome Count')
    plt.title(f'Nucleosomes within ±{flank_size//1000} kb of {len(ori_bed)} Origins')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_path / 'NumberNucs_bar.png', dpi=300)
    plt.close()
    
    # --- 5) Origin_Nuc_Density_violin_seaborn.png --------
    click.echo("Generating nucleosome density violin plot...")
    
    def count_df(bed):
        return (flank_bed.intersect(bed, c=True)
                         .to_dataframe(names=['chrom','start','end','nuc_count'])
                         .assign(density=lambda df: df['nuc_count'].replace(0, np.nan)/(2*flank_size/1000))
               )
    
    # Collect density data for all samples
    density_data = []
    condition_data = []
    
    for name, bed in nuc_beds.items():
        df = count_df(bed)
        density_data.extend(df['density'].tolist())
        condition_data.extend([name] * len(df))
    
    violin_df = pd.DataFrame({
        'density': density_data,
        'Condition': condition_data
    })
    
    plt.figure(figsize=(10,6))
    # Use different colors for each violin [[memory:8286407]]
    palette = dict(zip(counts.keys(), colors))
    ax = sns.violinplot(x='Condition', y='density', data=violin_df, 
                       inner='quartile', cut=0, scale='width', palette=palette)
    ax.set_ylabel('Nucleosome density\n(# nuc per kb)')
    ax.set_title(f'Per-origin nucleosome density (n={len(ori_bed)})')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_path / 'Origin_Nuc_Density_violin_seaborn.png', dpi=300)
    plt.close()
    
    # --- 6) Average signal around TSS using BigWig -------
    click.echo("Generating TSS profiles...")
    
    # Prepare gene TSS from annotation
    gene_df = gene_df.rename(columns={0:'chrom', 1:'start', 2:'end', 3:'strand'})
    gene_df['TSS'] = np.where(gene_df['strand']=='+', gene_df['start'], gene_df['end'])
    
    # Open BigWig files
    bw_objs = {name: pyBigWig.open(path) for name, path in bw_paths.items()}
    
    def profile_bw(bw, window, step=1):
        pos = np.arange(-window, window+1, step)
        sum_cov = np.zeros_like(pos, dtype=float)
        count = np.zeros_like(pos, dtype=float)
        for _, row in gene_df.iterrows():
            chrom, tss = row['chrom'], int(row['TSS'])
            start = max(0, tss-window); end = tss+window
            vals = bw.values(chrom, start, end+1, numpy=True)
            if vals is None or len(vals) < len(pos): continue
            if step != 1: vals = vals[::step]
            mask = ~np.isnan(vals)
            sum_cov[mask] += vals[mask]
            count[mask] += 1
        return pos, sum_cov / np.where(count==0, np.nan, count)
    
    # Generate TSS metaplots for ±1kb window, step=10bp
    window = 1000; step = 10
    plt.figure(figsize=(10,6))
    for i, (name, bw) in enumerate(bw_objs.items()):
        x, prof = profile_bw(bw, window, step)
        plt.plot(x, prof, label=name, color=colors[i], linewidth=2)
    
    plt.axvline(0, color='black', linestyle='--', alpha=0.7)
    plt.xlabel('Position relative to TSS (bp)')
    plt.ylabel('Mean coverage')
    plt.title('Nucleosome signal around TSS')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path / f'Profile_TSS_±{window}bp.png', dpi=300)
    plt.close()
    
    # --- 7) Average signal around Origins using BigWig -------
    click.echo(f"Generating origin profiles for windows: {origin_windows_list}")
    
    def profile_origin_bw(bw, origin_df, window, step=10):
        pos = np.arange(-window, window + 1, step)
        sum_cov = np.zeros_like(pos, dtype=float)
        count = np.zeros_like(pos, dtype=float)
        for _, row in origin_df.iterrows():
            chrom, center = row['chrom'], (row['start'] + row['end']) // 2
            start = max(0, center - window)
            end = center + window
            vals = bw.values(chrom, start, end + 1, numpy=True)
            if vals is None or len(vals) < len(pos): continue
            if step != 1: vals = vals[::step]
            mask = ~np.isnan(vals)
            sum_cov[mask] += vals[mask]
            count[mask] += 1
        return pos, sum_cov / np.where(count == 0, np.nan, count)
    
    # Generate origin profiles for each window size
    for w in origin_windows_list:
        plt.figure(figsize=(10,6))
        for i, (name, bw) in enumerate(bw_objs.items()):
            x, prof = profile_origin_bw(bw, ori_bed, w)
            plt.plot(x, prof, label=name, color=colors[i], linewidth=2)
        
        plt.axvline(0, color='black', linestyle='--', alpha=0.7)
        plt.xlabel('Position relative to origin (bp)')
        plt.ylabel('Mean coverage')
        plt.title(f'Nucleosome signal around origins (±{w//1000} kb)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_path / f'Profile_Origin_±{w//1000}kb.png', dpi=300)
        plt.close()
    
    # Close BigWig handles
    for bw in bw_objs.values(): 
        bw.close()
    
    click.echo(f"Analysis complete! Output files saved to: {output_path}")
    click.echo("Generated plots:")
    click.echo("  - NumberNucs_bar.png")
    click.echo("  - Origin_Nuc_Density_violin_seaborn.png")
    click.echo(f"  - Profile_TSS_±{window}bp.png")
    for w in origin_windows_list:
        click.echo(f"  - Profile_Origin_±{w//1000}kb.png")


if __name__ == '__main__':
    main()