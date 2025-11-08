#!/usr/bin/env python3
"""
NucDensity_cli_v2.py

Nucleosome density analysis CLI tool with improved exclusion handling:
- When deciding whether to exclude an origin, we check the *extended* region
  that will be profiled (max of flank_size and origin_windows), so origins are
  excluded if their profiling windows would overlap any excluded intervals.
"""
from pathlib import Path
import click
import numpy as np
import pandas as pd
import pyBigWig
import pybedtools as pybt
import matplotlib.pyplot as plt
import seaborn as sns

# Built-in default excluded regions (used when --exclude-region file is NOT provided)
_DEFAULT_EXCLUDE_REGIONS = [
    ("chr17", 39841222, 39851317),
    ("chr8", 19672781, 19938153),
    ("chr1", 24610569, 24617377),
]


@click.command()
@click.option('--nuc-files', '-n', multiple=True, required=True,
              help='Nucleosome TSV files (at least 2 required). Format: name:path')
@click.option('--bw-files', '-b', multiple=True, required=True,
              help='BigWig files corresponding to nucleosome files. Format: name:path')
@click.option('--origin-bed', '-o', required=True, type=click.Path(exists=True),
              help='Origins BED file (3 columns: chrom start end)')
@click.option('--gene-annotation', '-g', required=True, type=click.Path(exists=True),
              help='Gene annotation file (expects at least chrom,start,end,strand in cols 0-3)')
@click.option('--origin-windows', '-w', default='1000,10000,50000,100000',
              help='Comma-separated list of window sizes for origin analysis (bp). Default: 1000,10000,50000,100000')
@click.option('--flank-size', '-f', default=100000, type=int,
              help='Flanking region size around origins in bp (default: 100000)')
@click.option('--height-quantile', '-q', default=0.1, type=float,
              help='Height quantile cutoff for nucleosome filtering (default: 0.1)')
@click.option('--output-dir', '-d', default='.', type=click.Path(),
              help='Output directory for plots (default: current directory)')
@click.option('--exclude-region', '-e', 'exclude_region_file', type=click.Path(exists=True),
              help=('Optional TSV/BED-like file with at least 3 columns: chrom start end. '
                    'If omitted, built-in default regions (3 intervals) are used.'))
def main(nuc_files, bw_files, origin_bed, gene_annotation, origin_windows,
         flank_size, height_quantile, output_dir, exclude_region_file):
    """
    Nucleosome density analysis CLI tool.

    Notes:
      - Origins are excluded if their *profiling window* (center ± max_extend)
        overlaps any exclude interval. max_extend = max(flank_size, max(origin_windows_list)).
      - Genes are filtered by TSS ±1000bp as before.
    """
    # ------------- Input validation ----------------
    if len(nuc_files) < 2:
        raise click.BadParameter("At least 2 nucleosome files are required")

    if len(nuc_files) != len(bw_files):
        raise click.BadParameter("Number of nucleosome files must match number of BigWig files")

    # Parse file inputs (format: Name:/path/to/file)
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

    if set(nuc_data.keys()) != set(bw_paths.keys()):
        raise click.BadParameter("Sample names must match between nucleosome and BigWig files")

    # Parse origin windows list
    origin_windows_list = [int(w.strip()) for w in origin_windows.split(',') if w.strip()]
    if len(origin_windows_list) == 0:
        raise click.BadParameter("origin-windows must contain at least one integer")

    # Prepare output path
    outpath = Path(output_dir)
    outpath.mkdir(parents=True, exist_ok=True)

    click.echo(f"Loading data for {len(nuc_data)} samples: {', '.join(nuc_data.keys())}")
    if exclude_region_file:
        click.echo(f"Exclude-regions file provided: {exclude_region_file}")
    else:
        click.echo("No exclude-regions file provided: using built-in default 3 regions")

    # ----------------- 1) Load nucleosome TSVs -----------------------
    nucs = {}
    for name, path in nuc_data.items():
        click.echo(f"Loading nucleosome data: {name} -> {path}")
        nucs[name] = pd.read_table(path, header=None)

    # ----------------- Load origins and gene annotation ---------------
    click.echo("Loading origins and gene annotation...")
    ori_bed = pd.read_csv(origin_bed, sep=r'\s+', usecols=[0, 1, 2], header=None, names=['chrom', 'start', 'end'])
    gene_df = pd.read_table(gene_annotation, header=None)

    # ----------------- Build exclude BedTool -------------------------
    if exclude_region_file:
        # Load first 3 columns from provided file
        ex_df = pd.read_csv(exclude_region_file, sep=r'\s+', header=None, usecols=[0, 1, 2])
        ex_df.columns = ['chrom', 'start', 'end']
        click.echo(f"  Loaded {len(ex_df)} exclude intervals from file")
        exclude_bed = pybt.BedTool.from_dataframe(ex_df)
    else:
        # Use built-in defaults
        ex_df = pd.DataFrame(_DEFAULT_EXCLUDE_REGIONS, columns=['chrom', 'start', 'end'])
        click.echo(f"  Using {len(ex_df)} built-in exclude intervals")
        exclude_bed = pybt.BedTool.from_dataframe(ex_df)

    # ----------------- Determine max extension for origin exclusion ----
    # We want to exclude origins if any profiling window we will use overlaps exclude regions.
    max_origin_extend = max(flank_size, max(origin_windows_list))
    click.echo(f"Max origin extension used for exclusion checks: {max_origin_extend} bp")

    # Build an "extended origin" BedTool: origin center ± max_origin_extend
    ori_centers = ori_bed.copy()
    ori_centers['center'] = ((ori_centers['start'].astype(int) + ori_centers['end'].astype(int)) // 2).astype(int)
    ori_centers['ext_start'] = (ori_centers['center'] - max_origin_extend).clip(lower=0)
    ori_centers['ext_end'] = ori_centers['center'] + max_origin_extend
    ori_extended_bt = pybt.BedTool.from_dataframe(ori_centers[['chrom', 'ext_start', 'ext_end']].rename(columns={'ext_start':'start','ext_end':'end'}))

    # ----------------- Filter origins by overlap with exclude (using extended windows) -----------
    click.echo("Filtering origins whose extended profiling windows overlap excluded regions...")
    ori_filtered_bt = ori_extended_bt.intersect(exclude_bed, v=True)  # -v keeps those not overlapping
    ori_filtered_df = ori_filtered_bt.to_dataframe(names=['chrom','start','end'])
    # Note: ori_filtered_df currently contains the extended coordinates; map back to original origin centers to keep canonical start/end
    # Reconstruct a set of centers kept, then filter original ori_bed by center
    kept_centers = set(((ori_filtered_df['start'] + ori_filtered_df['end']) // 2).astype(int).tolist())
    # compute centers for original ori_bed and keep rows whose center in kept_centers
    ori_bed['center'] = ((ori_bed['start'].astype(int) + ori_bed['end'].astype(int)) // 2).astype(int)
    before_origins = len(ori_bed)
    ori_bed = ori_bed[ori_bed['center'].isin(kept_centers)].copy().reset_index(drop=True)
    # drop center column for later use
    ori_bed.drop(columns=['center'], inplace=True)
    click.echo(f"  Origins: {before_origins} -> {len(ori_bed)} after excluding by extended windows")

    if len(ori_bed) == 0:
        raise click.BadParameter("No origins remain after applying exclude regions based on extended profiling windows. Adjust exclusions or provide a different exclude file.")

    # ----------------- Filter genes for TSS profiling ----------------
    # The script uses a TSS profiling window of ±1000 bp; keep that consistent here.
    tss_profile_window = 1000
    click.echo(f"Filtering genes whose TSS ±{tss_profile_window} bp overlap excluded regions (if any).")
    try:
        gene_df = gene_df.rename(columns={0: 'chrom', 1: 'start', 2: 'end', 3: 'strand'})
        gene_df['TSS'] = np.where(gene_df['strand'] == '+', gene_df['start'], gene_df['end'])
        gene_df['tss_start'] = (gene_df['TSS'] - tss_profile_window).clip(lower=0)
        gene_df['tss_end'] = gene_df['TSS'] + tss_profile_window
        gene_bt = pybt.BedTool.from_dataframe(gene_df[['chrom', 'tss_start', 'tss_end']].rename(columns={'tss_start': 'start', 'tss_end': 'end'}))
        gene_filtered_bt = gene_bt.intersect(exclude_bed, v=True)
        gene_filtered_df = gene_filtered_bt.to_dataframe(names=['chrom', 'start', 'end'])
        gene_df['__key__'] = gene_df['chrom'].astype(str) + ':' + gene_df['tss_start'].astype(str) + '-' + gene_df['tss_end'].astype(str)
        gene_filtered_df['__key__'] = gene_filtered_df['chrom'].astype(str) + ':' + gene_filtered_df['start'].astype(str) + '-' + gene_filtered_df['end'].astype(str)
        keep_keys = set(gene_filtered_df['__key__'])
        before_genes = len(gene_df)
        gene_df = gene_df[gene_df['__key__'].isin(keep_keys)].copy().reset_index(drop=True)
        gene_df.drop(columns=['__key__'], inplace=True)
        click.echo(f"  Genes for TSS profiling: {before_genes} -> {len(gene_df)} after excluding overlaps")
    except Exception as e:
        click.echo("Warning: gene annotation parsing for TSS filtering failed: " + str(e))
        click.echo("Proceeding without filtering genes for TSS profiling.")
        # Try to ensure TSS exists even if parsing failed
        try:
            gene_df = gene_df.rename(columns={0: 'chrom', 1: 'start', 2: 'end', 3: 'strand'})
            gene_df['TSS'] = np.where(gene_df['strand'] == '+', gene_df['start'], gene_df['end'])
        except Exception:
            raise click.BadParameter("Gene annotation must contain chrom,start,end,strand in cols 0-3 for TSS profiling.")

    # ----------------- Filter unreliable nucleosomes ----------------
    first_sample = list(nucs.keys())[0]
    if 5 not in nucs[first_sample].columns and 5 >= nucs[first_sample].shape[1]:
        raise click.BadParameter("Nucleosome TSVs expected to have 'height' (score) in column index 5. Please ensure format.")
    height_cut = nucs[first_sample][5].quantile(height_quantile)
    click.echo(f"Filtering nucleosomes with height cutoff (quantile={height_quantile}): {height_cut:.3f}")
    for name in nucs:
        before_count = len(nucs[name])
        nucs[name] = nucs[name][nucs[name][5] > height_cut]
        after_count = len(nucs[name])
        click.echo(f"  {name}: {before_count} -> {after_count} nucleosomes")

    # ----------------- Build flanking windows for (filtered) origins ----
    click.echo(f"Building flanking regions (±{flank_size:,} bp) around filtered origins...")
    ori_bed['start_flank'] = (ori_bed['start'] - flank_size).clip(lower=0)
    ori_bed['end_flank'] = ori_bed['end'] + flank_size
    flank_bed = pybt.BedTool.from_dataframe(
        ori_bed[['chrom', 'start_flank', 'end_flank']].rename(columns={'start_flank': 'start', 'end_flank': 'end'})
    )

    # ----------------- Make BedTools for nucleosome calls ----------------
    nuc_beds = {}
    for name, data in nucs.items():
        # expect columns: 0=chrom,1=start,2=end, 5=score
        nuc_beds[name] = pybt.BedTool.from_dataframe(data[[0, 1, 2, 5]].rename(columns={5: 'score'}))

    # ----------------- NumberNucs_bar.png -------------------------------
    click.echo("Generating nucleosome count bar plot...")
    counts = {}
    for name, bed in nuc_beds.items():
        counts[name] = len(bed.intersect(flank_bed, u=True))

    plt.figure(figsize=(8, 6))
    colors = plt.cm.Set3(np.linspace(0, 1, len(counts)))
    bars = plt.bar(list(counts.keys()), list(counts.values()), color=colors)
    for bar in bars:
        y = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, y, f'{y:,}', ha='center', va='bottom')
    plt.ylabel('Nucleosome Count')
    plt.title(f'Nucleosomes within ±{flank_size // 1000} kb of {len(ori_bed)} Origins')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(outpath / 'NumberNucs_bar.png', dpi=300)
    plt.close()

    # ----------------- Origin_Nuc_Density_violin_seaborn.png -------------
    click.echo("Generating nucleosome density violin plot...")

    def count_df(bed):
        # flank_bed.intersect(bed, c=True) returns a BED with count in 4th column
        df = flank_bed.intersect(bed, c=True).to_dataframe(names=['chrom', 'start', 'end', 'nuc_count'])
        # density: number per kb across the total flank width (2*flank_size)
        df['density'] = df['nuc_count'].replace(0, np.nan) / (2 * flank_size / 1000.0)
        return df

    density_data = []
    condition_data = []
    for name, bed in nuc_beds.items():
        df = count_df(bed)
        density_data.extend(df['density'].tolist())
        condition_data.extend([name] * len(df))

    violin_df = pd.DataFrame({'density': density_data, 'Condition': condition_data})

    plt.figure(figsize=(10, 6))
    palette = dict(zip(counts.keys(), colors))
    ax = sns.violinplot(x='Condition', y='density', data=violin_df,
                        inner='quartile', cut=0, scale='width', palette=palette)
    ax.set_ylabel('Nucleosome density\n(# nuc per kb)')
    ax.set_title(f'Per-origin nucleosome density (n={len(ori_bed)})')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(outpath / 'Origin_Nuc_Density_violin_seaborn.png', dpi=300)
    plt.close()

    # ----------------- TSS profiling using BigWig -----------------------
    click.echo("Generating TSS profiles...")
    # Ensure TSS column exists in gene_df
    if 'TSS' not in gene_df.columns:
        try:
            gene_df = gene_df.rename(columns={0: 'chrom', 1: 'start', 2: 'end', 3: 'strand'})
            gene_df['TSS'] = np.where(gene_df['strand'] == '+', gene_df['start'], gene_df['end'])
        except Exception:
            raise click.BadParameter("Gene annotation must contain chrom,start,end,strand in cols 0-3 for TSS profiling.")

    # Open BigWig files
    bw_objs = {name: pyBigWig.open(path) for name, path in bw_paths.items()}

    def profile_bw(bw, window, step=1):
        pos = np.arange(-window, window + 1, step)
        sum_cov = np.zeros_like(pos, dtype=float)
        count = np.zeros_like(pos, dtype=float)
        for _, row in gene_df.iterrows():
            chrom, tss = row['chrom'], int(row['TSS'])
            start = max(0, tss - window)
            end = tss + window
            vals = bw.values(chrom, start, end + 1, numpy=True)
            if vals is None or len(vals) < len(pos):
                continue
            if step != 1:
                vals = vals[::step]
            mask = ~np.isnan(vals)
            sum_cov[mask] += vals[mask]
            count[mask] += 1
        return pos, sum_cov / np.where(count == 0, np.nan, count)

    # TSS metaplot ±1kb, step 10
    window = 1000
    step = 10
    plt.figure(figsize=(10, 6))
    for i, (name, bw) in enumerate(bw_objs.items()):
        x, prof = profile_bw(bw, window, step)
        plt.plot(x, prof, label=name, color=colors[i], linewidth=2)
    plt.axvline(0, color='black', linestyle='--', alpha=0.7)
    plt.xlabel('Position relative to TSS (bp)')
    plt.ylabel('Mean coverage')
    plt.title('Nucleosome signal around TSS')
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath / f'Profile_TSS_±{window}bp.png', dpi=300)
    plt.close()

    # ----------------- Origin profiling using BigWig --------------------
    click.echo(f"Generating origin profiles for windows: {origin_windows_list}")

    def profile_origin_bw(bw, origin_df, window, step=10):
        pos = np.arange(-window, window + 1, step)
        sum_cov = np.zeros_like(pos, dtype=float)
        count = np.zeros_like(pos, dtype=float)
        for _, row in origin_df.iterrows():
            chrom = row['chrom']
            center = (int(row['start']) + int(row['end'])) // 2
            start = max(0, center - window)
            end = center + window
            vals = bw.values(chrom, start, end + 1, numpy=True)
            if vals is None or len(vals) < len(pos):
                continue
            if step != 1:
                vals = vals[::step]
            mask = ~np.isnan(vals)
            sum_cov[mask] += vals[mask]
            count[mask] += 1
        return pos, sum_cov / np.where(count == 0, np.nan, count)

    for w in origin_windows_list:
        plt.figure(figsize=(10, 6))
        for i, (name, bw) in enumerate(bw_objs.items()):
            x, prof = profile_origin_bw(bw, ori_bed, w)
            plt.plot(x, prof, label=name, color=colors[i], linewidth=2)
        plt.axvline(0, color='black', linestyle='--', alpha=0.7)
        plt.xlabel('Position relative to origin (bp)')
        plt.ylabel('Mean coverage')
        plt.title(f'Nucleosome signal around origins (±{w // 1000} kb)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(outpath / f'Profile_Origin_±{w // 1000}kb.png', dpi=300)
        plt.close()

    # Close BigWig handles
    for bw in bw_objs.values():
        bw.close()

    # Final messages
    click.echo(f"Analysis complete! Output files saved to: {outpath}")
    click.echo("Generated plots:")
    click.echo("  - NumberNucs_bar.png")
    click.echo("  - Origin_Nuc_Density_violin_seaborn.png")
    click.echo(f"  - Profile_TSS_±{window}bp.png")
    for w in origin_windows_list:
        click.echo(f"  - Profile_Origin_±{w // 1000}kb.png")


if __name__ == '__main__':
    main()