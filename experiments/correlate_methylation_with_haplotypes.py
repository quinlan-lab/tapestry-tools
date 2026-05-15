import polars as pl
import matplotlib.pyplot as plt
import os
import plotly.express as px
import plotly.io as pio

# Set the font family to Arial
# https://g.co/gemini/share/3898a74b2d77
FONT_FAMILY = "Arial"

FONT_SIZE = 35

plt.rcParams["font.family"] = FONT_FAMILY
plt.rcParams["font.sans-serif"] = [FONT_FAMILY]

plt.rcParams.update({"font.size": FONT_SIZE})

pio.templates["custom"] = pio.templates["plotly"]
pio.templates["custom"].layout.font.family = FONT_FAMILY  # type:ignore
pio.templates["custom"].layout.font.size = FONT_SIZE  # type:ignore
pio.templates.default = "custom"

# Check if the environment variable in .env was actually set in the shell
print(f"PYTHONPATH: {os.environ.get('PYTHONPATH')}")

pl.Config.set_tbl_rows(10)


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_FIG_DIR = os.path.normpath(
    os.path.join(REPO_ROOT, "..", "tapestry", "manuscript", "fig5")
)


def save_methylation_vs_founder_by_allele(fig, snp_id, outdir=DEFAULT_FIG_DIR):
    """Save a methylation-vs-founder plot (colored by meQTL allele) to disk,
    with a filename that includes the rsID."""
    os.makedirs(outdir, exist_ok=True)
    safe_id = (snp_id or "unknown_rsID").replace("/", "_").replace(":", "_")
    path = os.path.join(outdir, f"methylation_vs_founder.by_allele.{safe_id}.pdf")
    fig.write_image(path, format="pdf")
    print(f"  saved {path}")
    return path


def get_parental_df(df_long, mode, parent_type, base_cols):
    meth_col_suffix = f"{mode}_based_meth"
    suffix = f"_{parent_type}"
    target_meth = f"{meth_col_suffix}{suffix}"
    target_founder = f"founder{suffix}"

    return (
        df_long.filter(pl.col("metric").is_in([target_meth, target_founder]))
        .pivot(index=base_cols + ["sample"], on="metric", values="value")
        .select(
            [
                pl.col(base_cols + ["sample"]),
                pl.col(target_founder).alias("founder"),
                pl.col(target_meth).alias("methylation").cast(pl.Float64, strict=False),
            ]
        )
    )


def correlate_methylation_with_haplotypes(mode="count"):
    """
    Plot methylation vs founder haplotype (coloring by parent)
    and methylation vs parent (coloring by founder haplotype)
    for each locus.

    Args:
        mode: "count" or "model" for methylation type
    """
    # Load Data
    df = pl.read_csv(
        "ASM-loci.meth.bed",
        separator="\t",
        comment_prefix="##",
        null_values=["", "null"],
    ).rename({"#chrom": "chrom"})

    # Add SNP column from original ASM-loci.bed
    SNPs = (
        pl.read_csv(
            "ASM-loci.bed",
            separator="\t",
            comment_prefix="##",
            null_values=["", "null"],
        )
        .rename({"#chrom": "chrom"})
        .select(["chrom", "start", "end", "SNP"])
    )
    df = df.join(SNPs, on=["chrom", "start", "end"], how="left")

    # Add gene column from original ASM-loci.bed
    genes = (
        pl.read_csv(
            "ASM-loci.bed",
            separator="\t",
            comment_prefix="##",
            null_values=["", "null"],
        )
        .rename({"#chrom": "chrom"})
        .select(["chrom", "start", "end", "gene"])
    )
    df = df.join(genes, on=["chrom", "start", "end"], how="left")

    base_cols = ["chrom", "start", "end", "SNP", "gene"]
    df_long = df.unpivot(
        index=base_cols, variable_name="raw_column", value_name="value"
    )

    # Parse Sample and Metric
    df_long = df_long.with_columns(
        [
            pl.col("raw_column").str.extract(r"^([^_]+)_", 1).alias("sample"),
            pl.col("raw_column").str.extract(r"^[^_]+_(.*)$", 1).alias("metric"),
        ]
    )

    # Tidy up
    df_long = df_long.drop("raw_column")

    # Per-sample, per-haplotype allele at each meQTL (produced by
    # extract_meqtl_alleles_from_vcf.py). Optional: if the file is missing,
    # skip the meQTL-allele plot rather than failing.
    alleles_path = "meqtl_alleles.tsv"
    if os.path.exists(alleles_path):
        alleles = pl.read_csv(alleles_path, separator="\t")
    else:
        alleles = None
        print(f"  {alleles_path} not found; fig3 (color by meQTL allele) will be skipped")

    def attach_allele(df, hap_col):
        if alleles is None:
            return df.with_columns(pl.lit(None).cast(pl.Utf8).alias("allele"))
        return df.join(
            alleles.select(
                ["chrom", "start", "end", "sample", pl.col(hap_col).alias("allele")]
            ),
            on=["chrom", "start", "end", "sample"],
            how="left",
        )

    # Combine Paternal and Maternal data
    final_plot_df = pl.concat(
        [
            attach_allele(
                get_parental_df(df_long, mode, "pat", base_cols), "meqtl_allele_pat"
            ).with_columns(pl.lit("Father").alias("parent")),
            attach_allele(
                get_parental_df(df_long, mode, "mat", base_cols), "meqtl_allele_mat"
            ).with_columns(pl.lit("Mother").alias("parent")),
        ]
    ).drop_nulls(subset=["methylation", "founder"])

    # Define consistent colors for haplotypes (alphabetically sorted)
    haplotypes = sorted(final_plot_df["founder"].unique())
    color_sequence = (
        px.colors.qualitative.Set2
        + px.colors.qualitative.Set1
        + px.colors.qualitative.Dark2
    )
    haplotype_color_map = {
        hap: color_sequence[i % len(color_sequence)] for i, hap in enumerate(haplotypes)
    }

    # Plotting Loop
    for row in df.iter_rows(named=True):
        snp_id = row.get("SNP", "")
        gene = row.get("gene", "")
        locus_name = f"{row['chrom']}:{row['start']}-{row['end']}"
        title = f"Locus to compute methylation: {locus_name}"
        if snp_id:
            title += f"<br>meQTL: {snp_id}"
        if gene:
            title += f"<br>Gene: {gene}"

        print(title.replace("<br>", "\n"))

        # Filter for the specific locus
        locus_df = final_plot_df.filter(
            (pl.col("chrom") == row["chrom"]) & (pl.col("start") == row["start"])
        ).to_pandas()

        if locus_df.empty:
            print(f"No valid data for locus {locus_name}")
            continue

        # Sort by parent (Father, Mother) then by founder for consistent ordering
        locus_df = locus_df.sort_values(["parent", "founder"])

        # Plot 1: Methylation vs Founder (color by Parent)
        fig1 = px.strip(
            locus_df,
            x="founder",
            y="methylation",
            color="parent",
            hover_data=["sample"],
            title=f"{title}",
            labels={
                "founder": "Haplotype label",
                "methylation": "Haplotype methylation",
            },
            category_orders={"founder": sorted(locus_df["founder"].unique())},
        )

        fig1.update_traces(marker_size=25)

        fig1.update_layout(
            yaxis_range=[0, 1],
            legend_title="Parent of origin",
            bargap=0.1,
            width=1000,
            height=800,
            title_y=0.95,
            title_font_size=30,
            margin=dict(t=150),
            plot_bgcolor="white",
            paper_bgcolor="white",
            xaxis=dict(linecolor="black", gridcolor="lightgray"),
            yaxis=dict(linecolor="black", gridcolor="lightgray"),
        )
        fig1.show()

        # Plot 2: Methylation vs Parent (color by Founder haplotype)
        fig2 = px.strip(
            locus_df,
            x="parent",
            y="methylation",
            color="founder",
            hover_data=["sample"],
            title=f"{title}",
            labels={
                "parent": "Parent of origin",
                "methylation": "Haplotype methylation",
                "founder": "Founder haplotype",
            },
            category_orders={"parent": ["Father", "Mother"], "founder": haplotypes},
            color_discrete_map=haplotype_color_map,
        )

        fig2.update_traces(marker_size=25, jitter=0.3, pointpos=0)

        fig2.update_layout(
            yaxis_range=[0, 1],
            legend_title="Founder haplotype",
            bargap=0.1,
            width=1000,
            height=800,
            title_y=0.95,
            title_font_size=30,
            margin=dict(t=150),
            plot_bgcolor="white",
            paper_bgcolor="white",
            xaxis=dict(linecolor="black", gridcolor="lightgray"),
            yaxis=dict(linecolor="black", gridcolor="lightgray"),
        )
        fig2.show()

        # Plot 3: Methylation vs Founder (color by allele at meQTL)
        if alleles is not None and locus_df["allele"].notna().any():
            allele_df = locus_df.dropna(subset=["allele"])
            # Assign red to the higher-methylation allele and blue to the
            # lower-methylation allele, matching the methylated/unmethylated
            # colors in Fig 1A/B. Higher-meth allele is listed first in the
            # legend.
            allele_means = (
                allele_df.groupby("allele")["methylation"].mean().sort_values(ascending=False)
            )
            allele_order = allele_means.index.tolist()
            allele_color_map = {
                allele_order[0]: "#dc2626",  # red — higher methylation
            }
            if len(allele_order) > 1:
                allele_color_map[allele_order[1]] = "#1d4ed8"  # blue — lower methylation
            fig3 = px.strip(
                allele_df,
                x="founder",
                y="methylation",
                color="allele",
                hover_data=["sample", "parent"],
                labels={
                    "founder": "Haplotype label",
                    "methylation": "Haplotype methylation",
                    "allele": "meQTL allele",
                },
                category_orders={
                    "founder": sorted(allele_df["founder"].unique()),
                    "allele": allele_order,
                },
                color_discrete_map=allele_color_map,
            )
            fig3.update_traces(marker_size=25, jitter=0.3, pointpos=0)
            fig3.update_layout(
                yaxis_range=[0, 1],
                legend_title="meQTL allele",
                bargap=0.1,
                width=1000,
                height=800,
                margin=dict(t=20),
                plot_bgcolor="white",
                paper_bgcolor="white",
                xaxis=dict(linecolor="black", gridcolor="lightgray"),
                yaxis=dict(linecolor="black", gridcolor="lightgray"),
            )
            fig3.show()
            save_methylation_vs_founder_by_allele(fig3, snp_id)
