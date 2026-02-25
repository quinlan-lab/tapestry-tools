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

    # Combine Paternal and Maternal data
    final_plot_df = pl.concat(
        [
            get_parental_df(df_long, mode, "pat", base_cols).with_columns(
                pl.lit("Father").alias("parent")
            ),
            get_parental_df(df_long, mode, "mat", base_cols).with_columns(
                pl.lit("Mother").alias("parent")
            ),
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
                "founder": "Founder haplotype",
                "methylation": f"{mode.capitalize()}-based methylation",
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
                "methylation": f"{mode.capitalize()}-based methylation",
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
