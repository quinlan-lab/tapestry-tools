from tqdm import tqdm 
from pathlib import Path 
import polars as pl 
import polars.selectors as cs

from .tile import get_tiles 
from .read_data import read_tapestry
from .methylation import compute_methylation
from .version_sort import version_sort
from .prefix_columns import prefix_columns
from .get_samples_and_paths import get_samples_and_paths

def compute_delta_methylation(df): 
    # List of the metric types you want to calculate a delta for
    metrics_to_diff = [
        "count_based_meth",
        "model_based_meth"
    ]
    delta_expressions = []
    for metric in metrics_to_diff:
        expr = (pl.col(f"{metric}_pat") - pl.col(f"{metric}_mat")).alias(f"delta_of_{metric}")
        delta_expressions.append(expr)
    df = df.with_columns(delta_expressions)

    df = df.drop([
        'count_based_meth',
        'model_based_meth'
    ])

    return df 

def compute_delta_methylation_all_samples(reference_genome, tile_size, sample_meth_beds, testing): 
    df_tiles = get_tiles(reference_genome, tile_size)

    sample_ids, meth_file_paths = get_samples_and_paths(sample_meth_beds, testing)

    df_all_samples = None
    for sample_id, bed_meth in tqdm(zip(sample_ids, meth_file_paths)):
        if Path(bed_meth).exists():
            df_meth = read_tapestry(bed_meth)
        else: 
            print(f"Could not read CpG sites at which count- and model-based methylation levels have been phased to founder haplotypes") 
            print(f"Required file does not exist: '{bed_meth}'")
            print(f"This may be because this sample is a founder and therefore cannot be inheritance-based phased")
            continue

        df_meth_free_from_allele_specific_cpgs = df_meth.filter(~pl.col('cpg_is_allele_specific'))
        df_tiles_with_meth = compute_methylation(df_tiles, df_meth_free_from_allele_specific_cpgs)        
        df_tiles_with_delta_meth = compute_delta_methylation(df_tiles_with_meth)
        join_keys = ['chrom', 'start', 'end']
        df_tiles_with_delta_meth = (
            df_tiles_with_delta_meth
            .select(join_keys + [
                'num_cpgs', 
                'num_cpgs_with_non_null_count_based_meth', 
                'num_cpgs_with_non_null_count_based_meth_pat', 
                'num_cpgs_with_non_null_count_based_meth_mat', 
                'delta_of_count_based_meth', 
                'delta_of_model_based_meth'
            ])
            .rename({
                'num_cpgs_with_non_null_count_based_meth': 'num_valid_cpgs',
                'num_cpgs_with_non_null_count_based_meth_pat': 'num_valid_cpgs_pat',
                'num_cpgs_with_non_null_count_based_meth_mat': 'num_valid_cpgs_mat',
            })
        )
        df_tiles_with_delta_meth = prefix_columns(df_tiles_with_delta_meth, prefix=sample_id, join_keys=join_keys)

        if df_all_samples is None:
            df_all_samples = df_tiles_with_delta_meth
        else:
            df_all_samples = df_all_samples.join(
                df_tiles_with_delta_meth,
                on=join_keys,
                # capture tiles in which at least one of the two dfs has a record: 
                how="full", 
                coalesce=True, 
            )   
    return version_sort(df_all_samples)

def format_locus(df):
    return df.with_columns(
        locus = pl.format(
            "{}:{}-{}", 
            pl.col("chrom"),
            pl.col("start"),
            pl.col("end")
        )
    )

def call_imprinted_loci(df, meth_mode, delta_meth_threshold, num_valid_cpgs_per_hap_threshold, valid_cpg_ratio_threshold):
    delta_meth_cols = [col for col in df.columns if col.endswith(f"_{meth_mode}")]

    # List to store expressions that return the 'prefix' (sample name) if 'condition' (see below) is met
    sample_match_exprs = []
    
    for delta_meth_col in delta_meth_cols:
        prefix = delta_meth_col.removesuffix(f"_{meth_mode}")
        num_cpgs_col = f"{prefix}_num_cpgs"
        num_valid_cpgs_pat_col = f"{prefix}_num_valid_cpgs_pat"
        num_valid_cpgs_mat_col = f"{prefix}_num_valid_cpgs_mat"

        condition = (
            (pl.col(delta_meth_col).abs() > delta_meth_threshold) & 
            (pl.col(num_valid_cpgs_pat_col) >= num_valid_cpgs_per_hap_threshold) & 
            (pl.col(num_valid_cpgs_mat_col) >= num_valid_cpgs_per_hap_threshold) &
            ((pl.col(num_valid_cpgs_pat_col) / pl.col(num_cpgs_col)) > valid_cpg_ratio_threshold) &
            ((pl.col(num_valid_cpgs_mat_col) / pl.col(num_cpgs_col)) > valid_cpg_ratio_threshold)
        )

        # If 'condition' is True, return the sample name (prefix); otherwise return null
        sample_match_exprs.append(
            pl.when(condition).then(pl.lit(prefix)).otherwise(None)
        )

    df= format_locus(
        df
        .with_columns(
            # Combine all sample checks into one list column and remove non-matches (nulls)
            pl.concat_list(sample_match_exprs).list
            .drop_nulls()
            .alias("imprinted_samples")
        )
        # Filter to keep only rows where at least one sample matched
        .filter(pl.col("imprinted_samples").list.len() > 0)
        .with_columns(
            pl.col("imprinted_samples").list.len().alias("num_imprinted_samples")
        )
        .select(
            'chrom',
            'start', 
            'end',
            'imprinted_samples',
            'num_imprinted_samples',
            # cs.contains("NA12885") # TESTING 
        )
    )

    print(f"Number of candidate imprinted loci: {len(df)}")
    return df
    
