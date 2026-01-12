# tapestry-tools

How to use this tool
You can run this from your terminal. It gives you full control over the filtering thresholds defined in your imprinting.py logic.

Basic Usage:

Bash

python call_imprinted_loci.py all_samples.delta_meth.bed candidates.tsv
Customizing Thresholds: If you want to be stricter about what counts as an imprinted locus (e.g., requiring a higher methylation difference or more CpGs), you can pass arguments:

Bash

python call_imprinted_loci.py \
    input_data.bed \
    strict_candidates.tsv \
    --meth-mode count \
    --delta-threshold 0.75 \
    --min-cpgs-per-hap 10 \
    --min-cpg-ratio 0.8

Wrangling the output
Since the output is a standard TSV (Tab Separated Values), you can immediately interrogate it using other tools as requested:

Bash: cut -f 1-3 candidates.tsv | head (Check coordinates)

R: df <- read.delim("candidates.tsv")

Python: pl.read_csv("candidates.tsv", separator='\t')
