# tapestry-tools

Install `uv`: 
```
curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env
```

Additionally, prior to installation, you may need to specify a directory for `TMPDIR` that you have write access to. 

Install `tapestry-tools`: 
```
uv tool install git+https://github.com/quinlan-lab/tapestry-tools.git
```

Check that `call-imprinted-loci` tool is available:
```
call-imprinted-loci
```

Given a set of methylation differences between haplotypes for a set of genomic tiles and a set of samples, call imprinted loci: 
```
call-imprinted-loci \
    --delta-meth-bed /scratch/ucgd/lustre-labs/quinlan/data-shared/tapestry-tools/delta_meth_all_samples.bed \
    --imprinted-bed imprinted_candidates.bed
```
