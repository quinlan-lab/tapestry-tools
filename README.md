# tapestry-tools

Install `uv`: 
```
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Install `tapestry-tools`: 
```
uv tool install git+https://github.com/quinlan-lab/tapestry-tools.git
```

Check that `call-imprinted-loci` tool is available:
```
call-imprinted-loci
```

Call imprinted loci: 
```
call-imprinted-loci \
    --delta-meth-bed /scratch/ucgd/lustre-labs/quinlan/data-shared/tapestry-tools/delta_meth_all_samples.bed \
    --imprinted-bed imprinted_candidates.bed
```
