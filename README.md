# tapestry-tools

## Installation 

Install `uv`: 
```
curl -LsSf https://astral.sh/uv/install.sh | sh
source $HOME/.local/bin/env
```

Note: prior to running the `uv` installation commands above, you may need to specify a directory for `TMPDIR` that you have write access to. 

Install `tapestry-tools`: 
```
uv tool install git+https://github.com/quinlan-lab/tapestry-tools.git
```

Check that `compute-delta-methylation` tool is available:
```
compute-delta-methylation
```

## Imprinting 

Test the computation of methylation difference between haplotypes for multiple samples:
```
# This will run a 2min test (using just 2 samples)
compute-delta-methylation \
    --sample-meth-beds /scratch/ucgd/lustre-labs/quinlan/data-shared/tapestry-tools/CEPH1463.GRCh38.hifi.founder-phased/all-cpgs.all-samples.tsv \
    --delta-meth-bed /scratch/ucgd/lustre-labs/quinlan/data-shared/tapestry-tools/CEPH1463.GRCh38.hifi.founder-phased/delta-meth.all-samples.test.bed \
    --testing
```

Compute methylation differences using all samples in the pedigree: 
```
# 25 mins to run for all samples
nohup compute-delta-methylation \
    --sample-meth-beds /scratch/ucgd/lustre-labs/quinlan/data-shared/tapestry-tools/CEPH1463.GRCh38.hifi.founder-phased/all-cpgs.all-samples.tsv \
    --delta-meth-bed /scratch/ucgd/lustre-labs/quinlan/data-shared/tapestry-tools/CEPH1463.GRCh38.hifi.founder-phased/delta-meth.all-samples.bed \
    > /scratch/ucgd/lustre-labs/quinlan/data-shared/tapestry-tools/CEPH1463.GRCh38.hifi.founder-phased/delta-meth.all-samples.log 2>&1 & 
```

Given the output of `compute-delta-methylation` (a set of methylation differences between haplotypes for a set of genomic tiles and a set of samples), call imprinted loci: 
```
call-imprinted-loci \
    --delta-meth-bed /scratch/ucgd/lustre-labs/quinlan/data-shared/tapestry-tools/CEPH1463.GRCh38.hifi.founder-phased/delta-meth.all-samples.bed \
    --imprinted-bed /scratch/ucgd/lustre-labs/quinlan/data-shared/tapestry-tools/CEPH1463.GRCh38.hifi.founder-phased/imprinted-candidates.all-samples.bed
```

Given a set of loci (e.g., imprinted loci), and a set of samples, compute founder-phased methylation averaged over CpGs in each locus: 
```
compute-methylation-loci \
    --loci-bed /scratch/ucgd/lustre-labs/quinlan/data-shared/tapestry-tools/CEPH1463.GRCh38.hifi.founder-phased/imprinted-candidates.all-samples.bed \
    --sample-meth-beds /scratch/ucgd/lustre-labs/quinlan/data-shared/tapestry-tools/CEPH1463.GRCh38.hifi.founder-phased/all-cpgs.all-samples.tsv \
    --loci-meth-bed /scratch/ucgd/lustre-labs/quinlan/data-shared/tapestry-tools/CEPH1463.GRCh38.hifi.founder-phased/imprinted-candidates.all-samples.meth.bed
```