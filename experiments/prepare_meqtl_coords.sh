#!/usr/bin/env bash
# Driver: derive meQTL coordinates for ASM-loci.bed and verify them against
# the phased multi-sample VCF.
#
# Step 1 (enrich):  ASM-loci.bed  ->  ASM-loci.with-meqtl-coords.bed
#                   via Ensembl REST (GRCh38). Cached in .ensembl_cache.json.
# Step 2 (verify):  cross-check REF/ALT at each meqtl_pos against the VCF.
#
# Step 1 needs internet; step 2 needs bcftools and read access to the VCF
# (typically run on the cluster). Set VCF to override the default path.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
PY="$REPO/.venv/bin/python"

VCF="${VCF:-/scratch/ucgd/lustre-labs/quinlan/data-shared/haplotype-maps/CEPH1463.GRCh38/CEPH1463.GRCh38.pass.sorted.vcf.gz}"

echo "[1/2] enriching ASM-loci.bed with meQTL coords from Ensembl"
"$PY" "$HERE/enrich_asm_loci_with_meqtl_coords.py"

echo
echo "[2/2] verifying meQTL alleles against $VCF"
if ! command -v bcftools >/dev/null 2>&1; then
    echo "  bcftools not on PATH; skipping verification." >&2
    echo "  Re-run this script on a host with bcftools and the VCF visible." >&2
    exit 0
fi
if [[ ! -r "$VCF" ]]; then
    echo "  VCF not readable at $VCF; skipping verification." >&2
    echo "  Set VCF=<path> and re-run on a host that can see it." >&2
    exit 0
fi
"$PY" "$HERE/verify_meqtl_alleles_against_vcf.py" --vcf "$VCF"
