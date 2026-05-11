"""
Verify that the Ensembl-derived meQTL coordinates in
experiments/ASM-loci.with-meqtl-coords.bed agree with a phased multi-sample
VCF (typically the CEPH1463 joint call on the cluster).

For each row with an rsID, query the VCF at meqtl_pos and report:
  - whether REF matches Ensembl
  - which ALT(s) are actually segregating in the cohort
  - whether the site is effectively bi-allelic

Run this on the host that can see the VCF (the cluster). Requires bcftools
on PATH.

Usage:
  .venv/bin/python experiments/verify_meqtl_alleles_against_vcf.py \\
      --vcf /scratch/.../CEPH1463.GRCh38.pass.sorted.vcf.gz
"""

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent
DEFAULT_BED = HERE / "ASM-loci.with-meqtl-coords.bed"


def query_site(vcf, chrom, pos):
    """Return list of (ref, alt) tuples at chrom:pos (one per VCF record)."""
    region = f"{chrom}:{pos}-{pos}"
    out = subprocess.run(
        ["bcftools", "view", "-H", "-r", region, vcf],
        check=True, capture_output=True, text=True,
    ).stdout
    records = []
    for line in out.splitlines():
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        # Only keep the record whose POS matches exactly (region can spill).
        if fields[1] != str(pos):
            continue
        records.append((fields[3], fields[4]))
    return records


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="path to phased multi-sample VCF")
    ap.add_argument("--bed", default=str(DEFAULT_BED), help="enriched ASM loci bed")
    args = ap.parse_args()

    rows = []
    with open(args.bed) as f:
        header = None
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))

    print(f"{'SNP':<14}{'chrom':<7}{'pos':>12}  "
          f"{'ens_ref':<8}{'ens_alts':<10}{'vcf_ref':<8}{'vcf_alt':<10}status")
    print("-" * 90)

    n_checked = n_ok = 0
    for r in rows:
        if not r.get("meqtl_pos"):
            continue
        n_checked += 1
        chrom = r["chrom"]
        pos = int(r["meqtl_pos"])
        ens_ref = r["meqtl_ref"]
        ens_alts = r["meqtl_alt"]
        try:
            records = query_site(args.vcf, chrom, pos)
        except subprocess.CalledProcessError as e:
            print(f"  bcftools failed for {chrom}:{pos}: {e.stderr}", file=sys.stderr)
            continue

        if not records:
            status = "MISSING"
            vcf_ref = vcf_alt = "-"
        elif len(records) > 1:
            status = f"MULTIPLE_RECORDS({len(records)})"
            vcf_ref, vcf_alt = records[0]
        else:
            vcf_ref, vcf_alt = records[0]
            ref_match = vcf_ref == ens_ref
            ens_alt_set = set(ens_alts.split(",")) if ens_alts else set()
            vcf_alt_set = set(vcf_alt.split(","))
            alt_subset = vcf_alt_set.issubset(ens_alt_set) if ens_alt_set else False
            biallelic = len(vcf_alt_set) == 1
            if ref_match and alt_subset and biallelic:
                status = "OK_BIALLELIC"
                n_ok += 1
            elif ref_match and alt_subset:
                status = "OK_MULTIALLELIC_IN_VCF"
            elif ref_match:
                status = "ALT_NOT_IN_ENSEMBL"
            else:
                status = "REF_MISMATCH"

        print(f"{r['SNP']:<14}{chrom:<7}{pos:>12}  "
              f"{ens_ref:<8}{ens_alts:<10}{vcf_ref:<8}{vcf_alt:<10}{status}")

    print("-" * 90)
    print(f"{n_ok}/{n_checked} sites are bi-allelic with matching REF and ALT⊆Ensembl")


if __name__ == "__main__":
    main()
