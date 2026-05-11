"""
Verify that the Ensembl-derived meQTL coordinates in
experiments/ASM-loci.with-meqtl-coords.bed agree with a phased multi-sample
VCF (typically the CEPH1463 joint call on the cluster).

For each row with an rsID, query the VCF at meqtl_pos and report:
  - whether REF matches Ensembl
  - which ALT(s) are actually segregating in the cohort
  - whether the site is effectively bi-allelic

Optionally writes an enriched bed file (--write-out) that adds the
cohort-observed ALT and a per-site status, and renames the SNP header to
meQTL_id. That file is the input downstream code should rely on.

Run this on a host with bcftools on PATH and read access to the VCF.

Usage:
  .venv/bin/python experiments/verify_meqtl_alleles_against_vcf.py \\
      --vcf /scratch/.../CEPH1463.GRCh38.pass.sorted.vcf.gz \\
      --write-out experiments/ASM-loci.meQTL.bed
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
        if fields[1] != str(pos):
            continue
        records.append((fields[3], fields[4]))
    return records


def classify(ens_ref, ens_alts, vcf_ref, vcf_alt):
    ref_match = vcf_ref == ens_ref
    ens_alt_set = set(ens_alts.split(",")) if ens_alts else set()
    vcf_alt_set = set(vcf_alt.split(","))
    alt_subset = vcf_alt_set.issubset(ens_alt_set) if ens_alt_set else False
    biallelic = len(vcf_alt_set) == 1
    if ref_match and alt_subset and biallelic:
        return "OK_BIALLELIC"
    if ref_match and alt_subset:
        return "OK_MULTIALLELIC_IN_VCF"
    if ref_match:
        return "ALT_NOT_IN_ENSEMBL"
    return "REF_MISMATCH"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="path to phased multi-sample VCF")
    ap.add_argument("--bed", default=str(DEFAULT_BED), help="enriched ASM loci bed")
    ap.add_argument(
        "--write-out",
        help="if set, write a new bed with cohort_alt + cohort_status columns "
             "and SNP renamed to meQTL_id (e.g. experiments/ASM-loci.meQTL.bed)",
    )
    args = ap.parse_args()

    # Parse input bed, preserving ## lines and header.
    pre_header_lines = []
    header = None
    rows = []
    with open(args.bed) as f:
        for line in f:
            if line.startswith("##"):
                pre_header_lines.append(line)
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))

    print(f"{'meQTL_id':<14}{'chrom':<7}{'pos':>12}  "
          f"{'ens_ref':<8}{'ens_alts':<10}{'vcf_ref':<8}{'vcf_alt':<10}status")
    print("-" * 90)

    n_checked = n_ok = 0
    for r in rows:
        r["cohort_alt"] = ""
        r["cohort_status"] = ""
        if not r.get("meqtl_pos"):
            r["cohort_status"] = "no_meqtl_pos"
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
            r["cohort_status"] = "bcftools_error"
            continue

        if not records:
            status = "MISSING"
            vcf_ref = vcf_alt = "-"
        elif len(records) > 1:
            status = f"MULTIPLE_RECORDS({len(records)})"
            vcf_ref, vcf_alt = records[0]
        else:
            vcf_ref, vcf_alt = records[0]
            status = classify(ens_ref, ens_alts, vcf_ref, vcf_alt)
            if status == "OK_BIALLELIC":
                n_ok += 1

        r["cohort_alt"] = vcf_alt if vcf_alt != "-" else ""
        r["cohort_status"] = status

        print(f"{r['SNP']:<14}{chrom:<7}{pos:>12}  "
              f"{ens_ref:<8}{ens_alts:<10}{vcf_ref:<8}{vcf_alt:<10}{status}")

    print("-" * 90)
    print(f"{n_ok}/{n_checked} sites are bi-allelic with matching REF and ALT⊆Ensembl")

    if args.write_out:
        out_path = Path(args.write_out)
        out_header = [("meQTL_id" if c == "SNP" else c) for c in header]
        out_header += ["cohort_alt", "cohort_status"]
        with out_path.open("w") as f:
            for line in pre_header_lines:
                f.write(line)
            f.write("#" + "\t".join(out_header) + "\n")
            for r in rows:
                f.write("\t".join(r.get(c if c != "meQTL_id" else "SNP", "")
                                  for c in out_header) + "\n")
        print(f"wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
