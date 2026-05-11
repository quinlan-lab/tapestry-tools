"""
Extract per-sample, per-haplotype alleles at each meQTL in
experiments/ASM-loci.meQTL.bed from a phased multi-sample VCF whose GTs are
phased as pat|mat (e.g. CEPH1463 joint call).

Output: a tidy TSV keyed by locus coordinates (so it joins directly into
the long methylation frame in correlate_methylation_with_haplotypes.py):

    chrom  start  end  meQTL_id  sample  meqtl_allele_pat  meqtl_allele_mat

meqtl_allele_pat/meqtl_allele_mat are the actual base strings (e.g. "C", "G"). Rows
where a sample has an unphased or missing GT at a site are dropped.

Run on a host with bcftools and the VCF visible.

Usage:
  .venv/bin/python experiments/extract_meqtl_alleles_from_vcf.py \\
      --vcf /scratch/.../CEPH1463.GRCh38.pass.sorted.vcf.gz \\
      --out experiments/meqtl_alleles.tsv
"""

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent
DEFAULT_BED = HERE / "ASM-loci.meQTL.bed"
DEFAULT_OUT = HERE / "meqtl_alleles.tsv"

# Cohort statuses we trust for downstream allele-vs-methylation analysis.
ACCEPT_STATUSES = {"OK_BIALLELIC", "OK_MULTIALLELIC_IN_VCF"}


def read_bed(path):
    header = None
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def query_site(vcf, chrom, pos):
    """Return (ref, alts_list, [(sample, gt), ...]) for the single record
    whose POS matches exactly. None if no such record."""
    region = f"{chrom}:{pos}-{pos}"
    fmt = "%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n"
    out = subprocess.run(
        ["bcftools", "query", "-r", region, "-f", fmt, vcf],
        check=True, capture_output=True, text=True,
    ).stdout
    for line in out.splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 5:
            continue
        if fields[1] != str(pos):
            continue
        ref = fields[2]
        alts = fields[3].split(",")
        sample_gts = []
        for sg in fields[4:]:
            if "=" not in sg:
                continue
            s, gt = sg.split("=", 1)
            sample_gts.append((s, gt))
        return ref, alts, sample_gts
    return None


def parse_phased_gt(gt, alleles):
    """alleles = [REF, ALT1, ALT2, ...]. Return (pat, mat) base strings, or
    None if GT is unphased or missing."""
    if "|" not in gt:
        return None
    a, b = gt.split("|", 1)
    if a == "." or b == ".":
        return None
    try:
        ai, bi = int(a), int(b)
    except ValueError:
        return None
    if ai >= len(alleles) or bi >= len(alleles):
        return None
    return alleles[ai], alleles[bi]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--bed", default=str(DEFAULT_BED))
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument(
        "--accept-status", action="append", default=None,
        help=f"cohort_status values to keep (default: {sorted(ACCEPT_STATUSES)})",
    )
    args = ap.parse_args()
    accept = set(args.accept_status) if args.accept_status else ACCEPT_STATUSES

    rows = read_bed(args.bed)
    out_path = Path(args.out)

    n_loci = n_loci_emitted = n_records = 0
    with out_path.open("w") as f:
        f.write("chrom\tstart\tend\tmeQTL_id\tsample\tmeqtl_allele_pat\tmeqtl_allele_mat\n")
        for r in rows:
            if not r.get("meqtl_pos"):
                continue
            status = r.get("cohort_status", "")
            if status not in accept:
                print(f"skipping {r.get('meQTL_id')} ({r['chrom']}:{r['meqtl_pos']}): "
                      f"status={status}", file=sys.stderr)
                continue
            n_loci += 1
            chrom = r["chrom"]
            pos = int(r["meqtl_pos"])
            result = query_site(args.vcf, chrom, pos)
            if result is None:
                print(f"  no VCF record at {chrom}:{pos}", file=sys.stderr)
                continue
            ref, alts, sample_gts = result
            alleles = [ref] + alts
            meqtl_id = r.get("meQTL_id", "")
            emitted_here = 0
            for sample, gt in sample_gts:
                parsed = parse_phased_gt(gt, alleles)
                if parsed is None:
                    continue
                ap_, am = parsed
                f.write(f"{chrom}\t{r['start']}\t{r['end']}\t"
                        f"{meqtl_id}\t{sample}\t{ap_}\t{am}\n")
                emitted_here += 1
            n_records += emitted_here
            if emitted_here:
                n_loci_emitted += 1
            print(f"  {meqtl_id} {chrom}:{pos}: {emitted_here} phased samples",
                  file=sys.stderr)

    print(f"wrote {out_path}: {n_records} rows across "
          f"{n_loci_emitted}/{n_loci} loci", file=sys.stderr)


if __name__ == "__main__":
    main()
