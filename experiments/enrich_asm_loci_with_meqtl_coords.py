"""
Enrich experiments/ASM-loci.bed with meQTL SNP coordinates (GRCh38) from
the Ensembl REST API.

For every row whose SNP column is an rsID, look it up at rest.ensembl.org
(GRCh38) and write meqtl_pos / meqtl_ref / meqtl_alt columns. Rows without
an rsID (e.g. "deletion", "unknown") are preserved with empty values.

Output: experiments/ASM-loci.with-meqtl-coords.bed
Cache:  experiments/.ensembl_cache.json  (re-runs are free)
"""

import json
import re
import sys
import time
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

HERE = Path(__file__).parent
IN_BED = HERE / "ASM-loci.bed"
OUT_BED = HERE / "ASM-loci.with-meqtl-coords.bed"
CACHE_FILE = HERE / ".ensembl_cache.json"

# Default REST endpoint serves GRCh38. For GRCh37 use grch37.rest.ensembl.org.
ENSEMBL_BASE = "https://rest.ensembl.org"
ASSEMBLY = "GRCh38"
MIN_INTERVAL_S = 0.1  # Ensembl asks for <=15 req/s; this gives 10/s headroom.

RSID_RE = re.compile(r"^rs\d+$")


def load_cache():
    if CACHE_FILE.exists():
        return json.loads(CACHE_FILE.read_text())
    return {}


def save_cache(cache):
    CACHE_FILE.write_text(json.dumps(cache, indent=2, sort_keys=True))


def fetch_variation(rsid, cache):
    if rsid in cache:
        return cache[rsid]
    url = f"{ENSEMBL_BASE}/variation/human/{rsid}?content-type=application/json"
    req = Request(url, headers={"Accept": "application/json"})
    for attempt in range(5):
        try:
            with urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read())
            cache[rsid] = data
            save_cache(cache)
            return data
        except HTTPError as e:
            if e.code == 429:
                wait = 2 ** attempt
                print(f"  429 from Ensembl, backing off {wait}s", file=sys.stderr)
                time.sleep(wait)
                continue
            if e.code == 404:
                cache[rsid] = None
                save_cache(cache)
                return None
            raise
        except URLError as e:
            wait = 2 ** attempt
            print(f"  network error {e}; retry in {wait}s", file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"failed to fetch {rsid} after retries")


def pick_mapping(variation, chrom):
    """Return (pos, ref, alt, note) for the GRCh38 mapping on `chrom`.

    `chrom` is the bed-style name (e.g. "chr2"); Ensembl returns bare names
    ("2"), so strip the "chr" prefix when comparing.
    """
    want = chrom[3:] if chrom.startswith("chr") else chrom
    candidates = [
        m for m in variation.get("mappings", [])
        if m.get("assembly_name") == ASSEMBLY and m.get("seq_region_name") == want
    ]
    if not candidates:
        return None, None, None, f"no {ASSEMBLY} mapping on {chrom}"

    note_bits = []
    if variation.get("merged"):
        note_bits.append(f"merged_into={variation.get('name')}")

    snv = next((m for m in candidates if m.get("start") == m.get("end")), None)
    m = snv or candidates[0]
    pos = m.get("start")
    allele_string = m.get("allele_string", "")  # e.g. "A/G" or "A/G/T"
    parts = allele_string.split("/") if allele_string else []
    if len(parts) < 2:
        return pos, "", "", "unparseable allele_string=" + allele_string
    ref = parts[0]
    alts = parts[1:]
    if len(alts) > 1:
        note_bits.append("multiallelic")
    alt = ",".join(alts)
    return pos, ref, alt, ";".join(note_bits)


def main():
    cache = load_cache()
    last_call = 0.0

    with IN_BED.open() as f:
        lines = f.readlines()

    out_lines = []
    for line in lines:
        if line.startswith("##"):
            out_lines.append(line)
            continue
        if line.startswith("#"):
            header = line.rstrip("\n").split("\t")
            header += ["meqtl_pos", "meqtl_ref", "meqtl_alt", "meqtl_note"]
            out_lines.append("\t".join(header) + "\n")
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 5:
            out_lines.append(line)
            continue
        chrom, start, end, snp, gene = fields[:5]

        meqtl_pos = ref = alt = note = ""
        if RSID_RE.match(snp):
            elapsed = time.time() - last_call
            if elapsed < MIN_INTERVAL_S:
                time.sleep(MIN_INTERVAL_S - elapsed)
            print(f"looking up {snp} ({chrom})", file=sys.stderr)
            variation = fetch_variation(snp, cache)
            last_call = time.time()
            if variation is None:
                note = "not_found_in_ensembl"
            else:
                pos, r, a, n = pick_mapping(variation, chrom)
                if pos is not None:
                    meqtl_pos = str(pos)
                ref, alt, note = r or "", a or "", n or ""
        else:
            note = "no_rsid"

        out_lines.append(
            "\t".join([chrom, start, end, snp, gene, meqtl_pos, ref, alt, note]) + "\n"
        )

    OUT_BED.write_text("".join(out_lines))
    print(f"wrote {OUT_BED}", file=sys.stderr)


if __name__ == "__main__":
    main()
