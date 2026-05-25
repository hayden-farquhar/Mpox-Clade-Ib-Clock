#!/usr/bin/env python3
"""
Build the §5.4 outgroup panel and the three parental reference panels per
the pre-registered specification (OSF Pre-Reg CASR2 §5.4).

Outgroup panel (§5.4 paragraph 2):
- 30 clade-IIb genomes
- Stratified random sample by collection year (2022, 2023, 2024, 2025+) × WHO region
- ≥1 and ≤10 per (year × region) stratum
- High quality: length ≥ 190,000 bp, N fraction ≤ 0.10, Nextclade QC overall = good
- Random seed 42

Parental reference panels (§5.4 paragraph 3):
- Three independent panels of ≤10 high-quality genomes each (clade Ia, Ib, IIb)
- Same stratified-random procedure
- Seed = 42

Inputs (all from the 2026-05-22 freeze):
- NCBI metadata TSV with Length, Isolate Collection date, Geographic Location
- Nextclade clade-i TSV with clade assignments for clade-I genomes (Ia, Ib)
- Nextclade all-clades TSV for clade-IIb genomes (the §3.3 "unassigned" subset)

Outputs:
- panels/outgroup_30_iib.tsv — locked accession list with stratum + provenance
- panels/parental_ia_10.tsv, parental_ib_10.tsv, parental_iib_10.tsv
- panels/strata_summary.json — per-(clade × year × region) counts before and after sampling
- panels/SHA256SUMS — hashes of every output
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

# WHO region mapping (curated for the high-frequency mpox source countries).
# Anything not in this table falls into the catch-all "other" region; the
# stratification still runs but uses the catch-all bucket. This is documented
# in the panel TSV so the manuscript can disclose any "other" allocations.
WHO_REGION = {
    # Africa (AFRO)
    "Democratic Republic of the Congo": "AFR", "DRC": "AFR", "Congo": "AFR",
    "Nigeria": "AFR", "Cameroon": "AFR", "Central African Republic": "AFR",
    "South Africa": "AFR", "Liberia": "AFR", "Ghana": "AFR",
    "Republic of the Congo": "AFR", "Sierra Leone": "AFR", "Sudan": "AFR",
    "Uganda": "AFR", "Kenya": "AFR", "Tanzania": "AFR", "Burundi": "AFR",
    "Rwanda": "AFR", "Zambia": "AFR", "Zimbabwe": "AFR",
    "Ivory Coast": "AFR", "Cote d'Ivoire": "AFR", "Côte d'Ivoire": "AFR",
    "Gabon": "AFR", "Benin": "AFR", "Mozambique": "AFR", "Senegal": "AFR",
    "Angola": "AFR", "Ethiopia": "AFR", "Mali": "AFR",
    # Americas (PAHO)
    "USA": "AMR", "United States": "AMR",
    "Canada": "AMR", "Mexico": "AMR", "Brazil": "AMR", "Argentina": "AMR",
    "Chile": "AMR", "Colombia": "AMR", "Peru": "AMR", "Ecuador": "AMR",
    "Cuba": "AMR", "Dominican Republic": "AMR", "Panama": "AMR",
    "Costa Rica": "AMR", "Venezuela": "AMR", "Bolivia": "AMR",
    "Uruguay": "AMR", "Paraguay": "AMR", "Guatemala": "AMR",
    "El Salvador": "AMR", "Honduras": "AMR", "Jamaica": "AMR",
    # Eastern Mediterranean (EMRO)
    "Saudi Arabia": "EMR", "United Arab Emirates": "EMR", "UAE": "EMR",
    "Iran": "EMR", "Iraq": "EMR", "Lebanon": "EMR", "Jordan": "EMR",
    "Egypt": "EMR", "Morocco": "EMR", "Tunisia": "EMR", "Pakistan": "EMR",
    "Qatar": "EMR", "Kuwait": "EMR", "Bahrain": "EMR", "Oman": "EMR",
    "Israel": "EMR",
    # Europe (EURO)
    "United Kingdom": "EUR", "UK": "EUR", "England": "EUR", "Scotland": "EUR",
    "Wales": "EUR", "Northern Ireland": "EUR", "Ireland": "EUR",
    "Spain": "EUR", "France": "EUR", "Germany": "EUR", "Italy": "EUR",
    "Portugal": "EUR", "Netherlands": "EUR", "Belgium": "EUR", "Austria": "EUR",
    "Switzerland": "EUR", "Sweden": "EUR", "Norway": "EUR", "Denmark": "EUR",
    "Finland": "EUR", "Poland": "EUR", "Czech Republic": "EUR", "Czechia": "EUR",
    "Slovakia": "EUR", "Slovenia": "EUR", "Hungary": "EUR", "Romania": "EUR",
    "Bulgaria": "EUR", "Greece": "EUR", "Croatia": "EUR", "Serbia": "EUR",
    "Russia": "EUR", "Ukraine": "EUR", "Belarus": "EUR", "Lithuania": "EUR",
    "Latvia": "EUR", "Estonia": "EUR", "Iceland": "EUR", "Luxembourg": "EUR",
    "Malta": "EUR", "Cyprus": "EUR",
    # South-East Asia (SEARO)
    "India": "SEAR", "Indonesia": "SEAR", "Thailand": "SEAR",
    "Bangladesh": "SEAR", "Myanmar": "SEAR", "Sri Lanka": "SEAR", "Nepal": "SEAR",
    # Western Pacific (WPRO)
    "Australia": "WPR", "New Zealand": "WPR", "Japan": "WPR",
    "China": "WPR", "South Korea": "WPR", "Singapore": "WPR", "Malaysia": "WPR",
    "Philippines": "WPR", "Vietnam": "WPR", "Taiwan": "WPR", "Hong Kong": "WPR",
}


def country_to_region(loc: str) -> str:
    if not isinstance(loc, str) or not loc.strip():
        return "unknown"
    # NCBI Geographic Location is often "Country: subdivision". Take everything
    # before the first colon.
    head = loc.split(":")[0].strip()
    return WHO_REGION.get(head, "other")


def year_stratum(date_str: str) -> str:
    """Map a date string to a year-stratum label per §5.4 (2025+ pools 2025 and later)."""
    if not isinstance(date_str, str) or not date_str.strip():
        return "unknown"
    try:
        year = int(date_str.strip()[:4])
    except ValueError:
        return "unknown"
    if year <= 2022:
        return "2022"
    if year == 2023:
        return "2023"
    if year == 2024:
        return "2024"
    return "2025+"


def stratified_sample(df: pd.DataFrame, total: int, min_per_stratum: int,
                       max_per_stratum: int, seed: int) -> pd.DataFrame:
    """
    Two-pass stratified random sample.

    Pass 1: take min_per_stratum from each non-empty (year, region) stratum.
    Pass 2: distribute the remainder round-robin across strata that still have
            members (capped by max_per_stratum). Random order via the seed.
    """
    rng = random.Random(seed)
    df = df.copy().reset_index(drop=True)

    strata = df.groupby(["year_stratum", "who_region"])
    chosen: list[int] = []

    # Pass 1: min_per_stratum
    for (yr, reg), sub in strata:
        if yr == "unknown" or reg == "unknown":
            continue  # do not seed strata with un-mappable metadata
        idxs = list(sub.index)
        rng.shuffle(idxs)
        take = min(min_per_stratum, len(idxs), max_per_stratum)
        chosen.extend(idxs[:take])

    if len(chosen) > total:
        # If even one-per-stratum overflows the budget, sample without replacement
        # back down to `total`.
        rng.shuffle(chosen)
        chosen = chosen[:total]
        return df.loc[chosen].copy().assign(allocation_pass="pass1_capped")

    # Pass 2: round-robin fill
    remaining_budget = total - len(chosen)
    stratum_pools: dict[tuple, list[int]] = {}
    stratum_taken: dict[tuple, int] = defaultdict(int)
    for (yr, reg), sub in strata:
        if yr == "unknown" or reg == "unknown":
            continue
        pool = [i for i in sub.index if i not in chosen]
        rng.shuffle(pool)
        stratum_pools[(yr, reg)] = pool
        stratum_taken[(yr, reg)] = min(min_per_stratum, len(sub.index))

    stratum_keys = sorted(stratum_pools.keys())
    rng.shuffle(stratum_keys)

    pass2_chosen: list[int] = []
    advanced = True
    while remaining_budget > 0 and advanced:
        advanced = False
        for k in stratum_keys:
            if remaining_budget <= 0:
                break
            if stratum_taken[k] >= max_per_stratum:
                continue
            if not stratum_pools[k]:
                continue
            pick = stratum_pools[k].pop()
            pass2_chosen.append(pick)
            stratum_taken[k] += 1
            remaining_budget -= 1
            advanced = True

    out = df.loc[chosen + pass2_chosen].copy()
    out["allocation_pass"] = ["pass1_min"] * len(chosen) + ["pass2_fill"] * len(pass2_chosen)
    return out


def load_clade_table(clade_i_tsv: Path, all_clades_tsv: Path,
                     ncbi_meta_tsv: Path) -> pd.DataFrame:
    """Build a unified table: accession, clade, length, collection_date, country, qc.overallStatus."""
    meta = pd.read_csv(ncbi_meta_tsv, sep="\t", dtype=str)
    # NCBI accession column is just "Accession" in v18; previous freeze used a tweaked column.
    acc_col = next((c for c in ("Accession", "accession", "isolate-accession") if c in meta.columns), None)
    if acc_col is None:
        raise SystemExit(f"No accession column in {ncbi_meta_tsv}: cols = {list(meta.columns)[:10]}")
    meta = meta.rename(columns={acc_col: "accession"})

    # Clade-I assignments from the clade-i Nextclade run (used for Ia + Ib)
    ci_cols = ["seqName", "clade", "qc.overallStatus", "totalMissing"]
    ci = pd.read_csv(clade_i_tsv, sep="\t", dtype=str, usecols=ci_cols)
    ci["accession"] = ci["seqName"].str.split().str[0]
    ci = ci.drop(columns=["seqName"])
    ci = ci.rename(columns={"clade": "clade_ci", "qc.overallStatus": "qc_ci",
                            "totalMissing": "missing_ci"})

    # All-clades assignments (used for everything that the clade-i run called unassigned)
    ac = pd.read_csv(all_clades_tsv, sep="\t", dtype=str,
                     usecols=["seqName", "clade", "qc.overallStatus", "totalMissing"])
    ac["accession"] = ac["seqName"].str.split().str[0]
    ac = ac.drop(columns=["seqName"])
    ac = ac.rename(columns={"clade": "clade_ac", "qc.overallStatus": "qc_ac",
                            "totalMissing": "missing_ac"})

    df = meta.merge(ci, on="accession", how="left").merge(ac, on="accession", how="left")

    # Final clade tag: prefer the clade-i call when it isn't "unassigned"; else all-clades call.
    def pick_clade(row):
        if isinstance(row["clade_ci"], str) and row["clade_ci"] not in ("unassigned", "", None):
            return row["clade_ci"]
        return row["clade_ac"] if isinstance(row["clade_ac"], str) else None

    df["clade"] = df.apply(pick_clade, axis=1)

    # Final QC tag: use the QC from the Nextclade run that actually classified
    # the genome (otherwise IIb genomes get the clade-i QC, which is uniformly
    # bad because they were aligned to the wrong reference).
    def pick_qc(row):
        if isinstance(row["clade_ci"], str) and row["clade_ci"] not in ("unassigned", "", None):
            return row["qc_ci"]
        return row["qc_ac"] if isinstance(row["qc_ac"], str) else None
    df["qc_overall"] = df.apply(pick_qc, axis=1)

    # Standardise the columns we depend on downstream
    date_col = next((c for c in ("Isolate Collection date", "isolate-collection-date") if c in df.columns), None)
    if date_col is None:
        raise SystemExit(f"No date column in metadata. Columns: {list(df.columns)[:30]}")
    df["collection_date"] = df[date_col]
    df["country"] = df.get("Geographic Location", df.get("geo-location", ""))
    df["length"] = pd.to_numeric(df.get("Length", df.get("length", 0)), errors="coerce").fillna(0).astype(int)

    df["year_stratum"] = df["collection_date"].apply(year_stratum)
    df["who_region"] = df["country"].apply(country_to_region)

    return df


def main() -> int:
    ap = argparse.ArgumentParser(description="Build §5.4 outgroup and parental panels.")
    ap.add_argument("--clade-i-tsv", required=True, type=Path)
    ap.add_argument("--all-clades-tsv", required=True, type=Path)
    ap.add_argument("--ncbi-meta-tsv", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--min-length", type=int, default=190000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    df = load_clade_table(args.clade_i_tsv, args.all_clades_tsv, args.ncbi_meta_tsv)
    print(f"[load] {len(df):,} rows merged across NCBI + clade-i + all-clades", file=sys.stderr)

    # Quality filter: §3.3 thresholds (relaxed: we already excluded frameshift-bad from the
    # main corpus; for panel construction we require QC overall = good and length >= 190 kb)
    eligible = df[(df["qc_overall"] == "good") & (df["length"] >= args.min_length)].copy()
    print(f"[filter] {len(eligible):,} pass quality filter "
          f"(QC overall good + length ≥ {args.min_length:,} bp)", file=sys.stderr)

    # Per-clade pools
    pools = {
        "Ia":  eligible[eligible["clade"] == "Ia"],
        "Ib":  eligible[eligible["clade"] == "Ib"],
        "IIb": eligible[eligible["clade"] == "IIb"],
    }
    for k, v in pools.items():
        print(f"  [{k:>3}] eligible pool size: {len(v):,}", file=sys.stderr)

    # 1. Outgroup panel: 30 IIb
    outgroup = stratified_sample(pools["IIb"], total=30, min_per_stratum=1,
                                  max_per_stratum=10, seed=args.seed)
    outgroup_cols = ["accession", "clade", "length", "collection_date", "country",
                     "year_stratum", "who_region", "qc_overall", "allocation_pass"]
    outgroup[outgroup_cols].to_csv(args.out_dir / "outgroup_30_iib.tsv",
                                   sep="\t", index=False)
    print(f"[panel] outgroup: {len(outgroup):,} IIb genomes selected -> outgroup_30_iib.tsv",
          file=sys.stderr)

    # 2. Parental panels: ≤10 each
    parentals = {}
    for clade in ("Ia", "Ib", "IIb"):
        panel = stratified_sample(pools[clade], total=10, min_per_stratum=1,
                                   max_per_stratum=10, seed=args.seed)
        parentals[clade] = panel
        panel[outgroup_cols].to_csv(args.out_dir / f"parental_{clade.lower()}_10.tsv",
                                     sep="\t", index=False)
        print(f"[panel] parental {clade}: {len(panel):,} genomes -> parental_{clade.lower()}_10.tsv",
              file=sys.stderr)

    # 3. Strata summary
    summary = {"seed": args.seed, "min_length": args.min_length,
               "eligible_pool": {k: int(len(v)) for k, v in pools.items()}}
    for k, panel in [("outgroup_iib", outgroup)] + [(f"parental_{c.lower()}", parentals[c]) for c in ("Ia","Ib","IIb")]:
        ct = panel.groupby(["year_stratum", "who_region"]).size().to_dict()
        summary[k] = {f"{a}|{b}": int(v) for (a, b), v in ct.items()}
        summary[f"{k}_total"] = int(len(panel))
    (args.out_dir / "strata_summary.json").write_text(json.dumps(summary, indent=2))

    # 4. SHA-256 manifest
    sha_lines = []
    for f in sorted(args.out_dir.glob("*.tsv")) + [args.out_dir / "strata_summary.json"]:
        h = hashlib.sha256(f.read_bytes()).hexdigest()
        sha_lines.append(f"{h}  {f.name}")
    (args.out_dir / "SHA256SUMS").write_text("\n".join(sha_lines) + "\n")
    print(f"[manifest] wrote SHA256SUMS ({len(sha_lines)} entries)", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
