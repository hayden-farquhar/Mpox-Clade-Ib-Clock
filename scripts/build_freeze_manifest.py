#!/usr/bin/env python3
"""
Build the freeze manifest — provenance, SHA-256 hashes, and tool versions.

This script is the reproducibility anchor for the data freeze. The manifest
JSON it produces (and the SHA256SUMS file alongside it) backs every claim
in the pre-registration about exact data identity at preprint time.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def file_record(path: Path, repo_root: Path) -> dict:
    abs_path = path.resolve()
    try:
        rel = abs_path.relative_to(repo_root.resolve())
        path_str = str(rel)
    except ValueError:
        path_str = str(abs_path)
    return {
        "path": path_str,
        "size_bytes": abs_path.stat().st_size,
        "sha256": sha256_of(abs_path),
    }


def directory_record(path: Path, repo_root: Path) -> dict:
    """Hash every file in a directory tree (used for the Nextclade dataset dir)."""
    files = []
    abs_path = path.resolve()
    for p in sorted(abs_path.rglob("*")):
        if p.is_file():
            files.append(file_record(p, repo_root))
    try:
        rel = abs_path.relative_to(repo_root.resolve())
        path_str = str(rel)
    except ValueError:
        path_str = str(abs_path)
    return {
        "path": path_str,
        "type": "directory",
        "file_count": len(files),
        "files": files,
    }


def tool_version(cmd: list[str]) -> str:
    """Run a tool's --version (or equivalent) and return the first line."""
    try:
        out = subprocess.run(
            cmd, capture_output=True, text=True, timeout=10
        ).stdout.strip()
        return out.splitlines()[0] if out else ""
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return "NOT FOUND"


def gather_tool_versions() -> dict:
    tools = {
        "datasets": ["datasets", "--version"],
        "dataformat": ["dataformat", "--version"],
        "nextclade": ["nextclade", "--version"],
        "mafft": ["mafft", "--version"],
        "iqtree": ["iqtree", "--version"],
        "treetime": ["treetime", "--version"],
        "hyphy": ["hyphy", "--version"],
        "python": ["python3", "--version"],
        "snakemake": ["snakemake", "--version"],
    }
    return {name: tool_version(cmd) for name, cmd in tools.items()}


def gather_python_packages() -> dict:
    """Capture installed versions of pinned Python packages."""
    targets = [
        "biopython", "dendropy", "ete3", "numpy", "pandas", "scipy",
        "statsmodels", "matplotlib", "seaborn", "snakemake", "requests",
        "pyyaml", "tqdm",
    ]
    versions = {}
    for pkg in targets:
        try:
            # importlib.metadata works without invoking pip
            from importlib import metadata
            versions[pkg] = metadata.version(pkg)
        except Exception:
            versions[pkg] = "NOT INSTALLED"
    return versions


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description="Build the freeze manifest.")
    ap.add_argument("--freeze-date", required=True)
    ap.add_argument("--raw-dir", required=True, type=Path)
    ap.add_argument("--interim-dir", required=True, type=Path)
    ap.add_argument("--processed-dir", required=True, type=Path)
    ap.add_argument("--out-manifest", required=True, type=Path)
    ap.add_argument("--out-sha256", required=True, type=Path)
    args = ap.parse_args()

    repo_root = Path.cwd()
    now = datetime.now(timezone.utc).isoformat()

    # --- Collect file records ---
    raw_files: list[dict] = []
    for entry in sorted(args.raw_dir.iterdir()):
        if entry.name.startswith("."):
            continue  # skip log files
        if entry.is_dir():
            raw_files.append(directory_record(entry, repo_root))
        else:
            raw_files.append(file_record(entry, repo_root))

    interim_files: list[dict] = [
        file_record(p, repo_root)
        for p in sorted(args.interim_dir.iterdir())
        if p.is_file() and not p.name.startswith(".")
    ]

    processed_files: list[dict] = [
        file_record(p, repo_root)
        for p in sorted(args.processed_dir.iterdir())
        if p.is_file() and not p.name.startswith(".") and p.name not in {"freeze_manifest.json", "SHA256SUMS"}
    ]

    # --- Read filter summary if present (for headline stats) ---
    filter_summary_path = args.processed_dir / "filter_summary.json"
    filter_summary = None
    if filter_summary_path.exists():
        filter_summary = json.loads(filter_summary_path.read_text())

    # --- Compose manifest ---
    manifest = {
        "schema_version": "1.0",
        "freeze_date": args.freeze_date,
        "manifest_built_at_utc": now,
        "machine": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "hostname": platform.node(),
        },
        "tools": gather_tool_versions(),
        "python_packages": gather_python_packages(),
        "data_sources": {
            "ncbi_virus_datasets_cli": {
                "endpoint": "https://api.ncbi.nlm.nih.gov/datasets/v2/",
                "command": "datasets download virus genome taxon 'Monkeypox virus' --include genome,annotation",
                "fetched_at_utc": now,
            },
            "nextstrain_clade_i": {
                "url": "https://data.nextstrain.org/mpox_clade-I.json",
                "fetched_at_utc": now,
            },
            "nextstrain_all_clades": {
                "url": "https://data.nextstrain.org/mpox_all-clades.json",
                "fetched_at_utc": now,
            },
            "nextclade_dataset_primary": {
                "name": "nextstrain/mpox/clade-i",
                "fetched_at_utc": now,
            },
        },
        "raw_files": raw_files,
        "interim_files": interim_files,
        "processed_files": processed_files,
        "filter_summary": filter_summary,
        "preregistration_doi": "10.17605/OSF.IO/CASR2",
        "osf_project_url": "https://osf.io/gt3vx/",
        "reference_accession": "NC_063383",
    }

    args.out_manifest.write_text(json.dumps(manifest, indent=2, default=str))
    print(f"[manifest] written to {args.out_manifest}")

    # --- SHA256SUMS plain-text manifest (for easy diff and download verification) ---
    lines: list[str] = []

    def flatten(recs: list[dict]):
        for r in recs:
            if r.get("type") == "directory":
                for sub in r["files"]:
                    lines.append(f"{sub['sha256']}  {sub['path']}")
            else:
                lines.append(f"{r['sha256']}  {r['path']}")

    flatten(raw_files)
    flatten(interim_files)
    flatten(processed_files)

    args.out_sha256.write_text("\n".join(lines) + "\n")
    print(f"[sha256] {len(lines)} files hashed, written to {args.out_sha256}")

    # --- §3.3 exit-criterion check ---
    if filter_summary:
        passed = filter_summary.get("passed", 0)
        if passed < 30:
            print(
                f"\n!! ATTENTION: surviving corpus is {passed} genomes (<30).\n"
                "   Per pre-registration §3.4, RQ1 / H1 / H2 are dropped if the freeze\n"
                "   contains fewer than 30 date-resolved clade-Ib genomes after filters.\n"
                "   The analysis proceeds as a recombination-scan-only study.\n",
                file=sys.stderr,
            )
        print(
            f"[exit-criterion-check] surviving corpus size: {passed}",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
