"""Input files of the model-based builder, and a manifest that catches when they change. Owner: answer_set.

The answer set depends on its inputs as much as on the model: on 2026-10-01 the k2/107 machine replaced its HG002
NIST VCF between two runs, and the same model gave slightly different answer sets for every GIAB round (every round
reads all seven donor lists). The manifest records sha256, size and mtime of every input:

  python -m answer_set.k2.inputs record          # after installing or deliberately updating inputs
  python -m answer_set.k2.inputs verify [--full] # size + mtime (instant); --full re-hashes everything

Every build checks size + mtime of the inputs it reads (check()) and reports any difference in report.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from answer_set.k2 import kg, masks
from answer_set.k2.util import CACHE, DONORS, donor_bed, donor_vcf, ref_path
from common.settings import SETTINGS

MANIFEST = CACHE / "inputs.json"            # per clone: it records this machine's file times


def files(chrom: Optional[str] = None) -> List[Path]:
    """Every input a build reads from the data dir: for one chromosome, or (None) for all of kg.CHROMS. (The models
    are versioned in git with the code.)"""
    out: List[Path] = []
    for d in DONORS:
        out += [donor_vcf(d), Path(f"{donor_vcf(d)}.tbi"), donor_bed(d)]
    out += [masks.MASK_DIR / name for name in masks.MASKS.values()]
    for c in ([chrom] if chrom else kg.CHROMS):
        fa = ref_path(c)
        out += [fa, Path(f"{fa}.fai"), kg.panel(c), Path(f"{kg.panel(c)}.tbi"), kg.index_path(c),
                kg.GATK_DIR / f"{c}.sites.vcf.gz", kg.GATK_DIR / f"{c}.sites.vcf.gz.tbi"]
    return out


def _key(p: Path, root: Path) -> str:
    return os.path.relpath(p, root)


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def record(paths: Iterable[Path], manifest: Path = MANIFEST, root: Optional[Path] = None) -> int:
    """Write the manifest for `paths` (missing files are skipped). Returns the number of files recorded."""
    root = root or Path(SETTINGS.data_dir)
    entries: Dict[str, dict] = {}
    for p in paths:
        if p.exists():
            st = p.stat()
            entries[_key(p, root)] = {"sha256": _sha256(p), "size": st.st_size, "mtime_ns": st.st_mtime_ns}
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(entries, indent=0, sort_keys=True))
    return len(entries)


def check(paths: Iterable[Path], manifest: Path = MANIFEST, root: Optional[Path] = None, full: bool = False) -> List[str]:
    """Differences between `paths` and the manifest ([] = all as recorded). Size + mtime only, unless full."""
    root = root or Path(SETTINGS.data_dir)
    if not manifest.exists():
        return [f"no input manifest {manifest} (python -m answer_set.k2.inputs record)"]
    entries = json.loads(manifest.read_text())
    problems = []
    for p in paths:
        k, e = _key(p, root), entries.get(_key(p, root))
        if not p.exists():
            problems.append(f"missing: {k}")
        elif e is None:
            problems.append(f"not in the manifest: {k}")
        else:
            st = p.stat()
            if st.st_size != e["size"] or (full and _sha256(p) != e["sha256"]):
                problems.append(f"changed: {k}")
            elif not full and st.st_mtime_ns != e["mtime_ns"]:
                problems.append(f"modified since recorded: {k}")
    return problems


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("action", choices=["record", "verify"])
    ap.add_argument("--full", action="store_true", help="verify: re-hash every file (minutes)")
    args = ap.parse_args()
    if args.action == "record":
        print(f"recorded {record(files())} files in {MANIFEST}")
    else:
        problems = check(files(), full=args.full)
        print("\n".join(problems) or "all inputs as recorded")
        raise SystemExit(1 if problems else 0)


if __name__ == "__main__":
    main()
