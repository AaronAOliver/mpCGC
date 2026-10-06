#!/usr/bin/env python3
"""Put original sequence IDs back after dbCAN ran on placeholders.

Reads the map written by mpcgc_rename_ids.py and rewrites the given files, or
every file under the given directories, in place. Placeholders look like
mpcgcP000000001 / mpcgcC000000001 (fixed width), so a single pass of one
regular expression finds every occurrence, including inside derived IDs such
as mpcgcC000000001_17 (a Pyrodigal protein on a renamed contig).

With an empty map (IDs were short enough) the script does nothing.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

TOKEN = re.compile(r"mpcgc[PC]\d{9}")


def load_map(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(path) as fh:
        next(fh, None)
        for ln in fh:
            ph, _, orig = ln.rstrip("\n").partition("\t")
            if ph:
                out[ph] = orig
    return out


def restore_file(path: Path, mapping: dict[str, str]) -> int:
    try:
        text = path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return 0                                   # binary or unreadable: leave as is
    if "mpcgc" not in text:
        return 0
    n = 0

    def sub(m: re.Match) -> str:
        nonlocal n
        orig = mapping.get(m.group(0))
        if orig is None:
            return m.group(0)
        n += 1
        return orig

    new = TOKEN.sub(sub, text)
    if n:
        path.write_text(new, encoding="utf-8")
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--map", required=True, type=Path)
    ap.add_argument("paths", nargs="+", type=Path, help="files or directories to rewrite in place")
    args = ap.parse_args()

    mapping = load_map(args.map)
    if not mapping:
        return 0
    files = []
    for p in args.paths:
        files.extend(sorted(q for q in p.rglob("*") if q.is_file() and not q.is_symlink()) if p.is_dir() else [p])
    total = sum(restore_file(f, mapping) for f in files)

    # anything left over would mean a placeholder escaped the map
    left = [f for f in files if f.is_file() and TOKEN.search(f.read_text(encoding="utf-8", errors="ignore"))]
    if left:
        sys.exit(f"[restore] placeholders remain in: {', '.join(str(f) for f in left[:5])}")
    print(f"[restore] {total:,} placeholder IDs restored across {len(files)} files", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
