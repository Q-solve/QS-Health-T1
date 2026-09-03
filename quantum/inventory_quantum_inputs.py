#!/usr/bin/env python3
"""Inventory structure and docking inputs without modifying them."""

import argparse
import json
import zipfile
from pathlib import Path

EXTENSIONS = {
    ".pdb": "protein/complex structure (direct)",
    ".pdbqt": "docking pose (direct ATOM records)",
    ".sdf": "ligand structure (convert to PDB/XYZ)",
    ".mol": "ligand structure (convert to PDB/XYZ)",
    ".mol2": "ligand structure (convert to PDB/XYZ)",
    ".cif": "structure (convert/parse)",
    ".mmcif": "structure (convert/parse)",
    ".xyz": "quantum geometry (direct)",
    ".csv": "docking score table",
    ".tsv": "docking score table",
    ".xlsx": "docking score table",
}


def xlsx_preview(path):
    try:
        with zipfile.ZipFile(path) as z:
            names = z.namelist()
            sheets = [
                n
                for n in names
                if n.startswith("xl/worksheets/") and n.endswith(".xml")
            ]
            return {
                "worksheets": len(sheets),
                "note": "Use LibreOffice/openpyxl for full table parsing",
            }
    except zipfile.BadZipFile:
        return {"error": "not a valid XLSX archive"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data")
    args = ap.parse_args()
    root = Path(args.data_dir)
    records = []
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        kind = EXTENSIONS.get(p.suffix.lower())
        if kind:
            record = {
                "file": str(p),
                "format": p.suffix.lower(),
                "usable_as": kind,
                "bytes": p.stat().st_size,
            }
            if p.suffix.lower() == ".xlsx":
                record.update(xlsx_preview(p))
            records.append(record)
    print(
        json.dumps(
            {"data_dir": str(root), "files": records, "count": len(records)}, indent=2
        )
    )


if __name__ == "__main__":
    main()
