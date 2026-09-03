#!/usr/bin/env python3
import argparse
import json
import math
from pathlib import Path


def parse_atom(line):
    try:
        element = line[76:78].strip() or line[12:14].strip()[0]
        xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    except (IndexError, ValueError) as exc:
        raise ValueError(f"Invalid coordinate record: {line!r}") from exc
    return {
        "element": element.capitalize(),
        "xyz": xyz,
        "name": line[12:16].strip(),
        "res": line[17:20].strip(),
        "resid": line[22:26].strip(),
    }


def unit_difference(first, second):
    vector = [first[i] - second[i] for i in range(3)]
    norm = math.sqrt(sum(value * value for value in vector))
    if norm == 0:
        raise ValueError("Cannot cap a zero-length boundary bond")
    return vector, norm


def main():
    parser = argparse.ArgumentParser(
        description="Build a capped ligand/contact fragment from a reduced PDB"
    )
    parser.add_argument("--template", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--residues", nargs="+", default=["MET", "PHE", "ILE"])
    args = parser.parse_args()
    source = Path(args.template)
    if not source.is_file():
        parser.error(f"template does not exist: {source}")
    ligand = []
    residues = []
    for line in source.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            current = parse_atom(line)
            if current["res"] in ("UNL", "LIG", "MOL"):
                ligand.append(current)
            elif current["res"] in ("MET", "PHE", "ILE", "LEU"):
                residues.append(current)
    if not ligand:
        parser.error("no ligand records with residue UNL/LIG/MOL found")
    requested = {name.upper() for name in args.residues}
    valid = {"MET", "PHE", "ILE", "LEU"}
    if requested - valid:
        parser.error(f"unsupported residues: {', '.join(sorted(requested - valid))}")
    selected = ligand + [item for item in residues if item["res"] in requested]
    for residue in sorted({item["res"] for item in selected if item["res"] != "UNL"}):
        cb = next(
            (
                item
                for item in selected
                if item["res"] == residue and item["name"] == "CB"
            ),
            None,
        )
        cg = next(
            (
                item
                for item in selected
                if item["res"] == residue and item["name"] in ("CG", "CG1")
            ),
            None,
        )
        if cb and cg:
            vector, norm = unit_difference(cb["xyz"], cg["xyz"])
            cap_xyz = [cb["xyz"][i] + 1.09 * vector[i] / norm for i in range(3)]
            selected.append(
                {
                    "element": "H",
                    "xyz": cap_xyz,
                    "name": "Hcap",
                    "res": residue,
                    "resid": "cap",
                }
            )
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    model = "ligand + " + " + ".join(sorted(requested))
    with (outdir / "fragment.xyz").open("w") as handle:
        handle.write(f"{len(selected)}\n{model}\n")
        for item in selected:
            x, y, z = item["xyz"]
            handle.write(f"{item['element']} {x:.6f} {y:.6f} {z:.6f}\n")
    provenance = {
        "source": str(source),
        "real_data": True,
        "model": model,
        "atom_count": len(selected),
        "excluded_residues": sorted(valid - requested),
    }
    (outdir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
