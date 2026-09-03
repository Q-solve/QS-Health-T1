#!/usr/bin/env python3
"""Compute reduced-fragment interaction energies for completed model templates."""

import argparse
import csv
import subprocess
from pathlib import Path


def parse_atoms(path):
    groups = {"ligand": [], "protein": []}
    selected = {"MET", "PHE", "ILE", "LEU"}
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        element = (line[76:78].strip() or line[12:16].strip()[0]).capitalize()
        try:
            xyz = tuple(float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))
        except ValueError as exc:
            raise ValueError(f"Invalid coordinates in {path}: {line!r}") from exc
        residue = line[17:20].strip()
        if residue in {"UNL", "LIG", "MOL"}:
            groups["ligand"].append((element, *xyz))
        elif residue in selected:
            groups["protein"].append((element, *xyz))
    if not groups["ligand"] or not groups["protein"]:
        raise ValueError(f"Missing ligand or protein contacts in {path}")
    return groups


def write_xyz(path, atoms, label):
    path.write_text(
        f"{len(atoms)}\n{label}\n"
        + "\n".join(f"{e} {x:.6f} {y:.6f} {z:.6f}" for e, x, y, z in atoms)
        + "\n"
    )


def energy(directory, root, maxiter):
    numbers = {
        "H": 1,
        "C": 6,
        "N": 7,
        "O": 8,
        "S": 16,
        "P": 15,
        "F": 9,
        "Cl": 17,
        "Br": 35,
        "I": 53,
    }
    atom_lines = (directory / "fragment.xyz").read_text().splitlines()[2:]
    electrons = sum(
        numbers.get(line.split()[0], 0) for line in atom_lines if line.split()
    )
    spin = electrons % 2
    run = subprocess.run(
        [
            "python",
            "scripts/run_quantum.py",
            "--xyz",
            str(directory / "fragment.xyz"),
            "--outdir",
            str(directory),
            "--active-electrons",
            "4",
            "--active-orbitals",
            "4",
            "--charge",
            "0",
            "--spin",
            str(spin),
        ],
        cwd=root,
        capture_output=True,
        text=True,
    )
    if run.returncode:
        return None, run.stderr[-500:]
    vqe = subprocess.run(
        [
            "python",
            "scripts/run_vqe.py",
            "--indir",
            str(directory),
            "--outdir",
            str(directory),
            "--maxiter",
            str(maxiter),
        ],
        cwd=root,
        capture_output=True,
        text=True,
    )
    if vqe.returncode:
        return None, vqe.stderr[-500:]
    import json

    try:
        summary = json.loads((directory / "vqe_summary.json").read_text())
        return summary["energy_hartree"], ""
    except (OSError, ValueError, KeyError) as exc:
        return None, str(exc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument(
        "--template-root",
        default="outputs/top10_three_residues/templates/reduced_complexes",
    )
    ap.add_argument("--outdir", default="outputs/interaction_energy")
    ap.add_argument("--ligands", nargs="+", required=True)
    ap.add_argument("--maxiter", type=int, default=20)
    args = ap.parse_args()
    root = Path(args.root).resolve()
    templates = root / args.template_root
    out = root / args.outdir
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for ligand in args.ligands:
        groups = parse_atoms(templates / f"{ligand}_reduced_template.pdb")
        dirs = {}
        for kind, atoms in groups.items():
            d = out / ligand / kind
            d.mkdir(parents=True, exist_ok=True)
            write_xyz(d / "fragment.xyz", atoms, f"{ligand} {kind}")
            dirs[kind] = d
        complex_dir = root / "outputs/top10_three_residues" / ligand
        energies = {}
        for kind, directory in [
            ("ligand", dirs["ligand"]),
            ("protein", dirs["protein"]),
            ("complex", complex_dir),
        ]:
            if kind == "complex":
                energies[kind] = None
                try:
                    import json

                    energies[kind] = json.loads(
                        (directory / "vqe_summary.json").read_text()
                    )["energy_hartree"]
                except (OSError, ValueError, KeyError):
                    pass
            else:
                energies[kind], error = energy(directory, root, args.maxiter)
        interaction = None
        if all(energies[k] is not None for k in ("complex", "ligand", "protein")):
            interaction = energies["complex"] - energies["ligand"] - energies["protein"]
        rows.append(
            {
                "ligand": ligand,
                **{
                    f"{k}_energy_hartree": energies[k]
                    for k in ("complex", "ligand", "protein")
                },
                "interaction_energy_hartree": interaction,
                "interaction_energy_kcal_mol": interaction * 627.5095
                if interaction is not None
                else None,
            }
        )
    with (out / "interaction_energy_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
