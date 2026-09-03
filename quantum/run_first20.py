import json
import subprocess
from pathlib import Path

LIGANDS = [
    "165417647",
    "121834274",
    "121833591",
    "10172517",
    "121834393",
    "120998631",
    "23650961",
    "138057076",
    "121834427",
    "121833811",
    "121833786",
    "166261488",
    "167834490",
    "121833380",
    "121834240",
    "121834408",
    "121833141",
    "121833120",
    "121834149",
    "121834331",
]
ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data/reduced_complexes"
OUTPUT = ROOT / "outputs/first20"


def atoms_from_pdb(path):
    atoms = []
    for record in path.read_text().splitlines():
        if not record.startswith(("ATOM", "HETATM")):
            continue
        element = (record[76:78].strip() or record[12:16].strip()[0]).capitalize()
        try:
            coords = tuple(
                float(record[start:end])
                for start, end in ((30, 38), (38, 46), (46, 54))
            )
        except ValueError as exc:
            raise ValueError(f"Invalid coordinates in {path}: {record!r}") from exc
        atoms.append((element, *coords))
    return atoms


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    results = []
    atomic_numbers = {
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
    for ligand in LIGANDS:
        pdb = DATA / f"{ligand}_reduced_template.pdb"
        atoms = atoms_from_pdb(pdb)
        electron_count = sum(atomic_numbers.get(atom[0], 0) for atom in atoms)
        shell = "closed" if electron_count % 2 == 0 else "open"
        result = {
            "ligand": ligand,
            "shell": shell,
            "electron_count": electron_count,
            "pdbqt_source": None,
            "structure_source": str(pdb.relative_to(ROOT)),
        }
        if shell == "closed":
            directory = OUTPUT / ligand
            directory.mkdir(exist_ok=True)
            xyz = directory / "fragment.xyz"
            coordinates = [
                f"{element} {x:.6f} {y:.6f} {z:.6f}" for element, x, y, z in atoms
            ]
            xyz.write_text(
                f"{len(atoms)}\n{ligand} reduced template\n"
                + "\n".join(coordinates)
                + "\n"
            )
            quantum = subprocess.run(
                [
                    "python",
                    "scripts/run_quantum.py",
                    "--xyz",
                    str(xyz),
                    "--outdir",
                    str(directory),
                    "--active-electrons",
                    "4",
                    "--active-orbitals",
                    "4",
                    "--charge",
                    "0",
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
            )
            result["quantum_returncode"] = quantum.returncode
            result["quantum_log"] = quantum.stdout[-1000:] + quantum.stderr[-1000:]
            if quantum.returncode == 0:
                vqe = subprocess.run(
                    [
                        "python",
                        "scripts/run_vqe.py",
                        "--indir",
                        str(directory),
                        "--outdir",
                        str(directory),
                        "--maxiter",
                        "60",
                    ],
                    cwd=ROOT,
                    text=True,
                    capture_output=True,
                )
                result["vqe_returncode"] = vqe.returncode
                result["vqe_log"] = vqe.stdout[-1500:] + vqe.stderr[-1500:]
                summary_file = directory / "vqe_summary.json"
                if summary_file.exists():
                    try:
                        result["vqe"] = json.loads(summary_file.read_text())
                    except (
                        OSError,
                        TypeError,
                        ValueError,
                        json.JSONDecodeError,
                    ) as exc:
                        result["vqe_parse_error"] = str(exc)
        results.append(result)
        print(
            ligand,
            shell,
            result.get("vqe", {}).get("energy_hartree"),
            result.get("quantum_returncode"),
            result.get("vqe_returncode"),
        )
    (OUTPUT / "summary.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
