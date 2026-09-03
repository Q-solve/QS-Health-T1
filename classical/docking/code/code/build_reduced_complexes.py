"""
build_reduced_complexes.py

Takes a "template" PDB that contains:
  1) A ligand block with residue name UNL (placeholder coordinates), and
  2) A fixed set of protein side-chain atoms (the pocket residues)

For every ligand's top-pose file in POSES_DIR (PDB or PDBQT), this:
  - strips out the old UNL block from the template
  - converts that ligand's own docked atoms into UNL-labeled ATOM lines
  - keeps the protein side-chain atoms from the template unchanged
  - writes one combined "reduced complex" PDB per ligand

Pure Python -- no external dependencies.
"""

import os
import glob

# ---------- CONFIG ----------
TEMPLATE_FILE = "template.pdb"          # protein side chains + one UNL ligand block
POSES_DIR = "results/top_poses"         # folder with per-ligand top-pose files
POSES_GLOB = "*_top_pose.pdbqt"         # pattern matching your top-pose files (adjust if .pdb)
OUTPUT_DIR = "results/reduced_complexes"
LIGAND_RESNAME_IN_TEMPLATE = "UNL"
# -----------------------------

# Map common AutoDock/Vina atom types (as seen in the last PDBQT column)
# to a standard one/two-letter element symbol.
ATOM_TYPE_TO_ELEMENT = {
    "A": "C", "C": "C",
    "N": "N", "NA": "N", "N1+": "N", "NS": "N",
    "OA": "O", "O": "O", "OS": "O",
    "SA": "S", "S": "S",
    "HD": "H", "H": "H", "HS": "H",
    "F": "F", "Cl": "Cl", "Br": "Br", "I": "I",
    "P": "P",
    "Mg": "Mg", "Ca": "Ca", "Mn": "Mn", "Fe": "Fe", "Zn": "Zn",
}


def guess_element(atom_name, atom_type=None):
    """Best-effort element guess from an AutoDock atom type or PDB atom name."""
    if atom_type:
        t = atom_type.strip()
        if t in ATOM_TYPE_TO_ELEMENT:
            return ATOM_TYPE_TO_ELEMENT[t]
    # fall back to first alphabetic character(s) of the atom name
    name = atom_name.strip().lstrip("0123456789")
    for elem in ("Cl", "Br", "Mg", "Ca", "Mn", "Fe", "Zn"):
        if name.upper().startswith(elem.upper()):
            return elem
    return name[0].upper() if name else "C"


def parse_template(template_file):
    """
    Splits the template into (protein_lines, ligand_atom_count) --
    protein_lines are the raw ATOM lines for everything that is NOT the
    UNL ligand block, preserved exactly as-is.
    """
    protein_lines = []
    with open(template_file, "r") as f:
        for line in f:
            record = line[0:6].strip()
            if record not in ("ATOM", "HETATM"):
                continue
            resname = line[17:20].strip()
            if resname == LIGAND_RESNAME_IN_TEMPLATE:
                continue  # skip old ligand block
            protein_lines.append(line.rstrip("\n"))
    return protein_lines


def parse_ligand_pose(pose_file):
    """
    Parses ATOM/HETATM lines from a ligand pose file (PDB or PDBQT).
    Returns a list of dicts with atom name, coord, and element.
    """
    atoms = []
    with open(pose_file, "r") as f:
        for line in f:
            record = line[0:6].strip()
            if record not in ("ATOM", "HETATM"):
                continue
            try:
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            # PDBQT files carry the AutoDock atom type in the last
            # whitespace-separated column; PDB files may have an element
            # symbol in columns 77-78. Try both, fall back to atom name.
            atom_type = None
            tail = line[70:].split()
            if tail:
                atom_type = tail[-1]

            element = guess_element(atom_name, atom_type)

            atoms.append({"name": atom_name, "coord": (x, y, z), "element": element})

    return atoms


def format_atom_line(serial, atom_name, resname, chain, resnum, coord, element):
    """Formats a single ATOM line in fixed-column PDB format."""
    x, y, z = coord
    name_field = atom_name if len(atom_name) >= 4 else f" {atom_name:<3}"
    return (
        f"ATOM  {serial:>5} {name_field:<4} {resname:<3} {chain:>1}{resnum:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{0.00:>6.2f}"
        f"          {element:>2}"
    )


def build_complex(ligand_name, ligand_atoms, protein_lines, output_path):
    lines_out = []
    serial = 1

    for atom in ligand_atoms:
        lines_out.append(
            format_atom_line(serial, atom["name"], LIGAND_RESNAME_IN_TEMPLATE, "",
                              1, atom["coord"], atom["element"])
        )
        serial += 1

    lines_out.append("TER")

    # Re-serial the protein lines too, continuing the numbering, but keep
    # their original atom name/resname/chain/resnum fields untouched.
    for line in protein_lines:
        new_serial = f"{serial:>5}"
        new_line = line[0:6] + new_serial + line[11:]
        lines_out.append(new_line)
        serial += 1

    lines_out.append("TER")
    lines_out.append("END")

    with open(output_path, "w") as f:
        f.write("\n".join(lines_out) + "\n")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    protein_lines = parse_template(TEMPLATE_FILE)
    print(f"Loaded {len(protein_lines)} protein side-chain atoms from template.")

    pose_files = sorted(glob.glob(os.path.join(POSES_DIR, POSES_GLOB)))
    if not pose_files:
        print(f"No pose files found in '{POSES_DIR}' matching '{POSES_GLOB}'.")
        return

    for pose_file in pose_files:
        filename = os.path.basename(pose_file)
        ligand_name = filename.replace("_top_pose.pdbqt", "").replace("_top_pose.pdb", "")

        ligand_atoms = parse_ligand_pose(pose_file)
        if not ligand_atoms:
            print(f"  Warning: no atoms parsed from {filename}, skipping.")
            continue

        output_path = os.path.join(OUTPUT_DIR, f"{ligand_name}_reduced_template.pdb")
        build_complex(ligand_name, ligand_atoms, protein_lines, output_path)

        print(f"  {ligand_name}: {len(ligand_atoms)} ligand atoms + "
              f"{len(protein_lines)} protein atoms -> {output_path}")

    print(f"\nAll reduced complexes written to: {os.path.abspath(OUTPUT_DIR)}")


if __name__ == "__main__":
    main()
