"""
Docking.py
Prepares ligands from an SDF file (splits + converts to PDBQT with Meeko),
runs AutoDock Vina docking against a receptor PDBQT, then extracts the
best binding affinity per ligand from the Vina logs and writes a ranked
Excel workbook (Sheet1 = top 3 compounds, Sheet2 = top 1 compound).
"""

import os
import glob
import subprocess
from rdkit import Chem
import pandas as pd

# ---------- CONFIG ----------
SDF_FILE = "compounds.sdf"
LIGAND_DIR = "ligands"
RESULTS_DIR = "results"
RECEPTOR = "protein.pdbqt"

CENTER = (26.238, 6.716, 60.058)   # x, y, z of binding site center
SIZE = (22, 17, 20)                # box size in Angstroms
EXHAUSTIVENESS = 16
NUM_MODES = 10
CPU = 8
# -----------------------------

os.makedirs(LIGAND_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# 1. Split SDF into individual .mol files
suppl = Chem.SDMolSupplier(SDF_FILE, removeHs=False)
mol_files = []

for i, mol in enumerate(suppl):
    if mol is None:
        print(f"Warning: could not parse molecule at index {i}, skipping.")
        continue
    name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"lig_{i}"
    name = name.strip().replace(" ", "_") or f"lig_{i}"
    mol_path = os.path.join(LIGAND_DIR, f"{name}.mol")
    Chem.MolToMolFile(mol, mol_path)
    mol_files.append(mol_path)

print(f"Split {len(mol_files)} molecules from {SDF_FILE}")

# 2. Convert each .mol to .pdbqt using Meeko
pdbqt_files = []
for mol_path in mol_files:
    pdbqt_path = mol_path.replace(".mol", ".pdbqt")
    result = subprocess.run(
        ["mk_prepare_ligand.py", "-i", mol_path, "-o", pdbqt_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"Failed to convert {mol_path}:\n{result.stderr}")
        continue
    pdbqt_files.append(pdbqt_path)

print(f"Converted {len(pdbqt_files)} ligands to PDBQT")

# 3. Run Vina on each ligand
for lig in pdbqt_files:
    name = os.path.splitext(os.path.basename(lig))[0]
    out_file = os.path.join(RESULTS_DIR, f"{name}_out.pdbqt")
    log_file = os.path.join(RESULTS_DIR, f"{name}.log")

    cmd = [
        "vina",
        "--receptor", RECEPTOR,
        "--ligand", lig,
        "--center_x", str(CENTER[0]),
        "--center_y", str(CENTER[1]),
        "--center_z", str(CENTER[2]),
        "--size_x", str(SIZE[0]),
        "--size_y", str(SIZE[1]),
        "--size_z", str(SIZE[2]),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(NUM_MODES),
        "--cpu", str(CPU),
        "--out", out_file,
    ]

    print(f"Docking {name} ...")
    with open(log_file, "w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)

    if result.returncode != 0:
        print(f"  Vina failed for {name}, check {log_file}")
    else:
        print(f"  Done -> {out_file}")

print("\nAll docking runs complete. Check the 'results' folder.")

# 4. Extract top-3 poses per ligand, and full energy breakdown for each
#    ligand's single best pose
def extract_top_modes_from_log(file_path, top_n=3):
    """
    Parses a Vina .log file and returns a list of the top N modes
    (rank, affinity) found in the scoring table, best first.
    """
    modes = []
    found_table_separator = False

    with open(file_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("-----+") and "+" in stripped:
                found_table_separator = True
                continue
            if found_table_separator:
                parts = line.split()
                if len(parts) >= 2 and parts[0].isdigit():
                    try:
                        modes.append({
                            "Mode": int(parts[0]),
                            "Affinity (kcal/mol)": float(parts[1])
                        })
                        if len(modes) >= top_n:
                            break
                    except ValueError:
                        pass
    return modes


def extract_pdbqt_energy_terms(pdbqt_path, model_num=1):
    """
    Parses a specific MODEL block (default: MODEL 1, the top pose) of a
    Vina output PDBQT file and extracts the INTER, INTRA, and
    INTER+INTRA energy terms from the REMARK lines.
    """
    terms = {"INTER+INTRA (kcal/mol)": None, "INTER (kcal/mol)": None, "INTRA (kcal/mol)": None}
    in_target_model = False

    with open(pdbqt_path, "r") as f:
        for line in f:
            if line.startswith("MODEL"):
                in_target_model = line.split()[-1] == str(model_num)
                continue
            if not in_target_model:
                continue
            if "INTER + INTRA:" in line:
                terms["INTER+INTRA (kcal/mol)"] = float(line.split(":")[1].strip())
            elif "INTER:" in line:
                terms["INTER (kcal/mol)"] = float(line.split(":")[1].strip())
            elif "INTRA:" in line:
                terms["INTRA (kcal/mol)"] = float(line.split(":")[1].strip())
            elif line.startswith("ENDMDL"):
                break

    return terms


log_files = sorted(glob.glob(os.path.join(RESULTS_DIR, "*.log")))

if not log_files:
    print("No log files found — skipping Excel output.")
else:
    top3_rows = []       # Sheet 1: 3 rows per ligand
    best_pose_rows = []  # Sheet 2: 1 row per ligand, with energy breakdown

    for log_path in log_files:
        ligand_name = os.path.splitext(os.path.basename(log_path))[0]
        top_modes = extract_top_modes_from_log(log_path, top_n=3)

        if not top_modes:
            print(f"  Warning: no scoring table found in {os.path.basename(log_path)}, skipping.")
            continue

        # Sheet 1 rows: top 3 poses for this ligand
        for mode in top_modes:
            top3_rows.append({
                "Ligand": ligand_name,
                "Mode": mode["Mode"],
                "Affinity (kcal/mol)": mode["Affinity (kcal/mol)"]
            })

        # Sheet 2 row: this ligand's single best pose + energy breakdown
        best_mode = top_modes[0]
        pdbqt_path = os.path.join(RESULTS_DIR, f"{ligand_name}_out.pdbqt")
        row = {
            "Ligand": ligand_name,
            "Best Affinity (kcal/mol)": best_mode["Affinity (kcal/mol)"]
        }
        if os.path.exists(pdbqt_path):
            row.update(extract_pdbqt_energy_terms(pdbqt_path, model_num=1))
        else:
            print(f"  Warning: could not find {pdbqt_path} for energy breakdown.")
        best_pose_rows.append(row)

    if not top3_rows:
        print("No valid Vina scoring tables found — skipping Excel output.")
    else:
        # Sort ligands within Sheet 2 by best affinity (strongest binder first)
        best_pose_rows.sort(key=lambda x: x["Best Affinity (kcal/mol)"])

        df_top3 = pd.DataFrame(top3_rows)
        df_best_pose = pd.DataFrame(best_pose_rows)

        excel_path = os.path.join(RESULTS_DIR, "docking_ranked_results.xlsx")
        with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
            df_top3.to_excel(writer, sheet_name="Top 3 Poses", index=False)
            df_best_pose.to_excel(writer, sheet_name="Best Pose Summary", index=False)

        print(f"\nRanked results written to: {os.path.abspath(excel_path)}")
        print(f"  Sheet 'Top 3 Poses': {len(top3_rows)} row(s) across {len(log_files)} ligand(s)")
        print(f"  Sheet 'Best Pose Summary': {len(best_pose_rows)} ligand(s)")
