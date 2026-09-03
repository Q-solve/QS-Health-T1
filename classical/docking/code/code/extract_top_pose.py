"""
extract_top_pose.py

Extracts just the top pose (MODEL 1 / best-scoring conformation) from each
multi-conformation AutoDock Vina output PDBQT file, and saves it as its own
clean single-model PDBQT file. Useful as input for downstream steps (e.g.
QM/semi-empirical rescoring) that expect a single structure, not an
ensemble of poses.

Also prints the affinity + INTER/INTRA energy terms for the extracted pose
so you can sanity-check which pose was pulled out.
"""

import os
import glob

RESULTS_DIR = "results"
OUTPUT_DIR = "results/top_poses"
MODEL_NUM = 1  # 1 = best-scoring pose


def extract_model(pdbqt_path, model_num=1):
    """
    Returns the full text block (including MODEL/ENDMDL lines) for the
    requested model number from a multi-model PDBQT file, along with any
    REMARK VINA RESULT / INTER / INTRA lines found inside it.
    """
    lines_out = []
    in_target_model = False
    info = {"affinity": None, "inter_intra": None, "inter": None, "intra": None}

    with open(pdbqt_path, "r") as f:
        for line in f:
            if line.startswith("MODEL"):
                current_model = int(line.split()[-1])
                in_target_model = (current_model == model_num)
                if in_target_model:
                    lines_out.append(line)
                continue

            if not in_target_model:
                continue

            lines_out.append(line)

            if "VINA RESULT:" in line:
                parts = line.split(":")[1].split()
                info["affinity"] = float(parts[0])
            elif "INTER + INTRA:" in line:
                info["inter_intra"] = float(line.split(":")[1].strip())
            elif "INTER:" in line:
                info["inter"] = float(line.split(":")[1].strip())
            elif "INTRA:" in line:
                info["intra"] = float(line.split(":")[1].strip())

            if line.startswith("ENDMDL"):
                break

    if not lines_out:
        return None, info

    return "".join(lines_out), info


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    pdbqt_files = sorted(glob.glob(os.path.join(RESULTS_DIR, "*_out.pdbqt")))

    if not pdbqt_files:
        print(f"No '_out.pdbqt' files found in '{RESULTS_DIR}'.")
        return

    print(f"{'Ligand':<20} {'Affinity':>10} {'INTER+INTRA':>13} {'INTER':>10} {'INTRA':>10}")
    print("-" * 68)

    for pdbqt_path in pdbqt_files:
        filename = os.path.basename(pdbqt_path)
        ligand_name = filename.replace("_out.pdbqt", "")

        model_text, info = extract_model(pdbqt_path, model_num=MODEL_NUM)

        if model_text is None:
            print(f"  Warning: MODEL {MODEL_NUM} not found in {filename}, skipping.")
            continue

        out_path = os.path.join(OUTPUT_DIR, f"{ligand_name}_top_pose.pdbqt")
        with open(out_path, "w") as f:
            f.write(model_text)
            if not model_text.rstrip().endswith("END"):
                f.write("END\n")

        aff = info["affinity"] if info["affinity"] is not None else float("nan")
        inter_intra = info["inter_intra"] if info["inter_intra"] is not None else float("nan")
        inter = info["inter"] if info["inter"] is not None else float("nan")
        intra = info["intra"] if info["intra"] is not None else float("nan")

        print(f"{ligand_name:<20} {aff:>10.3f} {inter_intra:>13.3f} {inter:>10.3f} {intra:>10.3f}")

    print(f"\nTop-pose PDBQT files written to: {os.path.abspath(OUTPUT_DIR)}")


if __name__ == "__main__":
    main()
