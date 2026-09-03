#!/usr/bin/env python3

from collections import Counter, defaultdict
import sys
import numpy as np

# Allow PDB file to be supplied as an argument
PDB_FILE = sys.argv[1] if len(sys.argv) > 1 else "6A2L.pdb"


def get_xyz(line):
    return np.array([
        float(line[30:38]),
        float(line[38:46]),
        float(line[46:54])
    ])


# --------------------------------------------------
# Read PDB
# --------------------------------------------------

atoms = []
hetatoms = []

with open(PDB_FILE) as f:
    for line in f:
        record = line[0:6].strip()

        if record == "ATOM":
            atoms.append(line)

        elif record == "HETATM":
            hetatoms.append(line)


print("=" * 65)
print("             6A2L PfDHFR-TS STRUCTURE REPORT")
print("=" * 65)

# --------------------------------------------------
# 1. Protein chains
# --------------------------------------------------

chains = Counter(line[21].strip() for line in atoms)

print("\n[1] PROTEIN CHAINS")

for chain, count in sorted(chains.items()):
    print(f"    Chain {chain}: {count} atoms")


# --------------------------------------------------
# 2. Hetero compounds
# --------------------------------------------------

hetero_residues = Counter(
    line[17:20].strip()
    for line in hetatoms
)

print("\n[2] HETERO COMPOUNDS")

for residue, count in hetero_residues.items():
    print(f"    {residue}: {count} atoms")


# --------------------------------------------------
# 3. 9QO by chain
# --------------------------------------------------

print("\n[3] 9QO LIGAND ANALYSIS")

ligand_by_chain = defaultdict(list)

for line in hetatoms:

    residue = line[17:20].strip()

    if residue == "9QO":

        chain = line[21].strip()

        ligand_by_chain[chain].append(
            get_xyz(line)
        )


for chain, coords_list in sorted(
    ligand_by_chain.items()
):

    coords = np.array(coords_list)

    centroid = coords.mean(axis=0)
    minimum = coords.min(axis=0)
    maximum = coords.max(axis=0)

    print(f"\n    Chain {chain} 9QO")

    print(f"    Atoms: {len(coords)}")

    print(
        f"    Centroid: "
        f"X={centroid[0]:.3f}, "
        f"Y={centroid[1]:.3f}, "
        f"Z={centroid[2]:.3f}"
    )

    print(
        f"    Minimum: "
        f"X={minimum[0]:.3f}, "
        f"Y={minimum[1]:.3f}, "
        f"Z={minimum[2]:.3f}"
    )

    print(
        f"    Maximum: "
        f"X={maximum[0]:.3f}, "
        f"Y={maximum[1]:.3f}, "
        f"Z={maximum[2]:.3f}"
    )


# --------------------------------------------------
# 4. All ligand/cofactor centroids
# --------------------------------------------------

print("\n[4] HETERO-MOLECULE CENTROIDS")

for residue in ["9QO", "NAP", "UMP"]:

    groups = defaultdict(list)

    for line in hetatoms:

        if line[17:20].strip() == residue:

            chain = line[21].strip()

            groups[chain].append(
                get_xyz(line)
            )

    print(f"\n    {residue}")

    for chain, coords_list in sorted(
        groups.items()
    ):

        coords = np.array(coords_list)

        centroid = coords.mean(axis=0)

        print(
            f"    Chain {chain}: "
            f"{len(coords)} atoms | "
            f"centroid = "
            f"({centroid[0]:.3f}, "
            f"{centroid[1]:.3f}, "
            f"{centroid[2]:.3f})"
        )


# --------------------------------------------------
# 5. Nearby residues around Chain A 9QO
# --------------------------------------------------

print("\n[5] RESIDUES WITHIN 5 Å OF CHAIN A 9QO")

if "A" in ligand_by_chain:

    ligand_centroid = np.array(
        ligand_by_chain["A"]
    ).mean(axis=0)

    nearby = {}

    for line in atoms:

        if line[21].strip() != "A":
            continue

        xyz = get_xyz(line)

        distance = np.linalg.norm(
            xyz - ligand_centroid
        )

        if distance <= 5.0:

            residue_number = line[22:26].strip()
            residue_name = line[17:20].strip()

            key = (
                residue_name,
                residue_number
            )

            if key not in nearby:

                nearby[key] = distance

            else:

                nearby[key] = min(
                    nearby[key],
                    distance
                )

    for (
        residue_name,
        residue_number
    ), distance in sorted(
        nearby.items(),
        key=lambda x: x[1]
    ):

        print(
            f"    {residue_name} "
            f"{residue_number}: "
            f"{distance:.2f} Å"
        )


# --------------------------------------------------
# 6. PDB metadata
# --------------------------------------------------

print("\n[6] STRUCTURE METADATA")

with open(PDB_FILE) as f:

    for line in f:

        if line.startswith("TITLE"):

            print("    " + line.strip())

        elif line.startswith("EXPDTA"):

            print("    " + line.strip())

        elif line.startswith("COMPND"):

            print("    " + line.strip())


# --------------------------------------------------
# 7. Water count
# --------------------------------------------------

water_count = sum(
    1
    for line in hetatoms
    if line[17:20].strip() == "HOH"
)

print("\n[7] WATER")

print(f"    Water atoms: {water_count}")


# --------------------------------------------------
# 8. Summary
# --------------------------------------------------

print("\n" + "=" * 65)
print("SUMMARY")
print("=" * 65)

print(f"Protein atoms : {len(atoms)}")
print(f"HETATM atoms  : {len(hetatoms)}")
print(
    f"Chains        : "
    f"{', '.join(sorted(chains))}"
)
print(
    f"9QO atoms     : "
    f"{hetero_residues.get('9QO', 0)}"
)
print(
    f"NAP atoms     : "
    f"{hetero_residues.get('NAP', 0)}"
)
print(
    f"UMP atoms     : "
    f"{hetero_residues.get('UMP', 0)}"
)
print(
    f"HOH atoms     : "
    f"{hetero_residues.get('HOH', 0)}"
)

print("=" * 65)
print("Inspection complete.")
print("=" * 65)
