#!/usr/bin/env python3

import numpy as np

PDB_FILE = "6A2L.pdb"
LIGAND = "9QO"

coords = []

with open(PDB_FILE) as f:
    for line in f:
        if line.startswith("HETATM") and line[17:20].strip() == LIGAND:
            coords.append([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ])

coords = np.array(coords)

if len(coords) == 0:
    raise ValueError(f"No {LIGAND} atoms found in {PDB_FILE}")

centroid = coords.mean(axis=0)
minimum = coords.min(axis=0)
maximum = coords.max(axis=0)

print("=" * 60)
print("        9QO LIGAND BINDING-SITE INSPECTION")
print("=" * 60)

print(f"\nPDB file: {PDB_FILE}")
print(f"Ligand: {LIGAND}")
print(f"Number of 9QO atoms: {len(coords)}")

print("\n9QO CENTROID")
print(f"X = {centroid[0]:.3f}")
print(f"Y = {centroid[1]:.3f}")
print(f"Z = {centroid[2]:.3f}")

print("\n9QO COORDINATE RANGE")

print(
    f"Minimum: "
    f"X={minimum[0]:.3f}, "
    f"Y={minimum[1]:.3f}, "
    f"Z={minimum[2]:.3f}"
)

print(
    f"Maximum: "
    f"X={maximum[0]:.3f}, "
    f"Y={maximum[1]:.3f}, "
    f"Z={maximum[2]:.3f}"
)

print("\n" + "=" * 60)
print("Inspection complete.")
print("=" * 60)
