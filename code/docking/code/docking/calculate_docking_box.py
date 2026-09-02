import numpy as np

coords = []

with open("6A2L.pdb") as f:
    for line in f:
        if (
            line.startswith("HETATM")
            and line[17:20].strip() == "9QO"
            and line[21].strip() == "A"
        ):
            coords.append([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54])
            ])

coords = np.array(coords)

minimum = coords.min(axis=0)
maximum = coords.max(axis=0)
size = maximum - minimum

print("9QO Chain A")
print("Number of atoms:", len(coords))
print("Minimum:", minimum)
print("Maximum:", maximum)
print("Ligand dimensions:", size)
print("Recommended box with 5 A padding:", size + 10)
