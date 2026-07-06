import numpy as np
from pathlib import Path


def read_atoms_by_name(path):
    atoms = {}

    with open(path, "r") as file:
        for line in file:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    if name in atoms:
                        raise ValueError(f"Duplicate atom name found: {name}")

                    atoms[name] = np.array([x, y, z], dtype=float)

                except ValueError as e:
                    print(f"Skipping or error in line: {line.strip()}")
                    print(e)

    return atoms


def raw_rmsd(P, Q):
    diff = P - Q
    return np.sqrt((diff * diff).sum() / len(P))


def kabsch_rmsd(P, Q):
    P_centered = P - P.mean(axis=0)
    Q_centered = Q - Q.mean(axis=0)

    C = np.dot(P_centered.T, Q_centered)
    V, S, Wt = np.linalg.svd(C)

    d = np.sign(np.linalg.det(np.dot(V, Wt)))
    D = np.diag([1.0, 1.0, d])

    U = np.dot(np.dot(V, D), Wt)
    P_rotated = np.dot(P_centered, U)

    diff = P_rotated - Q_centered
    return np.sqrt((diff * diff).sum() / len(P))


original_path = Path("prep/ligand_E20.pdb")
docked_path = Path("res/E20_docked_pose1.pdbqt")

original_atoms = read_atoms_by_name(original_path)
docked_atoms = read_atoms_by_name(docked_path)

common_names = sorted(set(original_atoms) & set(docked_atoms))

missing_in_docked = sorted(set(original_atoms) - set(docked_atoms))
extra_in_docked = sorted(set(docked_atoms) - set(original_atoms))

print(f"Original atoms: {len(original_atoms)}")
print(f"Docked atoms: {len(docked_atoms)}")
print(f"Common atoms: {len(common_names)}")

if missing_in_docked:
    print("\nMissing in docked:")
    print(missing_in_docked)

if extra_in_docked:
    print("\nExtra in docked:")
    print(extra_in_docked)

P = np.array([original_atoms[name] for name in common_names])
Q = np.array([docked_atoms[name] for name in common_names])

print("\nCentroids using matched atoms:")
print(f"Original centroid: {P.mean(axis=0)}")
print(f"Docked centroid:   {Q.mean(axis=0)}")
print(f"Centroid distance: {np.linalg.norm(P.mean(axis=0) - Q.mean(axis=0)):.3f} Å")

print("\nRMSD using atom-name matching:")
print(f"Raw RMSD = {raw_rmsd(P, Q):.3f} Å")
print(f"Kabsch RMSD = {kabsch_rmsd(P, Q):.3f} Å")