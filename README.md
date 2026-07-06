# Molecular Docking Redocking Project

This repository documents a complete beginner molecular docking workflow using AutoDock Vina.

The project performs a redocking experiment using a crystallographic protein-ligand complex from the Protein Data Bank. The goal was to test whether AutoDock Vina could reproduce the original binding pose of the ligand.

## Project summary

- Protein structure: `1EVE`
- Target: Acetylcholinesterase
- Ligand: `E20`
- Docking software: AutoDock Vina
- Final redocking RMSD: `0.572 Å`

A redocking RMSD below 2 Å is generally considered a good result for reproducing a crystallographic binding pose. In this project, the best docked pose closely matched the original ligand position.

---

## Objective

The objective of this project was to learn and document a complete molecular docking workflow from scratch:

1. Download a protein-ligand structure from the Protein Data Bank.
2. Identify the ligand and non-protein molecules.
3. Extract the receptor and ligand into separate files.
4. Convert files to `.pdbqt` format.
5. Define a docking search box.
6. Run AutoDock Vina.
7. Visualize the docked pose.
8. Compare the docked pose with the crystallographic ligand.
9. Calculate redocking RMSD.

---

## Repository structure

```text
.
├── README.md
├── data/
│   └── 1EVE.pdb
├── prep/
│   ├── ligand_E20.pdb
│   ├── ligand_E20.pdbqt
│   ├── receptor_1eve.pdb
│   └── receptor_1eve.pdbqt
├── config/
│   └── config.txt
├── res/
│   ├── E20_docked.pdbqt
│   ├── E20_docked_pose1.pdbqt
│   ├── vina_log.txt
│   └── redocking_comparison.png
├── scr/
│   └── calculate_rmsd.py
└── notes/
    ├── 001.md
    ├── 002.md
    ├── 003.md
    ├── 004.md
    ├── 005.md
    ├── 006.md
    ├── 007.md
    ├── 008.md
    ├── 009.md
    └── 010.md
```

---

## Tools used

- AutoDock Vina
- Open Babel
- UCSF ChimeraX
- Python
- Git and GitHub
- macOS Terminal

---

## Input structure

The starting file was:

```text
data/1EVE.pdb
```

This file contains:

- protein atoms
- ligand atoms
- water molecules
- non-standard residues
- metadata

The ligand of interest was identified as:

```text
E20
```

Other non-protein molecules included:

```text
HOH = water
NAG = N-acetylglucosamine
```

---

## Ligand and receptor extraction

The ligand was extracted from the original PDB file using:

```bash
grep "^HETATM" data/1EVE.pdb | grep "E20" > prep/ligand_E20.pdb
```

The receptor was extracted by keeping only standard protein atom records:

```bash
grep "^ATOM" data/1EVE.pdb > prep/receptor_1eve.pdb
```

This created two separate files:

```text
prep/ligand_E20.pdb
prep/receptor_1eve.pdb
```

---

## Ligand center calculation

The ligand center was calculated from the fixed-width coordinate columns of the PDB file.

Final ligand center:

```text
center_x = 2.78114
center_y = 64.3831
center_z = 67.9714
```

These coordinates were used as the center of the AutoDock Vina search box.

Important note: PDB files use fixed-width columns, so coordinate extraction should use column positions instead of simple whitespace splitting.

Example command:

```bash
awk '/^HETATM/ {
  x += substr($0,31,8);
  y += substr($0,39,8);
  z += substr($0,47,8);
  n++
}
END {
  print "center_x =", x/n;
  print "center_y =", y/n;
  print "center_z =", z/n
}' prep/ligand_E20.pdb
```

---

## PDBQT conversion

AutoDock Vina requires `.pdbqt` files.

The ligand was converted using Open Babel:

```bash
obabel prep/ligand_E20.pdb -O prep/ligand_E20.pdbqt
```

The receptor was converted as a rigid receptor using:

```bash
obabel prep/receptor_1eve.pdb -O prep/receptor_1eve.pdbqt -xr
```

The `-xr` option was used to keep the receptor rigid, which is standard for a basic AutoDock Vina workflow.

Open Babel produced an aromatic bond warning during receptor conversion, but the file was created successfully and contained 4254 `ATOM` records.

---

## AutoDock Vina configuration

The docking configuration file was:

```text
config/config.txt
```

Configuration:

```text
receptor = prep/receptor_1eve.pdbqt
ligand = prep/ligand_E20.pdbqt

center_x = 2.78114
center_y = 64.3831
center_z = 67.9714

size_x = 20
size_y = 20
size_z = 20

out = res/E20_docked.pdbqt
```

The docking box was centered around the crystallographic ligand position.

---

## Running AutoDock Vina

Docking was run with:

```bash
vina --config config/config.txt | tee res/vina_log.txt
```

This generated:

```text
res/E20_docked.pdbqt
res/vina_log.txt
```

---

## Docking results

AutoDock Vina produced the following docking scores:

```text
mode | affinity | rmsd l.b. | rmsd u.b.
-----+----------+-----------+----------
1    | -11.07   | 0         | 0
2    | -10.86   | 3.936     | 9.024
3    | -10.43   | 3.618     | 8.795
4    | -10.07   | 4.075     | 9.464
5    | -9.815   | 2.244     | 3.188
6    | -9.696   | 1.601     | 2.760
7    | -9.431   | 2.857     | 4.156
8    | -9.332   | 3.504     | 10.51
9    | -9.313   | 2.382     | 4.249
```

The best predicted pose was mode 1:

```text
Affinity = -11.07 kcal/mol
```

Docking scores are approximate and should not be interpreted as experimental binding affinities.

---

## Pose visualization

The docked pose was visualized in UCSF ChimeraX together with:

```text
prep/receptor_1eve.pdb
prep/ligand_E20.pdb
res/E20_docked.pdbqt
```

The docked ligand was found in the same binding pocket as the original crystallographic ligand.

The comparison image was saved as:

```text
res/redocking_comparison.png
```

---

## RMSD calculation

The best docking pose was extracted from the Vina output:

```bash
awk '/^MODEL 1/{flag=1} flag{print} /^ENDMDL/{if(flag){exit}}' res/E20_docked.pdbqt > res/E20_docked_pose1.pdbqt
```

RMSD was calculated using a Python script:

```text
scr/calculate_rmsd.py
```

The first RMSD attempt gave an incorrect high value because atom order changed during file conversion.

Example mismatch:

```text
original C10 | docked O24
```

To fix this, the script was updated to match atoms by atom name instead of relying on file order.

Final RMSD result:

```text
Original atoms: 28
Docked atoms: 28
Common atoms: 28

Centroid distance: 0.453 Å

Raw RMSD = 0.768 Å
Kabsch RMSD = 0.572 Å
```

---

## Final interpretation

The redocking was successful.

The best docked pose reproduced the crystallographic ligand position with a final Kabsch RMSD of:

```text
0.572 Å
```

Since RMSD values below 2 Å are commonly used as a practical threshold for successful redocking, this result suggests that the docking setup was reasonable for this protein-ligand system.

---

## Key learning points

This project covered several important molecular docking and structural bioinformatics concepts:

- PDB files contain atomic coordinates and metadata.
- `ATOM` records usually correspond to standard protein atoms.
- `HETATM` records often correspond to ligands, water molecules, ions, cofactors, or other non-standard residues.
- PDB files use fixed-width columns.
- Receptor and ligand files must be separated before docking.
- AutoDock Vina requires `.pdbqt` input files.
- The docking box can be centered around a known crystallographic ligand.
- Receptor preparation should usually keep the receptor rigid for basic Vina docking.
- Docking scores are approximate.
- Visual inspection is useful but not enough.
- RMSD calculation requires correct atom mapping.
- Atom order can change during file conversion and produce misleading RMSD values.

---

## Limitations

This was a beginner educational workflow.

Some parts of the workflow were simplified, especially receptor and ligand preparation. A more rigorous docking study could include:

- improved protonation state assignment
- better charge preparation
- validation of ligand chemistry
- comparison with multiple ligands
- testing different docking box sizes
- using specialized receptor preparation tools
- checking protein residues in the binding site
- analyzing molecular interactions
- comparing results with experimental data

---
