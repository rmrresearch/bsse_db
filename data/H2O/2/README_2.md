# `2/` — H₂O Dimer Interaction Energy Dataset

This directory contains the full set of quantum-chemistry calculations used to study
**water (H₂O) dimers**, with a focus on analyzing **Basis Set Superposition Error (BSSE)**
as a function of intermolecular geometry. Calculations are performed at the
CCSD(T)/aug-cc-pVDZ level of theory. Both monomers are treated as **rigid** throughout.

---

## 1. Configuration and Geometry Parameters

The geometry of the dimer is described in terms of the relative position and
orientation of two rigid H₂O monomers. Each configuration corresponds to a distinct
intermolecular arrangement, indexed by an integer label.

### Relative Position Vector

The translational displacement between monomer 1 and monomer 2 is defined as

$$\mathbf{r} = \mathbf{r}_2 - \mathbf{r}_1$$

where $\mathbf{r}_1$ and $\mathbf{r}_2$ are the Cartesian coordinates of the
**oxygen atoms** in monomer 1 and monomer 2, respectively. All coordinates are
given in **Ångströms**, which is the default unit in the NWChem input files used here.

Varying $\mathbf{r}$ directly controls the **oxygen–oxygen (O–O) separation**, which
implicitly also modifies the center-of-mass distance between the monomers.

---

## 2. Translational Configuration Dataset

In the translational dataset the **relative orientation of the two monomers is held
fixed**, and the O–O displacement vector $\mathbf{r} = (x, y, z)$ is varied.

Each configuration is labeled as:

```
n_x_y_z
```

where:

- `n` is the configuration index
- `(x, y, z)` are the Cartesian components (in Ångströms) of the O–O displacement
  vector $\mathbf{r}$

### Grid structure

The goal of the translational sampling is to construct **2D contour grids** of the
BSSE as a function of two displacement coordinates, with the third held fixed. This
gives rise to three families of slices:

- **Fixed x** → contour over (y, z)
- **Fixed y** → contour over (x, z)
- **Fixed z** → contour over (x, y)

### Grid completeness

The translational dataset was originally generated without guaranteeing full grid
coverage for every slice. As a result, most 2D slices are incomplete — not due to
convergence failures, but because the original sampling did not cover all required
grid points for a given fixed coordinate. The previous analysis had been producing
contour plots from these partial grids without flagging them as incomplete.

The only slices with full 7×7 grid coverage (49 points each), and therefore the
only ones considered for contour analysis, are the **fixed-y** slices at:

- y = 2.0, 2.25, 2.5, 2.75 Å

All other slices — fixed x and fixed z — have incomplete grids and are excluded
from the contour analysis. A full breakdown of coverage per slice is available by
running `analyze_trnl_grid_coverage.py` in `data_frames_plots/`.

---

## 3. Rotational–Radial Configuration Dataset

In the rotational–radial dataset both the **relative orientation** and the
**intermolecular separation** of the two monomers are varied. Euler angles are
sampled randomly and the O–O distance R is varied independently.

Each configuration is labeled as:

```
n_α_β_γ_R
```

where:

- `n` is the configuration index
- `(α, β, γ)` are Euler angles in degrees describing the rotation of monomer 2
  relative to monomer 1
- `R` is the O–O separation distance in Ångströms

### R-point coverage

Similarly to the translational case, not all rotational orientations have sufficient
R values to construct a meaningful BSSE-vs-R curve. Orientations with fewer than 5
distinct R values are excluded from the analysis. Of the 64 unique orientations
sampled, 32 meet this criterion and are used for plotting. A full breakdown is
available by running `analyze_rot_coverage.py` in `data_frames_plots/`.

---

## 4. Directory Structure and Organization

```
2/
├── translational/            # NWChem input/output files — translational scan
├── rotational_radial/        # NWChem input/output files — rotational/radial scan
├── translational_pymol_xyz/  # XYZ files extracted from translational inputs
├── rotational_pymol_xyz/     # XYZ files extracted from rotational inputs
├── data_frames_plots/        # Parsed CSV data frames, analysis scripts, and plots
├── load_translational.pml    # PyMOL script to load all translational configurations
├── load_rotational.pml       # PyMOL script to load all rotational configurations
└── make_pml.py               # Utility to regenerate the .pml loader files
```

### Configuration subdirectories

Within `translational/` and `rotational_radial/`, each configuration folder is
organized as:

```
<configuration>/<method>/<basis>/
```

where:

- `<configuration>` is one of the translational or rotational–radial labels above
- `<method>` reflects the energy contribution stored by the parser. All calculations
  are run under a single `CCSD_T/` directory; the NWChem output reports total SCF,
  MP2, and CCSD(T) energies, which are decomposed during parsing into:
  - `SCF` — total SCF energy
  - `MP2` — MP2 correlation energy (relative to SCF)
  - `CCSD_T` — CCSD(T) correlation component (relative to MP2)
- `<basis>` is the basis set used (e.g., `aug-cc-pvdz`)

Each `<basis>` directory contains the NWChem input (`.nw`) and output (`.out`) files
for the BSSE components (e.g., `A_A`, `A_AB`, `AB_AB`). Only configurations for which
**all required energy components converge** are retained for BSSE analysis.

### `data_frames_plots/`

Contains the parsed output data in CSV format, the parser scripts used to extract
energies from NWChem output files, and potential energy curve plots. See
`data_frames_plots/README.md` for full usage instructions.

### PyMOL visualization

The `translational_pymol_xyz/` and `rotational_pymol_xyz/` folders contain one XYZ
file per configuration, extracted from the corresponding `AB_AB.nw` input files.
To visualize all configurations in PyMOL, run from this directory:

```bash
pymol load_translational.pml
# or
pymol load_rotational.pml
```

To regenerate the `.pml` loader files if the XYZ folders are updated:

```bash
python make_pml.py
```
